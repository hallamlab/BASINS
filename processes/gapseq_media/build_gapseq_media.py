#!/usr/bin/env python3
"""Build data-trained, gapseq-compatible media for BASINS environmental groups."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import shutil
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


KEY_CANDIDATES = ("cruise_year_month_depth", "sample_id")
CRUISE_CANDIDATES = ("Cruise", "cruise")
DEPTH_CANDIDATES = ("Depth", "depth", "Depth_m", "depth_m")
RESP_RE = re.compile(r"^resp_(\d+)$")
SAFE_RE = re.compile(r"[^A-Za-z0-9._-]+")
DEFAULT_BASAL_MEDIUM = Path(__file__).with_name("default_basal_medium.csv")


@dataclass(frozen=True)
class Family:
    name: str
    labels: tuple[str, ...]
    weights: pd.DataFrame


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", required=True, type=Path)
    parser.add_argument("--cruise-groups", required=True, type=Path)
    parser.add_argument("--o2-assignments", required=True, type=Path)
    parser.add_argument("--gmm-assignments", required=True, type=Path)
    parser.add_argument("--hybrid-assignments", required=True, type=Path)
    parser.add_argument("--nutrients", required=True, type=Path)
    parser.add_argument(
        "--template",
        type=Path,
        help=(
            "Optional shared basal-medium CSV. When omitted, BASINS uses its "
            "bundled minimal glucose basal medium."
        ),
    )
    parser.add_argument("--mapping", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--method", choices=("ordinal", "binary"), default="ordinal")
    parser.add_argument("--min-effective-n", type=float, default=1.0)
    parser.add_argument("--detection-limit", type=float, default=0.0)
    parser.add_argument(
        "--depth-baseline-m",
        default="",
        help=(
            "Optional comma-separated anchored depths for conventional depth "
            "baseline media, for example 100,150,200."
        ),
    )
    return parser.parse_args()


def read_table(path: Path) -> pd.DataFrame:
    sep = "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","
    return pd.read_csv(path, sep=sep)


def find_col(frame: pd.DataFrame, candidates: tuple[str, ...], source: str) -> str:
    for col in candidates:
        if col in frame.columns:
            return col
    raise ValueError(f"{source} lacks any of the required columns: {', '.join(candidates)}")


def norm_text(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    return re.sub(r"\.0$", "", text)


def sanitize(value: str) -> str:
    cleaned = SAFE_RE.sub("_", value.strip()).strip("._-")
    return cleaned or "unlabelled"


def parse_depths(value: str) -> tuple[float, ...]:
    depths = []
    for item in value.split(","):
        item = item.strip()
        if not item:
            continue
        depth = float(item)
        if not np.isfinite(depth) or depth < 0:
            raise ValueError("--depth-baseline-m values must be finite and nonnegative")
        depths.append(depth)
    return tuple(dict.fromkeys(depths))


def depth_label(depth: float) -> str:
    text = f"{depth:g}".replace(".", "p")
    return f"depth_{text}m"


def load_nutrients(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", dtype=str)
    required = {"id", "name"}
    if not required.issubset(frame.columns):
        raise ValueError(f"Nutrient vocabulary lacks columns {sorted(required)}: {path}")
    if frame["id"].duplicated().any():
        raise ValueError(f"Nutrient vocabulary contains duplicate ids: {path}")
    return frame


def load_mapping(path: Path, matrix: pd.DataFrame, nutrients: pd.DataFrame) -> pd.DataFrame:
    mapping = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    required = {"chemistry_column", "compound_id", "compound_name"}
    if not required.issubset(mapping.columns):
        raise ValueError(f"Mapping lacks columns: {sorted(required - set(mapping.columns))}")
    if mapping["chemistry_column"].duplicated().any():
        raise ValueError("Mapping chemistry_column values must be unique")
    known = set(nutrients["id"])
    invalid = sorted(set(mapping["compound_id"]) - known)
    if invalid:
        raise ValueError(f"Mapping contains compound ids absent from nutrients.tsv: {invalid}")
    mapping["available_in_matrix"] = mapping["chemistry_column"].isin(matrix.columns)
    return mapping


def load_template(
    path: Path, nutrients: pd.DataFrame, source: str
) -> tuple[pd.DataFrame, dict[str, object]]:
    if not path.is_file():
        raise ValueError(f"Basal-medium template not found: {path}")
    known = set(nutrients["id"])
    frame = pd.read_csv(path)
    if list(frame.columns) != ["compounds", "name", "maxFlux"]:
        raise ValueError(f"Template must have compounds,name,maxFlux columns in order: {path}")
    frame["compounds"] = frame["compounds"].astype(str)
    frame["maxFlux"] = pd.to_numeric(frame["maxFlux"], errors="raise")
    if (~np.isfinite(frame["maxFlux"]) | (frame["maxFlux"] < 0)).any():
        raise ValueError(f"Template maxFlux values must be finite and nonnegative: {path}")
    invalid = sorted(set(frame["compounds"]) - known)
    if invalid:
        raise ValueError(f"Template {path.name} has ids absent from nutrients.tsv: {invalid}")
    if frame["compounds"].duplicated().any():
        raise ValueError(f"Template has duplicate compound ids: {path}")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    audit = {
        "template": path.stem,
        "source": source,
        "path": str(path.resolve()),
        "sha256": digest,
        "n_compounds": len(frame),
        "n_positive_flux": int((frame["maxFlux"] > 0).sum()),
    }
    return frame, audit


def align_by_key(matrix: pd.DataFrame, assignment: pd.DataFrame, source: str) -> pd.DataFrame:
    key_m = find_col(matrix, KEY_CANDIDATES, "matrix")
    key_a = find_col(assignment, KEY_CANDIDATES, source)
    if assignment[key_a].map(norm_text).duplicated().any():
        raise ValueError(f"{source} contains duplicate sample keys")
    lookup = assignment.copy()
    lookup["_join_key"] = lookup[key_a].map(norm_text)
    keys = matrix[key_m].map(norm_text)
    result = pd.DataFrame({"_join_key": keys}).merge(
        lookup.drop(columns=[key_a] if key_a != "_join_key" else []),
        on="_join_key",
        how="left",
        validate="many_to_one",
    )
    result.index = matrix.index
    return result


def assignment_family(matrix: pd.DataFrame, assignment: pd.DataFrame, name: str) -> Family:
    aligned = align_by_key(matrix, assignment, name)
    resp_cols = sorted(
        (col for col in aligned.columns if RESP_RE.match(col)),
        key=lambda col: int(RESP_RE.match(col).group(1)),
    )
    if not resp_cols:
        raise ValueError(f"{name} assignments contain no resp_N columns")
    weights = aligned[resp_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0).clip(lower=0.0)
    row_sums = weights.sum(axis=1)
    weights = weights.div(row_sums.where(row_sums > 0, 1.0), axis=0)
    if name == "legacy_o2":
        canonical = ("oxic", "dysoxic", "suboxic", "anoxic")
        labels = canonical[: len(resp_cols)] if len(resp_cols) <= len(canonical) else tuple(f"o2_{i}" for i in range(len(resp_cols)))
    elif name == "hybrid":
        # Hybrid builder orders the outer product C-major then G. Deriving the
        # full grid avoids losing names for cells that are never the hard argmax.
        if len(resp_cols) % 4 == 0:
            gmm_k = len(resp_cols) // 4
            labels = tuple(
                f"hybrid_c{component // gmm_k}_g{component % gmm_k}"
                for component in range(len(resp_cols))
            )
        else:
            labels = tuple(f"hybrid_{i}" for i in range(len(resp_cols)))
    else:
        labels = tuple(f"gmm_{i}" for i in range(len(resp_cols)))
    weights.columns = labels
    return Family(name, labels, weights)


def cruise_group_family(matrix: pd.DataFrame, assignments: pd.DataFrame) -> Family:
    cruise_m = find_col(matrix, CRUISE_CANDIDATES, "matrix")
    cruise_a = find_col(assignments, CRUISE_CANDIDATES, "cruise groups")
    resp_cols = sorted(
        (col for col in assignments.columns if RESP_RE.match(col)),
        key=lambda col: int(RESP_RE.match(col).group(1)),
    )
    if not resp_cols:
        raise ValueError("Cruise-group assignments contain no resp_N columns")
    if assignments[cruise_a].map(norm_text).duplicated().any():
        raise ValueError("Cruise-group assignments contain duplicate cruises")
    lookup = assignments.copy()
    lookup["_cruise_key"] = lookup[cruise_a].map(norm_text)
    aligned = pd.DataFrame({"_cruise_key": matrix[cruise_m].map(norm_text)}).merge(
        lookup[["_cruise_key", *resp_cols]], on="_cruise_key", how="left", validate="many_to_one"
    )
    weights = aligned[resp_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0).clip(lower=0.0)
    row_sums = weights.sum(axis=1)
    weights = weights.div(row_sums.where(row_sums > 0, 1.0), axis=0)
    labels = tuple(f"cruise_group_{i + 1}" for i in range(len(resp_cols)))
    weights.columns = labels
    weights.index = matrix.index
    return Family("cruise_group", labels, weights)


def weighted_values(values: pd.Series, weights: pd.Series) -> tuple[np.ndarray, np.ndarray, int]:
    numeric = pd.to_numeric(values, errors="coerce")
    negative = int((numeric < 0).sum())
    valid = numeric.notna() & (numeric >= 0) & weights.notna() & (weights > 0)
    return numeric[valid].to_numpy(float), weights[valid].to_numpy(float), negative


def weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    return float(np.average(values, weights=weights)) if len(values) and weights.sum() > 0 else math.nan


def transform(value: float, positives: np.ndarray, method: str, detection_limit: float) -> tuple[float, str]:
    if not np.isfinite(value) or value <= detection_limit:
        return 0.0, "absent_or_below_detection"
    if method == "binary":
        return 10.0, "detected"
    positives = positives[np.isfinite(positives) & (positives > detection_limit)]
    if not len(positives) or np.allclose(positives, positives[0]):
        return 1.0, "detected_no_variation"
    q25, q75 = np.quantile(positives, [0.25, 0.75])
    if value <= q25:
        return 0.1, "low_le_q25"
    if value <= q75:
        return 1.0, "moderate_q25_q75"
    return 10.0, "high_gt_q75"


def add_medium(
    *,
    matrix: pd.DataFrame,
    mapping: pd.DataFrame,
    weights: pd.Series,
    scope: str,
    cruise_group: str,
    family: str,
    compartment: str,
    basal_medium: pd.DataFrame,
    template_name: str,
    template_source: str,
    outdir: Path,
    method: str,
    min_effective_n: float,
    detection_limit: float,
) -> tuple[list[dict], dict]:
    medium = basal_medium.copy()
    medium["maxFlux"] = medium["maxFlux"].astype(float)
    medium = medium.set_index("compounds", drop=False)
    effective_n = float(weights.sum())
    provenance = []
    for row in mapping.itertuples(index=False):
        column = row.chemistry_column
        if column not in matrix.columns:
            provenance.append({
                "scope": scope, "cruise_group": cruise_group, "family": family,
                "compartment": compartment, "compound_id": row.compound_id,
                "compound_name": row.compound_name, "source_column": column,
                "measurement_unit": getattr(row, "measurement_unit", ""),
                "weighted_mean": "", "effective_n": effective_n,
                "measured_observations_n": 0,
                "measured_membership_weight": 0.0,
                "measured_membership_fraction": 0.0,
                "negative_values_rejected": 0,
                "transform_method": method, "transform_bin": "unmapped_matrix_column",
                "maxFlux": "", "value_source": "template_or_absent",
                "template": template_name, "template_source": template_source,
            })
            continue
        values, valid_weights, n_negative = weighted_values(matrix[column], weights)
        estimate = weighted_mean(values, valid_weights)
        source = "weighted_observations"
        if not np.isfinite(estimate):
            global_numeric = pd.to_numeric(matrix[column], errors="coerce")
            global_numeric = global_numeric[(global_numeric >= 0) & global_numeric.notna()]
            estimate = float(global_numeric.median()) if len(global_numeric) else math.nan
            source = "global_median_fallback"
        positives = pd.to_numeric(matrix[column], errors="coerce").to_numpy(float)
        flux, bin_name = transform(estimate, positives, method, detection_limit)
        nutrient_name = row.compound_name
        if row.compound_id in medium.index:
            nutrient_name = medium.at[row.compound_id, "name"]
            medium.at[row.compound_id, "maxFlux"] = flux
        else:
            medium.loc[row.compound_id] = [row.compound_id, nutrient_name, flux]
        provenance.append({
            "scope": scope, "cruise_group": cruise_group, "family": family,
            "compartment": compartment, "compound_id": row.compound_id,
            "compound_name": nutrient_name, "source_column": column,
            "measurement_unit": getattr(row, "measurement_unit", ""),
            "weighted_mean": estimate, "effective_n": effective_n,
            "measured_observations_n": len(values),
            "measured_membership_weight": float(valid_weights.sum()),
            "measured_membership_fraction": (
                float(valid_weights.sum() / effective_n)
                if effective_n > 0 else math.nan
            ),
            "negative_values_rejected": n_negative, "transform_method": method,
            "transform_bin": bin_name, "maxFlux": flux, "value_source": source,
            "template": template_name, "template_source": template_source,
        })
    medium["maxFlux"] = pd.to_numeric(medium["maxFlux"], errors="raise")
    recipe_text = "\n".join(
        f"{row.compounds}\t{float(row.maxFlux):.12g}"
        for row in medium.sort_index().itertuples(index=False)
    )
    recipe_id = hashlib.sha256(recipe_text.encode("utf-8")).hexdigest()[:12]
    relative = Path(scope) / sanitize(cruise_group or "all_cruises") / family / f"{sanitize(compartment or cruise_group)}.csv"
    destination = outdir / "media" / relative
    destination.parent.mkdir(parents=True, exist_ok=True)
    medium.reset_index(drop=True)[["compounds", "name", "maxFlux"]].to_csv(
        destination, index=False, quoting=csv.QUOTE_MINIMAL
    )
    manifest = {
        "scope": scope, "cruise_group": cruise_group, "family": family,
        "compartment": compartment, "medium_file": str(Path("media") / relative),
        "template": template_name, "template_source": template_source,
        "effective_n": effective_n,
        "used_global_fallback": any(row["value_source"] == "global_median_fallback" for row in provenance),
        "n_compounds": len(medium),
        "recipe_id": recipe_id,
    }
    return provenance, manifest


def main() -> None:
    args = parse_args()
    if args.min_effective_n < 0 or args.detection_limit < 0:
        raise ValueError("--min-effective-n and --detection-limit must be nonnegative")
    depth_baseline_m = parse_depths(args.depth_baseline_m)
    matrix = read_table(args.matrix)
    nutrients = load_nutrients(args.nutrients)
    mapping = load_mapping(args.mapping, matrix, nutrients)
    template_path = args.template if args.template is not None else DEFAULT_BASAL_MEDIUM
    template_source = "custom" if args.template is not None else "basin_builtin_default"
    basal_medium, template_audit = load_template(
        template_path, nutrients, template_source
    )
    cruise = cruise_group_family(matrix, read_table(args.cruise_groups))
    families = (
        assignment_family(matrix, read_table(args.o2_assignments), "legacy_o2"),
        assignment_family(matrix, read_table(args.gmm_assignments), "gmm"),
        assignment_family(matrix, read_table(args.hybrid_assignments), "hybrid"),
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    media_root = args.outdir / "media"
    if media_root.exists():
        shutil.rmtree(media_root)
    all_provenance: list[dict] = []
    manifest: list[dict] = []
    exclusions: list[dict] = []

    def supported(weights: pd.Series, scope: str, group: str, family: str, compartment: str) -> bool:
        effective_n = float(weights.sum())
        if effective_n >= args.min_effective_n:
            return True
        exclusions.append({
            "scope": scope,
            "cruise_group": group,
            "family": family,
            "compartment": compartment,
            "effective_n": effective_n,
            "min_effective_n": args.min_effective_n,
            "reason": "insufficient_soft_membership_support",
        })
        return False

    # Cruise-only profiles.
    for group in cruise.labels:
        if not supported(cruise.weights[group], "cruise_groups", group, "cruise_group", ""):
            continue
        prov, item = add_medium(
            matrix=matrix, mapping=mapping, weights=cruise.weights[group],
            scope="cruise_groups", cruise_group=group, family="cruise_group",
            compartment="", basal_medium=basal_medium,
            template_name=str(template_audit["template"]),
            template_source=template_source, outdir=args.outdir, method=args.method,
            min_effective_n=args.min_effective_n, detection_limit=args.detection_limit,
        )
        all_provenance.extend(prov); manifest.append(item)

    # Conventional anchored-depth profiles. These provide a compact physical
    # baseline for the compartment media without averaging across depth.
    if depth_baseline_m:
        depth_col = find_col(matrix, DEPTH_CANDIDATES, "matrix")
        observed_depth = pd.to_numeric(matrix[depth_col], errors="coerce")
        for depth in depth_baseline_m:
            label = depth_label(depth)
            weights = pd.Series(
                np.isclose(observed_depth, depth, rtol=0.0, atol=1e-6),
                index=matrix.index,
                dtype=float,
            )
            if not supported(weights, "depth_baseline", "", "depth", label):
                continue
            prov, item = add_medium(
                matrix=matrix, mapping=mapping, weights=weights,
                scope="depth_baseline", cruise_group="", family="depth",
                compartment=label, basal_medium=basal_medium,
                template_name=str(template_audit["template"]),
                template_source=template_source, outdir=args.outdir,
                method=args.method, min_effective_n=args.min_effective_n,
                detection_limit=args.detection_limit,
            )
            all_provenance.extend(prov); manifest.append(item)

    # Legacy compartment profiles and the complete cruise-group x compartment grid.
    for family in families:
        for label in family.labels:
            if not supported(
                family.weights[label], "compartments", "", family.name, label
            ):
                continue
            prov, item = add_medium(
                matrix=matrix, mapping=mapping, weights=family.weights[label],
                scope="compartments", cruise_group="", family=family.name,
                compartment=label, basal_medium=basal_medium,
                template_name=str(template_audit["template"]),
                template_source=template_source, outdir=args.outdir, method=args.method,
                min_effective_n=args.min_effective_n, detection_limit=args.detection_limit,
            )
            all_provenance.extend(prov); manifest.append(item)
        for group in cruise.labels:
            for label in family.labels:
                weights = cruise.weights[group] * family.weights[label]
                if not supported(
                    weights, "cruise_group_x_compartment", group, family.name, label
                ):
                    continue
                prov, item = add_medium(
                    matrix=matrix, mapping=mapping, weights=weights,
                    scope="cruise_group_x_compartment", cruise_group=group,
                    family=family.name, compartment=label, basal_medium=basal_medium,
                    template_name=str(template_audit["template"]),
                    template_source=template_source,
                    outdir=args.outdir, method=args.method,
                    min_effective_n=args.min_effective_n, detection_limit=args.detection_limit,
                )
                all_provenance.extend(prov); manifest.append(item)

    tables = args.outdir / "tables"
    tables.mkdir(exist_ok=True)
    pd.DataFrame(manifest).to_csv(tables / "gapseq_media_manifest.tsv", sep="\t", index=False)
    manifest_frame = pd.DataFrame(manifest)
    recipe_rows = []
    if not manifest_frame.empty:
        for recipe_id, frame in manifest_frame.groupby("recipe_id", sort=True):
            equivalent = frame["medium_file"].astype(str).tolist()
            for row in frame.itertuples(index=False):
                recipe_rows.append({
                    "recipe_id": recipe_id,
                    "scope": row.scope,
                    "cruise_group": row.cruise_group,
                    "family": row.family,
                    "compartment": row.compartment,
                    "medium_file": row.medium_file,
                    "equivalent_media_n": len(equivalent),
                    "equivalent_media": ";".join(equivalent),
                })
    pd.DataFrame(recipe_rows).to_csv(
        tables / "gapseq_media_recipe_audit.tsv", sep="\t", index=False
    )
    pd.DataFrame(all_provenance).to_csv(tables / "gapseq_media_provenance.tsv", sep="\t", index=False)
    pd.DataFrame(
        exclusions,
        columns=[
            "scope", "cruise_group", "family", "compartment", "effective_n",
            "min_effective_n", "reason",
        ],
    ).to_csv(tables / "gapseq_media_exclusions.tsv", sep="\t", index=False)
    mapping.to_csv(tables / "compound_mapping_audit.tsv", sep="\t", index=False)
    basal_medium.to_csv(tables / "basal_medium_used.csv", index=False)
    pd.DataFrame([template_audit]).to_csv(
        tables / "template_audit.tsv", sep="\t", index=False
    )
    summary = {
        "method": args.method,
        "detection_limit": args.detection_limit,
        "min_effective_n": args.min_effective_n,
        "n_matrix_rows": len(matrix),
        "n_cruise_groups": len(cruise.labels),
        "compartment_family_sizes": {
            family.name: len(family.labels) for family in families
        },
        "n_media": len(manifest),
        "n_unique_recipes": int(manifest_frame["recipe_id"].nunique())
        if not manifest_frame.empty else 0,
        "depth_baseline_m": list(depth_baseline_m),
        "n_media_using_global_fallback": sum(bool(item["used_global_fallback"]) for item in manifest),
        "n_excluded_unsupported_media": len(exclusions),
        "n_mapped_compounds": int(mapping["available_in_matrix"].sum()),
        "n_mapping_columns_absent": int((~mapping["available_in_matrix"]).sum()),
        "basal_medium": template_audit,
    }
    (tables / "gapseq_media_summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
