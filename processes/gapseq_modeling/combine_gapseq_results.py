#!/usr/bin/env python3
"""Combine per-genome gapseq environmental-comparison tables."""

from __future__ import annotations

import argparse
import csv
import glob
import itertools
import re
from pathlib import Path
from xml.etree import ElementTree

import numpy as np
import pandas as pd


TABLES = (
    "growth_comparison.tsv",
    "genome_performance_summary.tsv",
    "family_performance_summary.tsv",
    "all_exchange_fluxes.tsv",
    "counterfactual_growth.tsv",
    "nutrient_importance_summary.tsv",
)
TAXONOMY_RANKS = (
    "domain", "phylum", "class", "order", "family", "genus", "species"
)
TAXONOMY_COLUMNS = tuple(f"tax_{rank}" for rank in TAXONOMY_RANKS)


def count_sbml_elements(path: Path, element_name: str) -> int:
    """Count an SBML element without loading the complete XML tree."""
    count = 0
    for _, element in ElementTree.iterparse(path, events=("end",)):
        if element.tag.rsplit("}", 1)[-1] == element_name:
            count += 1
        element.clear()
    return count


def single_matching_file(
    directory: Path,
    pattern: str,
    label: str,
    preferred_name: str | None = None,
) -> Path:
    """Select one artifact, resolving stale publication copies by exact name."""
    matches = sorted(directory.glob(pattern))
    if len(matches) == 1:
        return matches[0]
    if preferred_name:
        preferred = [path for path in matches if path.name == preferred_name]
        if len(preferred) == 1:
            return preferred[0]
    candidate_names = ", ".join(path.name for path in matches) or "none"
    preferred_text = (
        f"; preferred exact filename was {preferred_name!r}"
        if preferred_name else ""
    )
    raise ValueError(
        f"Expected one resolvable {label} matching {pattern!r} in "
        f"{directory}, found {len(matches)}{preferred_text}. "
        f"Candidates: {candidate_names}"
    )


def gapseq_input_stem(input_file: object) -> str:
    """Return the filename stem used by gapseq for reconstruction artifacts."""
    name = Path(str(input_file)).name
    if name.lower().endswith(".gz"):
        name = name[:-3]
    for suffix in (".fasta", ".faa", ".fna", ".fa"):
        if name.lower().endswith(suffix):
            return name[:-len(suffix)]
    return Path(name).stem


def count_predicted_medium_compounds(path: Path) -> int:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or "compounds" not in reader.fieldnames:
            raise ValueError(
                f"Predicted medium must contain a compounds column: {path}"
            )
        return sum(
            1 for row in reader
            if str(row.get("compounds", "")).strip()
        )


def write_reconstruction_qc(
    reconstructions_root: Path,
    results_root: Path,
    source_metadata: pd.DataFrame,
    outdir: Path,
    expected_genome_ids: set[str] | None = None,
) -> pd.DataFrame:
    """Write one auditable reconstruction-QC row per genome."""
    rows = []
    if expected_genome_ids:
        reconstruction_dirs = [
            reconstructions_root / genome_id
            for genome_id in sorted(expected_genome_ids)
        ]
        missing = [path.name for path in reconstruction_dirs if not path.is_dir()]
        if missing:
            raise ValueError(
                "Missing reconstruction directories for current manifest: "
                + ", ".join(missing)
            )
    else:
        reconstruction_dirs = sorted(
            path for path in reconstructions_root.iterdir() if path.is_dir()
        )
    for reconstruction_dir in reconstruction_dirs:
        genome_id = reconstruction_dir.name
        provenance_file = reconstruction_dir / "reconstruction_provenance.tsv"
        if not provenance_file.is_file():
            raise ValueError(
                f"Missing reconstruction provenance: {provenance_file}"
            )
        provenance = pd.read_csv(provenance_file, sep="\t", dtype=str)
        if len(provenance) != 1:
            raise ValueError(
                f"Expected one provenance row in {provenance_file}, "
                f"found {len(provenance)}"
            )
        provenance_row = provenance.iloc[0]
        artifact_stem = gapseq_input_stem(provenance_row.get("input_file", ""))
        draft_candidates = sorted(reconstruction_dir.glob("*-draft.xml"))
        draft_model = single_matching_file(
            reconstruction_dir,
            "*-draft.xml",
            "draft SBML model",
            preferred_name=f"{artifact_stem}-draft.xml",
        )
        final_model = reconstruction_dir / f"{genome_id}.xml"
        if not final_model.is_file():
            raise ValueError(f"Missing final SBML model: {final_model}")
        medium_candidates = sorted(reconstruction_dir.glob("*-medium.csv"))
        predicted_medium = single_matching_file(
            reconstruction_dir,
            "*-medium.csv",
            "gapseq-predicted medium",
            preferred_name=f"{artifact_stem}-medium.csv",
        )
        model_summary_file = results_root / genome_id / "model_summary.tsv"
        if not model_summary_file.is_file():
            raise ValueError(
                f"Missing model summary for reconstructed genome {genome_id}: "
                f"{model_summary_file}"
            )
        model_summary = pd.read_csv(model_summary_file, sep="\t")
        if len(model_summary) != 1:
            raise ValueError(
                f"Expected one model-summary row in {model_summary_file}, "
                f"found {len(model_summary)}"
            )

        retry_file = reconstruction_dir / "reconstruction_retry_audit.tsv"
        if retry_file.is_file():
            retry_audit = pd.read_csv(retry_file, sep="\t", dtype=str)
            if len(retry_audit) != 1:
                raise ValueError(
                    f"Expected one retry-audit row in {retry_file}, "
                    f"found {len(retry_audit)}"
                )
            retry_row = retry_audit.iloc[0]
        else:
            retry_row = pd.Series(dtype=object)
        model_row = model_summary.iloc[0]
        draft_reactions = count_sbml_elements(draft_model, "reaction")
        final_reactions = int(model_row["reactions"])
        rows.append({
            "genome_id": genome_id,
            "reconstruction_status": "complete",
            "input_file": provenance_row.get("input_file", ""),
            "taxonomy": provenance_row.get("taxonomy", ""),
            "aligner": provenance_row.get("aligner", ""),
            "threads": provenance_row.get("threads", ""),
            "gapfill_medium_source": provenance_row.get(
                "gapfill_medium_source", ""
            ),
            "database_source": provenance_row.get("database_source", ""),
            "reconstruction_mode": provenance_row.get(
                "reconstruction_mode", "gapseq_doall_default"
            ),
            "default_doall_exit_status": provenance_row.get(
                "default_doall_exit_status", 0
            ),
            "effective_gapfill_minimum_growth": provenance_row.get(
                "effective_gapfill_minimum_growth", 0.01
            ),
            "effective_gapfill_solver": provenance_row.get(
                "effective_gapfill_solver", "auto"
            ),
            "gapfill_numeric_retry_attempted": retry_row.get(
                "retry_attempted", "false"
            ),
            "gapfill_numeric_retry_exit_status": retry_row.get(
                "retry_exit_status", ""
            ),
            "draft_model_file": str(draft_model),
            "draft_model_candidates_n": len(draft_candidates),
            "ignored_draft_model_files": "|".join(
                path.name for path in draft_candidates if path != draft_model
            ),
            "final_model_file": str(final_model),
            "predicted_medium_file": str(predicted_medium),
            "predicted_medium_candidates_n": len(medium_candidates),
            "ignored_predicted_medium_files": "|".join(
                path.name for path in medium_candidates
                if path != predicted_medium
            ),
            "draft_reactions": draft_reactions,
            "final_reactions": final_reactions,
            "net_reactions_added": final_reactions - draft_reactions,
            "final_metabolites": int(model_row["metabolites"]),
            "final_gene_products": int(model_row["genes"]),
            "final_exchanges": int(model_row["exchanges"]),
            "predicted_medium_compounds":
                count_predicted_medium_compounds(predicted_medium),
            "solver": model_row.get("solver", ""),
            "initial_solver_status": model_row.get("initial_status", ""),
            "initial_objective_value": model_row.get(
                "initial_objective_value", np.nan
            ),
        })

    result = pd.DataFrame(rows)
    if not result.empty and not source_metadata.empty:
        result = result.merge(
            source_metadata,
            on="genome_id",
            how="left",
            validate="one_to_one",
        )
    result.to_csv(
        outdir / "combined_reconstruction_qc.tsv", sep="\t", index=False
    )
    return result


def write_reconstruction_summary(
    reconstruction_qc: pd.DataFrame,
    outdir: Path,
) -> pd.DataFrame:
    """Write publication-ready descriptive reconstruction statistics."""
    rows = [
        {
            "metric": "reconstructed_genomes",
            "units": "genomes",
            "genomes_n": len(reconstruction_qc),
            "count": len(reconstruction_qc),
        },
        {
            "metric": "complete_reconstructions",
            "units": "genomes",
            "genomes_n": len(reconstruction_qc),
            "count": int(
                reconstruction_qc["reconstruction_status"].eq("complete").sum()
            ),
        },
        {
            "metric": "failed_reconstructions",
            "units": "genomes",
            "genomes_n": len(reconstruction_qc),
            "count": int(
                reconstruction_qc["reconstruction_status"].ne("complete").sum()
            ),
        },
        {
            "metric": "numeric_retry_reconstructions",
            "units": "genomes",
            "genomes_n": len(reconstruction_qc),
            "count": int(
                reconstruction_qc["reconstruction_mode"]
                .eq("gapseq_gapfill_reduced_minimum_growth_retry").sum()
            ),
        },
    ]
    numeric_metrics = (
        ("final_reactions", "reactions"),
        ("final_metabolites", "metabolites"),
        ("final_gene_products", "gene products"),
        ("final_exchanges", "exchange reactions"),
        ("net_reactions_added", "reactions"),
        ("predicted_medium_compounds", "compounds"),
    )
    for metric, units in numeric_metrics:
        values = pd.to_numeric(
            reconstruction_qc[metric], errors="coerce"
        ).dropna()
        rows.append({
            "metric": metric,
            "units": units,
            "genomes_n": len(values),
            "median": values.median(),
            "q25": values.quantile(0.25, interpolation="lower"),
            "q75": values.quantile(0.75, interpolation="lower"),
            "minimum": values.min(),
            "maximum": values.max(),
            "quantile_interpolation": "lower",
        })
    summary = pd.DataFrame(rows)
    ordered_columns = [
        "metric", "units", "genomes_n", "count", "median", "q25", "q75",
        "minimum", "maximum", "quantile_interpolation",
    ]
    summary = summary.reindex(columns=ordered_columns)
    summary.to_csv(
        outdir / "combined_reconstruction_summary.tsv", sep="\t", index=False
    )
    return summary


def infer_source_depth(path: str) -> float:
    matches = re.findall(r"(?i)(?:^|[/_-])(\d+(?:\.\d+)?)m(?:[/_.-]|$)", path)
    return float(matches[-1]) if matches else np.nan


def normalize_genome_id(value: object) -> str:
    text = str(value).strip().lower()
    match = re.search(
        r"ab-\d+_[a-z]\d+_(?:ab-\d+|[a-z]\d+)(?:\.sag(?:-xpg)?)?", text,
        flags=re.IGNORECASE,
    )
    if match:
        return match.group(0).lower()
    text = re.sub(r"(?i)\.(fna|fa|fasta|faa)(\.gz)?$", "", Path(text).name)
    return text


def load_genome_taxonomy(path: Path | None) -> pd.DataFrame:
    columns = ["genome_id", *TAXONOMY_COLUMNS]
    if path is None:
        return pd.DataFrame(columns=columns)
    first_line = path.read_text().splitlines()[0]
    delimiter = "\t" if "\t" in first_line else ","
    metadata = pd.read_csv(path, sep=delimiter, dtype=str)
    id_candidates = (
        "genome_id", "Genome_Id", "mag_genome_uid",
        "mag_native_genome_id_raw", "mag_genome_id",
    )
    id_column = next(
        (candidate for candidate in id_candidates if candidate in metadata),
        None,
    )
    if id_column is None:
        raise ValueError(
            "Genome metadata requires one identifier column from: "
            + ", ".join(id_candidates)
        )
    lower_columns = {column.lower(): column for column in metadata.columns}
    result = pd.DataFrame({
        "genome_id": metadata[id_column].map(normalize_genome_id)
    })
    for rank in TAXONOMY_RANKS:
        source = lower_columns.get(rank) or lower_columns.get(f"mag_{rank}")
        result[f"tax_{rank}"] = metadata[source] if source else pd.NA
    result = result.drop_duplicates()
    if result["genome_id"].duplicated().any():
        duplicates = sorted(
            result.loc[result["genome_id"].duplicated(False), "genome_id"]
            .unique()
        )
        raise ValueError(
            "Genome metadata contains conflicting duplicate identifiers: "
            + ", ".join(duplicates[:10])
        )
    return result.reindex(columns=columns)


def load_genome_source_metadata(path: Path | None) -> pd.DataFrame:
    columns = ["genome_id", "source_depth_m", "source_genome_file"]
    if path is None:
        return pd.DataFrame(columns=columns)
    delimiter = "\t" if "\t" in path.read_text().splitlines()[0] else ","
    manifest = pd.read_csv(path, sep=delimiter, dtype=str)
    if list(manifest.columns) != ["label", "filepath"]:
        raise ValueError("Genome manifest must have exactly label,filepath columns")
    rows = []
    for row in manifest.itertuples(index=False):
        source = Path(row.filepath)
        pattern = str(source if source.is_absolute() else path.parent / source)
        has_glob = any(character in pattern for character in "*?[]{}")
        matches = sorted(glob.glob(pattern)) if has_glob else [pattern]
        for match in matches:
            label = row.label
            if label.lower() == "auto" or (has_glob and len(matches) > 1):
                basename = re.sub(
                    r"(?i)\.(fna|fa|fasta|faa)(\.gz)?$", "", Path(match).name
                )
                label = basename if label.lower() == "auto" else f"{label}_{basename}"
            rows.append({
                "genome_id": label,
                "source_depth_m": infer_source_depth(str(Path(match).resolve())),
                "source_genome_file": str(Path(match).resolve()),
            })
    result = pd.DataFrame(rows, columns=columns)
    if result["genome_id"].duplicated().any():
        raise ValueError("Genome manifest resolves to duplicate genome identifiers")
    return result


def lineage_counts(frame: pd.DataFrame, rank: str) -> str:
    if frame.empty or rank not in frame:
        return ""
    values = frame[["genome_id", rank]].copy()
    values[rank] = values[rank].fillna("").astype(str).str.strip()
    values = values.loc[values[rank].ne("")]
    if values.empty:
        return ""
    counts = (
        values.drop_duplicates()
        .groupby(rank, sort=True)["genome_id"]
        .nunique()
        .sort_values(ascending=False)
    )
    return "; ".join(
        f"{lineage} (n={int(count)})"
        for lineage, count in counts.items()
    )


def write_counterfactual_taxonomy_summaries(
    nutrient_importance: pd.DataFrame,
    genome_metadata: pd.DataFrame,
    outdir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    required = {
        "genome_id", "compound_id", "compound_name",
        "media_with_numeric_effect", "importance_score",
        "maximum_relative_growth_loss",
        "fraction_media_growth_lost_completely",
    }
    if not required.issubset(nutrient_importance.columns):
        return pd.DataFrame(), pd.DataFrame()
    frame = nutrient_importance.copy()
    for column in (
        "media_with_numeric_effect", "importance_score",
        "mean_relative_growth_loss", "maximum_relative_growth_loss",
        "fraction_media_growth_lost_completely",
    ):
        if column in frame:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame["importance_gt_1pct"] = frame["importance_score"].gt(0.01)
    frame["growth_lost_any_medium"] = frame[
        "fraction_media_growth_lost_completely"
    ].gt(0)
    frame["growth_lost_all_tested_media"] = frame[
        "fraction_media_growth_lost_completely"
    ].ge(1 - 1e-12)

    metadata_columns = ["genome_id", *TAXONOMY_COLUMNS]
    if not set(TAXONOMY_COLUMNS).issubset(frame.columns):
        frame = frame.merge(
            genome_metadata[metadata_columns],
            on="genome_id",
            how="left",
            validate="many_to_one",
        )

    compound_rows = []
    lineage_rows = []
    group_columns = ["compound_id", "compound_name"]
    if "nutrient_source" in frame:
        group_columns.append("nutrient_source")
    for keys, compound_frame in frame.groupby(
        group_columns, dropna=False, sort=True
    ):
        if not isinstance(keys, tuple):
            keys = (keys,)
        identity = dict(zip(group_columns, keys))
        tested = compound_frame.loc[
            compound_frame["media_with_numeric_effect"].gt(0)
        ].copy()
        affected = tested.loc[tested["importance_gt_1pct"]]
        complete_loss = tested.loc[tested["growth_lost_any_medium"]]
        row = {
            **identity,
            "models_tested": tested["genome_id"].nunique(),
            "models_with_importance_gt_1pct":
                affected["genome_id"].nunique(),
            "models_with_complete_growth_loss":
                complete_loss["genome_id"].nunique(),
            "models_with_complete_loss_in_all_tested_media": tested.loc[
                tested["growth_lost_all_tested_media"], "genome_id"
            ].nunique(),
            "median_maximum_relative_growth_loss": tested[
                "maximum_relative_growth_loss"
            ].median(),
            "maximum_relative_growth_loss": tested[
                "maximum_relative_growth_loss"
            ].max(),
        }
        for rank in ("phylum", "class", "order", "family", "genus", "species"):
            tax_column = f"tax_{rank}"
            row[f"importance_gt_1pct_{rank}_counts"] = lineage_counts(
                affected, tax_column
            )
            row[f"complete_loss_{rank}_counts"] = lineage_counts(
                complete_loss, tax_column
            )
        compound_rows.append(row)

        for rank in TAXONOMY_RANKS:
            tax_column = f"tax_{rank}"
            rank_frame = tested.copy()
            rank_frame[tax_column] = (
                rank_frame[tax_column].fillna("").astype(str).str.strip()
            )
            rank_frame = rank_frame.loc[rank_frame[tax_column].ne("")]
            for lineage, lineage_frame in rank_frame.groupby(
                tax_column, sort=True
            ):
                lineage_rows.append({
                    **identity,
                    "taxonomic_rank": rank,
                    "lineage": lineage,
                    "models_tested": lineage_frame["genome_id"].nunique(),
                    "models_with_importance_gt_1pct": lineage_frame.loc[
                        lineage_frame["importance_gt_1pct"], "genome_id"
                    ].nunique(),
                    "models_with_complete_growth_loss": lineage_frame.loc[
                        lineage_frame["growth_lost_any_medium"], "genome_id"
                    ].nunique(),
                    "fraction_models_with_importance_gt_1pct": (
                        lineage_frame.groupby("genome_id")[
                            "importance_gt_1pct"
                        ].max().mean()
                    ),
                    "fraction_models_with_complete_growth_loss": (
                        lineage_frame.groupby("genome_id")[
                            "growth_lost_any_medium"
                        ].max().mean()
                    ),
                    "median_maximum_relative_growth_loss": lineage_frame[
                        "maximum_relative_growth_loss"
                    ].median(),
                    "maximum_relative_growth_loss": lineage_frame[
                        "maximum_relative_growth_loss"
                    ].max(),
                })
    compound_summary = pd.DataFrame(compound_rows)
    lineage_summary = pd.DataFrame(lineage_rows)
    compound_summary.to_csv(
        outdir / "combined_counterfactual_taxonomic_summary.tsv",
        sep="\t",
        index=False,
    )
    lineage_summary.to_csv(
        outdir / "combined_counterfactual_lineage_summary.tsv",
        sep="\t",
        index=False,
    )
    return compound_summary, lineage_summary


def mean_pairwise_absolute(values: pd.Series) -> float:
    numeric = pd.to_numeric(values, errors="coerce").dropna().to_numpy(float)
    if len(numeric) < 2:
        return np.nan
    return float(np.mean([abs(a - b) for a, b in itertools.combinations(numeric, 2)]))


def comparison_subsets(growth: pd.DataFrame):
    yield "all_genomes", growth
    if "source_depth_m" not in growth:
        return
    for depth, frame in growth.dropna(subset=["source_depth_m"]).groupby(
        "source_depth_m", sort=True
    ):
        yield f"source_depth_{depth:g}m", frame


def write_compartment_comparison(growth: pd.DataFrame, outdir: Path) -> None:
    required = {
        "genome_id", "scope", "family", "supplemented_biomass_flux",
        "within_genome_percent_of_max",
    }
    if not required.issubset(growth.columns):
        return
    primary = growth.loc[
        growth["scope"].isin(["depth_baseline", "compartments"])
        & growth["family"].isin(["depth", "legacy_o2", "gmm", "hybrid"])
    ].copy()
    if primary.empty:
        return
    primary["recipe_id"] = primary.get(
        "recipe_id", primary["medium"]
    ).fillna(primary["medium"])
    rows = []
    hybrid_rows = []
    for subset_name, subset in comparison_subsets(primary):
        for family, frame in subset.groupby("family", sort=True):
            nominal_media_n = frame["medium"].nunique()
            unique_recipe_n = frame["recipe_id"].nunique()
            per_genome = []
            for genome_id, genome_frame in frame.groupby("genome_id", sort=True):
                unique = genome_frame.drop_duplicates("recipe_id")
                percent = pd.to_numeric(
                    unique["within_genome_percent_of_max"], errors="coerce"
                )
                per_genome.append({
                    "genome_id": genome_id,
                    "distinct_growth_states": int(
                        pd.to_numeric(
                            unique["supplemented_biomass_flux"], errors="coerce"
                        ).round(9).nunique()
                    ),
                    "growth_range_pct_max": float(percent.max() - percent.min()),
                    "mean_pairwise_growth_separation_pct_max":
                        mean_pairwise_absolute(percent),
                })
            metrics = pd.DataFrame(per_genome)
            rows.append({
                "analysis_subset": subset_name,
                "family": family,
                "genomes": metrics["genome_id"].nunique(),
                "nominal_media_n": nominal_media_n,
                "unique_recipe_n": unique_recipe_n,
                "recipe_retention_fraction": unique_recipe_n / nominal_media_n,
                "median_distinct_growth_states": metrics[
                    "distinct_growth_states"
                ].median(),
                "fraction_genomes_with_multiple_growth_states": (
                    metrics["distinct_growth_states"] > 1
                ).mean(),
                "median_growth_range_pct_max": metrics[
                    "growth_range_pct_max"
                ].median(),
                "growth_range_pct_max_q25": metrics[
                    "growth_range_pct_max"
                ].quantile(0.25),
                "growth_range_pct_max_q75": metrics[
                    "growth_range_pct_max"
                ].quantile(0.75),
                "median_pairwise_growth_separation_pct_max": metrics[
                    "mean_pairwise_growth_separation_pct_max"
                ].median(),
            })

        hybrid = subset.loc[subset["family"].eq("hybrid")].copy()
        if hybrid.empty:
            continue
        hybrid["oxygen_component"] = hybrid["compartment"].astype(str).str.extract(
            r"hybrid_c(\d+)_", expand=False
        )
        oxygen_names = {
            "0": "oxic", "1": "dysoxic", "2": "suboxic", "3": "anoxic"
        }
        for component, frame in hybrid.dropna(
            subset=["oxygen_component"]
        ).groupby("oxygen_component", sort=True):
            per_genome = []
            for genome_id, genome_frame in frame.groupby("genome_id", sort=True):
                unique = genome_frame.drop_duplicates("recipe_id")
                percent = pd.to_numeric(
                    unique["within_genome_percent_of_max"], errors="coerce"
                )
                per_genome.append({
                    "genome_id": genome_id,
                    "distinct_growth_states": int(
                        pd.to_numeric(
                            unique["supplemented_biomass_flux"], errors="coerce"
                        ).round(9).nunique()
                    ),
                    "growth_range_pct_max": float(percent.max() - percent.min()),
                })
            metrics = pd.DataFrame(per_genome)
            hybrid_rows.append({
                "analysis_subset": subset_name,
                "oxygen_component": int(component),
                "oxygen_compartment": oxygen_names.get(component, component),
                "genomes": metrics["genome_id"].nunique(),
                "nominal_hybrid_media_n": frame["medium"].nunique(),
                "unique_recipe_n": frame["recipe_id"].nunique(),
                "genomes_with_added_resolution": int(
                    (metrics["distinct_growth_states"] > 1).sum()
                ),
                "fraction_genomes_with_added_resolution": (
                    metrics["distinct_growth_states"] > 1
                ).mean(),
                "median_distinct_growth_states": metrics[
                    "distinct_growth_states"
                ].median(),
                "median_within_oxygen_growth_range_pct_max": metrics[
                    "growth_range_pct_max"
                ].median(),
            })
    pd.DataFrame(rows).to_csv(
        outdir / "combined_compartment_system_comparison.tsv",
        sep="\t", index=False,
    )
    pd.DataFrame(hybrid_rows).to_csv(
        outdir / "combined_hybrid_incremental_resolution.tsv",
        sep="\t", index=False,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-root", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--genomes-manifest", type=Path)
    parser.add_argument("--genome-metadata", type=Path)
    parser.add_argument("--reconstructions-root", type=Path)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    source_metadata = load_genome_source_metadata(args.genomes_manifest)
    taxonomy = load_genome_taxonomy(args.genome_metadata)
    if not taxonomy.empty:
        source_metadata = source_metadata.merge(
            taxonomy,
            on="genome_id",
            how="left",
            validate="one_to_one",
        )
    source_metadata.to_csv(
        args.outdir / "combined_genome_source_metadata.tsv", sep="\t", index=False
    )
    expected_genome_ids = (
        set(source_metadata["genome_id"].dropna().astype(str))
        if not source_metadata.empty and "genome_id" in source_metadata
        else set()
    )
    reconstruction_ids = (
        {
            path.name for path in args.reconstructions_root.iterdir()
            if path.is_dir()
        }
        if args.reconstructions_root is not None
        and args.reconstructions_root.is_dir()
        else set()
    )
    comparison_ids = {
        path.name for path in args.results_root.iterdir() if path.is_dir()
    }
    directory_ids = sorted(
        expected_genome_ids | reconstruction_ids | comparison_ids
    )
    directory_audit = pd.DataFrame([
        {
            "genome_id": genome_id,
            "in_current_manifest": (
                genome_id in expected_genome_ids
                if expected_genome_ids else pd.NA
            ),
            "reconstruction_directory_present": genome_id in reconstruction_ids,
            "comparison_directory_present": genome_id in comparison_ids,
            "model_summary_present": (
                args.results_root / genome_id / "model_summary.tsv"
            ).is_file(),
            "selection_status": (
                "current_manifest"
                if not expected_genome_ids or genome_id in expected_genome_ids
                else "ignored_stale_directory"
            ),
        }
        for genome_id in directory_ids
    ])
    directory_audit.to_csv(
        args.outdir / "combined_genome_directory_selection_audit.tsv",
        sep="\t",
        index=False,
    )
    inventory = []
    if args.reconstructions_root is not None:
        reconstruction_qc = write_reconstruction_qc(
            args.reconstructions_root,
            args.results_root,
            source_metadata,
            args.outdir,
            expected_genome_ids=expected_genome_ids,
        )
        reconstruction_summary = write_reconstruction_summary(
            reconstruction_qc, args.outdir
        )
        inventory.append({
            "source_table": "derived_reconstruction_qc",
            "genomes": reconstruction_qc["genome_id"].nunique(),
            "rows": len(reconstruction_qc),
            "output": "combined_reconstruction_qc.tsv",
        })
        inventory.append({
            "source_table": "derived_reconstruction_summary",
            "genomes": reconstruction_qc["genome_id"].nunique(),
            "rows": len(reconstruction_summary),
            "output": "combined_reconstruction_summary.tsv",
        })
    combined_tables = {}
    if expected_genome_ids:
        result_directories = [
            args.results_root / genome_id
            for genome_id in sorted(expected_genome_ids)
        ]
        missing_results = [
            path.name for path in result_directories if not path.is_dir()
        ]
        if missing_results:
            raise ValueError(
                "Missing comparison directories for current manifest: "
                + ", ".join(missing_results)
            )
    else:
        result_directories = sorted(
            path for path in args.results_root.iterdir() if path.is_dir()
        )
    for table in TABLES:
        frames = []
        for genome_dir in result_directories:
            source = genome_dir / table
            if not source.is_file() or source.stat().st_size == 0:
                continue
            frame = pd.read_csv(source, sep="\t")
            frame.insert(0, "genome_id", genome_dir.name)
            frames.append(frame)
        combined = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
        if not combined.empty and not source_metadata.empty:
            annotation_columns = [
                column for column in (
                    "genome_id", "source_depth_m", *TAXONOMY_COLUMNS
                )
                if column in source_metadata
            ]
            combined = combined.merge(
                source_metadata[annotation_columns],
                on="genome_id", how="left", validate="many_to_one"
            )
        destination = args.outdir / f"combined_{table}"
        combined.to_csv(destination, sep="\t", index=False)
        combined_tables[table] = combined
        inventory.append({"source_table": table, "genomes": len(frames), "rows": len(combined), "output": destination.name})
    growth = combined_tables.get("growth_comparison.tsv", pd.DataFrame())
    if not growth.empty:
        write_compartment_comparison(growth, args.outdir)
        inventory.extend([
            {
                "source_table": "derived_compartment_comparison",
                "genomes": growth["genome_id"].nunique(),
                "rows": len(pd.read_csv(
                    args.outdir / "combined_compartment_system_comparison.tsv",
                    sep="\t",
                )),
                "output": "combined_compartment_system_comparison.tsv",
            },
            {
                "source_table": "derived_hybrid_incremental_resolution",
                "genomes": growth["genome_id"].nunique(),
                "rows": len(pd.read_csv(
                    args.outdir / "combined_hybrid_incremental_resolution.tsv",
                    sep="\t",
                )),
                "output": "combined_hybrid_incremental_resolution.tsv",
            },
        ])
    nutrient_importance = combined_tables.get(
        "nutrient_importance_summary.tsv", pd.DataFrame()
    )
    if not nutrient_importance.empty and not taxonomy.empty:
        compound_summary, lineage_summary = (
            write_counterfactual_taxonomy_summaries(
                nutrient_importance,
                source_metadata,
                args.outdir,
            )
        )
        inventory.extend([
            {
                "source_table": "derived_counterfactual_taxonomic_summary",
                "genomes": nutrient_importance["genome_id"].nunique(),
                "rows": len(compound_summary),
                "output":
                    "combined_counterfactual_taxonomic_summary.tsv",
            },
            {
                "source_table": "derived_counterfactual_lineage_summary",
                "genomes": nutrient_importance["genome_id"].nunique(),
                "rows": len(lineage_summary),
                "output":
                    "combined_counterfactual_lineage_summary.tsv",
            },
        ])
    pd.DataFrame(inventory).to_csv(args.outdir / "combined_results_inventory.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
