#!/usr/bin/env python3
"""Plot model-derived genome niche positions over the environmental PCA.

Genome coordinates are not direct PCA projections.  Each coordinate is the
barycenter of sample-supported hybrid-compartment centroids, weighted by the
genome model's positive growth contrast among the corresponding unique media
recipes.  The accompanying tables expose every intermediate used by the plot.
"""

from __future__ import annotations

import argparse
import colorsys
import math
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgb
from matplotlib.lines import Line2D

PROCESSES_DIR = Path(__file__).resolve().parents[1]
if str(PROCESSES_DIR) not in sys.path:
    sys.path.insert(0, str(PROCESSES_DIR))

from shared_plot_export import save_figure_all_formats
from shared_plot_style import (  # noqa: E402
    install_publication_style,
)


OXYGEN_NAMES = ("oxic", "dysoxic", "suboxic", "anoxic")
O2_COLORS = {
    "oxic": "red",
    "dysoxic": "green",
    "suboxic": "lightblue",
    "anoxic": "purple",
}
# ASPIRE's colorblind-friendly categorical colors, assigned deterministically
# to the exact phylum strings reported by the counterfactual summary tables.
PHYLUM_COLORS = {
    "Pseudomonadota": "#0072B2",
    "Campylobacterota": "#D55E00",
    "Bacteroidota": "#009E73",
    "SAR324": "#CC79A7",
    "Marinisomatota": "#E69F00",
    "Patescibacteria": "#6A3D9A",
    "Patescibacteriota": "#6A3D9A",
    "Thermoproteota": "#8C564B",
    "Thermoplasmatota": "#56B4E9",
    "Unclassified": "#7F7F7F",
}
RESPONSE_MARKERS = {
    "Neither": "o",
    "O2-sensitive only": "^",
    "H2S-dependent only": "s",
    "O2-sensitive + H2S-dependent": "D",
    "Not assessed": "P",
}


def read_table(path: Path) -> pd.DataFrame:
    try:
        return pd.read_csv(
            path,
            sep="\t"
            if path.suffix.lower() in {".tsv", ".tab", ".txt"}
            else ",",
        )
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def parse_seqkit_paired_libraries(
    seqkit: pd.DataFrame,
    *,
    file_col: str = "file",
    count_col: str = "num_seqs",
) -> tuple[pd.Series, pd.DataFrame]:
    """Parse paired FASTQ SeqKit rows into sample-level input fragments.

    Recruitment counts use normalized read-fragment identifiers, so a paired
    library containing N R1 records and N R2 records has N input fragments.
    Mate counts must agree exactly; summing the two rows would double-count the
    denominator.
    """
    required = {file_col, count_col}
    if not required.issubset(seqkit):
        raise ValueError(
            "SeqKit table lacks required columns: "
            + ", ".join(sorted(required - set(seqkit)))
        )
    rows = []
    for source_row, row in seqkit[[file_col, count_col]].iterrows():
        filename = Path(str(row[file_col]).strip()).name
        stem = re.sub(
            r"(?i)\.(?:fastq|fq)(?:\.(?:gz|bz2|xz|zst))?$", "", filename
        )
        match = re.fullmatch(
            r"(?i)(?P<sample>.+?)(?:_pe[._-]?|[._-]R?)(?P<mate>[12])"
            r"(?:_001)?",
            stem,
        )
        if match is None:
            raise ValueError(
                "Could not parse paired-end sample and mate from SeqKit file "
                f"'{row[file_col]}'. Expected suffixes such as _pe.1/_pe.2, "
                "_R1/_R2, or _1/_2."
            )
        count = pd.to_numeric(pd.Series([row[count_col]]), errors="coerce").iloc[0]
        if not np.isfinite(count) or count <= 0 or float(count) != int(count):
            raise ValueError(
                f"SeqKit {count_col} must be a positive integer for '{filename}'"
            )
        rows.append({
            "seqkit_source_row": int(source_row),
            "seqkit_file": str(row[file_col]),
            "seqkit_filename": filename,
            "sample": match.group("sample"),
            "mate": int(match.group("mate")),
            "num_seqs": int(count),
        })
    if not rows:
        raise ValueError("SeqKit table contains no paired-end FASTQ rows")
    audit = pd.DataFrame(rows)
    audit["duplicate_source_rows"] = audit.groupby(
        ["sample", "mate"]
    )["seqkit_source_row"].transform("size")
    duplicates = audit["duplicate_source_rows"].gt(1)
    conflicting = (
        audit.loc[duplicates]
        .groupby(["sample", "mate"])
        .agg(
            filenames=("seqkit_filename", "nunique"),
            counts=("num_seqs", "nunique"),
        )
    )
    conflicting = conflicting.loc[
        conflicting["filenames"].gt(1) | conflicting["counts"].gt(1)
    ]
    if not conflicting.empty:
        raise ValueError(
            "SeqKit table contains conflicting duplicate sample/mate rows: "
            + ", ".join(
                f"{sample}/R{mate}" for sample, mate in conflicting.index
            )
        )
    audit = audit.drop_duplicates(["sample", "mate"], keep="first").copy()
    pairs = audit.pivot(index="sample", columns="mate", values="num_seqs")
    missing = pairs.index[pairs.reindex(columns=[1, 2]).isna().any(axis=1)]
    if len(missing):
        raise ValueError(
            "SeqKit paired libraries lack R1 or R2 for: "
            + ", ".join(map(str, missing[:10]))
        )
    mismatch = pairs.index[pairs[1].ne(pairs[2])]
    if len(mismatch):
        details = ", ".join(
            f"{sample} (R1={int(pairs.loc[sample, 1])}, "
            f"R2={int(pairs.loc[sample, 2])})"
            for sample in mismatch[:10]
        )
        raise ValueError(f"SeqKit R1/R2 sequence counts disagree: {details}")
    library_sizes = pairs[1].astype(float).rename("input_fragments")
    audit["input_fragments"] = audit["sample"].map(library_sizes)
    audit["pair_validation"] = np.where(
        audit["duplicate_source_rows"].gt(1),
        "matched_duplicate_reference_collapsed",
        "matched",
    )
    return library_sizes, audit.sort_values(["sample", "mate"]).reset_index(
        drop=True
    )


def normalize_recruitment_counts(
    table: pd.DataFrame,
    *,
    method: str,
    input_fragments: pd.Series | None = None,
) -> tuple[pd.DataFrame, pd.Series, str]:
    """Normalize counts or validate an explicitly pre-normalized matrix."""
    coerced = table.apply(pd.to_numeric, errors="coerce")
    selected = method.lower()
    if selected in {"provided_fpkm", "provided_tpm"}:
        invalid = coerced.isna() & table.notna()
        if invalid.any().any():
            raise ValueError(
                f"{selected} input contains nonnumeric abundance values"
            )
        numeric = coerced.fillna(0.0)
        if (~np.isfinite(numeric.to_numpy(dtype=float))).any():
            raise ValueError(f"{selected} input contains nonfinite values")
        if numeric.lt(0).any().any():
            raise ValueError(f"{selected} input contains negative values")
        return numeric, pd.Series(dtype=float), selected
    numeric = coerced.fillna(0.0).clip(lower=0)
    if selected == "auto":
        if input_fragments is None:
            raise ValueError(
                "auto abundance normalization requires a paired-end SeqKit "
                "library table; select median_ratio explicitly to use the "
                "legacy normalization"
            )
        selected = "input_fragment_fpm"
    if selected == "input_fragment_fpm":
        if input_fragments is None:
            raise ValueError(
                "input_fragment_fpm normalization requires a SeqKit library table"
            )
        denominators = pd.to_numeric(
            input_fragments.reindex(numeric.columns), errors="coerce"
        )
        missing = denominators.index[
            denominators.isna() | ~np.isfinite(denominators) | denominators.le(0)
        ]
        if len(missing):
            raise ValueError(
                "SeqKit input-fragment counts are missing or invalid for "
                "recruitment samples: " + ", ".join(map(str, missing[:10]))
            )
        normalized = numeric.div(denominators, axis=1).mul(1_000_000.0)
        return normalized, denominators, "input_fragment_fpm"
    if selected == "median_ratio":
        normalized, factors = median_ratio_normalize(numeric)
        return normalized, factors, "median_ratio_positive_counts"
    if selected == "raw":
        return numeric, pd.Series(dtype=float), "raw_counts"
    if selected == "relative":
        denominators = numeric.sum(axis=0).replace(0, np.nan)
        return (
            numeric.div(denominators, axis=1).fillna(0.0),
            denominators,
            "selected_catalog_relative_abundance",
        )
    raise ValueError(
        "Abundance normalization must be auto, input_fragment_fpm, "
        "provided_fpkm, provided_tpm, median_ratio, raw, or relative"
    )


def normalized_hybrid_label(value: object) -> str:
    text = str(value)
    match = re.fullmatch(
        r"(?:hyb_C|hybrid_c)([0-3])_(?:G|g)([0-9]+)", text
    )
    if not match:
        raise ValueError(f"Unrecognized hybrid-compartment label: {text}")
    return f"hybrid_c{int(match.group(1))}_g{int(match.group(2))}"


def responsibility_column(compartment: str) -> str:
    match = re.fullmatch(r"hybrid_c([0-3])_g([0-9]+)", compartment)
    if not match:
        raise ValueError(f"Unrecognized hybrid compartment: {compartment}")
    return f"hyb_C{match.group(1)}_G{match.group(2)}"


def display_hybrid_label(compartment: str) -> str:
    match = re.fullmatch(r"hybrid_c([0-3])_g([0-9]+)", compartment)
    if not match:
        return compartment
    return f"{OXYGEN_NAMES[int(match.group(1))]}-GMM{int(match.group(2))}"


def hybrid_sort_key(compartment: str) -> tuple[int, int]:
    match = re.fullmatch(r"hybrid_c([0-3])_g([0-9]+)", compartment)
    return (
        (int(match.group(1)), int(match.group(2)))
        if match
        else (99, 99)
    )


def hybrid_palette(compartments: list[str]) -> dict[str, tuple[float, float, float]]:
    parsed = [(*hybrid_sort_key(value), value) for value in compartments]
    gmm_values = sorted({gmm for oxygen, gmm, _ in parsed if oxygen < 99})
    if len(gmm_values) == 1:
        lightness = {gmm_values[0]: 0.52}
    else:
        lightness = dict(
            zip(gmm_values, np.linspace(0.30, 0.76, len(gmm_values)))
        )
    result = {}
    for oxygen, gmm, raw in parsed:
        hue, _, saturation = colorsys.rgb_to_hls(
            *to_rgb(O2_COLORS[OXYGEN_NAMES[oxygen]])
        )
        result[raw] = colorsys.hls_to_rgb(
            hue, float(lightness[gmm]), max(0.58, saturation)
        )
    return result


def lighten(color, fraction: float = 0.68):
    rgb = np.asarray(to_rgb(color), dtype=float)
    return tuple(rgb + (1.0 - rgb) * fraction)


def merge_sample_inputs(
    scores: pd.DataFrame, hybrid: pd.DataFrame
) -> pd.DataFrame:
    if "cruise_year_month_depth" in scores and "cruise_year_month_depth" in hybrid:
        keys = ["cruise_year_month_depth"]
    else:
        keys = ["Cruise", "date", "Depth_anchored"]
    required_scores = set(keys) | {"PC1", "PC2"}
    if not required_scores.issubset(scores.columns):
        raise ValueError(
            f"PCA scores lack required columns: {sorted(required_scores - set(scores))}"
        )
    missing_keys = set(keys) - set(hybrid.columns)
    if missing_keys:
        raise ValueError(
            f"Hybrid assignments lack join columns: {sorted(missing_keys)}"
        )
    score_columns = keys + ["PC1", "PC2"]
    merged = hybrid.merge(
        scores[score_columns],
        on=keys,
        how="inner",
        validate="one_to_one",
    )
    if merged.empty:
        raise ValueError("PCA scores and hybrid assignments have no matched samples")
    return merged


def build_hybrid_centroids(
    samples: pd.DataFrame, manifest: pd.DataFrame
) -> pd.DataFrame:
    hybrid_media = manifest.loc[
        manifest["scope"].eq("compartments")
        & manifest["family"].eq("hybrid")
    ].copy()
    if hybrid_media.empty:
        raise ValueError("Media manifest contains no hybrid-compartment media")
    hybrid_media["compartment"] = hybrid_media["compartment"].map(
        normalized_hybrid_label
    )
    hybrid_media["effective_n"] = pd.to_numeric(
        hybrid_media["effective_n"], errors="coerce"
    )

    rows = []
    for row in hybrid_media.itertuples(index=False):
        column = responsibility_column(row.compartment)
        if column not in samples:
            raise ValueError(
                f"Hybrid assignments lack responsibility column {column}"
            )
        weights = pd.to_numeric(samples[column], errors="coerce").fillna(0.0)
        valid = (
            pd.to_numeric(samples["PC1"], errors="coerce").notna()
            & pd.to_numeric(samples["PC2"], errors="coerce").notna()
            & weights.gt(0)
        )
        weights = weights.loc[valid]
        if weights.sum() <= 0:
            continue
        rows.append({
            "compartment": row.compartment,
            "display_label": display_hybrid_label(row.compartment),
            "responsibility_column": column,
            "recipe_id": str(row.recipe_id),
            "effective_n_manifest": row.effective_n,
            "sample_weight_sum": float(weights.sum()),
            "PC1_centroid": float(
                np.average(samples.loc[valid, "PC1"], weights=weights)
            ),
            "PC2_centroid": float(
                np.average(samples.loc[valid, "PC2"], weights=weights)
            ),
        })
    centroids = pd.DataFrame(rows)
    if centroids.empty:
        raise ValueError("No hybrid centroids could be calculated")
    return centroids.sort_values(
        "compartment", key=lambda x: x.map(hybrid_sort_key)
    ).reset_index(drop=True)


def build_recipe_centroids(hybrid_centroids: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for recipe_id, frame in hybrid_centroids.groupby("recipe_id", sort=True):
        weights = pd.to_numeric(
            frame["sample_weight_sum"], errors="coerce"
        ).fillna(0.0)
        if weights.sum() <= 0:
            weights = pd.Series(np.ones(len(frame)), index=frame.index)
        rows.append({
            "recipe_id": recipe_id,
            "hybrid_compartments": ";".join(frame["compartment"]),
            "hybrid_labels": ";".join(frame["display_label"]),
            "nominal_compartments_n": len(frame),
            "sample_weight_sum": float(weights.sum()),
            "PC1_recipe_centroid": float(
                np.average(frame["PC1_centroid"], weights=weights)
            ),
            "PC2_recipe_centroid": float(
                np.average(frame["PC2_centroid"], weights=weights)
            ),
        })
    return pd.DataFrame(rows)


def counterfactual_flags(
    nutrient_importance: pd.DataFrame,
) -> pd.DataFrame:
    output_columns = [
        "genome_id", "o2_sensitive_gt_1pct",
        "h2s_dependent_complete_growth_loss",
        "counterfactual_response_class",
    ]
    required = {
        "genome_id", "compound_name", "importance_score",
        "fraction_media_growth_lost_completely",
    }
    if nutrient_importance.empty or not required.issubset(
        nutrient_importance.columns
    ):
        return pd.DataFrame(columns=output_columns)
    frame = nutrient_importance.copy()
    frame["compound_name_key"] = (
        frame["compound_name"].fillna("").astype(str).str.strip().str.lower()
    )
    frame["importance_score"] = pd.to_numeric(
        frame["importance_score"], errors="coerce"
    )
    frame["fraction_media_growth_lost_completely"] = pd.to_numeric(
        frame["fraction_media_growth_lost_completely"], errors="coerce"
    )
    rows = []
    for genome_id, genome in frame.groupby("genome_id", sort=True):
        o2 = genome.loc[genome["compound_name_key"].eq("o2")]
        h2s = genome.loc[genome["compound_name_key"].eq("h2s")]
        o2_sensitive = bool(o2["importance_score"].gt(0.01).any())
        h2s_dependent = bool(
            h2s["fraction_media_growth_lost_completely"].gt(0).any()
        )
        if o2_sensitive and h2s_dependent:
            response = "O2-sensitive + H2S-dependent"
        elif o2_sensitive:
            response = "O2-sensitive only"
        elif h2s_dependent:
            response = "H2S-dependent only"
        else:
            response = "Neither"
        rows.append({
            "genome_id": genome_id,
            "o2_sensitive_gt_1pct": o2_sensitive,
            "h2s_dependent_complete_growth_loss": h2s_dependent,
            "counterfactual_response_class": response,
        })
    return pd.DataFrame(rows, columns=output_columns)


def build_genome_positions(
    growth: pd.DataFrame,
    recipe_centroids: pd.DataFrame,
    nutrient_importance: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    hybrid_growth = growth.loc[
        growth["scope"].eq("compartments")
        & growth["family"].eq("hybrid")
    ].copy()
    hybrid_growth["supplemented_biomass_flux"] = pd.to_numeric(
        hybrid_growth["supplemented_biomass_flux"], errors="coerce"
    )
    taxonomy_columns = [
        column for column in (
            "tax_domain", "tax_phylum", "tax_class", "tax_order",
            "tax_family", "tax_genus", "tax_species",
        )
        if column in hybrid_growth
    ]
    # Duplicate nominal media with the same recipe are deliberately collapsed.
    per_recipe = (
        hybrid_growth.groupby(
            ["genome_id", "recipe_id"], as_index=False, dropna=False
        )
        .agg({
            "supplemented_biomass_flux": "median",
            **{column: "first" for column in taxonomy_columns},
        })
        .merge(
            recipe_centroids,
            on="recipe_id",
            how="inner",
            validate="many_to_one",
        )
    )
    if per_recipe.empty:
        raise ValueError("No genome growth rows matched hybrid media recipes")

    weight_rows = []
    position_rows = []
    flags = counterfactual_flags(nutrient_importance)
    for genome_id, frame in per_recipe.groupby("genome_id", sort=True):
        frame = frame.copy()
        values = frame["supplemented_biomass_flux"].to_numpy(float)
        finite = np.isfinite(values)
        minimum = float(np.nanmin(values)) if finite.any() else np.nan
        maximum = float(np.nanmax(values)) if finite.any() else np.nan
        contrast = np.where(finite, np.maximum(values - minimum, 0.0), 0.0)
        if np.sum(contrast) > 1e-12:
            weights = contrast / np.sum(contrast)
            weighting_status = "positive_growth_contrast"
        else:
            available = finite.astype(float)
            weights = available / available.sum()
            weighting_status = "uniform_no_growth_contrast"
        x = frame["PC1_recipe_centroid"].to_numpy(float)
        y = frame["PC2_recipe_centroid"].to_numpy(float)
        pc1 = float(np.sum(weights * x))
        pc2 = float(np.sum(weights * y))
        spread = float(
            np.sqrt(np.sum(weights * ((x - pc1) ** 2 + (y - pc2) ** 2)))
        )
        frame["supplemented_growth_contrast"] = contrast
        frame["genome_recipe_weight"] = weights
        frame["weighting_status"] = weighting_status
        weight_rows.append(frame)
        top = frame.iloc[int(np.argmax(weights))]
        metadata = {
            column: frame[column].dropna().iloc[0]
            if frame[column].notna().any() else pd.NA
            for column in taxonomy_columns
        }
        position_rows.append({
            "genome_id": genome_id,
            **metadata,
            "predicted_PC1": pc1,
            "predicted_PC2": pc2,
            "weighted_recipe_spread": spread,
            "hybrid_recipes_n": int(len(frame)),
            "distinct_hybrid_growth_states": int(
                pd.Series(values[finite]).round(9).nunique()
            ),
            "minimum_supplemented_biomass_flux": minimum,
            "maximum_supplemented_biomass_flux": maximum,
            "growth_range": maximum - minimum,
            "weighting_status": weighting_status,
            "highest_weight_recipe_id": top["recipe_id"],
            "highest_weight_hybrid_labels": top["hybrid_labels"],
            "highest_recipe_weight": float(np.max(weights)),
        })
    positions = pd.DataFrame(position_rows).merge(
        flags, on="genome_id", how="left", validate="one_to_one"
    )
    positions["counterfactual_response_class"] = positions[
        "counterfactual_response_class"
    ].fillna("Not assessed")
    positions["tax_phylum"] = (
        positions.get("tax_phylum", pd.Series(index=positions.index, dtype=object))
        .fillna("Unclassified").astype(str).replace("", "Unclassified")
    )
    return positions, pd.concat(weight_rows, ignore_index=True)


def build_species_positions(
    genome_positions: pd.DataFrame,
    genome_metadata: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Select one quality-ranked model per species and audit all-model medians.

    The representative hierarchy is independent of metabolic simulation:
    MIMAG high precedes medium, followed by greater assembly length, greater
    completeness, lower contamination, greater N50, and genome identifier.
    """
    required_metadata = {
        "Genome_Id", "mimag_tier", "sum_len", "Completeness",
        "Contamination", "N50",
    }
    missing_metadata = sorted(required_metadata.difference(genome_metadata.columns))
    if missing_metadata:
        raise ValueError(
            "Genome metadata lacks fields required for species representative "
            f"selection: {', '.join(missing_metadata)}"
        )

    metadata = genome_metadata[list(required_metadata)].copy()
    metadata["_genome_key"] = metadata["Genome_Id"].map(standardized_genome_id)
    if metadata["_genome_key"].duplicated().any():
        duplicated = sorted(
            metadata.loc[metadata["_genome_key"].duplicated(False), "Genome_Id"]
            .astype(str).unique()
        )
        raise ValueError(
            "Genome metadata contains duplicate normalized identifiers: "
            + ", ".join(duplicated[:10])
        )
    for column in ("sum_len", "Completeness", "Contamination", "N50"):
        metadata[column] = pd.to_numeric(metadata[column], errors="coerce")

    positioned = genome_positions.copy()
    positioned["_genome_key"] = positioned["genome_id"].map(
        standardized_genome_id
    )
    positioned = positioned.merge(
        metadata,
        on="_genome_key",
        how="left",
        validate="one_to_one",
    )
    unmatched = positioned.loc[positioned["Genome_Id"].isna(), "genome_id"]
    if not unmatched.empty:
        raise ValueError(
            "No genome-quality metadata matched modeled genomes: "
            + ", ".join(sorted(unmatched.astype(str).unique())[:10])
        )
    incomplete = positioned[list(required_metadata - {"Genome_Id"})].isna().any(axis=1)
    if incomplete.any():
        affected = positioned.loc[incomplete, "genome_id"].astype(str).tolist()
        raise ValueError(
            "Representative-selection metadata contains missing values for: "
            + ", ".join(sorted(affected)[:10])
        )

    tier_order = {"high": 2, "medium": 1}
    positioned["mimag_tier"] = positioned["mimag_tier"].astype(str).str.lower()
    positioned["representative_eligible"] = positioned["mimag_tier"].isin(
        tier_order
    )
    positioned["_mimag_tier_rank"] = positioned["mimag_tier"].map(tier_order)
    ineligible_species = sorted(
        set(positioned["tax_species"])
        - set(positioned.loc[positioned["representative_eligible"], "tax_species"])
    )
    if ineligible_species:
        raise ValueError(
            "No high- or medium-tier representative candidate was available for: "
            + ", ".join(map(str, ineligible_species))
        )

    positioned = positioned.sort_values(
        [
            "tax_phylum", "tax_species", "_mimag_tier_rank", "sum_len",
            "Completeness", "Contamination", "N50", "genome_id",
        ],
        ascending=[True, True, False, False, False, True, False, True],
        kind="mergesort",
    ).reset_index(drop=True)
    positioned["selection_rank_within_species"] = (
        positioned.groupby(["tax_phylum", "tax_species"], dropna=False)
        .cumcount()
        .add(1)
    )
    positioned["selected_species_representative"] = (
        positioned["representative_eligible"]
        & positioned["selection_rank_within_species"].eq(1)
    )
    positioned["selection_hierarchy"] = (
        "mimag_tier(high>medium)>sum_len>Completeness>"
        "Contamination(lower)>N50>genome_id"
    )

    rows = []
    sensitivity_rows = []
    group_columns = ["tax_phylum", "tax_species"]
    for (phylum, species), frame in positioned.groupby(
        group_columns, dropna=False, sort=True
    ):
        selected = frame.loc[frame["selected_species_representative"]]
        if len(selected) != 1:
            raise ValueError(
                f"Expected one selected representative for {species}; "
                f"found {len(selected)}"
            )
        representative = selected.iloc[0]
        o2_fraction = frame["o2_sensitive_gt_1pct"].fillna(False).astype(bool).mean()
        h2s_fraction = (
            frame["h2s_dependent_complete_growth_loss"]
            .fillna(False).astype(bool).mean()
        )
        o2_consensus = bool(o2_fraction >= 0.5)
        h2s_consensus = bool(h2s_fraction >= 0.5)
        if o2_consensus and h2s_consensus:
            response = "O2-sensitive + H2S-dependent"
        elif o2_consensus:
            response = "O2-sensitive only"
        elif h2s_consensus:
            response = "H2S-dependent only"
        else:
            response = "Neither"
        rows.append({
            "tax_phylum": phylum,
            "tax_species": species,
            "genomes_n": frame["genome_id"].nunique(),
            "representative_genome_id": representative["genome_id"],
            "representative_mimag_tier": representative["mimag_tier"],
            "representative_assembly_length_bp": representative["sum_len"],
            "representative_completeness": representative["Completeness"],
            "representative_contamination": representative["Contamination"],
            "representative_N50": representative["N50"],
            "predicted_PC1": representative["predicted_PC1"],
            "predicted_PC1_q25": frame["predicted_PC1"].quantile(0.25),
            "predicted_PC1_q75": frame["predicted_PC1"].quantile(0.75),
            "predicted_PC2": representative["predicted_PC2"],
            "predicted_PC2_q25": frame["predicted_PC2"].quantile(0.25),
            "predicted_PC2_q75": frame["predicted_PC2"].quantile(0.75),
            "o2_sensitive_genomes_n": int(
                frame["o2_sensitive_gt_1pct"].fillna(False).astype(bool).sum()
            ),
            "o2_sensitive_fraction": o2_fraction,
            "h2s_dependent_genomes_n": int(
                frame["h2s_dependent_complete_growth_loss"]
                .fillna(False).astype(bool).sum()
            ),
            "h2s_dependent_fraction": h2s_fraction,
            "counterfactual_classes_observed_n": (
                frame["counterfactual_response_class"].nunique()
            ),
            "within_species_response_heterogeneous": (
                frame["counterfactual_response_class"].nunique() > 1
            ),
            "all_model_consensus_counterfactual_response_class": response,
            "counterfactual_response_class": representative[
                "counterfactual_response_class"
            ],
            "species_response_rule": "response of selected representative model",
            "position_summary": (
                "highest MIMAG tier then greatest assembly length representative"
            ),
        })
        median_pc1 = float(frame["predicted_PC1"].median())
        median_pc2 = float(frame["predicted_PC2"].median())
        sensitivity_rows.append({
            "tax_phylum": phylum,
            "tax_species": species,
            "genomes_n": frame["genome_id"].nunique(),
            "representative_genome_id": representative["genome_id"],
            "representative_PC1": representative["predicted_PC1"],
            "representative_PC2": representative["predicted_PC2"],
            "all_model_median_PC1": median_pc1,
            "all_model_median_PC2": median_pc2,
            "representative_minus_median_PC1": (
                representative["predicted_PC1"] - median_pc1
            ),
            "representative_minus_median_PC2": (
                representative["predicted_PC2"] - median_pc2
            ),
            "representative_to_median_distance": float(np.hypot(
                representative["predicted_PC1"] - median_pc1,
                representative["predicted_PC2"] - median_pc2,
            )),
        })

    audit_columns = [
        "tax_phylum", "tax_species", "genome_id", "Genome_Id", "mimag_tier",
        "sum_len", "Completeness", "Contamination", "N50",
        "representative_eligible", "selection_rank_within_species",
        "selected_species_representative", "selection_hierarchy",
        "predicted_PC1", "predicted_PC2", "counterfactual_response_class",
    ]
    return (
        pd.DataFrame(rows),
        positioned[audit_columns].copy(),
        pd.DataFrame(sensitivity_rows),
    )


def standardized_genome_id(value: object) -> str:
    """Normalize recruitment and modeling identifiers to a shared xPG key."""
    text = str(value).split("__")[-1]
    text = re.sub(r"\.full$", "", text, flags=re.IGNORECASE)
    return text.lower().replace("_", "-")


def median_ratio_normalize(table: pd.DataFrame) -> tuple[pd.DataFrame, pd.Series]:
    """Normalize nonexclusive genome-recruitment counts across samples."""
    numeric = table.apply(pd.to_numeric, errors="coerce").fillna(0.0).clip(lower=0)
    positive = numeric.where(numeric > 0)
    geometric_means = np.exp(np.log(positive).mean(axis=1, skipna=True))
    geometric_means = geometric_means.replace([np.inf, -np.inf, 0], np.nan)
    ratios = numeric.div(geometric_means, axis=0).where(numeric > 0)
    size_factors = ratios.median(axis=0, skipna=True)
    valid = size_factors.replace([np.inf, -np.inf], np.nan).dropna()
    valid = valid.loc[valid > 0]
    if valid.empty:
        raise ValueError(
            "Genome recruitment normalization produced no positive size factors"
        )
    center = float(np.exp(np.log(valid).mean()))
    size_factors = size_factors / center
    missing = ~np.isfinite(size_factors) | size_factors.le(0)
    if missing.any():
        library_sizes = numeric.sum(axis=0)
        positive_sizes = library_sizes.loc[library_sizes > 0]
        fallback_center = float(np.exp(np.log(positive_sizes).mean()))
        size_factors.loc[missing] = (
            library_sizes.loc[missing] / fallback_center
        )
    if (~np.isfinite(size_factors) | size_factors.le(0)).any():
        raise ValueError("Genome recruitment normalization failed")
    return numeric.div(size_factors, axis=1), size_factors


def build_observed_genome_positions(
    abundance: pd.DataFrame,
    scores: pd.DataFrame,
    genome_positions: pd.DataFrame,
    *,
    sample_col: str,
    genome_col: str,
    value_col: str,
    normalization: str = "auto",
    input_fragments: pd.Series | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Position modeled xPGs by abundance across exact BASINS cruise-depth sites."""
    required = {sample_col, genome_col, value_col}
    if not required.issubset(abundance):
        raise ValueError(
            "Genome abundance table lacks required columns: "
            + ", ".join(sorted(required - set(abundance)))
        )
    depth_col = "Depth_anchored" if "Depth_anchored" in scores else "Depth"
    required_scores = {"Cruise", depth_col, "PC1", "PC2"}
    if not required_scores.issubset(scores):
        raise ValueError(
            "PCA scores lack abundance join columns: "
            + ", ".join(sorted(required_scores - set(scores)))
        )

    work = abundance[[sample_col, genome_col, value_col]].copy()
    work["_genome_key"] = work[genome_col].map(standardized_genome_id)
    work[value_col] = pd.to_numeric(work[value_col], errors="coerce").fillna(0.0)
    matrix = work.pivot_table(
        index="_genome_key", columns=sample_col, values=value_col,
        aggfunc="sum", fill_value=0.0,
    )
    normalized, normalization_factors, normalization_used = (
        normalize_recruitment_counts(
            matrix,
            method=normalization,
            input_fragments=input_fragments,
        )
    )

    samples = pd.DataFrame({sample_col: normalized.columns.astype(str)})
    extracted = samples[sample_col].str.extract(
        r"(?i)SI0*(\d+)_([0-9]+(?:\.[0-9]+)?)m"
    )
    samples["_cruise_key"] = pd.to_numeric(extracted[0], errors="coerce")
    samples["_depth_key"] = pd.to_numeric(extracted[1], errors="coerce")
    coords = scores[["Cruise", depth_col, "PC1", "PC2"]].copy()
    coords["_cruise_key"] = pd.to_numeric(coords["Cruise"], errors="coerce")
    coords["_depth_key"] = pd.to_numeric(coords[depth_col], errors="coerce")
    coords = coords.dropna(
        subset=["_cruise_key", "_depth_key", "PC1", "PC2"]
    ).drop_duplicates(["_cruise_key", "_depth_key"])
    matched = samples.merge(
        coords[["_cruise_key", "_depth_key", "PC1", "PC2"]],
        on=["_cruise_key", "_depth_key"], how="left", validate="many_to_one",
    )
    matched_samples = matched.loc[matched["PC1"].notna(), sample_col].tolist()
    if not matched_samples:
        raise ValueError(
            "No genome-abundance samples had exact BASINS cruise-depth coordinates"
        )
    coordinate_lookup = matched.set_index(sample_col).loc[matched_samples]

    modeled = genome_positions.copy()
    modeled["_genome_key"] = modeled["genome_id"].map(standardized_genome_id)
    if modeled["_genome_key"].duplicated().any():
        raise ValueError("Modeled genome identifiers are not unique after normalization")
    missing_genomes = sorted(set(modeled["_genome_key"]) - set(normalized.index))
    if missing_genomes:
        raise ValueError(
            "Modeled genomes absent from recruitment table: "
            + ", ".join(missing_genomes[:10])
        )

    rows = []
    for _, row in modeled.iterrows():
        key = row["_genome_key"]
        weights = pd.to_numeric(
            normalized.loc[key, matched_samples], errors="coerce"
        ).fillna(0.0)
        detected = weights.gt(0)
        if not detected.any() or weights.loc[detected].sum() <= 0:
            continue
        rows.append({
            "genome_id": row["genome_id"],
            "tax_phylum": row["tax_phylum"],
            "tax_species": row["tax_species"],
            "observed_PC1": float(np.average(
                coordinate_lookup.loc[detected, "PC1"],
                weights=weights.loc[detected],
            )),
            "observed_PC2": float(np.average(
                coordinate_lookup.loc[detected, "PC2"],
                weights=weights.loc[detected],
            )),
            "matched_samples_n": int(len(matched_samples)),
            "detected_matched_samples_n": int(detected.sum()),
            "normalized_abundance_weight_sum": float(weights.loc[detected].sum()),
            "abundance_normalization": normalization_used,
            "abundance_units": (
                "fragments_per_million_input_fragments"
                if normalization_used == "input_fragment_fpm"
                else {
                    "provided_fpkm": "fragments_per_kilobase_per_million",
                    "provided_tpm": "transcripts_per_million",
                }.get(normalization_used, normalization_used)
            ),
            "coordinate_matching": "exact cruise and depth",
            "position_interpretation": (
                "Observed recruitment-abundance-weighted PCA barycenter"
            ),
        })
    genome_observed = pd.DataFrame(rows)
    if genome_observed.empty:
        raise ValueError("No modeled genomes produced observed abundance positions")

    species_rows = []
    for (phylum, species), frame in genome_observed.groupby(
        ["tax_phylum", "tax_species"], dropna=False, sort=True
    ):
        species_rows.append({
            "tax_phylum": phylum,
            "tax_species": species,
            "genomes_n": int(frame["genome_id"].nunique()),
            "observed_PC1": float(frame["observed_PC1"].median()),
            "observed_PC1_q25": float(frame["observed_PC1"].quantile(0.25)),
            "observed_PC1_q75": float(frame["observed_PC1"].quantile(0.75)),
            "observed_PC2": float(frame["observed_PC2"].median()),
            "observed_PC2_q25": float(frame["observed_PC2"].quantile(0.25)),
            "observed_PC2_q75": float(frame["observed_PC2"].quantile(0.75)),
            "position_summary": (
                "median of genome-level observed recruitment barycenters"
            ),
        })
    sample_audit = pd.DataFrame([{
        "abundance_samples_total": int(len(samples)),
        "abundance_samples_with_parsed_coordinates": int(
            samples["_cruise_key"].notna().sum()
        ),
        "abundance_samples_exactly_matched_to_pca": int(len(matched_samples)),
        "modeled_genomes_total": int(modeled["genome_id"].nunique()),
        "modeled_genomes_positioned": int(genome_observed["genome_id"].nunique()),
        "modeled_species_positioned": int(
            genome_observed["tax_species"].nunique()
        ),
        "abundance_normalization": normalization_used,
        "normalization_factor_min": (
            float(normalization_factors.min())
            if not normalization_factors.empty else np.nan
        ),
        "normalization_factor_max": (
            float(normalization_factors.max())
            if not normalization_factors.empty else np.nan
        ),
        "normalization_factor_interpretation": (
            "input read-pair fragments"
            if normalization_used == "input_fragment_fpm"
            else normalization_used
        ),
    }])
    return genome_observed, pd.DataFrame(species_rows), sample_audit


def build_observed_hybrid_affinity(
    abundance: pd.DataFrame,
    scores: pd.DataFrame,
    genome_positions: pd.DataFrame,
    samples_with_hybrid: pd.DataFrame,
    hybrid_centroids: pd.DataFrame,
    *,
    sample_col: str,
    genome_col: str,
    value_col: str,
    normalization: str = "auto",
    input_fragments: pd.Series | None = None,
) -> pd.DataFrame:
    """Estimate each genome's observed affinity to supported hybrid media.

    Recruitment abundance at an exactly matched cruise-depth sample is
    multiplied by that sample's soft hybrid responsibility.  The accumulated
    mass is normalized over hybrid compartments for which a medium recipe was
    constructed.  This deliberately avoids converting either abundance or
    environmental membership to hard labels before defining the expected
    compartment.
    """
    required = {sample_col, genome_col, value_col}
    if not required.issubset(abundance):
        raise ValueError(
            "Genome abundance table lacks required columns: "
            + ", ".join(sorted(required - set(abundance)))
        )
    depth_col = (
        "Depth_anchored" if "Depth_anchored" in samples_with_hybrid else "Depth"
    )
    required_samples = {"Cruise", depth_col}
    if not required_samples.issubset(samples_with_hybrid):
        raise ValueError(
            "Hybrid sample table lacks abundance join columns: "
            + ", ".join(sorted(required_samples - set(samples_with_hybrid)))
        )

    work = abundance[[sample_col, genome_col, value_col]].copy()
    work["_genome_key"] = work[genome_col].map(standardized_genome_id)
    work[value_col] = pd.to_numeric(
        work[value_col], errors="coerce"
    ).fillna(0.0)
    matrix = work.pivot_table(
        index="_genome_key", columns=sample_col, values=value_col,
        aggfunc="sum", fill_value=0.0,
    )
    normalized, _, normalization_used = normalize_recruitment_counts(
        matrix, method=normalization, input_fragments=input_fragments
    )

    assay_samples = pd.DataFrame({sample_col: normalized.columns.astype(str)})
    extracted = assay_samples[sample_col].str.extract(
        r"(?i)SI0*(\d+)_([0-9]+(?:\.[0-9]+)?)m"
    )
    assay_samples["_cruise_key"] = pd.to_numeric(extracted[0], errors="coerce")
    assay_samples["_depth_key"] = pd.to_numeric(extracted[1], errors="coerce")

    responsibility_columns = hybrid_centroids[
        "responsibility_column"
    ].drop_duplicates().tolist()
    environmental = samples_with_hybrid[
        ["Cruise", depth_col] + responsibility_columns
    ].copy()
    environmental["_cruise_key"] = pd.to_numeric(
        environmental["Cruise"], errors="coerce"
    )
    environmental["_depth_key"] = pd.to_numeric(
        environmental[depth_col], errors="coerce"
    )
    environmental = environmental.drop_duplicates(
        ["_cruise_key", "_depth_key"]
    )
    matched = assay_samples.merge(
        environmental[["_cruise_key", "_depth_key"] + responsibility_columns],
        on=["_cruise_key", "_depth_key"], how="left", validate="many_to_one",
    )
    matched_samples = matched.loc[
        matched[responsibility_columns].notna().any(axis=1), sample_col
    ].tolist()
    if not matched_samples:
        raise ValueError(
            "No genome-abundance samples matched hybrid responsibilities"
        )
    responsibility_lookup = matched.set_index(sample_col).loc[matched_samples]

    modeled = genome_positions.copy()
    modeled["_genome_key"] = modeled["genome_id"].map(standardized_genome_id)
    taxonomy_columns = [
        column for column in (
            "tax_domain", "tax_phylum", "tax_class", "tax_order",
            "tax_family", "tax_genus", "tax_species",
        ) if column in modeled
    ]
    centroid_lookup = hybrid_centroids.set_index("responsibility_column")
    rows = []
    for _, model in modeled.iterrows():
        key = model["_genome_key"]
        if key not in normalized.index:
            continue
        abundance_weights = pd.to_numeric(
            normalized.loc[key, matched_samples], errors="coerce"
        ).fillna(0.0)
        if abundance_weights.sum() <= 0:
            continue
        masses = {}
        for responsibility in responsibility_columns:
            membership = pd.to_numeric(
                responsibility_lookup[responsibility], errors="coerce"
            ).fillna(0.0)
            masses[responsibility] = float(
                np.sum(abundance_weights.to_numpy(float)
                       * membership.to_numpy(float))
            )
        total_mass = float(sum(masses.values()))
        if total_mass <= 0:
            continue
        metadata = {
            column: model[column] for column in taxonomy_columns
        }
        for responsibility, mass in masses.items():
            centroid = centroid_lookup.loc[responsibility]
            rows.append({
                "genome_id": model["genome_id"],
                **metadata,
                "compartment": centroid["compartment"],
                "display_label": centroid["display_label"],
                "recipe_id": centroid["recipe_id"],
                "observed_abundance_responsibility_mass": mass,
                "observed_compartment_fraction": mass / total_mass,
                "observed_supported_hybrid_mass_total": total_mass,
                "matched_samples_n": len(matched_samples),
                "abundance_normalization": normalization_used,
                "expected_compartment_basis": (
                    "recruitment abundance x soft hybrid responsibility"
                ),
            })
    result = pd.DataFrame(rows)
    if result.empty:
        raise ValueError("No modeled genomes had observed hybrid affinity")
    return result


def build_species_hybrid_affinity(
    genome_affinity: pd.DataFrame,
) -> pd.DataFrame:
    """Collapse clonal genome records before expected-compartment testing."""
    taxonomy_columns = [
        column for column in (
            "tax_domain", "tax_phylum", "tax_class", "tax_order",
            "tax_family", "tax_genus", "tax_species",
        ) if column in genome_affinity
    ]
    keys = taxonomy_columns + ["compartment", "display_label", "recipe_id"]
    species = genome_affinity.groupby(keys, dropna=False, as_index=False).agg(
        observed_compartment_fraction=("observed_compartment_fraction", "median"),
        genomes_n=("genome_id", "nunique"),
    )
    normalizer = species.groupby(taxonomy_columns, dropna=False)[
        "observed_compartment_fraction"
    ].transform("sum")
    species["observed_compartment_fraction"] = (
        species["observed_compartment_fraction"] / normalizer
    )
    species["observed_affinity_summary"] = (
        "median genome-model fraction, renormalized within species"
    )
    return species


def build_species_recipe_growth(recipe_weights: pd.DataFrame) -> pd.DataFrame:
    taxonomy_columns = [
        column for column in (
            "tax_domain", "tax_phylum", "tax_class", "tax_order",
            "tax_family", "tax_genus", "tax_species",
        ) if column in recipe_weights
    ]
    keys = taxonomy_columns + [
        "recipe_id", "hybrid_compartments", "hybrid_labels"
    ]
    return recipe_weights.groupby(keys, dropna=False, as_index=False).agg(
        supplemented_biomass_flux=("supplemented_biomass_flux", "median"),
        genomes_n=("genome_id", "nunique"),
    )


def expected_compartment_performance(
    affinity: pd.DataFrame,
    recipe_growth: pd.DataFrame,
    *,
    unit_column: str | None,
    top_fraction: float,
) -> pd.DataFrame:
    """Rank the observed expected recipe against unique hybrid recipes."""
    if not 0 < top_fraction < 1:
        raise ValueError("top_fraction must be between zero and one")
    taxonomy_columns = [
        column for column in (
            "tax_domain", "tax_phylum", "tax_class", "tax_order",
            "tax_family", "tax_genus", "tax_species",
        ) if column in affinity and column in recipe_growth
    ]
    unit_keys = [unit_column] if unit_column else taxonomy_columns
    rows = []
    grouped = affinity.groupby(unit_keys, dropna=False, sort=True)
    for unit_value, observed in grouped:
        if unit_column:
            value = unit_value[0] if isinstance(unit_value, tuple) else unit_value
            growth = recipe_growth.loc[recipe_growth[unit_column].eq(value)].copy()
        else:
            values = unit_value if isinstance(unit_value, tuple) else (unit_value,)
            selector = np.ones(len(recipe_growth), dtype=bool)
            for column, value in zip(taxonomy_columns, values):
                if pd.isna(value):
                    selector &= recipe_growth[column].isna().to_numpy()
                else:
                    selector &= recipe_growth[column].eq(value).to_numpy()
            growth = recipe_growth.loc[selector].copy()
        if growth.empty:
            continue
        observed = observed.sort_values(
            ["observed_compartment_fraction", "compartment"],
            ascending=[False, True],
        )
        expected = observed.iloc[0]
        growth["supplemented_biomass_flux"] = pd.to_numeric(
            growth["supplemented_biomass_flux"], errors="coerce"
        )
        growth = growth.dropna(subset=["supplemented_biomass_flux"])
        growth = growth.drop_duplicates("recipe_id")
        expected_growth_rows = growth.loc[
            growth["recipe_id"].astype(str).eq(str(expected["recipe_id"]))
        ]
        if expected_growth_rows.empty or len(growth) < 2:
            continue
        expected_growth = float(
            expected_growth_rows.iloc[0]["supplemented_biomass_flux"]
        )
        ranks = growth["supplemented_biomass_flux"].rank(
            method="average", ascending=False
        )
        expected_index = expected_growth_rows.index[0]
        expected_rank = float(ranks.loc[expected_index])
        media_n = int(len(growth))
        top_n_target = int(math.ceil(top_fraction * media_n))
        top_mask = ranks.le(top_n_target)
        top_media_n = int(top_mask.sum())
        others = growth.loc[
            ~growth["recipe_id"].astype(str).eq(str(expected["recipe_id"])),
            "supplemented_biomass_flux",
        ]
        other_median = float(others.median())
        other_mean = float(others.mean())
        ratio = (
            expected_growth / other_median
            if np.isfinite(other_median) and other_median > 0 else np.nan
        )
        metadata = {
            column: expected[column] for column in taxonomy_columns
        }
        rows.append({
            **({unit_column: value} if unit_column else {}),
            **metadata,
            "analysis_unit": "genome" if unit_column else "species",
            "genomes_n": int(expected.get("genomes_n", 1)),
            "observed_expected_compartment": expected["compartment"],
            "observed_expected_compartment_label": expected["display_label"],
            "observed_expected_recipe_id": expected["recipe_id"],
            "observed_expected_compartment_fraction": float(
                expected["observed_compartment_fraction"]
            ),
            "unique_hybrid_recipes_n": media_n,
            "expected_compartment_growth_rank": expected_rank,
            "expected_compartment_growth_percentile": (
                1.0 - (expected_rank - 1.0) / (media_n - 1.0)
            ),
            "top_growth_fraction_configured": top_fraction,
            "top_growth_rank_target_n": top_n_target,
            "top_growth_recipes_observed_n": top_media_n,
            "expected_compartment_in_top_growth_fraction": bool(
                expected_rank <= top_n_target
            ),
            "expected_supplemented_biomass_flux": expected_growth,
            "nonexpected_median_supplemented_biomass_flux": other_median,
            "nonexpected_mean_supplemented_biomass_flux": other_mean,
            "expected_minus_nonexpected_median_flux": (
                expected_growth - other_median
            ),
            "expected_vs_nonexpected_median_growth_ratio": ratio,
            "top_growth_null_probability": top_media_n / media_n,
            "ranking_scope": "unique supported hybrid-compartment recipes",
            "tie_rule": "average rank; tied recipes retained together",
        })
    return pd.DataFrame(rows)


def expected_compartment_recipe_rank_audit(
    affinity: pd.DataFrame,
    recipe_growth: pd.DataFrame,
    *,
    unit_column: str | None,
    top_fraction: float,
) -> pd.DataFrame:
    """Return every recipe rank underlying expected-compartment summaries."""
    taxonomy_columns = [
        column for column in (
            "tax_domain", "tax_phylum", "tax_class", "tax_order",
            "tax_family", "tax_genus", "tax_species",
        ) if column in affinity and column in recipe_growth
    ]
    unit_keys = [unit_column] if unit_column else taxonomy_columns
    rows = []
    for unit_value, observed in affinity.groupby(
        unit_keys, dropna=False, sort=True
    ):
        if unit_column:
            value = unit_value[0] if isinstance(unit_value, tuple) else unit_value
            growth = recipe_growth.loc[recipe_growth[unit_column].eq(value)].copy()
        else:
            values = unit_value if isinstance(unit_value, tuple) else (unit_value,)
            selector = np.ones(len(recipe_growth), dtype=bool)
            for column, value in zip(taxonomy_columns, values):
                if pd.isna(value):
                    selector &= recipe_growth[column].isna().to_numpy()
                else:
                    selector &= recipe_growth[column].eq(value).to_numpy()
            growth = recipe_growth.loc[selector].copy()
        growth["supplemented_biomass_flux"] = pd.to_numeric(
            growth["supplemented_biomass_flux"], errors="coerce"
        )
        growth = growth.dropna(subset=["supplemented_biomass_flux"])
        growth = growth.drop_duplicates("recipe_id")
        if growth.empty:
            continue
        expected = observed.sort_values(
            ["observed_compartment_fraction", "compartment"],
            ascending=[False, True],
        ).iloc[0]
        growth_rank = growth["supplemented_biomass_flux"].rank(
            method="average", ascending=False
        )
        top_n_target = int(math.ceil(top_fraction * len(growth)))
        for index, recipe in growth.iterrows():
            rank = float(growth_rank.loc[index])
            rows.append({
                **({unit_column: value} if unit_column else {}),
                **{
                    column: expected[column] for column in taxonomy_columns
                },
                "analysis_unit": "genome" if unit_column else "species",
                "recipe_id": recipe["recipe_id"],
                "hybrid_compartments": recipe.get("hybrid_compartments", ""),
                "hybrid_labels": recipe.get("hybrid_labels", ""),
                "supplemented_biomass_flux": recipe[
                    "supplemented_biomass_flux"
                ],
                "growth_rank_average_ties": rank,
                "growth_percentile": (
                    1.0 - (rank - 1.0) / (len(growth) - 1.0)
                    if len(growth) > 1 else 1.0
                ),
                "in_top_growth_fraction": bool(rank <= top_n_target),
                "observed_expected_compartment": expected["compartment"],
                "observed_expected_compartment_label": expected["display_label"],
                "observed_expected_recipe_id": expected["recipe_id"],
                "is_observed_expected_recipe": bool(
                    str(recipe["recipe_id"]) == str(expected["recipe_id"])
                ),
                "top_growth_fraction_configured": top_fraction,
                "top_growth_rank_target_n": top_n_target,
                "ranking_scope": "unique supported hybrid-compartment recipes",
            })
    return pd.DataFrame(rows)


def poisson_binomial_upper_tail(probabilities: np.ndarray, observed: int) -> float:
    """Exact P(X >= observed) for independent, unequal Bernoulli nulls."""
    distribution = np.array([1.0])
    for probability in np.asarray(probabilities, dtype=float):
        distribution = np.convolve(
            distribution, np.array([1.0 - probability, probability])
        )
    return float(distribution[int(observed):].sum())


def benjamini_hochberg(values: pd.Series) -> pd.Series:
    result = pd.Series(np.nan, index=values.index, dtype=float)
    valid = pd.to_numeric(values, errors="coerce").dropna().sort_values()
    if valid.empty:
        return result
    m = len(valid)
    adjusted = valid.to_numpy(float) * m / np.arange(1, m + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result.loc[valid.index] = np.minimum(adjusted, 1.0)
    return result


def summarize_expected_compartment_performance(
    genome_performance: pd.DataFrame,
    species_performance: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    for unit, frame in (
        ("genome", genome_performance), ("species", species_performance)
    ):
        if frame.empty:
            continue
        success = frame["expected_compartment_in_top_growth_fraction"].astype(bool)
        ratios = pd.to_numeric(
            frame["expected_vs_nonexpected_median_growth_ratio"], errors="coerce"
        )
        rows.append({
            "analysis_unit": unit,
            "units_n": len(frame),
            "genome_models_represented_n": (
                frame["genome_id"].nunique() if "genome_id" in frame
                else int(frame.get("genomes_n", pd.Series(dtype=float)).sum())
            ),
            "expected_compartment_top_growth_n": int(success.sum()),
            "expected_compartment_top_growth_fraction": float(success.mean()),
            "median_expected_growth_rank": float(
                frame["expected_compartment_growth_rank"].median()
            ),
            "median_expected_growth_percentile": float(
                frame["expected_compartment_growth_percentile"].median()
            ),
            "median_expected_vs_nonexpected_growth_ratio": float(ratios.median()),
            "units_with_expected_growth_above_nonexpected_median_n": int(
                frame["expected_minus_nonexpected_median_flux"].gt(0).sum()
            ),
            "ranking_scope": "unique supported hybrid-compartment recipes",
        })
    return pd.DataFrame(rows)


def lineage_expected_compartment_enrichment(
    species_performance: pd.DataFrame,
) -> pd.DataFrame:
    """Test lineage/expected-state enrichment with species as replicates."""
    rows = []
    rank_columns = [
        column for column in (
            "tax_phylum", "tax_class", "tax_order", "tax_family", "tax_genus"
        ) if column in species_performance
    ]
    for rank_column in rank_columns:
        grouping = [rank_column, "observed_expected_compartment_label"]
        for (lineage, compartment), frame in species_performance.groupby(
            grouping, dropna=False, sort=True
        ):
            success = frame[
                "expected_compartment_in_top_growth_fraction"
            ].astype(bool)
            probabilities = pd.to_numeric(
                frame["top_growth_null_probability"], errors="coerce"
            ).to_numpy(float)
            observed_n = int(success.sum())
            expected_n = float(probabilities.sum())
            ratios = pd.to_numeric(
                frame["expected_vs_nonexpected_median_growth_ratio"],
                errors="coerce",
            )
            rows.append({
                "taxonomic_rank": rank_column.removeprefix("tax_"),
                "taxonomic_lineage": lineage,
                "observed_expected_compartment_label": compartment,
                "species_n": int(len(frame)),
                "genome_models_represented_n": int(
                    pd.to_numeric(
                        frame.get("genomes_n", pd.Series(1, index=frame.index)),
                        errors="coerce",
                    ).sum()
                ),
                "expected_compartment_top_growth_species_n": observed_n,
                "expected_compartment_top_growth_species_fraction": float(
                    success.mean()
                ),
                "null_expected_top_growth_species_n": expected_n,
                "top_growth_enrichment_ratio": (
                    observed_n / expected_n if expected_n > 0 else np.nan
                ),
                "poisson_binomial_enrichment_p": (
                    poisson_binomial_upper_tail(probabilities, observed_n)
                ),
                "median_expected_vs_nonexpected_growth_ratio": float(
                    ratios.median()
                ),
                "median_expected_minus_nonexpected_flux": float(
                    frame["expected_minus_nonexpected_median_flux"].median()
                ),
                "replicate_unit": "species",
                "null_model": (
                    "expected recipe exchangeable among assessable recipes "
                    "within each species"
                ),
            })
    result = pd.DataFrame(rows)
    if not result.empty:
        result["poisson_binomial_enrichment_q"] = benjamini_hochberg(
            result["poisson_binomial_enrichment_p"]
        )
    return result


def resolve_annotation_collisions(
    fig,
    ax,
    annotations: list,
    forbidden_xy: list[tuple[float, float]],
) -> None:
    """Move labels away from centroids/species points and other label boxes."""
    if not annotations:
        return
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    forbidden_display = [
        ax.transData.transform(point) for point in forbidden_xy
    ]
    shifts = [(0, 0)]
    for radius in (16, 28, 42, 60, 82, 108):
        shifts.extend([
            (0, radius), (0, -radius), (radius, 0), (-radius, 0),
            (radius, radius), (radius, -radius),
            (-radius, radius), (-radius, -radius),
        ])
    accepted = []
    for annotation in annotations:
        original_data = annotation.get_position()
        original_display = ax.transData.transform(original_data)
        initial_box = annotation.get_window_extent(renderer=renderer)
        chosen = (0, 0)
        for dx, dy in shifts:
            candidate = mtransforms.Bbox.from_extents(
                initial_box.x0 + dx,
                initial_box.y0 + dy,
                initial_box.x1 + dx,
                initial_box.y1 + dy,
            )
            padded = mtransforms.Bbox.from_extents(
                candidate.x0 - 6, candidate.y0 - 6,
                candidate.x1 + 6, candidate.y1 + 6,
            )
            covers_point = any(
                padded.x0 <= px <= padded.x1
                and padded.y0 <= py <= padded.y1
                for px, py in forbidden_display
            )
            overlaps_label = any(
                padded.overlaps(previous) for previous in accepted
            )
            if not covers_point and not overlaps_label:
                chosen = (dx, dy)
                accepted.append(padded)
                break
        new_display = (
            original_display[0] + chosen[0],
            original_display[1] + chosen[1],
        )
        annotation.set_position(ax.transData.inverted().transform(new_display))
    fig.canvas.draw()


def repel_marker_positions(
    fig,
    ax,
    x_values,
    y_values,
    *,
    minimum_distance_pixels: float = 8.2,
    maximum_displacement_pixels: float = 30.0,
    iterations: int = 500,
) -> tuple[np.ndarray, np.ndarray]:
    """Return deterministic display-repelled coordinates for dense markers."""
    fig.canvas.draw()
    original = ax.transData.transform(
        np.column_stack([
            np.asarray(x_values, dtype=float),
            np.asarray(y_values, dtype=float),
        ])
    )
    positions = original.copy()
    axes_box = ax.get_window_extent(renderer=fig.canvas.get_renderer())
    margin = minimum_distance_pixels / 2
    lower = np.array([axes_box.x0 + margin, axes_box.y0 + margin])
    upper = np.array([axes_box.x1 - margin, axes_box.y1 - margin])

    for iteration in range(iterations):
        displacement = np.zeros_like(positions)
        maximum_overlap = 0.0
        for left in range(len(positions)):
            for right in range(left + 1, len(positions)):
                delta = positions[right] - positions[left]
                distance = float(np.hypot(delta[0], delta[1]))
                overlap = minimum_distance_pixels - distance
                if overlap <= 0:
                    continue
                maximum_overlap = max(maximum_overlap, overlap)
                if distance < 1e-9:
                    angle = (
                        (left * 37 + right * 71) % 360
                    ) * np.pi / 180.0
                    direction = np.array([np.cos(angle), np.sin(angle)])
                else:
                    direction = delta / distance
                push = 0.52 * overlap * direction
                displacement[left] -= push
                displacement[right] += push

        # A weak spring limits visual distortion while collision forces dominate.
        displacement += 0.015 * (original - positions)
        positions += 0.72 * displacement

        offsets = positions - original
        offset_norm = np.linalg.norm(offsets, axis=1)
        too_far = offset_norm > maximum_displacement_pixels
        if np.any(too_far):
            offsets[too_far] *= (
                maximum_displacement_pixels / offset_norm[too_far]
            )[:, None]
            positions[too_far] = original[too_far] + offsets[too_far]
        positions = np.minimum(np.maximum(positions, lower), upper)
        if maximum_overlap < 0.15 and iteration > 20:
            break

    data_positions = ax.transData.inverted().transform(positions)
    return data_positions[:, 0], data_positions[:, 1]


def explained_axis_labels(explained: pd.DataFrame) -> tuple[str, str]:
    ratios = {}
    for row in explained.itertuples(index=False):
        pc = str(getattr(row, "PC"))
        ratio = float(getattr(row, "explained_variance_ratio"))
        ratios[pc] = ratio
    return (
        f"PC1 ({100 * ratios.get('PC1', np.nan):.1f}%)",
        f"PC2 ({100 * ratios.get('PC2', np.nan):.1f}%)",
    )


def species_position_zoom_limits(
    species_positions: pd.DataFrame,
) -> tuple[tuple[float, float], tuple[float, float]]:
    """Return the shared detail extent used by genome-position biplots."""
    x_values = pd.to_numeric(
        species_positions["predicted_PC1"], errors="coerce"
    ).dropna()
    y_values = pd.to_numeric(
        species_positions["predicted_PC2"], errors="coerce"
    ).dropna()
    if x_values.empty or y_values.empty:
        raise ValueError("Species positions contain no finite PC1/PC2 values")
    x_padding = max(0.12, 0.08 * (x_values.max() - x_values.min()))
    y_padding = max(0.06, 0.12 * (y_values.max() - y_values.min()))
    return (
        (float(x_values.min() - x_padding), float(x_values.max() + x_padding)),
        (float(y_values.min() - y_padding), float(y_values.max() + y_padding)),
    )


def paired_species_position_zoom_limits(
    species_paired: pd.DataFrame,
) -> tuple[tuple[float, float], tuple[float, float]]:
    """Return a padded detail extent containing both paired endpoints."""
    x_values = pd.concat(
        [
            pd.to_numeric(species_paired["predicted_PC1"], errors="coerce"),
            pd.to_numeric(species_paired["observed_PC1"], errors="coerce"),
        ],
        ignore_index=True,
    ).dropna()
    y_values = pd.concat(
        [
            pd.to_numeric(species_paired["predicted_PC2"], errors="coerce"),
            pd.to_numeric(species_paired["observed_PC2"], errors="coerce"),
        ],
        ignore_index=True,
    ).dropna()
    if x_values.empty or y_values.empty:
        raise ValueError("Paired species positions contain no finite PC1/PC2 values")
    x_padding = max(0.12, 0.08 * (x_values.max() - x_values.min()))
    y_padding = max(0.06, 0.12 * (y_values.max() - y_values.min()))
    return (
        (float(x_values.min() - x_padding), float(x_values.max() + x_padding)),
        (float(y_values.min() - y_padding), float(y_values.max() + y_padding)),
    )


def plot_biplot(
    samples: pd.DataFrame,
    centroids: pd.DataFrame,
    species_positions: pd.DataFrame,
    explained: pd.DataFrame,
    output_base: Path,
) -> None:
    install_publication_style()
    compartments = list(centroids["compartment"])
    palette = hybrid_palette(compartments)
    responsibility_columns = {
        row.compartment: row.responsibility_column
        for row in centroids.itertuples(index=False)
    }
    hard = pd.DataFrame({
        compartment: pd.to_numeric(
            samples[column], errors="coerce"
        ).fillna(0.0)
        for compartment, column in responsibility_columns.items()
    }).idxmax(axis=1)

    fig = plt.figure(figsize=(14.6, 8.4))
    grid = fig.add_gridspec(
        2, 2, width_ratios=(3.15, 1.25), height_ratios=(1.15, 1.0),
        wspace=0.30, hspace=0.28,
    )
    ax = fig.add_subplot(grid[:, 0])
    zoom_ax = fig.add_subplot(grid[0, 1])
    legend_ax = fig.add_subplot(grid[1, 1])
    legend_ax.axis("off")
    for compartment in sorted(compartments, key=hybrid_sort_key):
        mask = hard.eq(compartment)
        for target in (ax, zoom_ax):
            target.scatter(
                samples.loc[mask, "PC1"],
                samples.loc[mask, "PC2"],
                s=17,
                color=lighten(palette[compartment], 0.80),
                edgecolor="none",
                alpha=1.0,
                zorder=1,
            )
    ax.axhline(0, linewidth=0.8, color="0.84", zorder=0)
    ax.axvline(0, linewidth=0.8, color="0.84", zorder=0)
    zoom_ax.axhline(0, linewidth=0.8, color="0.84", zorder=0)
    zoom_ax.axvline(0, linewidth=0.8, color="0.84", zorder=0)

    # Freeze the full environmental extent before adding annotations.
    x_values = pd.to_numeric(samples["PC1"], errors="coerce").dropna()
    y_values = pd.to_numeric(samples["PC2"], errors="coerce").dropna()
    x_padding = max(0.35, 0.035 * (x_values.max() - x_values.min()))
    y_padding = max(0.35, 0.035 * (y_values.max() - y_values.min()))
    ax.set_xlim(x_values.min() - x_padding, x_values.max() + x_padding)
    ax.set_ylim(y_values.min() - y_padding, y_values.max() + y_padding)

    centroid_table = centroids.sort_values(
        "compartment", key=lambda values: values.map(hybrid_sort_key)
    ).reset_index(drop=True)
    centroid_table["centroid_number"] = np.arange(1, len(centroid_table) + 1)
    for row in centroid_table.itertuples(index=False):
        for target, size, font_size in (
            (ax, 230, 8.5), (zoom_ax, 175, 7.5)
        ):
            target.scatter(
                row.PC1_centroid, row.PC2_centroid,
                marker="o", s=size, facecolor="white",
                edgecolor="black", linewidth=1.5, zorder=20,
            )
            target.text(
                row.PC1_centroid, row.PC2_centroid,
                str(row.centroid_number), ha="center", va="center",
                fontsize=font_size, fontweight="bold",
                color="black", zorder=21, clip_on=True,
            )

    plot_data = species_positions.copy()
    phylum_totals = plot_data.groupby("tax_phylum")["tax_species"].nunique()
    plot_data["phylum_species_n"] = plot_data["tax_phylum"].map(phylum_totals)
    # Draw species-rich phyla first and rare phyla last.
    plot_data = plot_data.sort_values(
        ["phylum_species_n", "tax_phylum", "tax_species"],
        ascending=[False, True, True],
    ).reset_index(drop=True)
    zoom_xlim, zoom_ylim = species_position_zoom_limits(species_positions)
    zoom_ax.set_xlim(*zoom_xlim)
    zoom_ax.set_ylim(*zoom_ylim)
    main_plot_x, main_plot_y = repel_marker_positions(
        fig, ax, plot_data["predicted_PC1"], plot_data["predicted_PC2"]
    )
    zoom_plot_x, zoom_plot_y = repel_marker_positions(
        fig,
        zoom_ax,
        plot_data["predicted_PC1"],
        plot_data["predicted_PC2"],
        minimum_distance_pixels=8.7,
        maximum_displacement_pixels=28.0,
    )
    for draw_index, row in enumerate(plot_data.itertuples(index=False)):
        marker = RESPONSE_MARKERS[row.counterfactual_response_class]
        size = 48
        zorder = 8 + draw_index / max(1, len(plot_data))
        for target, plot_x, plot_y, multiplier in (
            (ax, main_plot_x[draw_index], main_plot_y[draw_index], 1.0),
            (
                zoom_ax,
                zoom_plot_x[draw_index],
                zoom_plot_y[draw_index],
                1.15,
            ),
        ):
            target.scatter(
                [plot_x], [plot_y],
                marker=marker, s=size * multiplier,
                color=PHYLUM_COLORS.get(
                    row.tax_phylum, PHYLUM_COLORS["Unclassified"]
                ),
                edgecolor="white", linewidth=0.35, alpha=1.0, zorder=zorder,
            )

    present_phyla = (
        species_positions.groupby("tax_phylum")["tax_species"].nunique()
        .sort_values(ascending=False)
    )
    phylum_handles = [
        Line2D(
            [], [], marker="o", linestyle="", markersize=7,
            markerfacecolor=PHYLUM_COLORS.get(phylum, "#7F7F7F"),
            markeredgecolor="white",
            label=f"{phylum} (n={count} species)",
        )
        for phylum, count in present_phyla.items()
    ]
    response_counts = (
        species_positions.groupby(
            "counterfactual_response_class"
        )["tax_species"].nunique()
    )
    response_handles = [
        Line2D(
            [], [], marker=marker, linestyle="", markersize=7,
            markerfacecolor="0.38", markeredgecolor="white",
            label=(
                f"{response.replace('O2', r'$O_2$').replace('H2S', r'$H_2S$')} "
                f"(n={int(response_counts.get(response, 0))})"
            ),
        )
        for response, marker in RESPONSE_MARKERS.items()
        if int(response_counts.get(response, 0)) > 0
    ]
    centroid_handles = [
        Line2D(
            [], [], marker="", linestyle="",
            label=f"{row.centroid_number}. {row.display_label}",
        )
        for row in centroid_table.itertuples(index=False)
    ]
    lineage_legend = legend_ax.legend(
        handles=phylum_handles, title="Genome phylum",
        loc="upper left", bbox_to_anchor=(0.0, 1.03), ncol=2,
        frameon=False, fontsize=6.8, title_fontsize=8.5,
        columnspacing=0.7, handletextpad=0.35,
    )
    legend_ax.add_artist(lineage_legend)
    response_legend = legend_ax.legend(
        handles=response_handles,
        title="Selected representative response",
        loc="upper left", bbox_to_anchor=(0.0, 0.66), frameon=False,
        fontsize=6.8, title_fontsize=8.5,
    )
    legend_ax.add_artist(response_legend)
    legend_ax.legend(
        handles=centroid_handles, title=r"Hybrid O$_2$-GMM centroid",
        loc="lower left", bbox_to_anchor=(0.0, -0.02), ncol=2,
        frameon=False, fontsize=5.9, title_fontsize=8.5,
        handlelength=0.0, handletextpad=0.0,
        labelspacing=0.3, columnspacing=0.8,
    )
    xlabel, ylabel = explained_axis_labels(explained)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(
        r"Quality-selected species metabolic niches across hybrid O$_2$-GMM compartments "
        "(display-repelled markers)"
    )
    ax.set_box_aspect(1.0)
    zoom_ax.set_title("Species-position detail (display-repelled)", fontsize=11)
    zoom_ax.set_xlabel("PC1", fontsize=9)
    zoom_ax.set_ylabel("PC2", fontsize=9)
    zoom_ax.tick_params(labelsize=8)
    zoom_ax.set_box_aspect(1.0)
    save_figure_all_formats(
        fig, output_base, dpi=300
    )
    plt.close(fig)


def plot_predicted_observed_biplot(
    samples: pd.DataFrame,
    centroids: pd.DataFrame,
    species_predicted: pd.DataFrame,
    species_observed: pd.DataFrame,
    explained: pd.DataFrame,
    output_base: Path,
    *,
    observed_label: str = "Observed abundance",
) -> None:
    """Link species metabolic predictions to observed recruitment positions."""
    install_publication_style()
    compartments = list(centroids["compartment"])
    palette = hybrid_palette(compartments)
    responsibility_columns = {
        row.compartment: row.responsibility_column
        for row in centroids.itertuples(index=False)
    }
    hard = pd.DataFrame({
        compartment: pd.to_numeric(samples[column], errors="coerce").fillna(0.0)
        for compartment, column in responsibility_columns.items()
    }).idxmax(axis=1)
    paired = species_predicted.merge(
        species_observed,
        on=["tax_phylum", "tax_species"],
        how="inner",
        validate="one_to_one",
        suffixes=("_predicted", "_observed"),
    )
    if paired.empty:
        raise ValueError(
            "Predicted and observed species positions had no shared species"
        )

    fig = plt.figure(figsize=(14.6, 8.4))
    grid = fig.add_gridspec(
        2, 2, width_ratios=(3.15, 1.25), height_ratios=(1.15, 1.0),
        wspace=0.30, hspace=0.28,
    )
    ax = fig.add_subplot(grid[:, 0])
    zoom_ax = fig.add_subplot(grid[0, 1])
    legend_ax = fig.add_subplot(grid[1, 1])
    legend_ax.axis("off")
    for compartment in sorted(compartments, key=hybrid_sort_key):
        mask = hard.eq(compartment)
        for target in (ax, zoom_ax):
            target.scatter(
                samples.loc[mask, "PC1"], samples.loc[mask, "PC2"],
                s=17, color=lighten(palette[compartment], 0.80),
                edgecolor="none", alpha=1.0, zorder=1,
            )
    for target in (ax, zoom_ax):
        target.axhline(0, linewidth=0.8, color="0.84", zorder=0)
        target.axvline(0, linewidth=0.8, color="0.84", zorder=0)

    centroid_table = centroids.sort_values(
        "compartment", key=lambda values: values.map(hybrid_sort_key)
    ).reset_index(drop=True)
    centroid_table["centroid_number"] = np.arange(1, len(centroid_table) + 1)
    for row in centroid_table.itertuples(index=False):
        for target, size, font_size in (
            (ax, 230, 8.5), (zoom_ax, 175, 7.5)
        ):
            target.scatter(
                row.PC1_centroid, row.PC2_centroid,
                marker="o", s=size, color=palette[row.compartment],
                edgecolor="black", linewidth=1.5, zorder=20,
            )
            target.text(
                row.PC1_centroid, row.PC2_centroid,
                str(row.centroid_number), ha="center", va="center",
                fontsize=font_size, fontweight="bold",
                color="black", zorder=21, clip_on=True,
            )

    phylum_counts = paired.groupby("tax_phylum")["tax_species"].nunique()
    paired["phylum_species_n"] = paired["tax_phylum"].map(phylum_counts)
    paired = paired.sort_values(
        ["phylum_species_n", "tax_phylum", "tax_species"],
        ascending=[False, True, True],
    )
    for draw_index, row in enumerate(paired.itertuples(index=False)):
        color = PHYLUM_COLORS.get(
            row.tax_phylum, PHYLUM_COLORS["Unclassified"]
        )
        zorder = 8 + draw_index / max(1, len(paired))
        for target, scale in ((ax, 1.0), (zoom_ax, 0.88)):
            target.plot(
                [row.predicted_PC1, row.observed_PC1],
                [row.predicted_PC2, row.observed_PC2],
                color=color, linewidth=1.1, alpha=0.78,
                zorder=zorder - 1,
            )
            target.scatter(
                row.predicted_PC1, row.predicted_PC2,
                marker="D", s=76 * scale, color=color,
                edgecolor="white", linewidth=0.6,
                alpha=1.0, zorder=zorder,
            )
            target.scatter(
                row.observed_PC1, row.observed_PC2,
                marker="o", s=112 * scale, color=color,
                edgecolor="black", linewidth=0.9,
                alpha=1.0, zorder=zorder + 0.5,
            )

    x_values = pd.to_numeric(samples["PC1"], errors="coerce").dropna()
    y_values = pd.to_numeric(samples["PC2"], errors="coerce").dropna()
    x_padding = max(0.35, 0.035 * (x_values.max() - x_values.min()))
    y_padding = max(0.35, 0.035 * (y_values.max() - y_values.min()))
    ax.set_xlim(x_values.min() - x_padding, x_values.max() + x_padding)
    ax.set_ylim(y_values.min() - y_padding, y_values.max() + y_padding)
    zoom_xlim, zoom_ylim = paired_species_position_zoom_limits(paired)
    zoom_ax.set_xlim(*zoom_xlim)
    zoom_ax.set_ylim(*zoom_ylim)
    xlabel, ylabel = explained_axis_labels(explained)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(
        f"Predicted metabolic and {observed_label.lower()} positions across hybrid "
        r"O$_2$-GMM environmental space"
    )
    ax.set_box_aspect(1.0)
    zoom_ax.set_title("Species-position detail", fontsize=11)
    zoom_ax.set_xlabel("PC1", fontsize=9)
    zoom_ax.set_ylabel("PC2", fontsize=9)
    zoom_ax.tick_params(labelsize=8)
    zoom_ax.set_box_aspect(1.0)

    phylum_handles = [
        Line2D(
            [], [], marker="o", linestyle="", markersize=7,
            markerfacecolor=PHYLUM_COLORS.get(phylum, "#7F7F7F"),
            markeredgecolor="white",
            label=f"{phylum} (n={count} species)",
        )
        for phylum, count in phylum_counts.sort_values(ascending=False).items()
    ]
    lineage_legend = legend_ax.legend(
        handles=phylum_handles, title="Genome phylum",
        loc="upper left", bbox_to_anchor=(0.0, 1.03), ncol=2,
        frameon=False, fontsize=6.8, title_fontsize=8.5,
        columnspacing=0.7, handletextpad=0.35,
    )
    legend_ax.add_artist(lineage_legend)
    endpoint_handles = [
        Line2D(
            [], [], marker="D", linestyle="", markersize=7,
            markerfacecolor="0.45", markeredgecolor="white",
            label="Model-predicted position",
        ),
        Line2D(
            [], [], marker="o", linestyle="", markersize=8,
            markerfacecolor="0.45", markeredgecolor="black",
            label=f"{observed_label} position",
        ),
        Line2D(
            [], [], linestyle="-", linewidth=1.1, color="0.45",
            label="Same species",
        ),
    ]
    endpoint_legend = legend_ax.legend(
        handles=endpoint_handles, title="Position type",
        loc="upper left", bbox_to_anchor=(0.0, 0.66),
        frameon=False, fontsize=6.8, title_fontsize=8.5,
    )
    legend_ax.add_artist(endpoint_legend)
    centroid_handles = [
        Line2D(
            [], [], marker="o", linestyle="", markersize=9,
            markerfacecolor=palette[row.compartment],
            markeredgecolor="black",
            label=f"{row.centroid_number}. {row.display_label}",
        )
        for row in centroid_table.itertuples(index=False)
    ]
    legend_ax.legend(
        handles=centroid_handles, title=r"Hybrid O$_2$-GMM centroid",
        loc="lower left", bbox_to_anchor=(0.0, -0.02), ncol=2,
        frameon=False, fontsize=5.9, title_fontsize=8.5,
        handletextpad=0.4, labelspacing=0.3, columnspacing=0.55,
    )
    save_figure_all_formats(fig, output_base, dpi=300)
    plt.close(fig)


def _lin_concordance_correlation(x: np.ndarray, y: np.ndarray) -> float:
    """Lin's concordance correlation coefficient for one coordinate axis."""
    if len(x) < 2:
        return np.nan
    variance_x = float(np.var(x, ddof=1))
    variance_y = float(np.var(y, ddof=1))
    covariance = float(np.cov(x, y, ddof=1)[0, 1])
    denominator = variance_x + variance_y + float((x.mean() - y.mean()) ** 2)
    return 2.0 * covariance / denominator if denominator > 0 else np.nan


def build_species_position_agreement(
    species_paired: pd.DataFrame,
    centroids: pd.DataFrame,
    *,
    modality: str,
    permutations: int,
    bootstrap_iterations: int,
    random_state: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Assess absolute predicted-observed agreement at the species level."""
    coordinate_columns = [
        "predicted_PC1", "predicted_PC2", "observed_PC1", "observed_PC2"
    ]
    paired = species_paired.dropna(subset=coordinate_columns).copy()
    if paired.empty:
        raise ValueError("No species had complete predicted and observed positions")
    predicted = paired[["predicted_PC1", "predicted_PC2"]].to_numpy(float)
    observed = paired[["observed_PC1", "observed_PC2"]].to_numpy(float)
    distances = np.linalg.norm(observed - predicted, axis=1)
    paired["predicted_observed_distance"] = distances

    centroid_coordinates = centroids[["PC1_centroid", "PC2_centroid"]].to_numpy(
        float
    )
    centroid_labels = centroids["display_label"].astype(str).to_numpy()
    predicted_nearest = np.argmin(
        np.square(predicted[:, None, :] - centroid_coordinates[None, :, :]).sum(
            axis=2
        ),
        axis=1,
    )
    observed_nearest = np.argmin(
        np.square(observed[:, None, :] - centroid_coordinates[None, :, :]).sum(
            axis=2
        ),
        axis=1,
    )
    paired["predicted_nearest_hybrid_compartment"] = centroid_labels[
        predicted_nearest
    ]
    paired["observed_nearest_hybrid_compartment"] = centroid_labels[
        observed_nearest
    ]
    paired["nearest_hybrid_compartment_agreement"] = (
        predicted_nearest == observed_nearest
    )

    rng = np.random.default_rng(random_state)
    null_mean_distances = np.empty(permutations, dtype=float)
    for iteration in range(permutations):
        permuted_observed = observed[rng.permutation(len(observed))]
        null_mean_distances[iteration] = np.linalg.norm(
            permuted_observed - predicted, axis=1
        ).mean()
    observed_mean_distance = float(distances.mean())
    permutation_p = float(
        (1 + np.count_nonzero(null_mean_distances <= observed_mean_distance))
        / (permutations + 1)
    )

    bootstrap_medians = np.empty(bootstrap_iterations, dtype=float)
    for iteration in range(bootstrap_iterations):
        indices = rng.integers(0, len(distances), size=len(distances))
        bootstrap_medians[iteration] = np.median(distances[indices])
    bootstrap_low, bootstrap_high = np.quantile(
        bootstrap_medians, [0.025, 0.975]
    )

    def correlation(left: np.ndarray, right: np.ndarray) -> float:
        if len(left) < 2 or np.std(left, ddof=1) == 0 or np.std(right, ddof=1) == 0:
            return np.nan
        return float(np.corrcoef(left, right)[0, 1])

    genome_count = pd.to_numeric(
        paired.get("genomes_n_predicted", pd.Series(dtype=float)),
        errors="coerce",
    ).sum()
    null_mean = float(null_mean_distances.mean())
    summary = pd.DataFrame([{
        "modality": modality,
        "analysis_unit": "species",
        "n_species": int(len(paired)),
        "n_genome_models_represented": int(genome_count),
        "mean_paired_distance_pc_units": observed_mean_distance,
        "median_paired_distance_pc_units": float(np.median(distances)),
        "median_distance_bootstrap_ci95_low": float(bootstrap_low),
        "median_distance_bootstrap_ci95_high": float(bootstrap_high),
        "rmse_paired_distance_pc_units": float(np.sqrt(np.mean(distances ** 2))),
        "minimum_paired_distance_pc_units": float(distances.min()),
        "maximum_paired_distance_pc_units": float(distances.max()),
        "null_mean_distance_pc_units": null_mean,
        "paired_distance_reduction_vs_null_fraction": (
            float(1.0 - observed_mean_distance / null_mean)
            if null_mean > 0 else np.nan
        ),
        "species_label_permutation_p": permutation_p,
        "permutation_alternative": "paired mean distance < shuffled mean distance",
        "permutation_iterations": int(permutations),
        "bootstrap_iterations": int(bootstrap_iterations),
        "random_state": int(random_state),
        "mean_observed_minus_predicted_PC1": float(
            (observed[:, 0] - predicted[:, 0]).mean()
        ),
        "mean_observed_minus_predicted_PC2": float(
            (observed[:, 1] - predicted[:, 1]).mean()
        ),
        "PC1_pearson_r": correlation(predicted[:, 0], observed[:, 0]),
        "PC2_pearson_r": correlation(predicted[:, 1], observed[:, 1]),
        "PC1_lin_concordance": _lin_concordance_correlation(
            predicted[:, 0], observed[:, 0]
        ),
        "PC2_lin_concordance": _lin_concordance_correlation(
            predicted[:, 1], observed[:, 1]
        ),
        "nearest_hybrid_compartment_agreement_n": int(
            (predicted_nearest == observed_nearest).sum()
        ),
        "nearest_hybrid_compartment_agreement_fraction": float(
            (predicted_nearest == observed_nearest).mean()
        ),
    }])
    return paired, summary


def process_observed_modality(
    *,
    abundance_path: Path,
    seqkit_path: Path | None,
    normalization: str,
    seqkit_file_col: str,
    seqkit_count_col: str,
    sample_col: str,
    genome_col: str,
    value_col: str,
    scores: pd.DataFrame,
    positions: pd.DataFrame,
    recipe_weights: pd.DataFrame,
    species_positions: pd.DataFrame,
    samples: pd.DataFrame,
    hybrid_centroids: pd.DataFrame,
    explained: pd.DataFrame,
    table_outdir: Path,
    plot_outdir: Path,
    modality: str,
    agreement_permutations: int,
    agreement_bootstrap_iterations: int,
    agreement_random_state: int,
    expected_compartment_top_fraction: float,
) -> list[dict[str, object]]:
    """Normalize and position one metagenome or metatranscriptome assay."""
    if modality == "metagenome":
        genome_name = "genome_observed_abundance_positions.tsv"
        species_name = "species_observed_abundance_positions.tsv"
        paired_name = "species_predicted_observed_positions.tsv"
        agreement_name = "species_predicted_observed_agreement.tsv"
        audit_name = "genome_observed_abundance_position_audit.tsv"
        library_name = "genome_abundance_library_sizes.tsv"
        seqkit_audit_name = "genome_abundance_seqkit_file_audit.tsv"
        plot_name = "genome_hybrid_compartment_predicted_observed_biplot"
        observed_label = "Observed abundance"
    elif modality == "metatranscriptome":
        genome_name = "genome_observed_expression_positions.tsv"
        species_name = "species_observed_expression_positions.tsv"
        paired_name = "species_predicted_observed_expression_positions.tsv"
        agreement_name = "species_predicted_observed_expression_agreement.tsv"
        audit_name = "genome_observed_expression_position_audit.tsv"
        library_name = "genome_expression_library_sizes.tsv"
        seqkit_audit_name = "genome_expression_seqkit_file_audit.tsv"
        plot_name = (
            "genome_hybrid_compartment_predicted_observed_expression_biplot"
        )
        observed_label = "Observed expression"
    else:
        raise ValueError(f"Unsupported observed modality: {modality}")

    input_fragments = None
    if seqkit_path is not None:
        input_fragments, seqkit_audit = parse_seqkit_paired_libraries(
            read_table(seqkit_path),
            file_col=seqkit_file_col,
            count_col=seqkit_count_col,
        )
        pd.DataFrame({
            "sample": input_fragments.index,
            "input_fragments": input_fragments.values,
            "library_layout": "paired_end",
            "denominator_rule": "R1 count after exact R1/R2 agreement",
        }).to_csv(table_outdir / library_name, sep="\t", index=False)
        seqkit_audit.to_csv(
            table_outdir / seqkit_audit_name, sep="\t", index=False
        )

    genome_observed, species_observed, abundance_audit = (
        build_observed_genome_positions(
            read_table(abundance_path),
            scores,
            positions,
            sample_col=sample_col,
            genome_col=genome_col,
            value_col=value_col,
            normalization=normalization,
            input_fragments=input_fragments,
        )
    )
    genome_observed.to_csv(table_outdir / genome_name, sep="\t", index=False)
    species_observed.to_csv(table_outdir / species_name, sep="\t", index=False)
    species_paired = species_positions.merge(
        species_observed,
        on=["tax_phylum", "tax_species"],
        how="inner",
        validate="one_to_one",
        suffixes=("_predicted", "_observed"),
    )
    species_paired["observed_minus_predicted_PC1"] = (
        species_paired["observed_PC1"] - species_paired["predicted_PC1"]
    )
    species_paired["observed_minus_predicted_PC2"] = (
        species_paired["observed_PC2"] - species_paired["predicted_PC2"]
    )
    species_paired, agreement_summary = build_species_position_agreement(
        species_paired,
        hybrid_centroids,
        modality=modality,
        permutations=agreement_permutations,
        bootstrap_iterations=agreement_bootstrap_iterations,
        random_state=agreement_random_state,
    )
    species_paired.to_csv(table_outdir / paired_name, sep="\t", index=False)
    agreement_summary.to_csv(
        table_outdir / agreement_name, sep="\t", index=False
    )
    abundance_audit["modality"] = modality
    abundance_audit["seqkit_library_table"] = (
        str(seqkit_path) if seqkit_path is not None else ""
    )
    abundance_audit.to_csv(table_outdir / audit_name, sep="\t", index=False)
    plot_predicted_observed_biplot(
        samples,
        hybrid_centroids,
        species_positions,
        species_observed,
        explained,
        plot_outdir / plot_name,
        observed_label=observed_label,
    )
    expected_outputs = []
    if modality == "metagenome":
        genome_affinity = build_observed_hybrid_affinity(
            read_table(abundance_path),
            scores,
            positions,
            samples,
            hybrid_centroids,
            sample_col=sample_col,
            genome_col=genome_col,
            value_col=value_col,
            normalization=normalization,
            input_fragments=input_fragments,
        )
        species_affinity = build_species_hybrid_affinity(genome_affinity)
        genome_performance = expected_compartment_performance(
            genome_affinity,
            recipe_weights,
            unit_column="genome_id",
            top_fraction=expected_compartment_top_fraction,
        )
        species_growth = build_species_recipe_growth(recipe_weights)
        species_performance = expected_compartment_performance(
            species_affinity,
            species_growth,
            unit_column=None,
            top_fraction=expected_compartment_top_fraction,
        )
        genome_recipe_ranks = expected_compartment_recipe_rank_audit(
            genome_affinity,
            recipe_weights,
            unit_column="genome_id",
            top_fraction=expected_compartment_top_fraction,
        )
        species_recipe_ranks = expected_compartment_recipe_rank_audit(
            species_affinity,
            species_growth,
            unit_column=None,
            top_fraction=expected_compartment_top_fraction,
        )
        performance_summary = summarize_expected_compartment_performance(
            genome_performance, species_performance
        )
        lineage_enrichment = lineage_expected_compartment_enrichment(
            species_performance
        )
        expected_tables = {
            "genome_observed_hybrid_affinity.tsv": genome_affinity,
            "species_observed_hybrid_affinity.tsv": species_affinity,
            "genome_expected_compartment_performance.tsv": genome_performance,
            "species_expected_compartment_performance.tsv": species_performance,
            "genome_hybrid_recipe_growth_ranks.tsv": genome_recipe_ranks,
            "species_hybrid_recipe_growth_ranks.tsv": species_recipe_ranks,
            "expected_compartment_performance_summary.tsv": performance_summary,
            "lineage_expected_compartment_enrichment.tsv": lineage_enrichment,
        }
        for filename, table in expected_tables.items():
            table.to_csv(table_outdir / filename, sep="\t", index=False)
            expected_outputs.append({
                "source_table": f"derived_{Path(filename).stem}",
                "genomes": genome_affinity["genome_id"].nunique(),
                "rows": len(table),
                "output": filename,
            })
    return [
        {
            "source_table": f"derived_{Path(genome_name).stem}",
            "genomes": genome_observed["genome_id"].nunique(),
            "rows": len(genome_observed),
            "output": genome_name,
        },
        {
            "source_table": f"derived_{Path(species_name).stem}",
            "genomes": genome_observed["genome_id"].nunique(),
            "rows": len(species_observed),
            "output": species_name,
        },
        {
            "source_table": f"derived_{Path(paired_name).stem}",
            "genomes": genome_observed["genome_id"].nunique(),
            "rows": len(species_paired),
            "output": paired_name,
        },
        {
            "source_table": f"derived_{Path(agreement_name).stem}",
            "genomes": genome_observed["genome_id"].nunique(),
            "rows": len(agreement_summary),
            "output": agreement_name,
        },
    ] + expected_outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", required=True, type=Path)
    parser.add_argument("--explained", required=True, type=Path)
    parser.add_argument("--hybrid-assignments", required=True, type=Path)
    parser.add_argument("--media-manifest", required=True, type=Path)
    parser.add_argument("--growth", required=True, type=Path)
    parser.add_argument("--nutrient-importance", required=True, type=Path)
    parser.add_argument(
        "--genome-metadata", required=True, type=Path,
        help=(
            "Genome-quality table used to select one species representative by "
            "MIMAG tier and assembly length"
        ),
    )
    parser.add_argument("--abundance", type=Path)
    parser.add_argument("--abundance-sample-col", default="sample")
    parser.add_argument("--abundance-genome-col", default="genome")
    parser.add_argument("--abundance-value-col", default="read_count")
    parser.add_argument(
        "--abundance-normalization",
        default="auto",
        choices=[
            "auto",
            "input_fragment_fpm",
            "provided_fpkm",
            "provided_tpm",
            "median_ratio",
            "raw",
            "relative",
        ],
    )
    parser.add_argument("--abundance-seqkit", type=Path)
    parser.add_argument("--abundance-seqkit-file-col", default="file")
    parser.add_argument("--abundance-seqkit-count-col", default="num_seqs")
    parser.add_argument("--transcript-abundance", type=Path)
    parser.add_argument("--transcript-abundance-sample-col", default="sample")
    parser.add_argument("--transcript-abundance-genome-col", default="genome")
    parser.add_argument("--transcript-abundance-value-col", default="read_count")
    parser.add_argument(
        "--transcript-abundance-normalization",
        default="auto",
        choices=[
            "auto",
            "input_fragment_fpm",
            "provided_fpkm",
            "provided_tpm",
            "median_ratio",
            "raw",
            "relative",
        ],
    )
    parser.add_argument("--transcript-abundance-seqkit", type=Path)
    parser.add_argument("--transcript-abundance-seqkit-file-col", default="file")
    parser.add_argument(
        "--transcript-abundance-seqkit-count-col", default="num_seqs"
    )
    parser.add_argument("--agreement-permutations", type=int, default=10_000)
    parser.add_argument(
        "--agreement-bootstrap-iterations", type=int, default=5_000
    )
    parser.add_argument("--agreement-random-state", type=int, default=42)
    parser.add_argument(
        "--expected-compartment-top-fraction", type=float, default=0.25,
        help=(
            "Fraction of unique supported hybrid recipes defining the "
            "top-growth set for observed-compartment validation"
        ),
    )
    parser.add_argument("--table-outdir", required=True, type=Path)
    parser.add_argument("--plot-outdir", required=True, type=Path)
    args = parser.parse_args()
    if args.agreement_permutations < 1:
        raise ValueError("--agreement-permutations must be positive")
    if args.agreement_bootstrap_iterations < 1:
        raise ValueError("--agreement-bootstrap-iterations must be positive")
    if not 0 < args.expected_compartment_top_fraction < 1:
        raise ValueError(
            "--expected-compartment-top-fraction must be between zero and one"
        )
    args.table_outdir.mkdir(parents=True, exist_ok=True)
    args.plot_outdir.mkdir(parents=True, exist_ok=True)

    scores = read_table(args.scores)
    explained = read_table(args.explained)
    hybrid = read_table(args.hybrid_assignments)
    manifest = read_table(args.media_manifest)
    growth = read_table(args.growth)
    nutrient_importance = read_table(args.nutrient_importance)
    genome_metadata = read_table(args.genome_metadata)

    samples = merge_sample_inputs(scores, hybrid)
    hybrid_centroids = build_hybrid_centroids(samples, manifest)
    recipe_centroids = build_recipe_centroids(hybrid_centroids)
    positions, weights = build_genome_positions(
        growth, recipe_centroids, nutrient_importance
    )
    species_positions, representative_audit, representative_sensitivity = (
        build_species_positions(positions, genome_metadata)
    )

    hybrid_centroids.to_csv(
        args.table_outdir / "genome_niche_hybrid_centroids.tsv",
        sep="\t", index=False,
    )
    recipe_centroids.to_csv(
        args.table_outdir / "genome_niche_recipe_centroids.tsv",
        sep="\t", index=False,
    )
    weights.to_csv(
        args.table_outdir / "genome_niche_recipe_weights.tsv",
        sep="\t", index=False,
    )
    positions.to_csv(
        args.table_outdir / "genome_predicted_niche_positions.tsv",
        sep="\t", index=False,
    )
    species_positions.to_csv(
        args.table_outdir / "species_predicted_niche_positions.tsv",
        sep="\t", index=False,
    )
    representative_audit.to_csv(
        args.table_outdir / "species_representative_selection_audit.tsv",
        sep="\t", index=False,
    )
    representative_sensitivity.to_csv(
        args.table_outdir / "species_representative_position_sensitivity.tsv",
        sep="\t", index=False,
    )
    plot_biplot(
        samples,
        hybrid_centroids,
        species_positions,
        explained,
        args.plot_outdir / "genome_hybrid_compartment_niche_biplot",
    )
    observed_inventory = []
    if args.abundance is not None:
        observed_inventory.extend(process_observed_modality(
            abundance_path=args.abundance,
            seqkit_path=args.abundance_seqkit,
            normalization=args.abundance_normalization,
            seqkit_file_col=args.abundance_seqkit_file_col,
            seqkit_count_col=args.abundance_seqkit_count_col,
            sample_col=args.abundance_sample_col,
            genome_col=args.abundance_genome_col,
            value_col=args.abundance_value_col,
            scores=scores,
            positions=positions,
            recipe_weights=weights,
            species_positions=species_positions,
            samples=samples,
            hybrid_centroids=hybrid_centroids,
            explained=explained,
            table_outdir=args.table_outdir,
            plot_outdir=args.plot_outdir,
            modality="metagenome",
            agreement_permutations=args.agreement_permutations,
            agreement_bootstrap_iterations=args.agreement_bootstrap_iterations,
            agreement_random_state=args.agreement_random_state,
            expected_compartment_top_fraction=(
                args.expected_compartment_top_fraction
            ),
        ))
    if args.transcript_abundance is not None:
        observed_inventory.extend(process_observed_modality(
            abundance_path=args.transcript_abundance,
            seqkit_path=args.transcript_abundance_seqkit,
            normalization=args.transcript_abundance_normalization,
            seqkit_file_col=args.transcript_abundance_seqkit_file_col,
            seqkit_count_col=args.transcript_abundance_seqkit_count_col,
            sample_col=args.transcript_abundance_sample_col,
            genome_col=args.transcript_abundance_genome_col,
            value_col=args.transcript_abundance_value_col,
            scores=scores,
            positions=positions,
            recipe_weights=weights,
            species_positions=species_positions,
            samples=samples,
            hybrid_centroids=hybrid_centroids,
            explained=explained,
            table_outdir=args.table_outdir,
            plot_outdir=args.plot_outdir,
            modality="metatranscriptome",
            agreement_permutations=args.agreement_permutations,
            agreement_bootstrap_iterations=args.agreement_bootstrap_iterations,
            agreement_random_state=args.agreement_random_state,
            expected_compartment_top_fraction=(
                args.expected_compartment_top_fraction
            ),
        ))
    inventory_path = args.table_outdir / "combined_results_inventory.tsv"
    if inventory_path.is_file():
        inventory = pd.read_csv(inventory_path, sep="\t")
        additions = pd.DataFrame([
            {
                "source_table": "derived_genome_niche_hybrid_centroids",
                "genomes": positions["genome_id"].nunique(),
                "rows": len(hybrid_centroids),
                "output": "genome_niche_hybrid_centroids.tsv",
            },
            {
                "source_table": "derived_genome_niche_recipe_centroids",
                "genomes": positions["genome_id"].nunique(),
                "rows": len(recipe_centroids),
                "output": "genome_niche_recipe_centroids.tsv",
            },
            {
                "source_table": "derived_genome_niche_recipe_weights",
                "genomes": positions["genome_id"].nunique(),
                "rows": len(weights),
                "output": "genome_niche_recipe_weights.tsv",
            },
            {
                "source_table": "derived_genome_predicted_niche_positions",
                "genomes": positions["genome_id"].nunique(),
                "rows": len(positions),
                "output": "genome_predicted_niche_positions.tsv",
            },
            {
                "source_table": "derived_species_predicted_niche_positions",
                "genomes": positions["genome_id"].nunique(),
                "rows": len(species_positions),
                "output": "species_predicted_niche_positions.tsv",
            },
            {
                "source_table": "derived_species_representative_selection_audit",
                "genomes": positions["genome_id"].nunique(),
                "rows": len(representative_audit),
                "output": "species_representative_selection_audit.tsv",
            },
            {
                "source_table": "derived_species_representative_position_sensitivity",
                "genomes": positions["genome_id"].nunique(),
                "rows": len(representative_sensitivity),
                "output": "species_representative_position_sensitivity.tsv",
            },
        ] + observed_inventory)
        inventory = inventory.loc[
            ~inventory["output"].isin(additions["output"])
        ]
        pd.concat([inventory, additions], ignore_index=True).to_csv(
            inventory_path, sep="\t", index=False
        )


if __name__ == "__main__":
    main()
