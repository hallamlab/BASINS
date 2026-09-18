#!/usr/bin/env python3
"""Build non-interpretive, integrated BASINS summary tables."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


MODULES = {
    "biochem_processing": "biochemical_processing",
    "env_missingness_sensitivity": "missingness_sensitivity",
    "env_pca": "environmental_pca",
    "env_compartments_selectk": "compartment_selection",
    "env_compartments_gmm": "gmm_compartments",
    "env_o2_soft_compartments": "oxygen_compartments",
    "env_hybrid_soft_compartments": "hybrid_compartments",
    "env_compare_compartments": "compartment_comparison",
    "env_o2_split_by_gmm": "oxygen_gmm_subcompartments",
    "env_stratification_index": "stratification",
    "env_state_transitions": "state_transitions",
    "env_succession_graphs": "succession_graphs",
    "env_compartment_feature_assoc": "feature_associations",
    "eof_pca": "eof_analysis",
    "eof_states": "eof_states",
    "eof_plots": "eof_modes",
    "gapseq_media": "gapseq_media",
    "genome_modeling": "genome_modeling",
}

KEY_OUTPUTS = (
    ("Best available oxygen and density", "biochem_processing/02_oxygen_best_available_density_RJM.tsv"),
    ("Selected missingness cutoff", "env_missingness_sensitivity/SELECTED_MISSINGNESS_CUTOFF.txt"),
    ("GMM assignments", "env_compartments_gmm/tables/compartments_assignments_smoothed.csv"),
    ("Oxygen assignments", "env_o2_soft_compartments/tables/o2_compartments_assignments_smoothed.csv"),
    ("Oxygen/GMM subcompartments", "env_o2_split_by_gmm/tables/merged_o2_split_by_gmm.csv"),
    ("gapseq media manifest", "gapseq_media/tables/gapseq_media_manifest.tsv"),
    ("Genome model comparison inventory", "genome_modeling/combined/combined_results_inventory.tsv"),
    ("Genome reconstruction QC", "genome_modeling/combined/combined_reconstruction_qc.tsv"),
    ("Genome reconstruction summary", "genome_modeling/combined/combined_reconstruction_summary.tsv"),
    ("Counterfactual taxonomic summary", "genome_modeling/combined/combined_counterfactual_taxonomic_summary.tsv"),
    ("Counterfactual lineage summary", "genome_modeling/combined/combined_counterfactual_lineage_summary.tsv"),
    ("Genome compartment-system comparison", "genome_modeling/combined/combined_compartment_system_comparison.tsv"),
    ("Hybrid incremental metabolic resolution", "genome_modeling/combined/combined_hybrid_incremental_resolution.tsv"),
    ("Predicted genome niche positions", "genome_modeling/combined/genome_predicted_niche_positions.tsv"),
    ("Species-level predicted niche positions", "genome_modeling/combined/species_predicted_niche_positions.tsv"),
    ("Observed genome abundance positions", "genome_modeling/combined/genome_observed_abundance_positions.tsv"),
    ("Species-level observed abundance positions", "genome_modeling/combined/species_observed_abundance_positions.tsv"),
    ("Paired predicted and observed species positions", "genome_modeling/combined/species_predicted_observed_positions.tsv"),
    ("Genome niche recipe weights", "genome_modeling/combined/genome_niche_recipe_weights.tsv"),
    ("Expected-compartment performance summary", "genome_modeling/combined/expected_compartment_performance_summary.tsv"),
    ("Species expected-compartment performance", "genome_modeling/combined/species_expected_compartment_performance.tsv"),
    ("Lineage expected-compartment enrichment", "genome_modeling/combined/lineage_expected_compartment_enrichment.tsv"),
    ("Stratification time series", "env_stratification_index/stratification_timeseries.tsv"),
)


def count_rows(path: Path) -> int | None:
    if not path.is_file():
        return None
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        sample = handle.read(8192)
        handle.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        except csv.Error:
            dialect = csv.excel_tab if path.suffix.lower() == ".tsv" else csv.excel
        rows = sum(1 for _ in csv.reader(handle, dialect))
    return max(0, rows - 1)


def write_tsv(path: Path, header: tuple[str, ...], rows: list[tuple[object, ...]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    root = args.input_root.resolve()
    output = args.output_dir.resolve()

    inventory = []
    for source_name, module_name in MODULES.items():
        source = root / source_name
        files = sorted(path for path in source.rglob("*") if path.is_file()) if source.is_dir() else []
        inventory.append((module_name, source_name, len(files), sum(path.stat().st_size for path in files)))
    write_tsv(
        output / "basin_module_inventory.tsv",
        ("module", "source_directory", "file_count", "total_size_bytes"),
        inventory,
    )

    important = []
    for label, relative in KEY_OUTPUTS:
        path = root / relative
        important.append((label, relative, path.is_file(), count_rows(path), path.stat().st_size if path.is_file() else None))
    write_tsv(
        output / "basin_key_outputs.tsv",
        ("description", "relative_path", "present", "data_rows", "size_bytes"),
        important,
    )

    selected_k = root / "env_compartments_selectk" / "SELECTED_K.txt"
    selected_missingness = (
        root / "env_missingness_sensitivity" / "SELECTED_MISSINGNESS_CUTOFF.txt"
    )
    overview = [
        ("modules_with_outputs", sum(row[2] > 0 for row in inventory)),
        ("total_module_files", sum(row[2] for row in inventory)),
        ("key_outputs_present", sum(bool(row[2]) for row in important)),
        ("selected_gmm_k", selected_k.read_text().strip() if selected_k.is_file() else ""),
        (
            "selected_missingness_cutoff",
            selected_missingness.read_text().strip()
            if selected_missingness.is_file()
            else "",
        ),
    ]
    write_tsv(output / "basin_run_overview.tsv", ("metric", "value"), overview)


if __name__ == "__main__":
    main()
