import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


SCRIPT = (
    Path(__file__).parents[1]
    / "processes"
    / "gapseq_modeling"
    / "plot_genome_niche_biplot.py"
)
SPEC = importlib.util.spec_from_file_location("plot_genome_niche_biplot", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def test_observed_genome_positions_use_exact_cruise_depth_coordinates():
    abundance = pd.DataFrame({
        "sample": ["SI001_10m", "SI001_20m", "SI001_10m", "SI001_20m"],
        "genome": ["10m__Genome_A", "10m__Genome_A", "20m__Genome_B", "20m__Genome_B"],
        "read_count": [10, 30, 20, 20],
    })
    scores = pd.DataFrame({
        "Cruise": [1, 1],
        "Depth_anchored": [10, 20],
        "PC1": [0.0, 2.0],
        "PC2": [0.0, 4.0],
    })
    modeled = pd.DataFrame({
        "genome_id": ["genome-a", "genome-b"],
        "tax_phylum": ["P1", "P1"],
        "tax_species": ["S1", "S1"],
    })

    genomes, species, audit = MODULE.build_observed_genome_positions(
        abundance,
        scores,
        modeled,
        sample_col="sample",
        genome_col="genome",
        value_col="read_count",
        normalization="median_ratio",
    )

    assert set(genomes["genome_id"]) == {"genome-a", "genome-b"}
    assert len(species) == 1
    assert audit.loc[0, "abundance_samples_exactly_matched_to_pca"] == 2


def test_seqkit_pairs_create_fragment_denominators_and_fpm():
    seqkit = pd.DataFrame({
        "file": [
            "fastq_renamed_metag/SI034_100m/SI034_100m_pe.1.fq.gz",
            "fastq_renamed_metag/SI034_100m/SI034_100m_pe.2.fq.gz",
            "fastq_renamed_metag/SI034_10m/SI034_10m_R1_001.fastq.gz",
            "fastq_renamed_metag/SI034_10m/SI034_10m_R2_001.fastq.gz",
        ],
        "num_seqs": [40_000_000, 40_000_000, 20_000_000, 20_000_000],
    })
    fragments, audit = MODULE.parse_seqkit_paired_libraries(seqkit)

    assert fragments.to_dict() == {
        "SI034_100m": 40_000_000.0,
        "SI034_10m": 20_000_000.0,
    }
    assert len(audit) == 4
    counts = pd.DataFrame(
        {"SI034_100m": [4_000, 2_000], "SI034_10m": [1_000, 500]},
        index=["g1", "g2"],
    )
    normalized, denominators, method = MODULE.normalize_recruitment_counts(
        counts, method="input_fragment_fpm", input_fragments=fragments
    )
    assert method == "input_fragment_fpm"
    assert denominators["SI034_100m"] == 40_000_000
    assert normalized.loc["g1", "SI034_100m"] == 100.0
    assert normalized.loc["g1", "SI034_10m"] == 50.0


def test_seqkit_pair_count_mismatch_is_rejected():
    seqkit = pd.DataFrame({
        "file": ["SI034_100m_pe.1.fq.gz", "SI034_100m_pe.2.fq.gz"],
        "num_seqs": [100, 99],
    })
    with pytest.raises(ValueError, match="R1/R2 sequence counts disagree"):
        MODULE.parse_seqkit_paired_libraries(seqkit)


def test_identical_duplicate_seqkit_references_are_collapsed():
    seqkit = pd.DataFrame({
        "file": [
            "old/SI060_200m_pe.1.fq.gz",
            "old/SI060_200m_pe.2.fq.gz",
            "current/SI060_200m_pe.1.fq.gz",
            "current/SI060_200m_pe.2.fq.gz",
        ],
        "num_seqs": [100, 100, 100, 100],
    })
    fragments, audit = MODULE.parse_seqkit_paired_libraries(seqkit)

    assert fragments.to_dict() == {"SI060_200m": 100.0}
    assert len(audit) == 2
    assert set(audit["duplicate_source_rows"]) == {2}
    assert set(audit["pair_validation"]) == {
        "matched_duplicate_reference_collapsed"
    }


def test_auto_normalization_does_not_silently_fallback_without_seqkit():
    counts = pd.DataFrame({"S1": [10]}, index=["g1"])
    with pytest.raises(ValueError, match="requires a paired-end SeqKit"):
        MODULE.normalize_recruitment_counts(counts, method="auto")


def test_provided_fpkm_is_not_renormalized():
    table = pd.DataFrame({"S1": [0.25, 2.0], "S2": [1.5, 0.0]}, index=["g1", "g2"])
    normalized, factors, method = MODULE.normalize_recruitment_counts(
        table, method="provided_fpkm"
    )
    pd.testing.assert_frame_equal(normalized, table)
    assert factors.empty
    assert method == "provided_fpkm"


def test_provided_tpm_rejects_negative_values():
    table = pd.DataFrame({"S1": [1.0, -0.1]}, index=["g1", "g2"])
    with pytest.raises(ValueError, match="negative"):
        MODULE.normalize_recruitment_counts(table, method="provided_tpm")


def test_expected_compartment_rank_and_ratio_use_unique_hybrid_recipes():
    affinity = pd.DataFrame({
        "genome_id": ["g1", "g1"],
        "tax_phylum": ["P1", "P1"],
        "tax_species": ["S1", "S1"],
        "compartment": ["hybrid_c0_g0", "hybrid_c2_g1"],
        "display_label": ["oxic-GMM0", "suboxic-GMM1"],
        "recipe_id": ["r1", "r2"],
        "observed_compartment_fraction": [0.8, 0.2],
    })
    growth = pd.DataFrame({
        "genome_id": ["g1", "g1", "g1", "g1"],
        "tax_phylum": ["P1"] * 4,
        "tax_species": ["S1"] * 4,
        "recipe_id": ["r1", "r2", "r3", "r4"],
        "hybrid_compartments": ["c1", "c2", "c3", "c4"],
        "hybrid_labels": ["C1", "C2", "C3", "C4"],
        "supplemented_biomass_flux": [4.0, 1.0, 2.0, 3.0],
    })
    result = MODULE.expected_compartment_performance(
        affinity, growth, unit_column="genome_id", top_fraction=0.25
    )

    assert len(result) == 1
    row = result.iloc[0]
    assert row["observed_expected_recipe_id"] == "r1"
    assert row["expected_compartment_growth_rank"] == 1.0
    assert bool(row["expected_compartment_in_top_growth_fraction"])
    assert row["expected_vs_nonexpected_median_growth_ratio"] == 2.0


def test_flat_growth_does_not_create_false_top_quartile_match():
    affinity = pd.DataFrame({
        "genome_id": ["g1"],
        "tax_phylum": ["P1"],
        "tax_species": ["S1"],
        "compartment": ["hybrid_c0_g0"],
        "display_label": ["oxic-GMM0"],
        "recipe_id": ["r1"],
        "observed_compartment_fraction": [1.0],
    })
    growth = pd.DataFrame({
        "genome_id": ["g1"] * 4,
        "tax_phylum": ["P1"] * 4,
        "tax_species": ["S1"] * 4,
        "recipe_id": ["r1", "r2", "r3", "r4"],
        "hybrid_compartments": ["c1", "c2", "c3", "c4"],
        "hybrid_labels": ["C1", "C2", "C3", "C4"],
        "supplemented_biomass_flux": [1.0] * 4,
    })
    result = MODULE.expected_compartment_performance(
        affinity, growth, unit_column="genome_id", top_fraction=0.25
    )

    assert result.iloc[0]["expected_compartment_growth_rank"] == 2.5
    assert not bool(
        result.iloc[0]["expected_compartment_in_top_growth_fraction"]
    )
    assert result.iloc[0]["top_growth_null_probability"] == 0.0


def test_poisson_binomial_enrichment_tail_is_exact():
    # With two fair Bernoulli trials, P(X >= 2) = 0.25.
    assert MODULE.poisson_binomial_upper_tail(np.array([0.5, 0.5]), 2) == 0.25


def test_observed_hybrid_affinity_combines_abundance_and_soft_membership():
    abundance = pd.DataFrame({
        "sample": ["SI001_10m", "SI001_20m"],
        "genome": ["Genome_A", "Genome_A"],
        "read_count": [3.0, 1.0],
    })
    scores = pd.DataFrame({
        "Cruise": [1, 1], "Depth_anchored": [10, 20],
        "PC1": [0.0, 1.0], "PC2": [0.0, 1.0],
    })
    samples = pd.DataFrame({
        "Cruise": [1, 1], "Depth_anchored": [10, 20],
        "hyb_C0_G0": [0.8, 0.2], "hyb_C2_G1": [0.2, 0.8],
    })
    modeled = pd.DataFrame({
        "genome_id": ["genome-a"], "tax_phylum": ["P1"],
        "tax_species": ["S1"],
    })
    centroids = pd.DataFrame({
        "compartment": ["hybrid_c0_g0", "hybrid_c2_g1"],
        "display_label": ["oxic-GMM0", "suboxic-GMM1"],
        "responsibility_column": ["hyb_C0_G0", "hyb_C2_G1"],
        "recipe_id": ["r1", "r2"],
    })

    affinity = MODULE.build_observed_hybrid_affinity(
        abundance, scores, modeled, samples, centroids,
        sample_col="sample", genome_col="genome", value_col="read_count",
        normalization="raw",
    )

    fractions = affinity.set_index("recipe_id")["observed_compartment_fraction"]
    # r1 mass = 3*0.8 + 1*0.2 = 2.6; total supported mass = 4.
    assert fractions["r1"] == pytest.approx(0.65)
    assert fractions["r2"] == pytest.approx(0.35)


def test_expected_compartment_outputs_are_wired_into_nextflow():
    pipeline = (Path(__file__).parents[1] / "basin_pipeline.nf").read_text()
    assert "--expected-compartment-top-fraction" in pipeline
    for filename in (
        "genome_expected_compartment_performance.tsv",
        "species_expected_compartment_performance.tsv",
        "expected_compartment_performance_summary.tsv",
        "lineage_expected_compartment_enrichment.tsv",
    ):
        assert filename in pipeline


def test_species_position_agreement_uses_species_as_analysis_units():
    paired = pd.DataFrame({
        "tax_phylum": ["P1", "P1", "P2"],
        "tax_species": ["S1", "S2", "S3"],
        "genomes_n_predicted": [2, 3, 1],
        "predicted_PC1": [0.0, 1.0, 0.0],
        "predicted_PC2": [0.0, 0.0, 1.0],
        "observed_PC1": [0.0, 1.0, 0.0],
        "observed_PC2": [0.0, 0.0, 1.0],
    })
    centroids = pd.DataFrame({
        "display_label": ["oxic-GMM0", "suboxic-GMM1", "anoxic-GMM2"],
        "PC1_centroid": [0.0, 1.0, 0.0],
        "PC2_centroid": [0.0, 0.0, 1.0],
    })

    audited, summary = MODULE.build_species_position_agreement(
        paired,
        centroids,
        modality="metagenome",
        permutations=100,
        bootstrap_iterations=100,
        random_state=42,
    )

    assert summary.loc[0, "analysis_unit"] == "species"
    assert summary.loc[0, "n_species"] == 3
    assert summary.loc[0, "n_genome_models_represented"] == 6
    assert summary.loc[0, "mean_paired_distance_pc_units"] == 0.0
    assert summary.loc[0, "nearest_hybrid_compartment_agreement_fraction"] == 1.0
    assert np.all(audited["nearest_hybrid_compartment_agreement"])


def test_paired_zoom_contains_every_predicted_and_observed_endpoint():
    paired = pd.DataFrame({
        "predicted_PC1": [0.0, 1.0],
        "predicted_PC2": [0.2, 0.4],
        "observed_PC1": [-2.0, 3.0],
        "observed_PC2": [-1.0, 2.0],
    })

    x_limits, y_limits = MODULE.paired_species_position_zoom_limits(paired)

    assert x_limits[0] < -2.0
    assert x_limits[1] > 3.0
    assert y_limits[0] < -1.0
    assert y_limits[1] > 2.0
