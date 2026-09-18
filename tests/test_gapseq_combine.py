import subprocess
import sys
from pathlib import Path
from xml.etree import ElementTree

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "processes/gapseq_modeling/combine_gapseq_results.py"


def write_sbml(path, reactions):
    root = ElementTree.Element(
        "{http://www.sbml.org/sbml/level3/version1/core}sbml"
    )
    model = ElementTree.SubElement(root, "model")
    reaction_list = ElementTree.SubElement(model, "listOfReactions")
    for index in range(reactions):
        ElementTree.SubElement(reaction_list, "reaction", id=f"R{index}")
    ElementTree.ElementTree(root).write(path, encoding="utf-8")


def test_combines_source_depth_and_compartment_comparisons(tmp_path):
    results = tmp_path / "results"
    rows = [
        ("depth100", "depth_baseline", "depth", "depth_100m", "r1", 1.0),
        ("depth150", "depth_baseline", "depth", "depth_150m", "r2", 2.0),
        ("oxic", "compartments", "legacy_o2", "oxic", "r3", 0.5),
        ("suboxic", "compartments", "legacy_o2", "suboxic", "r4", 2.0),
        ("gmm0", "compartments", "gmm", "gmm_0", "r3", 0.5),
        ("gmm1", "compartments", "gmm", "gmm_1", "r5", 1.5),
        ("h00", "compartments", "hybrid", "hybrid_c0_g0", "r3", 0.5),
        ("h01", "compartments", "hybrid", "hybrid_c0_g1", "r6", 1.0),
        ("h20", "compartments", "hybrid", "hybrid_c2_g0", "r4", 2.0),
        ("h21", "compartments", "hybrid", "hybrid_c2_g1", "r4", 2.0),
    ]
    for genome_id, scale in (("g100", 1.0), ("g200", 0.8)):
        genome_dir = results / genome_id
        genome_dir.mkdir(parents=True)
        frame = pd.DataFrame([
            {
                "medium": medium,
                "scope": scope,
                "cruise_group": "",
                "family": family,
                "compartment": compartment,
                "recipe_id": recipe,
                "supplemented_biomass_flux": growth * scale,
                "within_genome_percent_of_max": growth / 2.0 * 100,
            }
            for medium, scope, family, compartment, recipe, growth in rows
        ])
        frame.to_csv(genome_dir / "growth_comparison.tsv", sep="\t", index=False)
        pd.DataFrame([{
            "model_id": genome_id,
            "model_file": f"{genome_id}.xml",
            "solver": "glpk",
            "reactions": 12,
            "metabolites": 9,
            "genes": 7,
            "exchanges": 3,
            "objective": "biomass",
            "initial_status": "optimal",
            "initial_objective_value": 1.0,
        }]).to_csv(genome_dir / "model_summary.tsv", sep="\t", index=False)
        pd.DataFrame([{
            "compound_id": "cpd00009",
            "compound_name": "Phosphate",
            "nutrient_source": "basin_observation_trained",
            "media_tested": 2,
            "media_with_numeric_effect": 2,
            "importance_score": 1.0,
            "mean_relative_growth_loss": 1.0,
            "maximum_relative_growth_loss": 1.0,
            "fraction_media_growth_lost_completely": 1.0,
        }]).to_csv(
            genome_dir / "nutrient_importance_summary.tsv",
            sep="\t",
            index=False,
        )

    reconstructions = tmp_path / "reconstructions"
    for genome_id in ("g100", "g200"):
        genome_dir = reconstructions / genome_id
        genome_dir.mkdir(parents=True)
        write_sbml(genome_dir / f"{genome_id}-draft.xml", reactions=10)
        write_sbml(genome_dir / f"{genome_id}.xml", reactions=12)
        pd.DataFrame({
            "compounds": ["cpd1", "cpd2"],
            "name": ["one", "two"],
            "maxFlux": [1, 1],
        }).to_csv(genome_dir / f"{genome_id}-medium.csv", index=False)
        if genome_id == "g200":
            # Publication directories from older runs may contain artifacts
            # produced from a differently named input. Provenance must select
            # the current exact-stem files and audit the ignored candidates.
            write_sbml(genome_dir / "old-g200-draft.xml", reactions=99)
            pd.DataFrame({
                "compounds": ["stale"],
                "name": ["stale"],
                "maxFlux": [1],
            }).to_csv(genome_dir / "old-g200-medium.csv", index=False)
        provenance = {
            "genome_id": genome_id,
            "input_file": f"{genome_id}.faa",
            "taxonomy": "auto",
            "aligner": "diamond",
            "threads": 1,
            "gapfill_medium_source": "gapseq_doall_default",
            "database_source": "gapseq_environment_default",
            "model_file": f"{genome_id}.xml",
        }
        if genome_id == "g200":
            provenance.update({
                "reconstruction_mode": "gapseq_gapfill_reduced_minimum_growth_retry",
                "default_doall_exit_status": 1,
                "effective_gapfill_minimum_growth": 0.001,
                "effective_gapfill_solver": "glpk",
            })
            pd.DataFrame([{
                "genome_id": genome_id,
                "retry_attempted": "true",
                "default_exit_status": 1,
                "retry_exit_status": 0,
                "retry_minimum_growth": 0.001,
                "retry_solver": "glpk",
            }]).to_csv(
                genome_dir / "reconstruction_retry_audit.tsv",
                sep="\t", index=False,
            )
        pd.DataFrame([provenance]).to_csv(
            genome_dir / "reconstruction_provenance.tsv",
            sep="\t",
            index=False,
        )

    # A directory from a previous manifest must be audited but excluded from
    # every current combined result.
    stale_reconstruction = reconstructions / "stale_genome"
    stale_reconstruction.mkdir()
    write_sbml(stale_reconstruction / "stale_genome-draft.xml", reactions=99)
    stale_result = results / "stale_genome"
    stale_result.mkdir()
    pd.DataFrame([{
        "model_id": "stale_genome", "reactions": 99,
    }]).to_csv(stale_result / "model_summary.tsv", sep="\t", index=False)
    pd.DataFrame([{
        "medium": "stale", "scope": "compartments", "family": "hybrid",
        "compartment": "stale", "recipe_id": "stale",
        "supplemented_biomass_flux": 99,
        "within_genome_percent_of_max": 100,
    }]).to_csv(stale_result / "growth_comparison.tsv", sep="\t", index=False)

    manifest = tmp_path / "genomes.tsv"
    pd.DataFrame({
        "label": ["g100", "g200"],
        "filepath": [
            str(tmp_path / "inputs/100m/g100.faa"),
            str(tmp_path / "inputs/200m/g200.faa"),
        ],
    }).to_csv(manifest, sep="\t", index=False)
    genome_metadata = tmp_path / "genome_metadata.tsv"
    pd.DataFrame({
        "Genome_Id": ["g100", "g200"],
        "Domain": ["Bacteria", "Bacteria"],
        "Phylum": ["P1", "P2"],
        "Class": ["C1", "C2"],
        "Order": ["O1", "O2"],
        "Family": ["F1", "F2"],
        "Genus": ["G1", "G2"],
        "Species": ["S1", "S2"],
    }).to_csv(genome_metadata, sep="\t", index=False)
    outdir = tmp_path / "combined"
    subprocess.run([
        sys.executable,
        str(SCRIPT),
        "--results-root",
        str(results),
        "--outdir",
        str(outdir),
        "--reconstructions-root",
        str(reconstructions),
        "--genomes-manifest",
        str(manifest),
        "--genome-metadata",
        str(genome_metadata),
    ], check=True)

    metadata = pd.read_csv(
        outdir / "combined_genome_source_metadata.tsv", sep="\t"
    )
    assert dict(zip(metadata["genome_id"], metadata["source_depth_m"])) == {
        "g100": 100.0,
        "g200": 200.0,
    }
    reconstruction_qc = pd.read_csv(
        outdir / "combined_reconstruction_qc.tsv", sep="\t"
    )
    assert len(reconstruction_qc) == 2
    assert set(reconstruction_qc["reconstruction_status"]) == {"complete"}
    assert set(reconstruction_qc["draft_reactions"]) == {10}
    assert set(reconstruction_qc["final_reactions"]) == {12}
    assert set(reconstruction_qc["net_reactions_added"]) == {2}
    assert set(reconstruction_qc["final_metabolites"]) == {9}
    assert set(reconstruction_qc["final_gene_products"]) == {7}
    assert set(reconstruction_qc["predicted_medium_compounds"]) == {2}
    retried = reconstruction_qc.loc[
        reconstruction_qc["genome_id"].eq("g200")
    ].iloc[0]
    assert retried["reconstruction_mode"] == "gapseq_gapfill_reduced_minimum_growth_retry"
    assert retried["effective_gapfill_minimum_growth"] == 0.001
    assert bool(retried["gapfill_numeric_retry_attempted"])
    assert retried["draft_model_candidates_n"] == 2
    assert retried["ignored_draft_model_files"] == "old-g200-draft.xml"
    assert retried["predicted_medium_candidates_n"] == 2
    assert retried["ignored_predicted_medium_files"] == "old-g200-medium.csv"
    directory_audit = pd.read_csv(
        outdir / "combined_genome_directory_selection_audit.tsv", sep="\t"
    ).set_index("genome_id")
    assert directory_audit.loc[
        "stale_genome", "selection_status"
    ] == "ignored_stale_directory"
    reconstruction_summary = pd.read_csv(
        outdir / "combined_reconstruction_summary.tsv", sep="\t"
    ).set_index("metric")
    assert reconstruction_summary.loc[
        "complete_reconstructions", "count"
    ] == 2
    assert reconstruction_summary.loc[
        "numeric_retry_reconstructions", "count"
    ] == 1
    assert reconstruction_summary.loc["final_reactions", "median"] == 12
    assert reconstruction_summary.loc["net_reactions_added", "median"] == 2
    assert reconstruction_summary.loc[
        "predicted_medium_compounds", "median"
    ] == 2
    counterfactual = pd.read_csv(
        outdir / "combined_counterfactual_taxonomic_summary.tsv", sep="\t"
    )
    phosphate = counterfactual.loc[
        counterfactual["compound_name"].eq("Phosphate")
    ].iloc[0]
    assert phosphate["models_tested"] == 2
    assert phosphate["models_with_importance_gt_1pct"] == 2
    assert phosphate["models_with_complete_growth_loss"] == 2
    assert "P1 (n=1)" in phosphate["complete_loss_phylum_counts"]
    assert "P1 (n=1)" in phosphate[
        "importance_gt_1pct_phylum_counts"
    ]
    lineage = pd.read_csv(
        outdir / "combined_counterfactual_lineage_summary.tsv", sep="\t"
    )
    phyla = lineage.loc[lineage["taxonomic_rank"].eq("phylum")]
    assert set(phyla["lineage"]) == {"P1", "P2"}
    comparison = pd.read_csv(
        outdir / "combined_compartment_system_comparison.tsv", sep="\t"
    )
    assert set(comparison["analysis_subset"]) == {
        "all_genomes", "source_depth_100m", "source_depth_200m"
    }
    all_genomes = comparison.loc[
        comparison["analysis_subset"].eq("all_genomes")
    ]
    assert set(all_genomes["family"]) == {"depth", "legacy_o2", "gmm", "hybrid"}

    hybrid = pd.read_csv(
        outdir / "combined_hybrid_incremental_resolution.tsv", sep="\t"
    )
    all_hybrid = hybrid.loc[hybrid["analysis_subset"].eq("all_genomes")]
    oxic = all_hybrid.loc[all_hybrid["oxygen_compartment"].eq("oxic")].iloc[0]
    suboxic = all_hybrid.loc[
        all_hybrid["oxygen_compartment"].eq("suboxic")
    ].iloc[0]
    assert oxic["fraction_genomes_with_added_resolution"] == 1.0
    assert suboxic["fraction_genomes_with_added_resolution"] == 0.0
