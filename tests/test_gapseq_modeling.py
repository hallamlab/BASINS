import importlib.util
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest


cobra = pytest.importorskip("cobra")
from cobra import Model, Reaction, Metabolite
from cobra.io import write_sbml_model


SCRIPT = Path(__file__).parents[1] / "processes" / "gapseq_modeling" / "compare_gapseq_media.py"
SPEC = importlib.util.spec_from_file_location("compare_gapseq_media", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def exchange(reaction_id, metabolite_id, coefficient, bounds=(-1000, 1000)):
    metabolite = Metabolite(metabolite_id, compartment="e")
    reaction = Reaction(reaction_id, lower_bound=bounds[0], upper_bound=bounds[1])
    reaction.add_metabolites({metabolite: coefficient})
    return reaction


def test_modelseed_id_extraction_handles_annotations():
    assert MODULE.modelseed_ids({"seed.compound": ["cpd00007", "CPD00209"]}) == {
        "cpd00007", "cpd00209"
    }


def test_medium_validation_rejects_duplicates_and_negative_flux(tmp_path):
    path = tmp_path / "bad.csv"
    pd.DataFrame({
        "compounds": ["cpd00007", "cpd00007"],
        "name": ["O2", "O2"],
        "maxFlux": [1, -1],
    }).to_csv(path, index=False)
    with pytest.raises(ValueError):
        MODULE.validate_medium(path)


def test_negative_and_positive_uptake_orientation():
    model = Model("orientation")
    negative = exchange("EX_cpd00007_e", "cpd00007_e", -1)
    positive = exchange("EX_cpd00209_e", "cpd00209_e", 1)
    model.add_reactions([negative, positive])
    index, details = MODULE.exchange_index(model)
    MODULE.close_all_uptake(model, details)
    assert negative.lower_bound == 0
    assert positive.upper_bound == 0
    medium = pd.DataFrame({
        "compounds": ["cpd00007", "cpd00209"],
        "name": ["O2", "Nitrate"],
        "maxFlux": [3.0, 4.0],
    })
    mapping, supplied = MODULE.apply_medium(model, medium, index, details, "test")
    assert negative.lower_bound == -3.0
    assert positive.upper_bound == 4.0
    assert supplied == {"cpd00007", "cpd00209"}
    assert set(mapping.mapping_status) == {"mapped"}


def test_zero_flux_closes_uptake_without_closing_secretion():
    model = Model("zero")
    reaction = exchange("EX_cpd00007_e", "cpd00007_e", -1)
    model.add_reactions([reaction])
    index, details = MODULE.exchange_index(model)
    MODULE.close_all_uptake(model, details)
    medium = pd.DataFrame({"compounds": ["cpd00007"], "name": ["O2"], "maxFlux": [0.0]})
    MODULE.apply_medium(model, medium, index, details, "zero")
    assert reaction.lower_bound == 0
    assert reaction.upper_bound == 1000


def test_unmapped_compound_is_reported():
    model = Model("unmapped")
    reaction = exchange("EX_cpd00007_e", "cpd00007_e", -1)
    model.add_reactions([reaction])
    index, details = MODULE.exchange_index(model)
    medium = pd.DataFrame({"compounds": ["cpd00209"], "name": ["Nitrate"], "maxFlux": [1.0]})
    mapping, supplied = MODULE.apply_medium(model, medium, index, details, "test")
    assert mapping.iloc[0].mapping_status == "no_exchange_found"
    assert not supplied


def test_model_copy_prevents_bound_leakage():
    model = Model("copy")
    reaction = exchange("EX_cpd00007_e", "cpd00007_e", -1)
    model.add_reactions([reaction])
    copied = model.copy()
    copied.reactions.get_by_id(reaction.id).lower_bound = 0
    assert model.reactions.get_by_id(reaction.id).lower_bound == -1000


def test_end_to_end_fba_fva_pairwise_and_counterfactual(tmp_path):
    model = Model("toy")
    carbon = Metabolite("cpd00027_e", name="D-Glucose", compartment="e")
    uptake = Reaction("EX_cpd00027_e", lower_bound=-1000, upper_bound=1000)
    uptake.add_metabolites({carbon: -1})
    biomass = Reaction("BIOMASS", lower_bound=0, upper_bound=1000)
    biomass.add_metabolites({carbon: -1})
    model.add_reactions([uptake, biomass])
    model.objective = biomass
    model_file = tmp_path / "toy.xml"
    write_sbml_model(model, model_file)

    media_root = tmp_path / "generated"
    (media_root / "media").mkdir(parents=True)
    for name, flux in (("carbon", 10.0), ("no_carbon", 0.0)):
        pd.DataFrame({
            "compounds": ["cpd00027"], "name": ["D-Glucose"], "maxFlux": [flux]
        }).to_csv(media_root / "media" / f"{name}.csv", index=False)
    manifest = media_root / "manifest.tsv"
    pd.DataFrame({
        "medium_file": ["media/carbon.csv", "media/no_carbon.csv"]
    }).to_csv(manifest, sep="\t", index=False)
    outdir = tmp_path / "results"
    subprocess.run([
        sys.executable, str(SCRIPT), "--model", str(model_file),
        "--media-manifest", str(manifest), "--media-root", str(media_root),
        "--outdir", str(outdir), "--counterfactual-compounds", "cpd00027",
    ], check=True)
    growth = pd.read_csv(outdir / "growth_comparison.tsv", sep="\t")
    assert sorted(growth["biomass_flux"]) == [0.0, 10.0]
    assert (growth["strict_biomass_flux"] == growth["supplemented_biomass_flux"]).all()
    assert not growth["growth_rescued_by_supplements"].any()
    assert (growth["supplement_count"] == 0).all()
    assert growth["within_genome_percent_of_max"].max() == pytest.approx(100.0)
    assert growth["exact_maximum_medium"].sum() == 1
    assert (outdir / "genome_performance_summary.tsv").is_file()
    importance = pd.read_csv(outdir / "nutrient_importance_summary.tsv", sep="\t")
    carbon = importance.loc[importance["compound_id"].eq("cpd00027")].iloc[0]
    assert carbon["importance_score"] == pytest.approx(1.0)
    assert carbon["fraction_media_growth_lost_completely"] == pytest.approx(1.0)
    overview = pd.read_csv(outdir / "genome_performance_summary.tsv", sep="\t")
    assert overview.iloc[0]["nutrient_importance_cpd00027"] == pytest.approx(1.0)
    assert (outdir / "family_performance_summary.tsv").is_file()
    assert (outdir / "all_exchange_fva.tsv").is_file()
    assert len(pd.read_csv(outdir / "pairwise_growth_comparison.tsv", sep="\t")) == 1
    counterfactual = pd.read_csv(outdir / "counterfactual_growth.tsv", sep="\t")
    assert counterfactual.iloc[0]["growth_lost_completely"]
