#!/usr/bin/env python3
"""Simulate one fixed gapseq SBML model under BASINS-generated media."""

from __future__ import annotations

import argparse
import json
import math
import re
import warnings
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from cobra.flux_analysis import flux_variability_analysis
from cobra.io import read_sbml_model


CPD_RE = re.compile(r"cpd\d{5}", re.IGNORECASE)
DEFAULT_COUNTERFACTUALS = ("all",)
BASIN_TRAINED_COMPOUNDS = {
    "cpd00007", "cpd00009", "cpd00013", "cpd00075",
    "cpd00209", "cpd00239", "cpd00659", "cpd01024",
}
COUNTERFACTUAL_NAMES = {
    "cpd00007": ("oxygen", "O2"),
    "cpd00009": ("phosphate", "Phosphate"),
    "cpd00013": ("ammonium", "Ammonium"),
    "cpd00075": ("nitrite", "Nitrite"),
    "cpd00209": ("nitrate", "Nitrate"),
    "cpd00239": ("hydrogen_sulfide", "H2S"),
    "cpd00268": ("thiosulfate", "Thiosulfate"),
    "cpd00659": ("nitrous_oxide", "Nitrous oxide"),
    "cpd01024": ("methane", "Methane"),
}


def args_parser() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", required=True, type=Path)
    parser.add_argument("--supplement-medium", type=Path)
    parser.add_argument("--media-manifest", required=True, type=Path)
    parser.add_argument("--media-root", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--solver", default="glpk")
    parser.add_argument("--fraction-of-optimum", type=float, default=1.0)
    parser.add_argument("--flux-threshold", type=float, default=1e-9)
    parser.add_argument("--supplement-max-flux", type=float, default=10.0)
    parser.add_argument("--run-fva", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--run-counterfactuals", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--counterfactual-compounds", default=",".join(DEFAULT_COUNTERFACTUALS))
    return parser.parse_args()


def modelseed_ids(value: object) -> set[str]:
    if value is None:
        return set()
    if isinstance(value, dict):
        values = list(value.keys()) + list(value.values())
        return set().union(*(modelseed_ids(item) for item in values)) if values else set()
    if isinstance(value, (list, tuple, set)):
        return set().union(*(modelseed_ids(item) for item in value)) if value else set()
    return {match.lower() for match in CPD_RE.findall(str(value))}


def exchange_info(reaction) -> dict:
    metabolites = list(reaction.metabolites)
    ids = set()
    for metabolite in metabolites:
        ids |= modelseed_ids(metabolite.id)
        ids |= modelseed_ids(metabolite.name)
        ids |= modelseed_ids(metabolite.annotation)
    ids |= modelseed_ids(reaction.id)
    ids |= modelseed_ids(reaction.annotation)
    if len(metabolites) != 1:
        return {"valid": False, "coefficient": math.nan, "compound_ids": sorted(ids)}
    coefficient = float(reaction.metabolites[metabolites[0]])
    return {
        "valid": coefficient != 0,
        "coefficient": coefficient,
        "compound_ids": sorted(ids),
        "metabolite": metabolites[0],
    }


def exchange_index(model) -> tuple[dict[str, list], dict[str, dict]]:
    index: dict[str, list] = {}
    details = {}
    for reaction in model.exchanges:
        info = exchange_info(reaction)
        details[reaction.id] = info
        for compound in info["compound_ids"]:
            index.setdefault(compound, []).append(reaction)
    return index, details


def validate_medium(path: Path) -> pd.DataFrame:
    medium = pd.read_csv(path)
    required = ["compounds", "name", "maxFlux"]
    missing = [col for col in required if col not in medium]
    if missing:
        raise ValueError(f"{path}: missing columns {missing}")
    medium = medium[required].copy()
    medium["compounds"] = medium["compounds"].astype(str).str.lower()
    invalid = ~medium["compounds"].str.fullmatch(r"cpd\d{5}")
    if invalid.any():
        raise ValueError(f"{path}: invalid ModelSEED ids {medium.loc[invalid, 'compounds'].tolist()}")
    medium["maxFlux"] = pd.to_numeric(medium["maxFlux"], errors="raise")
    if (medium["maxFlux"] < 0).any():
        raise ValueError(f"{path}: negative maxFlux values are not allowed")
    if medium["compounds"].duplicated().any():
        duplicates = medium.loc[medium["compounds"].duplicated(False), "compounds"].tolist()
        raise ValueError(f"{path}: duplicate compounds {duplicates}")
    return medium


def close_all_uptake(model, details: dict[str, dict]) -> list[str]:
    ambiguous = []
    for reaction in model.exchanges:
        info = details[reaction.id]
        if not info["valid"]:
            ambiguous.append(reaction.id)
            continue
        if info["coefficient"] < 0:
            reaction.lower_bound = max(0.0, reaction.lower_bound)
        else:
            reaction.upper_bound = min(0.0, reaction.upper_bound)
    return ambiguous


def apply_medium(model, medium: pd.DataFrame, index, details, medium_name: str) -> tuple[pd.DataFrame, set[str]]:
    rows, supplied = [], set()
    for item in medium.itertuples(index=False):
        matches = index.get(item.compounds, [])
        if not matches:
            rows.append({
                "medium": medium_name, "compound_id": item.compounds, "compound_name": item.name,
                "requested_maxFlux": item.maxFlux, "exchange_reaction_id": "",
                "mapping_status": "no_exchange_found", "exchange_orientation": "",
                "lower_bound_after": "", "upper_bound_after": "",
            })
            continue
        status = "mapped" if len(matches) == 1 else "multiple_exchanges"
        for reaction in matches:
            info = details[reaction.id]
            if not info["valid"]:
                mapping_status, orientation = "invalid_exchange", "ambiguous"
            elif info["coefficient"] < 0:
                reaction.lower_bound = -float(item.maxFlux)
                mapping_status, orientation = status, "negative_flux_uptake"
            else:
                reaction.upper_bound = float(item.maxFlux)
                mapping_status, orientation = status, "positive_flux_uptake"
            if item.maxFlux > 0 and mapping_status in {"mapped", "multiple_exchanges"}:
                supplied.add(item.compounds)
            rows.append({
                "medium": medium_name, "compound_id": item.compounds, "compound_name": item.name,
                "requested_maxFlux": item.maxFlux, "exchange_reaction_id": reaction.id,
                "mapping_status": mapping_status, "exchange_orientation": orientation,
                "lower_bound_after": reaction.lower_bound, "upper_bound_after": reaction.upper_bound,
            })
    return pd.DataFrame(rows), supplied


def finite(value) -> float:
    try:
        result = float(value)
        return result if np.isfinite(result) else math.nan
    except (TypeError, ValueError):
        return math.nan


def reaction_table(model, solution, medium: str, threshold: float) -> pd.DataFrame:
    fluxes = solution.fluxes if solution.status == "optimal" else pd.Series(dtype=float)
    rows = []
    exchange_ids = {reaction.id for reaction in model.exchanges}
    for reaction in model.reactions:
        flux = finite(fluxes.get(reaction.id, math.nan))
        rows.append({
            "medium": medium, "reaction_id": reaction.id, "reaction_name": reaction.name,
            "subsystem": reaction.subsystem or "", "equation": reaction.reaction,
            "flux": flux, "absolute_flux": abs(flux) if np.isfinite(flux) else math.nan,
            "active": bool(np.isfinite(flux) and abs(flux) > threshold),
            "is_exchange": reaction.id in exchange_ids,
            "lower_bound": reaction.lower_bound, "upper_bound": reaction.upper_bound,
            "gene_reaction_rule": reaction.gene_reaction_rule,
            "gene_ids": ",".join(sorted(gene.id for gene in reaction.genes)),
        })
    return pd.DataFrame(rows)


def exchange_table(model, solution, medium: str, supplied: set[str], details, threshold: float) -> pd.DataFrame:
    fluxes = solution.fluxes if solution.status == "optimal" else pd.Series(dtype=float)
    rows = []
    for reaction in model.exchanges:
        info = details[reaction.id]
        raw = finite(fluxes.get(reaction.id, math.nan))
        environmental = info["coefficient"] * raw if info["valid"] and np.isfinite(raw) else math.nan
        if not info["valid"] or not np.isfinite(environmental):
            direction = "ambiguous"
        elif environmental > threshold:
            direction = "uptake"
        elif environmental < -threshold:
            direction = "secretion"
        else:
            direction = "inactive"
        compounds = info["compound_ids"]
        rows.append({
            "medium": medium, "reaction_id": reaction.id, "reaction_name": reaction.name,
            "compound_ids": ",".join(compounds),
            "compound_names": getattr(info.get("metabolite"), "name", ""),
            "raw_reaction_flux": raw, "environmental_flux": environmental,
            "direction": direction, "lower_bound": reaction.lower_bound,
            "upper_bound": reaction.upper_bound,
            "supplied_in_medium": any(compound in supplied for compound in compounds),
            "mapping_status": "mapped" if info["valid"] and compounds else "ambiguous",
        })
    return pd.DataFrame(rows)


def run_one(
    base_model,
    medium_path: Path,
    medium_name: str,
    args,
    per_medium: Path,
    background_supplements: pd.DataFrame,
    reference_growth: float,
) -> dict:
    medium = validate_medium(medium_path)

    # Strict scenario: only the BASINS environmental medium is available.
    strict_model = base_model.copy()
    strict_model.solver = args.solver
    strict_index, strict_details = exchange_index(strict_model)
    close_all_uptake(strict_model, strict_details)
    apply_medium(strict_model, medium, strict_index, strict_details, medium_name)
    strict_solution = strict_model.optimize()
    strict_biomass = (
        finite(strict_solution.objective_value)
        if strict_solution.status == "optimal" else math.nan
    )

    # Supplemented scenario: rescue nutrients are applied first and BASINS
    # environmental bounds take precedence for represented compounds.
    model = base_model.copy()
    model.solver = args.solver
    index, details = exchange_index(model)
    ambiguous = close_all_uptake(model, details)
    background_mapping = pd.DataFrame()
    background_supplied: set[str] = set()
    if not background_supplements.empty:
        background_mapping, background_supplied = apply_medium(
            model, background_supplements, index, details, medium_name
        )
        background_mapping.insert(1, "medium_role", "genome_growth_supplement")
    mapping, supplied = apply_medium(model, medium, index, details, medium_name)
    mapping.insert(1, "medium_role", "basin_environment")
    if not background_mapping.empty:
        mapping = pd.concat([background_mapping, mapping], ignore_index=True)
    supplied |= background_supplied
    solution = model.optimize()
    reactions = reaction_table(model, solution, medium_name, args.flux_threshold)
    exchanges = exchange_table(model, solution, medium_name, supplied, details, args.flux_threshold)
    mapping.to_csv(per_medium / f"{medium_name}.medium_mapping.tsv", sep="\t", index=False)
    reactions.to_csv(per_medium / f"{medium_name}.reaction_fluxes.tsv", sep="\t", index=False)
    exchanges.to_csv(per_medium / f"{medium_name}.exchange_fluxes.tsv", sep="\t", index=False)
    fva = pd.DataFrame()
    if args.run_fva and solution.status == "optimal":
        try:
            raw_fva = flux_variability_analysis(
                model, reaction_list=list(model.exchanges),
                fraction_of_optimum=args.fraction_of_optimum,
            ).reset_index(names="reaction_id")
            raw_fva.insert(0, "medium", medium_name)
            raw_fva["zero_feasible"] = (raw_fva.minimum <= 0) & (raw_fva.maximum >= 0)
            raw_fva["obligate_forward"] = raw_fva.minimum > args.flux_threshold
            raw_fva["obligate_reverse"] = raw_fva.maximum < -args.flux_threshold
            raw_fva["variable_flux"] = raw_fva.maximum - raw_fva.minimum
            fva = raw_fva
            fva.to_csv(per_medium / f"{medium_name}.exchange_fva.tsv", sep="\t", index=False)
        except Exception as exc:
            warnings.warn(f"FVA failed for {medium_name}: {exc}")
    biomass = finite(solution.objective_value) if solution.status == "optimal" else math.nan
    strict_feasible = bool(
        strict_solution.status == "optimal"
        and np.isfinite(strict_biomass)
        and strict_biomass > args.flux_threshold
    )
    supplemented_feasible = bool(
        solution.status == "optimal"
        and np.isfinite(biomass)
        and biomass > args.flux_threshold
    )
    rescue_effect = (
        biomass - strict_biomass
        if np.isfinite(biomass) and np.isfinite(strict_biomass) else math.nan
    )
    supplement_ids = sorted(
        background_supplements.loc[
            background_supplements["maxFlux"] > 0, "compounds"
        ].astype(str)
    )
    active = reactions.active.sum()
    active_ex = exchanges.direction.isin(["uptake", "secretion"]).sum()
    summary = {
        "model_id": model.id, "medium": medium_name,
        "strict_solver_status": strict_solution.status,
        "strict_growth_feasible": strict_feasible,
        "strict_biomass_flux": strict_biomass,
        "supplemented_solver_status": solution.status,
        "supplemented_growth_feasible": supplemented_feasible,
        "supplemented_biomass_flux": biomass,
        "growth_rescued_by_supplements": bool(not strict_feasible and supplemented_feasible),
        "absolute_rescue_effect": rescue_effect,
        "relative_rescue_effect": (
            rescue_effect / abs(strict_biomass)
            if np.isfinite(rescue_effect) and np.isfinite(strict_biomass)
            and abs(strict_biomass) > args.flux_threshold else math.nan
        ),
        "supplement_count": len(supplement_ids),
        "supplement_compounds": ",".join(supplement_ids),
        # Backward-compatible aliases describe the supplemented scenario.
        "solver_status": solution.status,
        "growth_feasible": supplemented_feasible,
        "biomass_flux": biomass,
        "active_reactions": int(active),
        "gapseq_default_biomass_flux": reference_growth,
        "absolute_growth_change_vs_gapseq_default": (
            biomass - reference_growth
            if np.isfinite(biomass) and np.isfinite(reference_growth) else math.nan
        ),
        "growth_fold_change_vs_gapseq_default": (
            biomass / reference_growth
            if np.isfinite(biomass) and np.isfinite(reference_growth) and reference_growth > args.flux_threshold
            else math.nan
        ),
        "active_internal_reactions": int((reactions.active & ~reactions.is_exchange).sum()),
        "active_exchanges": int(active_ex),
        "total_uptake_flux": exchanges.loc[exchanges.direction == "uptake", "environmental_flux"].sum(),
        "total_secretion_flux": -exchanges.loc[exchanges.direction == "secretion", "environmental_flux"].sum(),
        "supplied_compounds": len(supplied),
        "environmental_compounds_supplied": int((medium.maxFlux > 0).sum()),
        "background_supplements_supplied": int((background_supplements.maxFlux > 0).sum()),
        "mapped_compounds": int(mapping.mapping_status.isin(["mapped", "multiple_exchanges"]).groupby(mapping.compound_id).any().sum()),
        "unmapped_compounds": int((mapping.groupby("compound_id").mapping_status.apply(lambda x: (x == "no_exchange_found").all())).sum()),
        "ambiguous_exchanges": len(ambiguous),
    }
    pd.DataFrame([summary]).to_csv(per_medium / f"{medium_name}.summary.tsv", sep="\t", index=False)
    environmental_supplied = supplied - background_supplied
    return {
        "summary": summary, "reactions": reactions, "exchanges": exchanges, "fva": fva,
        "model": model, "medium": medium, "supplied": supplied, "details": details,
        "environmental_supplied": environmental_supplied,
        "background_supplied": background_supplied,
        "background_names": (
            background_supplements.set_index("compounds")["name"].to_dict()
            if not background_supplements.empty else {}
        ),
    }


def pairwise(results: dict[str, dict], outdir: Path, threshold: float) -> None:
    growth, reaction_rows, exchange_rows = [], [], []
    for name_a, name_b in combinations(results, 2):
        a, b = results[name_a], results[name_b]
        ga, gb = a["summary"]["biomass_flux"], b["summary"]["biomass_flux"]
        strict_a = a["summary"]["strict_biomass_flux"]
        strict_b = b["summary"]["strict_biomass_flux"]
        growth.append({
            "medium_A": name_a, "medium_B": name_b, "biomass_A": ga, "biomass_B": gb,
            "absolute_difference": gb - ga if np.isfinite(ga) and np.isfinite(gb) else math.nan,
            "relative_difference": (gb - ga) / abs(ga) if np.isfinite(ga) and ga != 0 and np.isfinite(gb) else math.nan,
            "fold_change": gb / ga if np.isfinite(ga) and ga != 0 and np.isfinite(gb) else math.nan,
            "strict_biomass_A": strict_a, "strict_biomass_B": strict_b,
            "strict_absolute_difference": (
                strict_b - strict_a
                if np.isfinite(strict_a) and np.isfinite(strict_b) else math.nan
            ),
            "supplemented_biomass_A": ga, "supplemented_biomass_B": gb,
            "supplemented_absolute_difference": (
                gb - ga if np.isfinite(ga) and np.isfinite(gb) else math.nan
            ),
        })
        for key, target in (("reactions", reaction_rows), ("exchanges", exchange_rows)):
            id_col = "reaction_id"
            flux_col = "flux" if key == "reactions" else "environmental_flux"
            merged = a[key][[id_col, flux_col]].merge(
                b[key][[id_col, flux_col]], on=id_col, how="outer", suffixes=("_A", "_B")
            ).fillna(0)
            for row in merged.itertuples(index=False):
                fa, fb = getattr(row, f"{flux_col}_A"), getattr(row, f"{flux_col}_B")
                active_a, active_b = abs(fa) > threshold, abs(fb) > threshold
                state = "active_both" if active_a and active_b else "A_only" if active_a else "B_only" if active_b else "inactive_both"
                target.append({
                    "medium_A": name_a, "medium_B": name_b, "reaction_id": row.reaction_id,
                    "flux_A": fa, "flux_B": fb, "delta_flux": fb - fa,
                    "absolute_delta_flux": abs(fb - fa), "active_A": active_a,
                    "active_B": active_b, "activity_state": state,
                })
    pd.DataFrame(growth).to_csv(outdir / "pairwise_growth_comparison.tsv", sep="\t", index=False)
    pd.DataFrame(reaction_rows).sort_values("absolute_delta_flux", ascending=False).to_csv(
        outdir / "pairwise_reaction_comparison.tsv", sep="\t", index=False
    )
    pd.DataFrame(exchange_rows).sort_values("absolute_delta_flux", ascending=False).to_csv(
        outdir / "pairwise_exchange_comparison.tsv", sep="\t", index=False
    )


def counterfactuals(
    results: dict[str, dict], compounds: set[str], args, outdir: Path
) -> pd.DataFrame:
    rows = []
    for name, result in results.items():
        for compound in sorted(compounds & result["supplied"]):
            environmental = compound in result["environmental_supplied"]
            if environmental and compound in BASIN_TRAINED_COMPOUNDS:
                nutrient_source = "basin_observation_trained"
            elif environmental:
                nutrient_source = "template_baseline"
            else:
                nutrient_source = "gapseq_auxotrophy_supplement"
            compound_names = result["medium"].set_index("compounds")["name"].to_dict()
            compound_names.update(result["background_names"])
            model = result["model"].copy()
            index, details = exchange_index(model)
            for reaction in index.get(compound, []):
                info = details[reaction.id]
                if info["coefficient"] < 0:
                    reaction.lower_bound = max(0.0, reaction.lower_bound)
                else:
                    reaction.upper_bound = min(0.0, reaction.upper_bound)
            solution = model.optimize()
            observed = result["summary"]["biomass_flux"]
            perturbed = finite(solution.objective_value) if solution.status == "optimal" else math.nan
            rows.append({
                "medium": name, "removed_compound": compound,
                "compound_name": compound_names.get(compound, compound),
                "nutrient_source": nutrient_source,
                "observed_biomass_flux": observed, "perturbed_biomass_flux": perturbed,
                "absolute_biomass_effect": observed - perturbed if np.isfinite(observed) and np.isfinite(perturbed) else math.nan,
                "relative_biomass_effect": (observed - perturbed) / abs(observed) if np.isfinite(observed) and observed != 0 and np.isfinite(perturbed) else math.nan,
                "growth_lost_completely": bool(np.isfinite(observed) and observed > args.flux_threshold and (not np.isfinite(perturbed) or perturbed <= args.flux_threshold)),
                "solver_status": solution.status,
            })
    columns = [
        "medium", "removed_compound", "compound_name", "nutrient_source",
        "observed_biomass_flux",
        "perturbed_biomass_flux", "absolute_biomass_effect",
        "relative_biomass_effect", "growth_lost_completely", "solver_status",
    ]
    frame = pd.DataFrame(rows, columns=columns)
    frame.to_csv(outdir / "counterfactual_growth.tsv", sep="\t", index=False)
    return frame


def summarize_nutrient_importance(
    counterfactual: pd.DataFrame, compounds: set[str], outdir: Path
) -> dict[str, float]:
    """Summarize proportional growth loss after removing each supplied nutrient."""
    rows = []
    wide = {}
    for compound in sorted(compounds):
        slug, name = COUNTERFACTUAL_NAMES.get(compound, (compound, compound))
        subset = counterfactual.loc[
            counterfactual["removed_compound"].eq(compound)
        ] if not counterfactual.empty else pd.DataFrame()
        if not subset.empty:
            name = str(subset["compound_name"].iloc[0])
            nutrient_source = str(subset["nutrient_source"].iloc[0])
        else:
            nutrient_source = ""
        effects = (
            pd.to_numeric(subset["relative_biomass_effect"], errors="coerce")
            .replace([np.inf, -np.inf], np.nan)
            .dropna()
            .clip(lower=0.0, upper=1.0)
            if not subset.empty else pd.Series(dtype=float)
        )
        score = finite(effects.median()) if not effects.empty else math.nan
        wide[f"nutrient_importance_{slug}"] = score
        rows.append({
            "compound_id": compound,
            "compound_name": name,
            "nutrient_source": nutrient_source,
            "media_tested": len(subset),
            "media_with_numeric_effect": len(effects),
            "importance_score": score,
            "mean_relative_growth_loss": (
                finite(effects.mean()) if not effects.empty else math.nan
            ),
            "maximum_relative_growth_loss": (
                finite(effects.max()) if not effects.empty else math.nan
            ),
            "fraction_media_growth_lost_completely": (
                float(subset["growth_lost_completely"].astype(bool).mean())
                if not subset.empty else math.nan
            ),
        })
    pd.DataFrame(rows).to_csv(
        outdir / "nutrient_importance_summary.tsv", sep="\t", index=False
    )
    return wide


def write_performance_summaries(
    growth: pd.DataFrame,
    outdir: Path,
    nutrient_importance: dict[str, float] | None = None,
) -> pd.DataFrame:
    """Add rankings and write one-genome overview and family summaries."""
    growth = growth.copy()
    for column, default in (
        ("scope", "unspecified"),
        ("cruise_group", ""),
        ("family", "unspecified"),
        ("compartment", ""),
        ("medium_file", ""),
    ):
        if column not in growth:
            growth[column] = default
    maximum = finite(growth["supplemented_biomass_flux"].max())
    minimum = finite(growth["supplemented_biomass_flux"].min())
    growth["within_genome_fraction_of_max"] = (
        growth["supplemented_biomass_flux"] / maximum
        if np.isfinite(maximum) and maximum > 0 else math.nan
    )
    growth["within_genome_percent_of_max"] = 100 * growth["within_genome_fraction_of_max"]
    growth["supplemented_growth_rank"] = growth["supplemented_biomass_flux"].rank(
        method="dense", ascending=False
    ).astype("Int64")
    growth["exact_maximum_medium"] = np.isclose(
        growth["supplemented_biomass_flux"], maximum, rtol=1e-9, atol=1e-12
    )
    growth["within_95_percent_of_maximum"] = growth["within_genome_fraction_of_max"] >= 0.95

    exact = growth.loc[growth["exact_maximum_medium"], "medium"].astype(str).tolist()
    top95 = growth.loc[growth["within_95_percent_of_maximum"], "medium"].astype(str).tolist()
    supplement_count = int(growth["supplement_count"].iloc[0])
    overview_row = {
        "model_id": growth["model_id"].iloc[0],
        "media_tested": len(growth),
        "strict_media_supporting_growth": int(growth["strict_growth_feasible"].sum()),
        "strict_growth_fraction": float(growth["strict_growth_feasible"].mean()),
        "supplemented_media_supporting_growth": int(growth["supplemented_growth_feasible"].sum()),
        "supplemented_growth_fraction": float(growth["supplemented_growth_feasible"].mean()),
        "media_rescued_by_supplements": int(growth["growth_rescued_by_supplements"].sum()),
        "supplement_count": supplement_count,
        "supplement_compounds": growth["supplement_compounds"].iloc[0],
        "supplemented_growth_min": minimum,
        "supplemented_growth_median": finite(growth["supplemented_biomass_flux"].median()),
        "supplemented_growth_max": maximum,
        "supplemented_max_to_min_fold": (
            maximum / minimum if np.isfinite(minimum) and minimum > 0 else math.nan
        ),
        "distinct_supplemented_growth_levels": int(
            growth["supplemented_biomass_flux"].round(9).nunique()
        ),
        "exact_maximum_media_n": len(exact),
        "exact_maximum_media": ";".join(exact),
        "top_95_percent_media_n": len(top95),
        "top_95_percent_media": ";".join(top95),
        "interpretation": (
            "supplement-dependent reduced-genome response"
            if supplement_count > 0 else "self-sufficient model response"
        ),
    }
    overview_row.update(nutrient_importance or {})
    overview = pd.DataFrame([overview_row])
    overview.to_csv(outdir / "genome_performance_summary.tsv", sep="\t", index=False)

    family_rows = []
    group_columns = ["scope", "cruise_group", "family"]
    for keys, frame in growth.groupby(group_columns, dropna=False, sort=True):
        family_max = finite(frame["supplemented_biomass_flux"].max())
        family_min = finite(frame["supplemented_biomass_flux"].min())
        exact_mask = np.isclose(
            frame["supplemented_biomass_flux"], family_max, rtol=1e-9, atol=1e-12
        )
        top_mask = frame["supplemented_biomass_flux"] >= 0.95 * family_max
        display = frame["compartment"].fillna(frame["cruise_group"]).fillna(frame["medium"])
        family_rows.append({
            "model_id": frame["model_id"].iloc[0],
            "scope": keys[0],
            "cruise_group": "" if pd.isna(keys[1]) else keys[1],
            "family": keys[2],
            "media_compared": len(frame),
            "family_growth_min": family_min,
            "family_growth_median": finite(frame["supplemented_biomass_flux"].median()),
            "family_growth_max": family_max,
            "family_max_to_min_fold": (
                family_max / family_min
                if np.isfinite(family_min) and family_min > 0 else math.nan
            ),
            "exact_winner_n": int(exact_mask.sum()),
            "exact_winners": ";".join(display.loc[exact_mask].astype(str)),
            "top_95_percent_n": int(top_mask.sum()),
            "top_95_percent_media": ";".join(display.loc[top_mask].astype(str)),
        })
    pd.DataFrame(family_rows).to_csv(
        outdir / "family_performance_summary.tsv", sep="\t", index=False
    )
    return growth


def main() -> None:
    args = args_parser()
    if not 0 < args.fraction_of_optimum <= 1:
        raise ValueError("--fraction-of-optimum must be in (0, 1]")
    if args.supplement_max_flux <= 0:
        raise ValueError("--supplement-max-flux must be positive")
    args.outdir.mkdir(parents=True, exist_ok=True)
    per_medium = args.outdir / "per_medium"
    per_medium.mkdir(exist_ok=True)
    captured = []
    with warnings.catch_warnings(record=True) as warning_records:
        warnings.simplefilter("always")
        model = read_sbml_model(str(args.model))
        model.solver = args.solver
        initial = model.optimize()
        reference_growth = finite(initial.objective_value) if initial.status == "optimal" else math.nan
        objective = str(model.objective.expression)
        pd.DataFrame([{
            "model_id": model.id, "model_file": str(args.model.resolve()), "solver": args.solver,
            "reactions": len(model.reactions), "metabolites": len(model.metabolites),
            "genes": len(model.genes), "exchanges": len(model.exchanges),
            "objective": objective, "initial_status": initial.status,
            "initial_objective_value": finite(initial.objective_value),
        }]).to_csv(args.outdir / "model_summary.tsv", sep="\t", index=False)
        manifest = pd.read_csv(args.media_manifest, sep="\t")
        if "medium_file" not in manifest:
            raise ValueError("Media manifest lacks medium_file")
        prepared_media = []
        basin_controlled_compounds: set[str] = set()
        for row in manifest.itertuples(index=False):
            medium_path = args.media_root / row.medium_file
            medium_name = re.sub(r"[^A-Za-z0-9._-]+", "_", Path(row.medium_file).with_suffix("").as_posix())
            medium = validate_medium(medium_path)
            basin_controlled_compounds |= set(medium["compounds"])
            metadata = {
                key: value for key, value in row._asdict().items()
                if key in {"scope", "cruise_group", "family", "compartment",
                           "medium_file", "template", "effective_n", "recipe_id"}
            }
            prepared_media.append((medium_path, medium_name, metadata))

        predicted_supplements = pd.DataFrame(columns=["compounds", "name", "maxFlux"])
        if args.supplement_medium:
            predicted_supplements = validate_medium(args.supplement_medium)
        supplement_audit = predicted_supplements.rename(
            columns={"maxFlux": "gapseq_predicted_maxFlux"}
        ).copy()
        supplement_audit["represented_in_basin_media"] = supplement_audit["compounds"].isin(
            basin_controlled_compounds
        )
        supplement_audit["applied_as_background_supplement"] = (
            ~supplement_audit["represented_in_basin_media"]
            & (supplement_audit["gapseq_predicted_maxFlux"] > 0)
        )
        supplement_audit["applied_supplement_maxFlux"] = np.where(
            supplement_audit["applied_as_background_supplement"],
            args.supplement_max_flux,
            0.0,
        )
        supplement_audit["policy_reason"] = np.where(
            supplement_audit["represented_in_basin_media"],
            "BASINS media control this compound",
            np.where(
                supplement_audit["gapseq_predicted_maxFlux"] > 0,
                "Absent from all BASINS media; supplied at non-limiting BASINS background bound",
                "Zero-flux prediction; not supplied",
            ),
        )
        supplement_audit.to_csv(args.outdir / "growth_supplement_audit.tsv", sep="\t", index=False)
        background_supplements = supplement_audit.loc[
            supplement_audit.get(
                "applied_as_background_supplement",
                pd.Series(False, index=supplement_audit.index),
            ),
            ["compounds", "name", "applied_supplement_maxFlux"],
        ].copy()
        background_supplements = background_supplements.rename(
            columns={"applied_supplement_maxFlux": "maxFlux"}
        )

        results = {}
        qc = []
        for medium_path, medium_name, metadata in prepared_media:
            try:
                results[medium_name] = run_one(
                    model,
                    medium_path,
                    medium_name,
                    args,
                    per_medium,
                    background_supplements,
                    reference_growth,
                )
                results[medium_name]["summary"].update(metadata)
                qc.append({"medium": medium_name, "status": "completed", "message": ""})
            except Exception as exc:
                qc.append({"medium": medium_name, "status": "failed", "message": str(exc)})
        if not results:
            raise RuntimeError("All medium simulations failed")
        requested = {
            item.strip().lower()
            for item in args.counterfactual_compounds.split(",")
            if item.strip()
        }
        if "all" in requested:
            requested = set().union(
                *(result["supplied"] for result in results.values())
            )
        counterfactual = pd.DataFrame()
        nutrient_importance = {}
        if args.run_counterfactuals:
            counterfactual = counterfactuals(results, requested, args, args.outdir)
            nutrient_importance = summarize_nutrient_importance(
                counterfactual, requested, args.outdir
            )
        growth = pd.DataFrame([item["summary"] for item in results.values()])
        growth = write_performance_summaries(
            growth, args.outdir, nutrient_importance
        )
        growth.to_csv(args.outdir / "growth_comparison.tsv", sep="\t", index=False)
        pd.concat([item["reactions"] for item in results.values()], ignore_index=True).to_csv(
            args.outdir / "all_reaction_fluxes.tsv", sep="\t", index=False
        )
        pd.concat([item["exchanges"] for item in results.values()], ignore_index=True).to_csv(
            args.outdir / "all_exchange_fluxes.tsv", sep="\t", index=False
        )
        fva_frames = [item["fva"] for item in results.values() if not item["fva"].empty]
        if fva_frames:
            pd.concat(fva_frames, ignore_index=True).to_csv(args.outdir / "all_exchange_fva.tsv", sep="\t", index=False)
        pairwise(results, args.outdir, args.flux_threshold)
        pd.DataFrame(qc).to_csv(args.outdir / "run_qc.tsv", sep="\t", index=False)
        captured = [str(item.message) for item in warning_records]
    (args.outdir / "warnings.log").write_text("\n".join(captured) + ("\n" if captured else ""))
    (args.outdir / "run_metadata.json").write_text(json.dumps({
        "model": str(args.model.resolve()), "solver": args.solver,
        "fraction_of_optimum": args.fraction_of_optimum,
        "flux_threshold": args.flux_threshold, "media_completed": len(results),
        "background_supplements": len(background_supplements),
        "supplement_max_flux": args.supplement_max_flux,
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
