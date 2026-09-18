from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pandas as pd


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "processes"
    / "missingness_sensitivity"
    / "run_missingness_sensitivity.py"
)
SPEC = importlib.util.spec_from_file_location("missingness_sensitivity", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_parse_thresholds_sorts_and_deduplicates():
    assert MODULE.parse_thresholds("0.25,0.10,0.25,0.20") == [0.10, 0.20, 0.25]
    assert MODULE.threshold_label(0.25) == "cutoff_0p25"
    try:
        MODULE.parse_thresholds("0.20,nan")
    except ValueError:
        pass
    else:
        raise AssertionError("A non-finite threshold must be rejected.")


def test_covariance_defaults_match_production_design():
    assert MODULE.DEFAULT_SELECTION_COVARIANCE_TYPE == "tied"
    assert MODULE.DEFAULT_FINAL_COVARIANCE_TYPE == "tied"


def test_cutoff_selection_maximizes_features_then_samples_and_ignores_pairwise_ari(
    tmp_path,
):
    summary = pd.DataFrame(
        {
            "missingness_cutoff": [0.20, 0.25, 0.30, 0.35],
            "status": ["complete"] * 4,
            "n_core_features": [8, 9, 9, 10],
            "n_rows_pca": [100, 98, 97, 80],
            "n_pcs_kept": [2, 2, 2, 2],
            "selected_stability_median_ari": [0.90, 0.88, 0.86, 0.90],
            "selected_min_cluster_frac": [0.05, 0.05, 0.05, 0.05],
        }
    )
    pairwise = pd.DataFrame(
        {
            "cutoff_left": [0.20, 0.25, 0.30],
            "cutoff_right": [0.25, 0.30, 0.35],
            # The 0.20-to-0.25 disagreement is deliberately complete. It must
            # remain diagnostic and cannot disqualify the otherwise feasible 0.25 run.
            "adjusted_rand_index": [0.0, 0.92, 0.95],
        }
    )

    tables_dir = tmp_path / "tables"
    tables_dir.mkdir()
    selected = MODULE.select_cutoff(
        summary,
        pairwise,
        tables_dir,
        min_cluster_frac=0.02,
        min_gmm_stability_ari=0.70,
    )

    assert selected == 0.35
    assert (tmp_path / "SELECTED_MISSINGNESS_CUTOFF.txt").read_text() == "0.35\n"
    decision = pd.read_csv(
        tables_dir / "missingness_cutoff_selection_decision.tsv", sep="\t"
    )
    assert decision.loc[decision["SELECTED"], "missingness_cutoff"].tolist() == [0.35]
    assert "passes_row_retention" not in decision
    rationale = json.loads(
        (tables_dir / "missingness_cutoff_selection_rationale.json").read_text()
    )
    assert rationale["selected_missingness_cutoff"] == 0.35
    assert rationale["sample_retention_use"] == "reported_and_used_only_as_a_tiebreaker"
    assert rationale["cross_cutoff_ari_use"] == "review_only_not_used_for_feasibility_or_selection"
