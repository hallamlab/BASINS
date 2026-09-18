from __future__ import annotations

import importlib.util
from argparse import Namespace
from pathlib import Path

import pandas as pd


SCRIPT = Path(__file__).parents[1] / "processes/eof_state_cluster/eof_state_clustering.py"
SPEC = importlib.util.spec_from_file_location("eof_state_clustering", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def test_select_k_requires_fraction_count_and_stability() -> None:
    metrics = pd.DataFrame([
        {"K": 2, "AIC": 20, "BIC": 22, "ICL": 23, "mean_resp_entropy": .1,
         "min_cluster_n": 4, "min_cluster_frac": .20, "CV_loglik_mean": -2,
         "CV_loglik_std": .1, "stability_median_ARI": .9,
         "stability_mean_ARI": .9, "stability_n_reps": 50},
        {"K": 3, "AIC": 18, "BIC": 20, "ICL": 21, "mean_resp_entropy": .2,
         "min_cluster_n": 8, "min_cluster_frac": .12, "CV_loglik_mean": -2.1,
         "CV_loglik_std": .1, "stability_median_ARI": .7,
         "stability_mean_ARI": .7, "stability_n_reps": 50},
    ])
    args = Namespace(min_cluster_frac=.10, min_cluster_n=5, stability_min_ari=.60,
                     allow_no_feasible_fallback=False, select_by="icl", select_delta=5,
                     cv_folds=5)
    chosen, decision = MODULE.select_k(metrics, args)
    assert chosen == 3
    assert not bool(decision.loc[decision.K.eq(2), "passes_min_cluster"].iloc[0])
    assert bool(decision.loc[decision.K.eq(3), "SELECTED"].iloc[0])
