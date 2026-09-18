from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    ROOT
    / "processes"
    / "compare_compartments"
    / "env_compare_compartments.py"
)


def load_module():
    spec = importlib.util.spec_from_file_location("env_compare_compartments", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_redundancy_statistics_identify_equivalent_labels():
    module = load_module()
    frame = pd.DataFrame(
        {
            "first": ["a", "a", "b", "b", "c", "c"],
            "equivalent": ["x", "x", "y", "y", "z", "z"],
            "other": ["u", "v", "u", "v", "u", "v"],
        }
    )

    result = module.grouping_redundancy_table(
        frame,
        {"First": "first", "Equivalent": "equivalent", "Other": "other"},
    )
    equivalent = result[
        (result["grouping_a"] == "First")
        & (result["grouping_b"] == "Equivalent")
    ].iloc[0]

    assert equivalent["normalized_mutual_information"] == 1.0
    assert equivalent["cramers_v_bias_corrected"] == 1.0


def test_sparse_benchmark_uses_matched_leave_cruise_out_rows():
    module = load_module()
    rows = []
    for cruise in range(8):
        for sample in range(8):
            gmm = "G1" if sample < 4 else "G2"
            depth = float(10 + sample * 10)
            rows.append(
                {
                    "Cruise": f"C{cruise}",
                    "Depth_anchored": depth,
                    "Season": "DJF" if cruise % 2 == 0 else "JJA",
                    "o2_compartment": "oxic" if sample < 5 else "low_o2",
                    "component": gmm,
                    "hybrid_compartment": f"{gmm}_{'upper' if sample < 4 else 'lower'}",
                    "Sparse feature": (
                        1.0
                        + 0.01 * depth
                        + (4.0 if gmm == "G2" else 0.0)
                        + 0.02 * cruise
                    ),
                }
            )
    frame = pd.DataFrame(rows)

    folds, audit = module.heldout_sparse_cv(
        frame=frame,
        feature="Sparse feature",
        depth_col="Depth_anchored",
        season_col="Season",
        cruise_col="Cruise",
        grouping_columns={
            "Legacy O2": "o2_compartment",
            "GMM": "component",
            "Hybrid": "hybrid_compartment",
        },
        minimum_group_n=2,
        minimum_cruises=2,
    )
    summary = module.summarize_sparse_cv(folds)

    assert audit.loc[0, "n_observations"] == len(frame)
    assert audit.loc[0, "n_cruises"] == 8
    assert folds["fold_cruise"].nunique() == 8
    assert set(folds["model"]) == {
        "Depth + season baseline",
        "Legacy O2",
        "GMM",
        "Hybrid",
    }
    assert (folds.groupby("model")["n_test"].sum() == len(frame)).all()
    gmm_r2 = summary.loc[summary["model"] == "GMM", "pooled_cv_r2"].iloc[0]
    baseline_r2 = summary.loc[
        summary["model"] == "Depth + season baseline", "pooled_cv_r2"
    ].iloc[0]
    assert np.isfinite(gmm_r2)
    assert gmm_r2 > baseline_r2


def test_confidence_weighted_silhouette_applies_weights():
    module = load_module()
    x = np.array([[0.0], [0.1], [0.2], [1.0], [2.0], [2.1], [2.2]])
    labels = pd.Series(["a", "a", "a", "a", "b", "b", "b"])
    weights = pd.Series([1.0, 1.0, 1.0, 0.01, 1.0, 1.0, 1.0])

    from sklearn.metrics import silhouette_samples

    samples = silhouette_samples(x, labels)
    expected = np.average(samples, weights=weights)
    observed = module.responsibility_weighted_silhouette(x, labels, weights)

    assert np.isclose(observed, expected)
    assert not np.isclose(observed, samples.mean())


def test_within_block_permutation_actually_changes_labels():
    module = load_module()
    x = np.concatenate(
        [np.linspace(-2.0, -1.0, 20), np.linspace(1.0, 2.0, 20)]
    )[:, None]
    labels = pd.Series(["a"] * 20 + ["b"] * 20)
    blocks = pd.Series(np.tile(np.arange(10), 4).astype(str))

    result = module.within_block_permutation_test(
        x,
        labels,
        blocks,
        n_perm=199,
        random_state=42,
    )

    assert result["observed"] > 0.5
    assert result["p_value"] < 0.05


def test_sparse_benchmark_reports_insufficient_cruise_support():
    module = load_module()
    frame = pd.DataFrame(
        {
            "Cruise": ["C1", "C1", "C2", "C2"],
            "Depth_anchored": [10.0, 20.0, 10.0, 20.0],
            "Season": ["DJF"] * 4,
            "component": ["G1", "G1", "G2", "G2"],
            "Sparse feature": [0.0, 1.0, 2.0, 3.0],
        }
    )

    folds, audit = module.heldout_sparse_cv(
        frame=frame,
        feature="Sparse feature",
        depth_col="Depth_anchored",
        season_col="Season",
        cruise_col="Cruise",
        grouping_columns={"GMM": "component"},
        minimum_group_n=1,
        minimum_cruises=3,
    )

    assert folds.empty
    assert audit.loc[0, "status"] == "insufficient_cruises"
