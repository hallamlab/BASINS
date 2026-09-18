from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    ROOT
    / "processes"
    / "eof_state_cluster"
    / "cruise_group_season_benchmark.py"
)


def load_module():
    spec = importlib.util.spec_from_file_location(
        "cruise_group_season_benchmark", SCRIPT
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_season_group_redundancy_detects_equivalent_partitions():
    module = load_module()
    frame = pd.DataFrame(
        {
            "Season": ["DJF", "MAM", "JJA", "SON"] * 8,
            "cruise_group": ["A", "B", "C", "D"] * 8,
            "Year": np.repeat(np.arange(2000, 2008), 4),
        }
    )
    result = module.season_group_redundancy(
        frame,
        season_col="Season",
        group_col="cruise_group",
        year_col="Year",
        permutations=99,
        seed=42,
    ).iloc[0]

    assert result["normalized_mutual_information"] == 1.0
    assert result["cramers_v_bias_corrected"] == 1.0
    assert result["permutations_completed"] == 99


def test_redundancy_plot_uses_ordered_seasons_grayscale_and_external_legend(
    monkeypatch, tmp_path
):
    module = load_module()
    captured = {}

    def capture(fig, *_args, **_kwargs):
        captured["fig"] = fig

    monkeypatch.setattr(module, "save_figure_all_formats", capture)
    counts = pd.DataFrame(
        {"cruise_group_1": [3, 4, 5, 6], "cruise_group_2": [6, 5, 4, 3]},
        index=["SON", "JJA", "DJF", "MAM"],
    )
    module.plot_redundancy(counts, tmp_path / "redundancy")

    left, right = captured["fig"].axes[:2]
    assert [tick.get_text() for tick in left.get_yticklabels()] == [
        "DJF", "MAM", "JJA", "SON"
    ]
    assert left.images[0].get_cmap().name == "Greys"
    assert [text.get_text() for text in left.texts] == [
        "5", "4", "6", "3", "4", "5", "3", "6"
    ]
    for patch in right.patches:
        red, green, blue, _ = patch.get_facecolor()
        assert red == green == blue
    # Group 2 is drawn first at the bottom; group 1 is stacked above it.
    bars_per_group = len(counts.index)
    assert all(
        patch.get_y() == 0 for patch in right.patches[:bars_per_group]
    )
    assert all(
        patch.get_y() > 0 for patch in right.patches[bars_per_group:]
    )
    assert [text.get_text() for text in right.get_legend().get_texts()] == [
        "cruise_group_1", "cruise_group_2"
    ]
    assert right.get_legend().get_bbox_to_anchor()._bbox.x0 > 1.0
    plt.close(captured["fig"])


def test_leave_year_out_comparison_can_add_cruise_group_beyond_season():
    module = load_module()
    rows = []
    for year in range(2000, 2008):
        for sample in range(8):
            season = ["DJF", "MAM", "JJA", "SON"][sample % 4]
            group = "G1" if sample < 4 else "G2"
            rows.append(
                {
                    "Year": year,
                    "Season": season,
                    "cruise_group": group,
                    "Sparse feature": (
                        1.0
                        + 0.2 * (sample % 4)
                        + (4.0 if group == "G2" else 0.0)
                        + 0.01 * (year - 2000)
                    ),
                }
            )
    frame = pd.DataFrame(rows)
    folds = module.leave_year_out_cv(
        frame,
        feature="Sparse feature",
        season_col="Season",
        group_col="cruise_group",
        year_col="Year",
        analysis_set="all_assignments",
    )
    summary = module.summarize_cv(folds)

    season_r2 = summary.loc[summary["model"] == "Season", "pooled_cv_r2"].iloc[0]
    combined_r2 = summary.loc[
        summary["model"] == "Season + cruise group", "pooled_cv_r2"
    ].iloc[0]
    assert folds["fold_year"].nunique() == 8
    assert (folds.groupby("model")["n_test"].sum() == len(frame)).all()
    assert combined_r2 > season_r2
