from __future__ import annotations

import importlib.util
from pathlib import Path

import matplotlib.figure
import numpy as np
import pandas as pd


ROOT = Path(__file__).parents[1]


def load_script(relative: str, name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def frame() -> pd.DataFrame:
    dates = pd.date_range("2010-01-15", periods=24, freq="30D")
    month = dates.month
    return pd.DataFrame({
        "date": dates,
        "cruise_group": np.where(np.arange(24) % 2, "cruise_group_1", "cruise_group_2"),
        "assignment_uncertain": np.arange(24) % 7 == 0,
        "oxygen_intrusion_class": np.where(np.arange(24) % 5 == 0, "intrusion", "baseline"),
        "depth_centroid_distance": 80 + 15 * np.sin(2 * np.pi * month / 12),
        "pea_J_m3": 250 + 50 * np.cos(2 * np.pi * month / 12),
        "pea_upper_J_m3": 30 + 10 * np.cos(2 * np.pi * month / 12),
        "pea_lower_J_m3": 80 + 15 * np.sin(2 * np.pi * month / 12),
    })


def assert_external_right_legend(fig, panel_prefixes: list[str]) -> None:
    assert len(fig.legends) == 1
    legend = fig.legends[0]
    assert legend._loc == 2  # upper left, anchored in reserved right margin
    assert legend.get_bbox_to_anchor()._bbox.x0 >= 0.80
    visible_axes = [axis for axis in fig.axes if axis.get_visible()]
    assert [axis.get_title(loc="left")[:1] for axis in visible_axes] == panel_prefixes


def test_renewal_phase_prefers_nitrate_qualified_classification() -> None:
    module = load_script(
        "processes/eof_state_cluster/eof_state_interpretation.py",
        "eof_state_interpretation_renewal_phase",
    )
    data = pd.DataFrame({
        "renewal_phase": [
            "renewal", "post-renewal", "unknown", "baseline", "baseline",
        ],
        "oxygen_intrusion_class": [
            "baseline", "baseline", "intrusion", "intrusion", "baseline",
        ],
        "oxygen_intrusion_onset": [False, False, True, True, False],
    })
    assert module.renewal_phase(data).tolist() == [
        "renewal", "post-renewal", "baseline", "baseline", "baseline",
    ]


def test_stratification_merge_retains_nitrate_qualified_renewal_fields(tmp_path) -> None:
    module = load_script(
        "processes/stratification_anomaly/env_stratification_anomaly_detection.py",
        "env_stratification_anomaly_renewal_fields",
    )
    source = pd.DataFrame({
        "profile_date": ["2009-05-13"],
        "Cruise": [33],
        "pea_J_m3": [102.35],
        "oxygen_intrusion_class": ["intrusion"],
        "renewal_phase": ["unknown"],
        "renewal_onset": [False],
        "renewal_phase_inferred": [False],
        "renewal_nitrate_status": ["insufficient_coverage"],
        "renewal_nitrate_bottom_depths_n": [1],
        "renewal_nitrate_bottom_depths_m": ["165,185,200"],
        "oxygen_anomaly_class": ["nitrate-unresolved"],
    })
    path = tmp_path / "stratification_summary.tsv"
    source.to_csv(path, sep="\t", index=False)

    retained = module._load_pea_timeseries(path, "profile_date")

    assert retained.loc[0, "renewal_phase"] == "unknown"
    assert not bool(retained.loc[0, "renewal_phase_inferred"])
    assert retained.loc[0, "renewal_nitrate_status"] == "insufficient_coverage"
    assert retained.loc[0, "oxygen_anomaly_class"] == "nitrate-unresolved"
    assert retained.loc[0, "oxygen_intrusion_class"] == "intrusion"


def test_cruise_group_monthly_figures_use_external_legends_and_panel_letters(monkeypatch) -> None:
    module = load_script(
        "processes/eof_state_cluster/eof_state_interpretation.py", "eof_state_interpretation"
    )
    captured = []
    monkeypatch.setattr(module, "save", lambda fig, base, formats: captured.append(fig))
    data = frame()
    colors = {"cruise_group_1": "#0072B2", "cruise_group_2": "#E69F00"}
    module.monthly_group_profile(
        data,
        [("depth_centroid_distance", "Depth-centroid distance (D)", "Depth centroid"),
         ("pea_J_m3", "PEA (J m$^{-3}$)", "PEA")],
        "cruise_group", colors, Path("unused"), ["png"],
    )
    module.maintext_monthly_profile(
        data, data, "depth_centroid_distance", "cruise_group", colors,
        Path("unused"), ["png"],
    )
    assert_external_right_legend(captured[0], ["A", "B"])
    assert len(captured[1].legends) == 0
    assert [axis.get_title(loc="left")[:1] for axis in captured[1].axes] == ["A", "B"]
    assert all(axis.get_legend() is not None for axis in captured[1].axes)
    assert [text.get_text() for text in captured[1].axes[0].get_legend().get_texts()] == [
        "All-cruise observation", "Monthly median", "Monthly IQR",
        "Annual maximum", "Annual minimum",
    ]
    assert [text.get_text() for text in captured[1].axes[1].get_legend().get_texts()] == [
        "cruise_group_1", "cruise_group_2", "Renewal", "Post-renewal",
        "Monthly median", "Monthly IQR",
    ]
    top_box, bottom_box = captured[1].axes[0].get_position(), captured[1].axes[1].get_position()
    assert top_box.y0 > bottom_box.y0
    assert abs(top_box.x0 - bottom_box.x0) < 0.02


def test_cruise_group_grayscale_matches_season_redundancy_mapping() -> None:
    module = load_script(
        "processes/eof_state_cluster/eof_state_interpretation.py",
        "eof_state_interpretation_grayscale",
    )
    colors = module.cruise_group_grayscale(
        ["cruise_group_2", "cruise_group_1"]
    )
    assert list(colors) == ["cruise_group_1", "cruise_group_2"]
    expected = module.plt.cm.Greys(np.linspace(0.35, 0.85, 2))
    for observed, target in zip(colors.values(), expected):
        assert np.allclose(observed, target)
        assert observed[0] == observed[1] == observed[2]
    assert module.cruise_group_marker_edge(colors["cruise_group_1"]) == "0.25"
    assert module.cruise_group_marker_edge(colors["cruise_group_1"], "renewal") == "0.25"


def test_stratification_monthly_figure_uses_external_legend_and_panel_letters(monkeypatch) -> None:
    module = load_script(
        "processes/stratification_anomaly/env_stratification_anomaly_detection.py",
        "env_stratification_anomaly_detection",
    )
    captured = []
    monkeypatch.setattr(matplotlib.figure.Figure, "savefig", lambda self, *a, **k: None)
    monkeypatch.setattr(module.plt, "close", lambda fig: captured.append(fig))
    data = frame().rename(columns={"depth_centroid_distance": "stratification_score"})
    data["coverage"] = 1.0
    module.plot_physical_biochem_monthly_profile(data, Path("unused"))
    assert_external_right_legend(captured[-1], ["A", "B", "C", "D"])
