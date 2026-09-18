import importlib.util
import sys
import types
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT = Path(__file__).parents[1] / "processes/stratification_metrics/env_stratification_metrics.py"


def _load_module():
    sys.modules.setdefault("gsw", types.SimpleNamespace())
    spec = importlib.util.spec_from_file_location("env_stratification_metrics", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_pea_kmeans_and_seasonal_intrusion_classification():
    module = _load_module()
    rng = np.random.default_rng(3)
    pea = np.concatenate(
        [rng.normal(100, 10, 40), rng.normal(250, 15, 40), rng.normal(420, 20, 40)]
    )
    dates = pd.date_range("2008-01-15", periods=len(pea), freq="45D")
    lower_sigma0 = np.asarray(
        24
        + 0.2 * np.sin(2 * np.pi * dates.dayofyear / 365.2425)
        + rng.normal(0, 0.025, len(pea))
    )
    lower_sigma0[[30, 31, 80, 81, 110]] += 0.18
    summary = pd.DataFrame(
        {
            "profile_date": dates,
            "pea_J_m3": pea,
            "sigma0_lower_mean_kg_m3": lower_sigma0,
        }
    )

    classified, model, bootstrap = module._classify_pea(
        summary,
        mode="global",
        bootstrap_iterations=10,
        random_state=42,
    )
    assert set(classified["pea_class"]) == {"mixed", "intermediate", "stratified"}
    assert classified["pea_class_method"].eq("kmeans_3class").all()
    assert len(model) == 3
    assert model["cluster_center_J_m3"].is_monotonic_increasing
    assert len(bootstrap) == 10

    classified, diagnostics = module._classify_deep_intrusion(
        classified,
        mode="global",
        quantile=0.90,
        value_col="sigma0_lower_mean_kg_m3",
        date_col="profile_date",
    )
    intrusions = classified["deep_intrusion_class"].eq("intrusion")
    assert intrusions.any()
    expected_intrusions = classified["deep_intrusion_density_anomaly_kg_m3"].ge(
        classified["deep_intrusion_threshold"]
    )
    assert intrusions.equals(expected_intrusions)
    assert diagnostics.loc[0, "residual_quantile"] == 0.90
    assert diagnostics["method"].eq("centered_3month_median_residual").all()
    assert not diagnostics["requires_positive_change"].any()


def test_oxygen_intrusion_uses_matched_depth_onset_and_persists_to_baseline():
    module = _load_module()
    dates = pd.date_range("2008-01-15", periods=5, freq="30D")
    summary = pd.DataFrame({
        "Cruise": range(1, 6),
        "profile_date": dates,
        "deep_intrusion_class": ["baseline", "intrusion", "baseline", "baseline", "baseline"],
    })
    values = {
        1: [20, 10, 5],
        2: [25, 15, 10],
        3: [23, 13, 8],
        4: [20, 10, 5],
        5: [100, 105, 110],
    }
    observations = pd.DataFrame([
        {"Cruise": cruise, "Depth": depth, "Oxygen": oxygen}
        for cruise, oxygen_values in values.items()
        for depth, oxygen in zip([60, 100, 150], oxygen_values)
    ])
    classified, events = module._classify_oxygen_intrusions(
        summary, observations, ["Cruise"], "Depth", "Oxygen", 40.0, 90.0,
        bottom_n=3, onset_threshold=4.5, persistence_threshold=4.5,
        end_consecutive=2,
    )
    assert classified["oxygen_intrusion_class"].tolist() == [
        "baseline", "intrusion", "baseline", "baseline", "baseline"
    ]
    assert classified["oxygen_intrusion_onset"].tolist() == [False, True, False, False, False]
    assert classified.loc[1, "oxygen_intrusion_density_supported"]
    assert classified.loc[4, "oxygen_compartment_breakdown"]
    assert classified.loc[1, "oxygen_intrusion_last_active"]
    assert not classified.loc[1, "oxygen_intrusion_end_confirmed"]
    assert classified.loc[3, "oxygen_intrusion_end_confirmed"]
    assert classified.loc[3, "oxygen_intrusion_class"] == "baseline"
    assert len(events) == 1
    assert events.loc[0, "n_cruises"] == 1
    assert str(events.loc[0, "end_confirmed_date"].date()) == "2008-04-14"


def test_nitrate_qualifies_renewal_and_controls_post_renewal():
    module = _load_module()
    dates = pd.date_range("2008-01-15", periods=8, freq="30D")
    summary = pd.DataFrame({
        "Cruise": range(1, 9),
        "profile_date": dates,
        "oxygen_intrusion_onset": [False, True, False, False, False, False, True, True],
        "oxygen_intrusion_change_median_um": [np.nan, 10, 0, 0, 0, 0, 8, 9],
        "oxygen_intrusion_bottom_median_um": [5, 15, 5, 4, 3, 2, 12, 13],
    })
    nitrate_values = {
        1: [0, 0, 0],
        2: [5, 10, 15],
        3: [4, 8, 12],       # O2 may decline; nitrate maintains post-renewal.
        4: [-1, -1, 7],      # Negative is unmeasured; phase is unresolved.
        5: [2, 4, 6],        # Detection resumes the same post-renewal event.
        6: [0, 0, 0],        # Adequate measured non-detect terminates it.
        7: [0, 0, 0],        # O2 event with measured nitrate absence.
        8: [-1, -1, 0.2],    # O2 event with insufficient nitrate coverage.
    }
    observations = pd.DataFrame([
        {"Cruise": cruise, "Depth": depth, "Nitrate": nitrate}
        for cruise, values in nitrate_values.items()
        for depth, nitrate in zip([60, 100, 150], values)
    ])

    classified, events, anomalies = module._classify_nitrate_renewal(
        summary, observations, ["Cruise"], "Depth", "Nitrate", 40.0,
        bottom_n=3, min_depths=2, detection_limit=0.0,
    )

    assert classified["renewal_phase"].tolist() == [
        "baseline", "renewal", "post-renewal", "unknown",
        "post-renewal", "baseline", "baseline", "unknown",
    ]
    assert classified["renewal_onset"].tolist() == [
        False, True, False, False, False, False, False, False,
    ]
    assert classified.loc[2, "renewal_phase"] == "post-renewal"
    assert classified.loc[3, "renewal_nitrate_status"] == "insufficient_coverage"
    assert classified.loc[5, "renewal_end_confirmed"]
    assert classified.loc[6, "oxygen_anomaly_class"] == "O2-only"
    assert classified.loc[7, "oxygen_anomaly_class"] == "nitrate-unresolved"
    assert len(events) == 1
    assert events.loc[0, "n_supported_cruises"] == 3
    assert events.loc[0, "n_post_renewal_cruises"] == 2
    assert anomalies["oxygen_anomaly_class"].tolist() == [
        "O2-only", "nitrate-unresolved"
    ]

    bridged, bridged_events, _ = module._classify_nitrate_renewal(
        summary, observations, ["Cruise"], "Depth", "Nitrate", 40.0,
        bottom_n=3, min_depths=2, detection_limit=0.0,
        bridge_enabled=True, bridge_max_cruises=2, bridge_max_days=150,
    )
    assert bridged.loc[3, "renewal_phase"] == "post-renewal"
    assert bridged.loc[3, "renewal_phase_inferred"]
    assert (
        bridged.loc[3, "renewal_inference_reason"]
        == "bracketed_insufficient_nitrate_coverage"
    )
    assert not bridged.loc[5, "renewal_phase_inferred"]
    assert bridged_events.loc[0, "n_inferred_post_renewal_cruises"] == 1


def test_four_depth_nitrate_window_can_use_depth_above_oxygen_window():
    module = _load_module()
    summary = pd.DataFrame({
        "Cruise": [1, 2],
        "profile_date": pd.to_datetime(["2008-01-15", "2008-02-15"]),
        "oxygen_intrusion_onset": [False, True],
        "oxygen_intrusion_change_median_um": [np.nan, 10],
        "oxygen_intrusion_bottom_median_um": [5, 15],
    })
    observations = pd.DataFrame([
        {"Cruise": cruise, "Depth": depth, "Nitrate": nitrate}
        for cruise, values in {
            1: [0, 0, 0, 0],
            # Only the upper two of the four-depth nitrate window are valid.
            2: [8, 0.2, -1, -1],
        }.items()
        for depth, nitrate in zip([50, 60, 100, 150], values)
    ])

    classified, events, anomalies = module._classify_nitrate_renewal(
        summary, observations, ["Cruise"], "Depth", "Nitrate", 40.0,
        bottom_n=4, min_depths=2, detection_limit=0.0,
    )

    assert classified.loc[1, "renewal_phase"] == "renewal"
    assert classified.loc[1, "renewal_nitrate_bottom_depths_m"] == "50,60,100,150"
    assert classified.loc[1, "renewal_nitrate_bottom_depths_n"] == 2
    assert len(events) == 1
    assert anomalies.empty


def test_pea_kmeans_boundaries_are_data_driven():
    module = _load_module()
    summary = pd.DataFrame({"pea_J_m3": np.linspace(100, 200, 120)})
    classified, model, bootstrap = module._classify_pea(
        summary,
        mode="global",
        bootstrap_iterations=5,
        random_state=42,
    )
    assert classified["pea_class_method"].eq("kmeans_3class").all()
    assert len(bootstrap) == 5
    assert not model.empty
    assert set(classified["pea_class"]) == {"mixed", "intermediate", "stratified"}


def test_multimetric_physical_regime_gmm_outputs_probabilities_and_diagnostics():
    module = _load_module()
    rng = np.random.default_rng(8)
    rows = []
    for group, pea_center in enumerate([80, 240, 460]):
        for index in range(20):
            rows.append(
                {
                    "profile_label": f"{group}-{index}",
                    "profile_date": pd.Timestamp("2010-01-01") + pd.Timedelta(days=len(rows) * 30),
                    "pea_J_m3": rng.normal(pea_center, 8),
                    "n2_max_s-2": rng.normal([2e-5, 8e-5, 1.8e-4][group], 2e-6),
                    "sigma0_upper_lower_diff_kg_m3": rng.normal([0.15, 0.5, 1.0][group], 0.03),
                    "mld_depth_m_dr0p03": rng.normal([150, 70, 25][group], 3),
                    "pycnocline_depth_m": rng.normal([155, 90, 45][group], 4),
                    "depth_min_m": 10,
                    "depth_max_m": 200,
                }
            )
    summary = pd.DataFrame(rows)
    classified, selection, clusters, projection = module._classify_physical_regimes(
        summary, random_state=42, k_max=5
    )
    assert set(classified["physical_regime_class"]) == {"mixed", "intermediate", "stratified"}
    probabilities = classified[
        [
            "physical_regime_probability_mixed",
            "physical_regime_probability_intermediate",
            "physical_regime_probability_stratified",
        ]
    ]
    assert np.allclose(probabilities.sum(axis=1), 1)
    assert selection["components"].tolist() == [1, 2, 3, 4, 5]
    assert clusters["physical_regime_class"].tolist() == ["mixed", "intermediate", "stratified"]
    assert len(projection) == len(summary)
