#!/usr/bin/env python3
"""
env_stratification_metrics.py

Compute standard physical-oceanography stratification metrics from CTD data:
  - In-situ density (rho) and potential density anomaly (sigma0)
  - Brunt-Vaisala frequency (N^2)
  - Mixed Layer Depth (MLD) by density threshold
  - Potential Energy Anomaly (PEA, Simpson-Hunter)
  - Pycnocline depth (max N^2) and layer-specific metrics
  - Optional global layer split via mld125

Inputs: salinity, temperature, depth, and either pressure or latitude (for depth->pressure).
Outputs: density profiles, N^2 profiles, per-profile summary metrics, upper/lower feature traces, and review plots.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shared_plot_style import install_publication_style

install_publication_style()

from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler

try:
    import gsw  # TEOS-10
except ImportError as exc:
    raise SystemExit(
        "Missing dependency 'gsw'. Install with: conda install -c conda-forge gsw"
    ) from exc


def _split_csv(s: str | None) -> list[str]:
    if not s:
        return []
    return [p.strip() for p in s.split(",") if p.strip()]


def _split_float_csv(s: str | None) -> list[float]:
    parts = _split_csv(s)
    if not parts:
        return []
    return [float(p) for p in parts]


def _require_cols(df: pd.DataFrame, cols: Iterable[str]) -> None:
    missing = [c for c in cols if c and c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")


def _infer_plot_feature_cols(df: pd.DataFrame, exclude: set[str]) -> list[str]:
    cols = []
    for col in df.columns:
        if col in exclude:
            continue
        if str(col).startswith("ASV"):
            continue
        vals = pd.to_numeric(df[col], errors="coerce")
        if np.isfinite(vals).any():
            cols.append(col)
    return cols


def _numeric(series: pd.Series) -> np.ndarray:
    return pd.to_numeric(series, errors="coerce").to_numpy()


def _resolve_vector(
    df: pd.DataFrame,
    col: str | None,
    value: float | None,
    name: str,
    n: int,
) -> np.ndarray | None:
    if col and value is not None:
        raise ValueError(f"Provide either --{name} or --{name}-col, not both.")
    if col:
        return _numeric(df[col])
    if value is not None:
        return np.full(n, float(value))
    return None


def _aggregate_by_depth(
    df: pd.DataFrame,
    depth_col: str,
    cols: list[str],
) -> pd.DataFrame:
    grouped = df.groupby(depth_col, sort=True)
    out = grouped[cols].mean().reset_index()
    out["__n_samples"] = grouped.size().to_numpy()
    return out


def _compute_pressure(
    depth: np.ndarray,
    pressure: np.ndarray | None,
    lat: np.ndarray | None,
) -> np.ndarray:
    if pressure is not None:
        return pressure
    if lat is None:
        raise ValueError("Latitude is required to compute pressure from depth.")
    return gsw.p_from_z(-depth, lat)


def _integrate_trapezoid(y: np.ndarray, x: np.ndarray) -> float:
    # np.trapezoid was added in newer NumPy; use trapz for compatibility.
    if hasattr(np, "trapezoid"):
        return float(np.trapezoid(y, x))
    return float(np.trapz(y, x))


def _compute_pea(depth: np.ndarray, sigma0: np.ndarray, lat: float | None) -> float:
    depth_rel = depth - np.nanmin(depth)
    h = np.nanmax(depth_rel)
    if not np.isfinite(h) or h <= 0:
        return np.nan
    rho_bar = _integrate_trapezoid(sigma0, depth_rel) / h
    g = float(gsw.grav(lat, 0)) if lat is not None else 9.81
    pea = g / h * _integrate_trapezoid((sigma0 - rho_bar) * depth_rel, depth_rel)
    return float(pea)


def _compute_mld(
    depth: np.ndarray,
    sigma0: np.ndarray,
    ref_depth: float,
    delta_rho: float,
) -> tuple[float, float, float]:
    if len(depth) == 0:
        return (np.nan, np.nan, np.nan)
    idx_ref = int(np.nanargmin(np.abs(depth - ref_depth)))
    rho_ref = sigma0[idx_ref]
    valid = np.isfinite(sigma0) & np.isfinite(depth)
    depth_v = depth[valid]
    sigma_v = sigma0[valid]
    if len(depth_v) == 0:
        return (np.nan, np.nan, np.nan)
    candidates = depth_v[(sigma_v - rho_ref) <= delta_rho]
    if len(candidates) == 0:
        return (np.nan, rho_ref, np.nan)
    mld = float(np.nanmax(candidates))
    deeper_mask = depth_v > mld
    if not deeper_mask.any():
        delta_below = np.nan
    else:
        idx_below = int(np.nanargmin(depth_v[deeper_mask]))
        sigma_below = sigma_v[deeper_mask][idx_below]
        delta_below = float(sigma_below - rho_ref)
    return (mld, rho_ref, delta_below)


def _format_threshold_tag(value: float) -> str:
    s = f"{value:g}"
    return s.replace(".", "p")


def _adaptive_delta_rho(
    depth: np.ndarray,
    sigma0: np.ndarray,
    layer_max_depth: float,
    percentile: float,
) -> float:
    if depth.size == 0:
        return np.nan
    min_depth = float(np.nanmin(depth))
    surface_mask = depth <= (min_depth + layer_max_depth)
    surface_sigma = sigma0[surface_mask & np.isfinite(sigma0)]
    if surface_sigma.size < 2:
        return np.nan
    low = float(np.nanmin(surface_sigma))
    high = float(np.nanpercentile(surface_sigma, percentile))
    delta = high - low
    return float(delta) if np.isfinite(delta) and delta > 0 else np.nan


def _split_layers(
    depth: np.ndarray,
    split_depth: float,
) -> tuple[np.ndarray, np.ndarray]:
    if not np.isfinite(split_depth):
        return (np.array([], dtype=bool), np.array([], dtype=bool))
    upper = depth <= split_depth
    lower = depth > split_depth
    return upper, lower


def _compute_global_pycnocline_depth(
    df: pd.DataFrame,
    profile_cols: list[str],
    sal_col: str,
    temp_col: str,
    depth_col: str,
    pressure_col: str | None,
    latitude_col: str | None,
    longitude_col: str | None,
    latitude: float | None,
    longitude: float | None,
) -> float:
    max_depth = np.nan
    warned_lat_zero = False

    for _, group in df.groupby(profile_cols):
        g = group.copy()
        g["__sal"] = _numeric(g[sal_col])
        g["__temp"] = _numeric(g[temp_col])
        g["__depth"] = _numeric(g[depth_col])
        g["__press"] = _numeric(g[pressure_col]) if pressure_col else np.nan
        lat_vec = _resolve_vector(g, latitude_col, latitude, "latitude", len(g))
        lon_vec = _resolve_vector(g, longitude_col, longitude, "longitude", len(g))
        g["__lat"] = lat_vec if lat_vec is not None else np.nan
        g["__lon"] = lon_vec if lon_vec is not None else np.nan

        req_mask = np.isfinite(g["__sal"]) & np.isfinite(g["__temp"]) & np.isfinite(g["__depth"])
        g = g.loc[req_mask]
        if g.empty:
            continue

        depth_df = _aggregate_by_depth(
            g,
            depth_col=depth_col,
            cols=["__sal", "__temp", "__press", "__lat", "__lon"],
        )
        depth = _numeric(depth_df[depth_col])
        sal = _numeric(depth_df["__sal"])
        temp = _numeric(depth_df["__temp"])
        lat = _numeric(depth_df["__lat"]) if "__lat" in depth_df else None
        lon = _numeric(depth_df["__lon"]) if "__lon" in depth_df else None
        if lat is not None and np.isnan(lat).all():
            lat = None
        if lon is not None and np.isnan(lon).all():
            lon = None
        pressure = _numeric(depth_df["__press"]) if pressure_col else None

        if lat is None:
            if pressure_col:
                lat = np.zeros(len(depth))
                if not warned_lat_zero:
                    print("[i] Latitude not provided; using 0 deg for SA conversion.")
                    warned_lat_zero = True
            else:
                raise ValueError("Latitude is required when pressure is not provided.")
        if lon is None:
            lon = np.zeros(len(depth))

        pressure = _compute_pressure(depth, pressure, lat)
        if len(depth) < 2:
            continue

        SA = gsw.SA_from_SP(sal, pressure, lon, lat)
        CT = gsw.CT_from_t(SA, temp, pressure)
        n2, p_mid = gsw.Nsquared(SA, CT, pressure, lat)
        if n2.size == 0:
            continue
        lat_mid = float(np.nanmean(lat)) if np.isfinite(np.nanmean(lat)) else 0.0
        depth_mid = -gsw.z_from_p(p_mid, lat_mid)
        idx_max = int(np.nanargmax(n2))
        depth_max = float(depth_mid[idx_max]) if depth_mid.size else np.nan
        if np.isfinite(depth_max):
            max_depth = depth_max if not np.isfinite(max_depth) else max(max_depth, depth_max)

    return max_depth


def _compute_global_mld_depth(
    df: pd.DataFrame,
    profile_cols: list[str],
    sal_col: str,
    temp_col: str,
    depth_col: str,
    pressure_col: str | None,
    latitude_col: str | None,
    longitude_col: str | None,
    latitude: float | None,
    longitude: float | None,
    ref_depth: float,
    delta_rho: float,
    stat: str,
) -> float:
    mld_vals: list[float] = []
    warned_lat_zero = False

    for _, group in df.groupby(profile_cols):
        g = group.copy()
        g["__sal"] = _numeric(g[sal_col])
        g["__temp"] = _numeric(g[temp_col])
        g["__depth"] = _numeric(g[depth_col])
        g["__press"] = _numeric(g[pressure_col]) if pressure_col else np.nan
        lat_vec = _resolve_vector(g, latitude_col, latitude, "latitude", len(g))
        lon_vec = _resolve_vector(g, longitude_col, longitude, "longitude", len(g))
        g["__lat"] = lat_vec if lat_vec is not None else np.nan
        g["__lon"] = lon_vec if lon_vec is not None else np.nan

        req_mask = np.isfinite(g["__sal"]) & np.isfinite(g["__temp"]) & np.isfinite(g["__depth"])
        g = g.loc[req_mask]
        if g.empty:
            continue

        depth_df = _aggregate_by_depth(
            g,
            depth_col=depth_col,
            cols=["__sal", "__temp", "__press", "__lat", "__lon"],
        )
        depth = _numeric(depth_df[depth_col])
        sal = _numeric(depth_df["__sal"])
        temp = _numeric(depth_df["__temp"])
        lat = _numeric(depth_df["__lat"]) if "__lat" in depth_df else None
        lon = _numeric(depth_df["__lon"]) if "__lon" in depth_df else None
        if lat is not None and np.isnan(lat).all():
            lat = None
        if lon is not None and np.isnan(lon).all():
            lon = None
        pressure = _numeric(depth_df["__press"]) if pressure_col else None

        if lat is None:
            if pressure_col:
                lat = np.zeros(len(depth))
                if not warned_lat_zero:
                    print("[i] Latitude not provided; using 0 deg for SA conversion.")
                    warned_lat_zero = True
            else:
                raise ValueError("Latitude is required when pressure is not provided.")
        if lon is None:
            lon = np.zeros(len(depth))

        pressure = _compute_pressure(depth, pressure, lat)
        if len(depth) < 2:
            continue

        SA = gsw.SA_from_SP(sal, pressure, lon, lat)
        CT = gsw.CT_from_t(SA, temp, pressure)
        sigma0 = gsw.sigma0(SA, CT)
        mld_depth, _, _ = _compute_mld(
            depth=depth,
            sigma0=sigma0,
            ref_depth=ref_depth,
            delta_rho=delta_rho,
        )
        if np.isfinite(mld_depth):
            mld_vals.append(float(mld_depth))

    if not mld_vals:
        return np.nan
    arr = np.array(mld_vals, dtype=float)
    if stat == "max":
        return float(np.nanmax(arr))
    if stat == "p90":
        return float(np.nanpercentile(arr, 90))
    if stat == "median":
        return float(np.nanmedian(arr))
    raise ValueError(f"Unknown --layer-split-stat '{stat}'.")


def _profile_label(keys: dict[str, object]) -> str:
    return "|".join([str(v) for v in keys.values()])


def _fit_ordered_pea_kmeans(values: np.ndarray, random_state: int) -> tuple:
    x = np.asarray(values, dtype=float).reshape(-1, 1)
    model = KMeans(
        n_clusters=3,
        n_init=10,
        random_state=random_state,
    ).fit(x)
    centers = np.sort(model.cluster_centers_.ravel())
    boundaries = [(centers[0] + centers[1]) / 2, (centers[1] + centers[2]) / 2]
    return model, centers, boundaries


def _classify_pea(
    summary_df: pd.DataFrame,
    mode: str,
    bootstrap_iterations: int,
    random_state: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if "pea_J_m3" not in summary_df.columns:
        return summary_df, pd.DataFrame(), pd.DataFrame()

    df = summary_df.copy()
    df["pea_class"] = "unknown"
    df["pea_threshold_low"] = np.nan
    df["pea_threshold_high"] = np.nan
    df["pea_class_method"] = "kmeans_3class"
    df["pea_threshold_low_ci_lower"] = np.nan
    df["pea_threshold_low_ci_upper"] = np.nan
    df["pea_threshold_high_ci_lower"] = np.nan
    df["pea_threshold_high_ci_upper"] = np.nan

    values = df["pea_J_m3"].to_numpy(dtype=float)
    valid = np.isfinite(values)
    if not valid.any():
        return df, pd.DataFrame(), pd.DataFrame()

    cluster_rows = []
    bootstrap_rows = []
    if mode != "global":
        raise ValueError("Three-class PEA k-means supports --pea-class-mode global only.")
    x = values[valid]
    if len(x) < 3:
        raise ValueError("Three-class PEA k-means requires at least three valid profiles.")

    model, centers, boundaries = _fit_ordered_pea_kmeans(x, random_state)
    low, high = boundaries
    df["pea_threshold_low"] = low
    df["pea_threshold_high"] = high
    df.loc[valid & (values <= low), "pea_class"] = "mixed"
    df.loc[valid & (values > low) & (values < high), "pea_class"] = "intermediate"
    df.loc[valid & (values >= high), "pea_class"] = "stratified"
    labels = ["mixed", "intermediate", "stratified"]
    score = silhouette_score(x.reshape(-1, 1), model.labels_) if len(x) > 3 else np.nan
    for label, center in zip(labels, centers):
        cluster_rows.append(
            {
                "method": "kmeans_3class",
                "pea_class": label,
                "cluster_center_J_m3": center,
                "n_profiles": int((df["pea_class"] == label).sum()),
                "pea_threshold_low": low,
                "pea_threshold_high": high,
                "inertia": model.inertia_,
                "silhouette_score": score,
            }
        )

    rng = np.random.default_rng(random_state)
    for iteration in range(bootstrap_iterations):
        sample = rng.choice(x, size=len(x), replace=True)
        try:
            _, _, boot_boundaries = _fit_ordered_pea_kmeans(sample, random_state + iteration + 1)
            bootstrap_rows.append(
                {
                    "iteration": iteration + 1,
                    "pea_threshold_low": boot_boundaries[0],
                    "pea_threshold_high": boot_boundaries[1],
                }
            )
        except (ValueError, np.linalg.LinAlgError):
            continue
    if bootstrap_rows:
        boot = pd.DataFrame(bootstrap_rows)
        df["pea_threshold_low_ci_lower"] = boot["pea_threshold_low"].quantile(0.025)
        df["pea_threshold_low_ci_upper"] = boot["pea_threshold_low"].quantile(0.975)
        df["pea_threshold_high_ci_lower"] = boot["pea_threshold_high"].quantile(0.025)
        df["pea_threshold_high_ci_upper"] = boot["pea_threshold_high"].quantile(0.975)
    return df, pd.DataFrame(cluster_rows), pd.DataFrame(
        bootstrap_rows,
        columns=["iteration", "pea_threshold_low", "pea_threshold_high"],
    )


PHYSICAL_REGIME_FEATURES = [
    "log1p_pea",
    "log1p_n2_max",
    "sigma0_upper_lower_diff",
    "relative_mld03",
    "relative_pycnocline_depth",
]


def _physical_regime_matrix(summary_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    required = [
        "pea_J_m3", "n2_max_s-2", "sigma0_upper_lower_diff_kg_m3",
        "mld_depth_m_dr0p03", "pycnocline_depth_m", "depth_min_m", "depth_max_m",
    ]
    missing = [column for column in required if column not in summary_df.columns]
    if missing:
        raise ValueError(f"Physical-regime classification missing columns: {missing}")
    raw = summary_df[required].apply(pd.to_numeric, errors="coerce")
    thickness = raw["depth_max_m"] - raw["depth_min_m"]
    matrix = pd.DataFrame(index=summary_df.index)
    matrix["log1p_pea"] = np.log1p(raw["pea_J_m3"].clip(lower=0))
    matrix["log1p_n2_max"] = np.log1p(raw["n2_max_s-2"].clip(lower=0) * 1e6)
    matrix["sigma0_upper_lower_diff"] = raw["sigma0_upper_lower_diff_kg_m3"]
    matrix["relative_mld03"] = (raw["mld_depth_m_dr0p03"] - raw["depth_min_m"]) / thickness
    matrix["relative_pycnocline_depth"] = (
        raw["pycnocline_depth_m"] - raw["depth_min_m"]
    ) / thickness
    matrix = matrix.replace([np.inf, -np.inf], np.nan)
    return matrix, raw


def _normalized_entropy(probabilities: np.ndarray) -> np.ndarray:
    return -(
        probabilities * np.log(np.clip(probabilities, 1e-12, 1.0))
    ).sum(axis=1) / np.log(probabilities.shape[1])


def _classify_physical_regimes(
    summary_df: pd.DataFrame,
    random_state: int,
    k_max: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    df = summary_df.copy()
    matrix, raw = _physical_regime_matrix(df)
    valid = matrix.notna().all(axis=1)
    if valid.sum() < 6:
        raise ValueError("Physical-regime GMM requires at least six complete cruise profiles.")
    x = matrix.loc[valid, PHYSICAL_REGIME_FEATURES].to_numpy(dtype=float)
    scaler = StandardScaler().fit(x)
    z = scaler.transform(x)

    selection_rows = []
    upper_k = min(max(3, k_max), len(z) - 1)
    for k in range(1, upper_k + 1):
        candidate = GaussianMixture(
            n_components=k, covariance_type="full", n_init=10,
            reg_covar=1e-6, random_state=random_state,
        ).fit(z)
        probabilities = candidate.predict_proba(z)
        labels = probabilities.argmax(axis=1)
        entropy_total = float(-np.sum(probabilities * np.log(np.clip(probabilities, 1e-12, 1))))
        counts = np.bincount(labels, minlength=k)
        selection_rows.append(
            {
                "components": k,
                "aic": candidate.aic(z),
                "bic": candidate.bic(z),
                "icl": candidate.bic(z) + 2 * entropy_total,
                "silhouette_score": (
                    silhouette_score(z, labels) if k > 1 and len(np.unique(labels)) > 1 else np.nan
                ),
                "mean_responsibility_entropy": float(_normalized_entropy(probabilities).mean()) if k > 1 else 0.0,
                "minimum_cluster_n": int(counts.min()),
                "minimum_cluster_fraction": float(counts.min() / len(z)),
                "used_for_classification": k == 3,
            }
        )

    model = GaussianMixture(
        n_components=3, covariance_type="full", n_init=20,
        reg_covar=1e-6, random_state=random_state,
    ).fit(z)
    centers_z = model.means_
    # PEA, N2, and density contrast increase with stratification; deeper MLD decreases it.
    ranking_score = centers_z @ np.asarray([1.0, 1.0, 1.0, -1.0, 0.0])
    order = np.argsort(ranking_score)
    probabilities = model.predict_proba(z)[:, order]
    labels = np.asarray(["mixed", "intermediate", "stratified"], dtype=object)
    assigned = labels[np.argmax(probabilities, axis=1)]

    df["physical_regime_class"] = "unknown"
    df["physical_regime_probability_mixed"] = np.nan
    df["physical_regime_probability_intermediate"] = np.nan
    df["physical_regime_probability_stratified"] = np.nan
    df["physical_regime_max_probability"] = np.nan
    df["physical_regime_entropy"] = np.nan
    df.loc[valid, "physical_regime_class"] = assigned
    df.loc[valid, "physical_regime_probability_mixed"] = probabilities[:, 0]
    df.loc[valid, "physical_regime_probability_intermediate"] = probabilities[:, 1]
    df.loc[valid, "physical_regime_probability_stratified"] = probabilities[:, 2]
    df.loc[valid, "physical_regime_max_probability"] = probabilities.max(axis=1)
    df.loc[valid, "physical_regime_entropy"] = _normalized_entropy(probabilities)

    projection = PCA(n_components=2, random_state=random_state).fit_transform(z)
    projection_df = df.loc[valid, [c for c in ["profile_label", "profile_date"] if c in df]].copy()
    projection_df["physical_pc1"] = projection[:, 0]
    projection_df["physical_pc2"] = projection[:, 1]
    projection_df["physical_regime_class"] = assigned
    projection_df["physical_regime_max_probability"] = probabilities.max(axis=1)

    cluster_rows = []
    for ordered_index, label in enumerate(labels):
        component = int(order[ordered_index])
        member_index = df.index[valid][assigned == label]
        row = {
            "physical_regime_class": label,
            "n_profiles": int(len(member_index)),
            "stratification_ranking_score": float(ranking_score[component]),
        }
        for feature_index, feature in enumerate(PHYSICAL_REGIME_FEATURES):
            row[f"standardized_center_{feature}"] = float(centers_z[component, feature_index])
        for column in raw.columns:
            row[f"median_{column}"] = float(raw.loc[member_index, column].median())
        cluster_rows.append(row)

    return df, pd.DataFrame(selection_rows), pd.DataFrame(cluster_rows), projection_df


def _classify_deep_intrusion(
    summary_df: pd.DataFrame,
    mode: str,
    quantile: float,
    value_col: str,
    date_col: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if value_col not in summary_df.columns:
        return summary_df, pd.DataFrame()

    df = summary_df.copy()
    df["deep_intrusion_class"] = "unknown"
    df["deep_intrusion_threshold"] = np.nan
    df["deep_intrusion_score"] = np.nan
    df["deep_intrusion_expected_sigma0_kg_m3"] = np.nan
    df["deep_intrusion_density_anomaly_kg_m3"] = np.nan
    df["deep_intrusion_density_change_kg_m3"] = np.nan

    if mode == "none":
        return df, pd.DataFrame()
    if mode != "global":
        raise ValueError(f"Unknown --deep-intrusion-mode '{mode}'.")

    values = df[value_col].to_numpy()
    valid = np.isfinite(values)
    if not valid.any():
        return df, pd.DataFrame()

    dates = pd.to_datetime(df[date_col], errors="coerce")
    seasonal_valid = valid & dates.notna().to_numpy()
    if seasonal_valid.sum() < 12:
        raise ValueError("Seasonal intrusion classification requires at least 12 dated profiles.")

    months = dates.dt.month.to_numpy(dtype=float)
    expected_all = np.full(len(df), np.nan)
    diagnostic_rows = []
    for target_month in range(1, 13):
        window_months = [((target_month - 2) % 12) + 1, target_month, (target_month % 12) + 1]
        window_mask = seasonal_valid & np.isin(months, window_months)
        target_mask = seasonal_valid & (months == target_month)
        expected_month = float(np.nanmedian(values[window_mask])) if window_mask.any() else np.nan
        expected_all[target_mask] = expected_month
        diagnostic_rows.append(
            {
                "target_month": target_month,
                "centered_window_months": ",".join(str(month) for month in window_months),
                "n_window_profiles": int(window_mask.sum()),
                "n_target_profiles": int(target_mask.sum()),
                "expected_lower_sigma0_kg_m3": expected_month,
            }
        )

    anomaly_all = values - expected_all
    anomaly = anomaly_all[seasonal_valid]
    threshold = float(np.nanquantile(anomaly, quantile))
    median = float(np.nanmedian(anomaly))
    mad = float(np.nanmedian(np.abs(anomaly - median)))
    robust_scale = 1.4826 * mad

    df.loc[seasonal_valid, "deep_intrusion_expected_sigma0_kg_m3"] = expected_all[seasonal_valid]
    df.loc[seasonal_valid, "deep_intrusion_density_anomaly_kg_m3"] = anomaly
    if robust_scale > 0:
        df.loc[seasonal_valid, "deep_intrusion_score"] = (anomaly - median) / robust_scale
    df["deep_intrusion_threshold"] = threshold

    ordered_index = df.loc[seasonal_valid].assign(__date=dates[seasonal_valid]).sort_values("__date").index
    changes = df.loc[ordered_index, value_col].astype(float).diff()
    df.loc[ordered_index, "deep_intrusion_density_change_kg_m3"] = changes.to_numpy()
    intrusion = df["deep_intrusion_density_anomaly_kg_m3"].ge(threshold)
    df.loc[seasonal_valid, "deep_intrusion_class"] = "baseline"
    df.loc[intrusion, "deep_intrusion_class"] = "intrusion"

    diagnostics = pd.DataFrame(diagnostic_rows)
    diagnostics.insert(0, "method", "centered_3month_median_residual")
    diagnostics["n_profiles"] = int(seasonal_valid.sum())
    diagnostics["residual_quantile"] = quantile
    diagnostics["residual_threshold_kg_m3"] = threshold
    diagnostics["residual_median_kg_m3"] = median
    diagnostics["residual_mad_kg_m3"] = mad
    diagnostics["n_intrusions"] = int(intrusion.sum())
    diagnostics["requires_positive_change"] = False
    return df, diagnostics


def _classify_oxygen_intrusions(
    summary_df: pd.DataFrame,
    observations: pd.DataFrame,
    profile_cols: list[str],
    depth_col: str,
    oxygen_col: str,
    split_depth: float,
    low_oxygen_max: float,
    bottom_n: int,
    onset_threshold: float,
    persistence_threshold: float,
    end_consecutive: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Detect bottom-water O2 intrusions with sensor-aware onset/end hysteresis."""
    df = summary_df.copy()
    output_defaults = {
        "oxygen_intrusion_class": "unknown",
        "oxygen_intrusion_onset": False,
        "oxygen_intrusion_event_id": pd.NA,
        "oxygen_low_compartment_present": False,
        "oxygen_compartment_breakdown": False,
        "oxygen_intrusion_bottom_depths_n": 0,
        "oxygen_intrusion_bottom_depths_m": "",
        "oxygen_intrusion_bottom_median_um": np.nan,
        "oxygen_intrusion_change_median_um": np.nan,
        "oxygen_intrusion_baseline_change_median_um": np.nan,
        "oxygen_intrusion_end_pending": False,
        "oxygen_intrusion_last_active": False,
        "oxygen_intrusion_end_confirmed": False,
        "oxygen_intrusion_density_supported": False,
        "oxygen_intrusion_onset_threshold_um": onset_threshold,
        "oxygen_intrusion_persistence_threshold_um": persistence_threshold,
        "oxygen_intrusion_end_consecutive_cruises": end_consecutive,
    }
    for column, default in output_defaults.items():
        df[column] = default

    if oxygen_col not in observations or not np.isfinite(split_depth):
        return df, pd.DataFrame()

    work = observations[profile_cols + [depth_col, oxygen_col]].copy()
    work[depth_col] = pd.to_numeric(work[depth_col], errors="coerce")
    work[oxygen_col] = pd.to_numeric(work[oxygen_col], errors="coerce")
    work = work.dropna(subset=[depth_col, oxygen_col])
    profiles = {
        tuple(keys) if isinstance(keys, tuple) else (keys,): group.groupby(depth_col)[oxygen_col].mean()
        for keys, group in work.groupby(profile_cols, sort=False)
    }
    ordered = df.assign(__date=pd.to_datetime(df["profile_date"], errors="coerce")).sort_values("__date")
    active = False
    event_id = 0
    baseline_profile: pd.Series | None = None
    previous_profile: pd.Series | None = None
    pending_end_indices: list[int] = []
    end_confirmed_dates: dict[int, pd.Timestamp] = {}
    event_rows: list[dict] = []

    def paired_change(current: pd.Series, reference: pd.Series) -> np.ndarray:
        common = current.index.intersection(reference.index)
        common = common[common > split_depth]
        common = common.sort_values()[-bottom_n:]
        return (current.loc[common] - reference.loc[common]).dropna().to_numpy(dtype=float)

    for idx, row in ordered.iterrows():
        key = tuple(row[col] for col in profile_cols)
        current = profiles.get(key, pd.Series(dtype=float))
        current_bottom_depths = sorted(float(x) for x in current.index if x > split_depth)[-bottom_n:]
        if len(current_bottom_depths) == bottom_n:
            df.at[idx, "oxygen_intrusion_bottom_median_um"] = float(
                np.median(current.loc[current_bottom_depths].to_numpy(dtype=float))
            )
        low_present = bool((current <= low_oxygen_max).any())
        df.at[idx, "oxygen_low_compartment_present"] = low_present
        df.at[idx, "oxygen_compartment_breakdown"] = not low_present
        df.at[idx, "oxygen_intrusion_class"] = "baseline"

        delta = paired_change(current, previous_profile) if previous_profile is not None else np.array([])
        common_depths = current.index.intersection(previous_profile.index) if previous_profile is not None else []
        common_depths = sorted(float(x) for x in common_depths if x > split_depth)[-bottom_n:]
        df.at[idx, "oxygen_intrusion_bottom_depths_n"] = int(delta.size)
        df.at[idx, "oxygen_intrusion_bottom_depths_m"] = ",".join(f"{x:g}" for x in common_depths)
        if delta.size:
            change = float(np.median(delta))
            df.at[idx, "oxygen_intrusion_change_median_um"] = change

        if active and (not low_present or baseline_profile is None):
            if pending_end_indices:
                for pending_idx in pending_end_indices:
                    df.at[pending_idx, "oxygen_intrusion_class"] = "baseline"
                    df.at[pending_idx, "oxygen_intrusion_event_id"] = pd.NA
                    df.at[pending_idx, "oxygen_intrusion_end_pending"] = False
            if not low_present:
                end_confirmed_dates[event_id] = pd.to_datetime(row["profile_date"], errors="coerce")
                df.at[idx, "oxygen_intrusion_end_confirmed"] = True
            pending_end_indices = []
            active = False
            baseline_profile = None

        if active and baseline_profile is not None:
            baseline_delta = paired_change(current, baseline_profile)
            baseline_change = (
                float(np.median(baseline_delta)) if baseline_delta.size == bottom_n else np.nan
            )
            df.at[idx, "oxygen_intrusion_baseline_change_median_um"] = baseline_change
            if np.isfinite(baseline_change) and baseline_change >= persistence_threshold and low_present:
                for pending_idx in pending_end_indices:
                    df.at[pending_idx, "oxygen_intrusion_class"] = "intrusion"
                    df.at[pending_idx, "oxygen_intrusion_event_id"] = event_id
                    df.at[pending_idx, "oxygen_intrusion_end_pending"] = False
                pending_end_indices = []
                df.at[idx, "oxygen_intrusion_class"] = "intrusion"
                df.at[idx, "oxygen_intrusion_event_id"] = event_id
            else:
                pending_end_indices.append(idx)
                df.at[idx, "oxygen_intrusion_class"] = "intrusion"
                df.at[idx, "oxygen_intrusion_event_id"] = event_id
                df.at[idx, "oxygen_intrusion_end_pending"] = True
                if len(pending_end_indices) >= end_consecutive:
                    for pending_idx in pending_end_indices:
                        df.at[pending_idx, "oxygen_intrusion_class"] = "baseline"
                        df.at[pending_idx, "oxygen_intrusion_event_id"] = pd.NA
                        df.at[pending_idx, "oxygen_intrusion_end_pending"] = False
                    end_confirmed_dates[event_id] = pd.to_datetime(row["profile_date"], errors="coerce")
                    df.at[idx, "oxygen_intrusion_end_confirmed"] = True
                    pending_end_indices = []
                    active = False
                    baseline_profile = None

        onset = (
            not active
            and low_present
            and delta.size == bottom_n
            and np.isfinite(df.at[idx, "oxygen_intrusion_change_median_um"])
            and df.at[idx, "oxygen_intrusion_change_median_um"] >= onset_threshold
            and previous_profile is not None
        )
        if onset:
            event_id += 1
            active = True
            baseline_profile = previous_profile.copy()
            baseline_delta = paired_change(current, baseline_profile)
            df.at[idx, "oxygen_intrusion_onset"] = True
            df.at[idx, "oxygen_intrusion_class"] = "intrusion"
            df.at[idx, "oxygen_intrusion_event_id"] = event_id
            df.at[idx, "oxygen_intrusion_baseline_change_median_um"] = float(np.median(baseline_delta))

        density_supported = (
            df.at[idx, "oxygen_intrusion_class"] == "intrusion"
            and row.get("deep_intrusion_class") == "intrusion"
        )
        df.at[idx, "oxygen_intrusion_density_supported"] = bool(density_supported)
        previous_profile = current.copy() if not current.empty else previous_profile

    df["oxygen_intrusion_density_supported"] = (
        df["oxygen_intrusion_class"].eq("intrusion")
        & df["deep_intrusion_class"].eq("intrusion")
    )
    for eid in sorted(pd.to_numeric(df["oxygen_intrusion_event_id"], errors="coerce").dropna().unique()):
        members = df[pd.to_numeric(df["oxygen_intrusion_event_id"], errors="coerce").eq(eid)]
        member_dates = pd.to_datetime(members["profile_date"], errors="coerce")
        if member_dates.notna().any():
            df.at[member_dates.idxmax(), "oxygen_intrusion_last_active"] = True
        event_rows.append({
            "oxygen_intrusion_event_id": int(eid),
            "start_date": member_dates.min(),
            "end_date": member_dates.max(),
            "n_cruises": int(len(members)),
            "max_baseline_change_um": float(members["oxygen_intrusion_baseline_change_median_um"].max()),
            "n_density_supported_cruises": int(members["oxygen_intrusion_density_supported"].sum()),
            "end_confirmed_date": end_confirmed_dates.get(int(eid), pd.NaT),
            "low_oxygen_max_um": low_oxygen_max,
            "bottom_n": bottom_n,
            "onset_threshold_um": onset_threshold,
            "persistence_threshold_um": persistence_threshold,
            "end_consecutive_cruises": end_consecutive,
        })
    return df, pd.DataFrame(event_rows)


def _classify_nitrate_renewal(
    summary_df: pd.DataFrame,
    observations: pd.DataFrame,
    profile_cols: list[str],
    depth_col: str,
    nitrate_col: str,
    split_depth: float,
    bottom_n: int,
    min_depths: int,
    detection_limit: float,
    bridge_enabled: bool = False,
    bridge_max_cruises: int = 2,
    bridge_max_days: float = 150.0,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Qualify O2-event onsets and track post-renewal using deep nitrate.

    Oxygen identifies candidate onsets. A candidate becomes a renewal only
    when nitrate is measured at ``min_depths`` of the deepest ``bottom_n``
    sampled depths and the median is above ``detection_limit``. Once a renewal
    begins, nitrate detection alone maintains the post-renewal phase. Negative
    nitrate values are invalid measurements; zero is a measured non-detect.
    """
    df = summary_df.copy()
    output_defaults = {
        "renewal_phase": "baseline",
        "renewal_onset": False,
        "renewal_event_id": pd.NA,
        "renewal_last_active": False,
        "renewal_end_confirmed": False,
        "renewal_phase_inferred": False,
        "renewal_inference_reason": "",
        "renewal_nitrate_status": "insufficient_coverage",
        "renewal_nitrate_supported": False,
        "renewal_nitrate_bottom_depths_n": 0,
        "renewal_nitrate_bottom_depths_m": "",
        "renewal_nitrate_bottom_median_um": np.nan,
        "renewal_nitrate_detection_limit_um": detection_limit,
        "renewal_nitrate_min_depths": min_depths,
        "oxygen_anomaly_class": "none",
    }
    for column, default in output_defaults.items():
        df[column] = default

    if (
        nitrate_col not in observations
        or not np.isfinite(split_depth)
        or "oxygen_intrusion_onset" not in df
    ):
        return df, pd.DataFrame(), pd.DataFrame()

    work = observations[profile_cols + [depth_col, nitrate_col]].copy()
    work[depth_col] = pd.to_numeric(work[depth_col], errors="coerce")
    work[nitrate_col] = pd.to_numeric(work[nitrate_col], errors="coerce")
    # Project QC convention: negative chemistry values are invalid/not
    # measured, whereas zero is a valid measured non-detect.
    work.loc[work[nitrate_col] < 0, nitrate_col] = np.nan
    work = work.dropna(subset=[depth_col])

    nitrate_profiles: dict[tuple, dict[str, object]] = {}
    for keys, group in work.groupby(profile_cols, sort=False):
        key = tuple(keys) if isinstance(keys, tuple) else (keys,)
        sampled_depths = sorted(
            float(value)
            for value in group[depth_col].dropna().unique()
            if float(value) > split_depth
        )[-bottom_n:]
        depth_values = (
            group.loc[group[depth_col].isin(sampled_depths)]
            .groupby(depth_col, sort=True)[nitrate_col]
            .mean()
            .reindex(sampled_depths)
        )
        valid = depth_values.dropna()
        median = float(valid.median()) if len(valid) >= min_depths else np.nan
        if len(valid) < min_depths:
            status = "insufficient_coverage"
        elif median > detection_limit:
            status = "detected"
        else:
            status = "non_detect"
        nitrate_profiles[key] = {
            "sampled_depths": sampled_depths,
            "valid_n": int(len(valid)),
            "median": median,
            "status": status,
        }

    ordered = (
        df.assign(__date=pd.to_datetime(df["profile_date"], errors="coerce"))
        .sort_values("__date")
    )
    active = False
    event_id = 0
    event_end_dates: dict[int, pd.Timestamp] = {}

    for idx, row in ordered.iterrows():
        key = tuple(row[col] for col in profile_cols)
        nitrate = nitrate_profiles.get(
            key,
            {
                "sampled_depths": [],
                "valid_n": 0,
                "median": np.nan,
                "status": "insufficient_coverage",
            },
        )
        status = str(nitrate["status"])
        detected = status == "detected"
        df.at[idx, "renewal_nitrate_status"] = status
        df.at[idx, "renewal_nitrate_supported"] = detected
        df.at[idx, "renewal_nitrate_bottom_depths_n"] = int(nitrate["valid_n"])
        df.at[idx, "renewal_nitrate_bottom_depths_m"] = ",".join(
            f"{depth:g}" for depth in nitrate["sampled_depths"]
        )
        df.at[idx, "renewal_nitrate_bottom_median_um"] = nitrate["median"]

        onset_raw = row.get("oxygen_intrusion_onset", False)
        if isinstance(onset_raw, str):
            oxygen_onset = onset_raw.strip().lower() in {"true", "1", "yes", "y"}
        else:
            oxygen_onset = bool(onset_raw) if pd.notna(onset_raw) else False

        # A new oxygen pulse takes precedence over an existing nitrate tail.
        if oxygen_onset:
            if detected:
                event_id += 1
                active = True
                df.at[idx, "renewal_phase"] = "renewal"
                df.at[idx, "renewal_onset"] = True
                df.at[idx, "renewal_event_id"] = event_id
            else:
                if active and status == "non_detect":
                    event_end_dates[event_id] = pd.to_datetime(
                        row["profile_date"], errors="coerce"
                    )
                    df.at[idx, "renewal_end_confirmed"] = True
                    active = False
                df.at[idx, "renewal_phase"] = (
                    "unknown" if status == "insufficient_coverage" else "baseline"
                )
                df.at[idx, "oxygen_anomaly_class"] = (
                    "nitrate-unresolved"
                    if status == "insufficient_coverage"
                    else "O2-only"
                )
            continue

        if active:
            df.at[idx, "renewal_event_id"] = event_id
            if detected:
                df.at[idx, "renewal_phase"] = "post-renewal"
            elif status == "non_detect":
                df.at[idx, "renewal_phase"] = "baseline"
                df.at[idx, "renewal_end_confirmed"] = True
                event_end_dates[event_id] = pd.to_datetime(
                    row["profile_date"], errors="coerce"
                )
                active = False
            else:
                # Missing nitrate cannot satisfy the post-renewal requirement
                # and cannot demonstrate termination.
                df.at[idx, "renewal_phase"] = "unknown"

    if bridge_enabled:
        chronological_indices = list(ordered.index)
        position = 0
        supported_phases = {"renewal", "post-renewal"}
        while position < len(chronological_indices):
            idx = chronological_indices[position]
            if df.at[idx, "renewal_phase"] != "unknown":
                position += 1
                continue
            block_start = position
            while (
                position < len(chronological_indices)
                and df.at[chronological_indices[position], "renewal_phase"] == "unknown"
            ):
                position += 1
            block_indices = chronological_indices[block_start:position]
            if (
                len(block_indices) > bridge_max_cruises
                or block_start == 0
                or position >= len(chronological_indices)
            ):
                continue
            previous_idx = chronological_indices[block_start - 1]
            next_idx = chronological_indices[position]
            if (
                df.at[previous_idx, "renewal_phase"] not in supported_phases
                or df.at[next_idx, "renewal_phase"] not in supported_phases
            ):
                continue
            previous_event = pd.to_numeric(
                pd.Series([df.at[previous_idx, "renewal_event_id"]]),
                errors="coerce",
            ).iloc[0]
            next_event = pd.to_numeric(
                pd.Series([df.at[next_idx, "renewal_event_id"]]),
                errors="coerce",
            ).iloc[0]
            if (
                not np.isfinite(previous_event)
                or not np.isfinite(next_event)
                or int(previous_event) != int(next_event)
            ):
                continue
            previous_date = pd.to_datetime(
                df.at[previous_idx, "profile_date"], errors="coerce"
            )
            next_date = pd.to_datetime(
                df.at[next_idx, "profile_date"], errors="coerce"
            )
            if (
                pd.isna(previous_date)
                or pd.isna(next_date)
                or (next_date - previous_date).days > bridge_max_days
            ):
                continue
            for bridge_idx in block_indices:
                # Candidate O2 anomalies remain separately unresolved and are
                # never promoted into primary renewal plots.
                if df.at[bridge_idx, "oxygen_anomaly_class"] != "none":
                    continue
                df.at[bridge_idx, "renewal_phase"] = "post-renewal"
                df.at[bridge_idx, "renewal_event_id"] = int(previous_event)
                df.at[bridge_idx, "renewal_phase_inferred"] = True
                df.at[bridge_idx, "renewal_inference_reason"] = (
                    "bracketed_insufficient_nitrate_coverage"
                )

    event_rows: list[dict] = []
    event_values = pd.to_numeric(df["renewal_event_id"], errors="coerce")
    for eid in sorted(event_values.dropna().unique()):
        members = df[
            event_values.eq(eid)
            & df["renewal_phase"].isin(["renewal", "post-renewal"])
        ]
        if members.empty:
            continue
        member_dates = pd.to_datetime(members["profile_date"], errors="coerce")
        if member_dates.notna().any():
            df.at[member_dates.idxmax(), "renewal_last_active"] = True
        event_rows.append({
            "renewal_event_id": int(eid),
            "start_date": member_dates.min(),
            "last_supported_date": member_dates.max(),
            "n_supported_cruises": int(len(members)),
            "n_post_renewal_cruises": int(
                members["renewal_phase"].eq("post-renewal").sum()
            ),
            "n_inferred_post_renewal_cruises": int(
                members["renewal_phase_inferred"].fillna(False).sum()
            ),
            "max_deep_nitrate_median_um": float(
                pd.to_numeric(
                    members["renewal_nitrate_bottom_median_um"], errors="coerce"
                ).max()
            ),
            "end_confirmed_date": event_end_dates.get(int(eid), pd.NaT),
            "bottom_n": bottom_n,
            "min_nitrate_depths": min_depths,
            "nitrate_detection_limit_um": detection_limit,
        })

    anomaly_columns = [
        *profile_cols,
        "profile_date",
        "oxygen_anomaly_class",
        "oxygen_intrusion_change_median_um",
        "oxygen_intrusion_bottom_median_um",
        "renewal_nitrate_status",
        "renewal_nitrate_bottom_depths_n",
        "renewal_nitrate_bottom_depths_m",
        "renewal_nitrate_bottom_median_um",
    ]
    anomaly_df = df.loc[
        df["oxygen_anomaly_class"].ne("none"),
        [column for column in anomaly_columns if column in df],
    ].copy()
    return df, pd.DataFrame(event_rows), anomaly_df


def _pea_class_validation(summary_df: pd.DataFrame) -> pd.DataFrame:
    metrics = [
        "pea_J_m3",
        "n2_mean_s-2",
        "n2_max_s-2",
        "sigma0_upper_lower_diff_kg_m3",
        "mld_depth_m_dr0p03",
        "mld_depth_m_dr0p125",
    ]
    rows = []
    for label in ["mixed", "intermediate", "stratified"]:
        sub = summary_df[summary_df["pea_class"] == label]
        row = {"pea_class": label, "n_profiles": len(sub)}
        for metric in metrics:
            if metric in sub.columns:
                row[f"median_{metric}"] = pd.to_numeric(sub[metric], errors="coerce").median()
        rows.append(row)
    return pd.DataFrame(rows)


def _save_figure_formats(fig: plt.Figure, base_path: Path) -> None:
    for suffix in [".pdf", ".png", ".svg"]:
        fig.savefig(base_path.with_suffix(suffix), dpi=300, bbox_inches="tight")
    plt.close(fig)


def _plot_pea_diagnostics(summary_df: pd.DataFrame, cluster_df: pd.DataFrame, out_dir: Path) -> None:
    values = pd.to_numeric(summary_df["pea_J_m3"], errors="coerce").dropna()
    if values.empty:
        return
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    colors = {"mixed": "tab:blue", "intermediate": "0.55", "stratified": "tab:green"}
    for label, color in colors.items():
        sub = pd.to_numeric(
            summary_df.loc[summary_df["pea_class"] == label, "pea_J_m3"], errors="coerce"
        ).dropna()
        if not sub.empty:
            axes[0].hist(sub, bins=18, alpha=0.55, color=color, label=label)
    for column, style in [("pea_threshold_low", "--"), ("pea_threshold_high", "--")]:
        threshold = pd.to_numeric(summary_df[column], errors="coerce").dropna()
        if not threshold.empty:
            axes[0].axvline(threshold.iloc[0], color="black", linestyle=style)
    axes[0].set_xlabel("PEA (J m$^{-3}$)")
    axes[0].set_ylabel("Profiles")
    axes[0].set_title("Empirical PEA classes")
    axes[0].legend()

    if not cluster_df.empty:
        bar_colors = [colors.get(label, "0.5") for label in cluster_df["pea_class"]]
        axes[1].bar(cluster_df["pea_class"], cluster_df["cluster_center_J_m3"], color=bar_colors)
        for i, row in cluster_df.reset_index(drop=True).iterrows():
            axes[1].text(i, row["cluster_center_J_m3"], f"n={int(row['n_profiles'])}", ha="center", va="bottom")
        axes[1].set_ylabel("Cluster center PEA (J m$^{-3}$)")
        axes[1].set_title("Three-class 1-D k-means")
    else:
        axes[1].axis("off")
        axes[1].text(0.5, 0.5, "Quantile classification", ha="center", va="center")
    fig.tight_layout()
    _save_figure_formats(fig, out_dir / "stratification_pea_classification")


def _plot_physical_regime_selection(selection_df: pd.DataFrame, out_dir: Path) -> None:
    if selection_df.empty:
        return
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    axes[0, 0].plot(selection_df["components"], selection_df["bic"], marker="o", label="BIC")
    axes[0, 0].plot(selection_df["components"], selection_df["aic"], marker="o", label="AIC")
    axes[0, 0].set_title("Information criteria")
    axes[0, 0].legend()
    axes[0, 1].plot(selection_df["components"], selection_df["icl"], marker="o", color="tab:purple")
    axes[0, 1].set_title("ICL (lower is better)")
    axes[1, 0].plot(
        selection_df["components"], selection_df["silhouette_score"], marker="o", color="tab:green"
    )
    axes[1, 0].set_title("Silhouette (higher is better)")
    axes[1, 1].plot(
        selection_df["components"], selection_df["minimum_cluster_fraction"],
        marker="o", label="Minimum cluster fraction",
    )
    axes[1, 1].plot(
        selection_df["components"], selection_df["mean_responsibility_entropy"],
        marker="o", label="Mean entropy",
    )
    axes[1, 1].set_title("Cluster size and uncertainty")
    axes[1, 1].legend()
    for ax in axes.ravel():
        ax.axvline(3, color="black", linestyle="--", alpha=0.5)
        ax.set_xlabel("GMM components")
        ax.grid(axis="y", linestyle="--", alpha=0.3)
    fig.suptitle("Physical-regime GMM selection diagnostics")
    fig.tight_layout()
    _save_figure_formats(fig, out_dir / "stratification_physical_regime_selection")


def _plot_physical_regime_clusters(
    projection_df: pd.DataFrame,
    cluster_df: pd.DataFrame,
    summary_df: pd.DataFrame,
    out_dir: Path,
) -> None:
    colors = {"mixed": "tab:blue", "intermediate": "0.55", "stratified": "tab:green"}
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
    for label, color in colors.items():
        sub = projection_df[projection_df["physical_regime_class"] == label]
        axes[0].scatter(
            sub["physical_pc1"], sub["physical_pc2"], c=color,
            s=35 + 45 * sub["physical_regime_max_probability"], alpha=0.75, label=label,
        )
    axes[0].set_xlabel("Physical PCA 1 (display only)")
    axes[0].set_ylabel("Physical PCA 2 (display only)")
    axes[0].set_title("Cruise physical regimes")
    axes[0].legend()

    center_columns = [f"standardized_center_{feature}" for feature in PHYSICAL_REGIME_FEATURES]
    heat = cluster_df.set_index("physical_regime_class")[center_columns]
    image = axes[1].imshow(heat.to_numpy(), aspect="auto", cmap="coolwarm", vmin=-2, vmax=2)
    axes[1].set_yticks(range(len(heat.index)), heat.index)
    axes[1].set_xticks(range(len(PHYSICAL_REGIME_FEATURES)), PHYSICAL_REGIME_FEATURES, rotation=45, ha="right")
    axes[1].set_title("Standardized component centers")
    fig.colorbar(image, ax=axes[1], label="Standard deviations")

    comparison = pd.crosstab(summary_df["pea_class"], summary_df["physical_regime_class"])
    comparison = comparison.reindex(
        index=["mixed", "intermediate", "stratified"],
        columns=["mixed", "intermediate", "stratified"], fill_value=0,
    )
    image2 = axes[2].imshow(comparison.to_numpy(), cmap="Blues")
    axes[2].set_xticks(range(3), comparison.columns, rotation=30, ha="right")
    axes[2].set_yticks(range(3), comparison.index)
    axes[2].set_xlabel("Multimetric physical regime")
    axes[2].set_ylabel("PEA-only class")
    axes[2].set_title("PEA versus multimetric classes")
    for i in range(3):
        for j in range(3):
            axes[2].text(j, i, int(comparison.iloc[i, j]), ha="center", va="center")
    fig.colorbar(image2, ax=axes[2], label="Cruises")
    fig.tight_layout()
    _save_figure_formats(fig, out_dir / "stratification_physical_regime_clusters")


def _plot_intrusion_diagnostics(summary_df: pd.DataFrame, out_dir: Path) -> None:
    required = {
        "profile_date",
        "sigma0_lower_mean_kg_m3",
        "deep_intrusion_expected_sigma0_kg_m3",
        "deep_intrusion_density_anomaly_kg_m3",
    }
    if not required.issubset(summary_df.columns):
        return
    plot_df = summary_df.copy()
    plot_df["profile_date"] = pd.to_datetime(plot_df["profile_date"], errors="coerce")
    plot_df = plot_df.dropna(subset=["profile_date"]).sort_values("profile_date")
    if plot_df.empty:
        return
    intrusion = plot_df["deep_intrusion_class"].eq("intrusion")
    fig, axes = plt.subplots(2, 1, figsize=(16, 8), sharex=True)
    axes[0].plot(
        plot_df["profile_date"], plot_df["sigma0_lower_mean_kg_m3"],
        color="black", marker="o", label="Observed lower-layer sigma0",
    )
    axes[0].plot(
        plot_df["profile_date"], plot_df["deep_intrusion_expected_sigma0_kg_m3"],
        color="tab:blue", linestyle="--", label="Seasonal expectation",
    )
    axes[0].scatter(
        plot_df.loc[intrusion, "profile_date"],
        plot_df.loc[intrusion, "sigma0_lower_mean_kg_m3"],
        marker="s", s=60, color="tab:red", edgecolor="black", label="Intrusion",
    )
    axes[0].set_ylabel("Lower-layer sigma0 (kg m$^{-3}$)")
    axes[0].legend()
    axes[0].grid(axis="y", linestyle="--", alpha=0.3)

    axes[1].plot(
        plot_df["profile_date"], plot_df["deep_intrusion_density_anomaly_kg_m3"],
        color="tab:purple", marker="o", label="Seasonally adjusted anomaly",
    )
    threshold = pd.to_numeric(plot_df["deep_intrusion_threshold"], errors="coerce").dropna()
    if not threshold.empty:
        axes[1].axhline(threshold.iloc[0], color="tab:red", linestyle="--", label="Upper threshold")
    axes[1].axhline(0, color="0.6", linewidth=1)
    axes[1].set_ylabel("Density anomaly (kg m$^{-3}$)")
    axes[1].set_xlabel("Date")
    axes[1].legend()
    axes[1].grid(axis="y", linestyle="--", alpha=0.3)
    fig.tight_layout()
    _save_figure_formats(fig, out_dir / "stratification_deep_intrusion")


def _plot_oxygen_intrusion_diagnostics(
    summary_df: pd.DataFrame, out_dir: Path
) -> None:
    required = {
        "profile_date",
        "oxygen_intrusion_bottom_median_um",
        "oxygen_intrusion_change_median_um",
        "oxygen_intrusion_class",
    }
    if not required.issubset(summary_df.columns):
        return
    plot_df = summary_df.copy()
    plot_df["profile_date"] = pd.to_datetime(plot_df["profile_date"], errors="coerce")
    plot_df = plot_df.dropna(subset=["profile_date"]).sort_values("profile_date")
    if plot_df.empty:
        return
    intrusion = plot_df["oxygen_intrusion_class"].eq("intrusion")
    onset = plot_df["oxygen_intrusion_onset"].fillna(False).astype(bool)
    fig, axes = plt.subplots(2, 1, figsize=(16, 8), sharex=True)
    axes[0].plot(
        plot_df["profile_date"], plot_df["oxygen_intrusion_bottom_median_um"],
        color="black", marker="o", label="Observed bottom-three median oxygen",
    )
    axes[0].scatter(
        plot_df.loc[intrusion, "profile_date"],
        plot_df.loc[intrusion, "oxygen_intrusion_bottom_median_um"],
        marker="s", s=60, color="tab:red", edgecolor="black",
        label="Oxygen intrusion active",
    )
    axes[0].scatter(
        plot_df.loc[onset, "profile_date"],
        plot_df.loc[onset, "oxygen_intrusion_bottom_median_um"],
        marker="*", s=140, color="gold", edgecolor="black", label="Onset",
    )
    axes[0].set_ylabel("Bottom-three median oxygen (µM)")
    axes[0].legend()
    axes[0].grid(axis="y", linestyle="--", alpha=0.3)

    change = pd.to_numeric(plot_df["oxygen_intrusion_change_median_um"], errors="coerce")
    valid_change = change.notna()
    axes[1].plot(
        plot_df.loc[valid_change, "profile_date"], change[valid_change],
        marker="o", color="tab:green", label="Bottom-three median change",
    )
    threshold = pd.to_numeric(
        plot_df.get("oxygen_intrusion_onset_threshold_um"), errors="coerce"
    ).dropna()
    if not threshold.empty:
        axes[1].axhline(
            threshold.iloc[0], color="tab:red", linestyle="--", label="Onset threshold"
        )
    axes[1].axhline(0, color="0.6", linewidth=1)
    axes[1].set_ylabel("Oxygen change (µM)")
    axes[1].set_xlabel("Date")
    axes[1].legend()
    axes[1].grid(axis="y", linestyle="--", alpha=0.3)
    fig.tight_layout()
    _save_figure_formats(fig, out_dir / "stratification_oxygen_intrusion_events")


def _plot_review(
    summary_df: pd.DataFrame,
    mld_long_df: pd.DataFrame,
    feature_df: pd.DataFrame | None,
    feature_cols: list[str] | None,
    out_path: Path,
    date_col: str | None,
    layer_split_depth: float | None,
) -> None:
    if summary_df.empty:
        print("[i] No summary rows to plot.")
        return

    plot_df = summary_df.copy()
    use_dates = False
    if date_col and "profile_date" in plot_df.columns:
        plot_df["plot_x"] = pd.to_datetime(plot_df["profile_date"], errors="coerce")
        use_dates = plot_df["plot_x"].notna().sum() >= 2
    if not use_dates:
        plot_df["plot_x"] = np.arange(len(plot_df))

    plot_df = plot_df.sort_values("plot_x")
    x_map = dict(zip(plot_df["profile_label"], plot_df["plot_x"]))

    has_pea_class = "pea_class" in plot_df.columns
    has_intrusion = "deep_intrusion_class" in plot_df.columns
    plot_features = feature_cols or []
    if feature_df is None:
        plot_features = []
    else:
        kept = []
        for feature in plot_features:
            upper_col = f"{feature}_upper"
            lower_col = f"{feature}_lower"
            if upper_col in feature_df.columns or lower_col in feature_df.columns:
                kept.append(feature)
        dropped = [f for f in plot_features if f not in kept]
        plot_features = kept
        if dropped:
            print(f"  [i] Dropping features with no upper/lower traces: {dropped}")
    if not plot_features:
        print("  [i] No feature trace panels to plot (empty feature set).")
    nrows = 4 + int(has_intrusion) + int(has_pea_class) + len(plot_features)
    fig, axes = plt.subplots(nrows, 1, figsize=(20, 3.0 * nrows), sharex=True)
    axes = np.atleast_1d(axes)

    def _place_legend(ax: plt.Axes) -> None:
        ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1), borderaxespad=0)

    # MLD + pycnocline panel
    ax = axes[0]
    if not mld_long_df.empty:
        for label, sub in mld_long_df.groupby("threshold_label"):
            sub = sub.copy()
            sub["plot_x"] = sub["profile_label"].map(x_map)
            sub = sub.dropna(subset=["plot_x", "mld_depth_m"])
            if sub.empty:
                continue
            sub = sub.sort_values("plot_x")
            ax.plot(sub["plot_x"], sub["mld_depth_m"], marker="o", label=label)
        ax.invert_yaxis()
        ax.set_ylabel("Depth (m)")
    if np.isfinite(layer_split_depth or np.nan):
        ax.plot(
            plot_df["plot_x"],
            np.full(len(plot_df), layer_split_depth),
            color="gray",
            linestyle=":",
            linewidth=1.5,
            label="layer split depth",
        )
    if "pycnocline_depth_m" in plot_df:
        ax.plot(
            plot_df["plot_x"],
            plot_df["pycnocline_depth_m"],
            color="black",
            linestyle="--",
            marker="o",
            label="pycnocline",
        )
    if ax.get_legend_handles_labels()[1]:
        _place_legend(ax)
    ax.set_title("Mixed Layer Depth + Pycnocline")

    # Adaptive delta panel
    ax = axes[1]
    adaptive = mld_long_df[mld_long_df["threshold_type"] == "adaptive"].copy()
    if not adaptive.empty:
        adaptive["plot_x"] = adaptive["profile_label"].map(x_map)
        adaptive = adaptive.dropna(subset=["plot_x", "delta_rho_kg_m3"])
        adaptive = adaptive.sort_values("plot_x")
        ax.plot(adaptive["plot_x"], adaptive["delta_rho_kg_m3"], marker="o", color="tab:orange")
    ax.set_ylabel("Adaptive delta_rho (kg/m3)")
    ax.set_title("Adaptive threshold used")

    # N2 panel
    ax = axes[2]
    ax.plot(plot_df["plot_x"], plot_df["n2_mean_s-2"], marker="o", label="N2 mean")
    if "n2_mean_upper_s-2" in plot_df:
        ax.plot(plot_df["plot_x"], plot_df["n2_mean_upper_s-2"], marker="o", label="N2 mean upper")
    if "n2_mean_lower_s-2" in plot_df:
        ax.plot(plot_df["plot_x"], plot_df["n2_mean_lower_s-2"], marker="o", label="N2 mean lower")
    ax.plot(plot_df["plot_x"], plot_df["n2_max_s-2"], marker="o", label="N2 max")
    ax.set_ylabel("N2 (s^-2)")
    ax.set_title("Brunt-Vaisala frequency")
    _place_legend(ax)

    # PEA panel
    ax = axes[3]
    ax.plot(plot_df["plot_x"], plot_df["pea_J_m3"], marker="o", color="tab:green", label="PEA total")
    if "pea_upper_J_m3" in plot_df:
        ax.plot(plot_df["plot_x"], plot_df["pea_upper_J_m3"], marker="o", label="PEA upper")
    if "pea_lower_J_m3" in plot_df:
        ax.plot(plot_df["plot_x"], plot_df["pea_lower_J_m3"], marker="o", label="PEA lower")
    ax.set_ylabel("PEA (J/m3)")
    ax.set_title("Potential Energy Anomaly")
    _place_legend(ax)

    panel_idx = 4
    if has_intrusion:
        ax = axes[panel_idx]
        colors = {
            "intrusion": "tab:red",
            "baseline": "tab:blue",
            "unknown": "tab:gray",
        }
        for cls in ["intrusion", "baseline", "unknown"]:
            sub = plot_df[plot_df["deep_intrusion_class"] == cls]
            if sub.empty:
                continue
            ax.scatter(
                sub["plot_x"],
                sub["sigma0_lower_mean_kg_m3"],
                label=cls,
                color=colors[cls],
                s=30,
            )
        if "deep_intrusion_threshold" in plot_df:
            ax.plot(
                plot_df["plot_x"],
                plot_df["deep_intrusion_threshold"],
                linestyle="--",
                color="black",
                label="intrusion threshold",
            )
        ax.set_ylabel("Sigma0 lower mean (kg/m3)")
        ax.set_title("Deep intrusion indicator")
        _place_legend(ax)
        panel_idx += 1

    # PEA classification panel
    if has_pea_class:
        ax = axes[panel_idx]
        colors = {
            "mixed": "tab:blue",
            "intermediate": "tab:gray",
            "stratified": "tab:red",
            "unknown": "tab:purple",
        }
        for cls in ["mixed", "intermediate", "stratified", "unknown"]:
            sub = plot_df[plot_df["pea_class"] == cls]
            if sub.empty:
                continue
            ax.scatter(sub["plot_x"], sub["pea_J_m3"], label=cls, color=colors[cls], s=30)
        if "pea_threshold_low" in plot_df and "pea_threshold_high" in plot_df:
            ax.plot(plot_df["plot_x"], plot_df["pea_threshold_low"], linestyle="--", color="black", label="PEA low")
            ax.plot(plot_df["plot_x"], plot_df["pea_threshold_high"], linestyle="--", color="black", label="PEA high")
        ax.set_ylabel("PEA (J/m3)")
        ax.set_title("PEA classification")
        _place_legend(ax)
        panel_idx += 1

    if plot_features:
        feature_plot = feature_df.copy()
        feature_plot["plot_x"] = feature_plot["profile_label"].map(x_map)
        for feature in plot_features:
            ax = axes[panel_idx]
            upper_col = f"{feature}_upper"
            lower_col = f"{feature}_lower"
            upper_n = 0
            lower_n = 0
            if upper_col in feature_plot:
                sub = feature_plot[["plot_x", upper_col]].dropna()
                if not sub.empty:
                    sub = sub.sort_values("plot_x")
                    ax.scatter(sub["plot_x"], sub[upper_col], color="tab:blue", s=25, label="upper")
                    upper_n = len(sub)
            if lower_col in feature_plot:
                sub = feature_plot[["plot_x", lower_col]].dropna()
                if not sub.empty:
                    sub = sub.sort_values("plot_x")
                    ax.scatter(sub["plot_x"], sub[lower_col], color="tab:orange", s=25, label="lower")
                    lower_n = len(sub)
            if upper_n == 0 and lower_n == 0:
                print(f"  [i] Feature '{feature}' has no upper/lower points to plot.")
            ax.set_ylabel(feature)
            ax.set_title(feature)
            if ax.get_legend_handles_labels()[1]:
                _place_legend(ax)
            panel_idx += 1

    if use_dates:
        axes[-1].set_xlabel("Date")
    else:
        axes[-1].set_xlabel("Profile index")

    fig.tight_layout(rect=[0, 0, 0.8, 1])
    _save_figure_formats(fig, Path(out_path).with_suffix(""))
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Compute stratification metrics (density, N^2, MLD, PEA) from CTD data."
    )
    ap.add_argument("--input", required=True, help="Input TSV file.")
    ap.add_argument("--output-dir", type=Path, required=True, help="Output directory.")
    ap.add_argument("--sep", default="\t", help="Input/output separator (default: tab).")

    ap.add_argument("--salinity-col", required=True, help="Salinity column (SP).")
    ap.add_argument("--temperature-col", required=True, help="Temperature column (deg C).")
    ap.add_argument("--depth-col", required=True, help="Depth column (m, positive down).")
    ap.add_argument("--pressure-col", default=None, help="Pressure column (dbar).")

    ap.add_argument("--latitude", type=float, default=None, help="Latitude (deg).")
    ap.add_argument("--latitude-col", default=None, help="Latitude column (deg).")
    ap.add_argument("--longitude", type=float, default=None, help="Longitude (deg).")
    ap.add_argument("--longitude-col", default=None, help="Longitude column (deg).")

    ap.add_argument(
        "--profile-cols",
        default=None,
        help="Comma-separated columns defining each profile (default: treat all rows as one profile).",
    )

    ap.add_argument("--depth-min", type=float, default=None, help="Minimum depth to include.")
    ap.add_argument("--depth-max", type=float, default=None, help="Maximum depth to include.")

    ap.add_argument(
        "--mld-delta-rho",
        default="0.03,0.125",
        help="Comma-separated MLD density thresholds (kg/m^3, default 0.03,0.125).",
    )
    ap.add_argument(
        "--mld-reference-depth",
        type=float,
        default=10.0,
        help="Reference depth for MLD (m, default 10).",
    )
    ap.add_argument(
        "--adaptive-layer-max-depth",
        type=float,
        default=10.0,
        help="Surface layer depth for adaptive MLD threshold (m, default 10).",
    )
    ap.add_argument(
        "--adaptive-percentile",
        type=float,
        default=90.0,
        help="Percentile for adaptive delta_rho (default 90).",
    )
    ap.add_argument(
        "--layer-split-mode",
        choices=["global_pycnocline", "mld125", "explicit"],
        default="global_pycnocline",
        help="Layer split depth mode (default global_pycnocline).",
    )
    ap.add_argument(
        "--layer-split-stat",
        choices=["max", "p90", "median"],
        default="max",
        help="Statistic for mld125 layer split depth (default max).",
    )
    ap.add_argument(
        "--layer-split-depth",
        type=float,
        default=None,
        help="Explicit depth for layer split (m). Use with --layer-split-mode explicit.",
    )
    ap.add_argument("--date-col", default=None, help="Optional date column for plotting.")
    ap.add_argument(
        "--pea-class-mode",
        choices=["global", "none"],
        default="global",
        help="PEA classification mode (default global).",
    )
    ap.add_argument(
        "--pea-bootstrap-iterations",
        type=int,
        default=500,
        help="Profile bootstrap iterations for k-means threshold intervals (default 500).",
    )
    ap.add_argument("--pea-random-state", type=int, default=42, help="PEA k-means random seed (default 42).")
    ap.add_argument(
        "--physical-regime-k-max", type=int, default=6,
        help="Largest GMM component count evaluated for physical-regime diagnostics (default 6).",
    )
    ap.add_argument(
        "--plot-feature-cols",
        default=None,
        help="Comma-separated feature columns to plot as traces (default: infer non-metadata numeric columns).",
    )
    ap.add_argument(
        "--deep-intrusion-mode",
        choices=["global", "none"],
        default="global",
        help="Deep intrusion classification mode (default global).",
    )
    ap.add_argument(
        "--deep-intrusion-quantile",
        type=float,
        default=0.90,
        help="Upper quantile of seasonally adjusted lower-density residuals (default 0.90).",
    )
    ap.add_argument(
        "--deep-intrusion-oxygen-col",
        default="Oxygen",
        help="Oxygen column used for matched-depth intrusion detection (default Oxygen).",
    )
    ap.add_argument(
        "--oxygen-low-compartment-max",
        type=float,
        default=90.0,
        help="Largest O2 value considered part of a low-O2 compartment, in µM (default 90).",
    )
    ap.add_argument(
        "--oxygen-intrusion-bottom-n",
        type=int,
        default=3,
        help="Number of deepest matched samples used for O2 intrusion detection (default 3).",
    )
    ap.add_argument(
        "--oxygen-intrusion-onset-threshold", type=float, default=4.5,
        help="Required bottom-layer median O2 increase at onset, in µM (default 4.5).",
    )
    ap.add_argument(
        "--oxygen-intrusion-persistence-threshold", type=float, default=4.5,
        help="Required O2 elevation above the pre-event baseline to remain active (default 4.5 µM).",
    )
    ap.add_argument(
        "--oxygen-intrusion-end-consecutive", type=int, default=2,
        help="Consecutive below-persistence cruises required to confirm event end (default 2).",
    )
    ap.add_argument(
        "--renewal-nitrate-col",
        default="Nitrate",
        help="Nitrate column used to qualify renewal and post-renewal (default Nitrate).",
    )
    ap.add_argument(
        "--renewal-nitrate-min-depths",
        type=int,
        default=2,
        help="Minimum measured depths among the nitrate qualification depths (default 2).",
    )
    ap.add_argument(
        "--renewal-nitrate-bottom-n",
        type=int,
        default=3,
        help="Number of deepest sampled depths considered for nitrate qualification (default 3).",
    )
    ap.add_argument(
        "--renewal-nitrate-detection-limit",
        type=float,
        default=0.0,
        help="Deep median nitrate must exceed this measured non-detect limit (default 0 µM).",
    )
    ap.add_argument(
        "--renewal-bridge-enabled",
        action="store_true",
        help="Infer post-renewal across short, bracketed insufficient-nitrate gaps.",
    )
    ap.add_argument(
        "--renewal-bridge-max-cruises",
        type=int,
        default=2,
        help="Largest consecutive unresolved block eligible for continuity inference (default 2).",
    )
    ap.add_argument(
        "--renewal-bridge-max-days",
        type=float,
        default=150.0,
        help="Largest elapsed time between supported flanks eligible for bridging (default 150 days).",
    )
    ap.add_argument(
        "--no-plots",
        action="store_true",
        help="Skip plot generation.",
    )
    return ap.parse_args()


def main() -> None:
    args = parse_args()

    if args.pea_bootstrap_iterations < 0:
        raise ValueError("--pea-bootstrap-iterations must be nonnegative.")
    if args.physical_regime_k_max < 3:
        raise ValueError("--physical-regime-k-max must be at least 3.")
    if not 0 < args.deep_intrusion_quantile < 1:
        raise ValueError("--deep-intrusion-quantile must be between 0 and 1.")
    if args.oxygen_intrusion_bottom_n < 1:
        raise ValueError("--oxygen-intrusion-bottom-n must be positive.")
    if args.oxygen_intrusion_onset_threshold < 0 or args.oxygen_intrusion_persistence_threshold < 0:
        raise ValueError("Oxygen intrusion thresholds must be nonnegative.")
    if args.oxygen_intrusion_end_consecutive < 1:
        raise ValueError("--oxygen-intrusion-end-consecutive must be positive.")
    if args.renewal_nitrate_bottom_n < 1:
        raise ValueError("--renewal-nitrate-bottom-n must be positive.")
    if not 1 <= args.renewal_nitrate_min_depths <= args.renewal_nitrate_bottom_n:
        raise ValueError(
            "--renewal-nitrate-min-depths must be between 1 and "
            "--renewal-nitrate-bottom-n."
        )
    if args.renewal_nitrate_detection_limit < 0:
        raise ValueError("--renewal-nitrate-detection-limit must be nonnegative.")
    if args.renewal_bridge_max_cruises < 1:
        raise ValueError("--renewal-bridge-max-cruises must be positive.")
    if args.renewal_bridge_max_days <= 0:
        raise ValueError("--renewal-bridge-max-days must be positive.")

    out_dir = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    for obsolete_name in [
        "stratification_oxygen_supported_deep_intrusion_model.tsv",
        "stratification_oxygen_supported_deep_intrusion.pdf",
        "stratification_oxygen_supported_deep_intrusion.png",
        "stratification_oxygen_supported_deep_intrusion.svg",
    ]:
        obsolete_path = out_dir / obsolete_name
        if obsolete_path.exists() or obsolete_path.is_symlink():
            obsolete_path.unlink()

    print("\n" + "=" * 70)
    print("STRATIFICATION METRICS (TEOS-10)")
    print("=" * 70)

    print("\n[1/9] Loading input table...")
    df = pd.read_csv(args.input, sep=args.sep)

    required = [
        args.salinity_col,
        args.temperature_col,
        args.depth_col,
    ]
    if args.pressure_col:
        required.append(args.pressure_col)
    if args.latitude_col:
        required.append(args.latitude_col)
    if args.longitude_col:
        required.append(args.longitude_col)
    if args.date_col:
        required.append(args.date_col)
    _require_cols(df, required)

    profile_cols = _split_csv(args.profile_cols)
    if not profile_cols:
        df = df.copy()
        df["profile_id"] = "all"
        profile_cols = ["profile_id"]
    _require_cols(df, profile_cols)

    if args.depth_min is not None:
        df = df[pd.to_numeric(df[args.depth_col], errors="coerce") >= args.depth_min]
    if args.depth_max is not None:
        df = df[pd.to_numeric(df[args.depth_col], errors="coerce") <= args.depth_max]

    print(f"  Profiles: {df[profile_cols].drop_duplicates().shape[0]}")

    exclude_cols = {
        args.salinity_col,
        args.temperature_col,
        args.depth_col,
        args.pressure_col,
        args.latitude_col,
        args.longitude_col,
        args.date_col,
    }
    exclude_cols.update(profile_cols)
    exclude_cols = {c for c in exclude_cols if c}

    if args.plot_feature_cols:
        plot_feature_cols = _split_csv(args.plot_feature_cols)
        missing = [c for c in plot_feature_cols if c not in df.columns]
        if missing:
            raise ValueError(f"--plot-feature-cols includes columns not found in input: {missing}")
    else:
        plot_feature_cols = _infer_plot_feature_cols(df, exclude_cols)
    plot_feature_cols = [c for c in plot_feature_cols if c not in profile_cols]
    if plot_feature_cols:
        print(f"  [i] Plotting {len(plot_feature_cols)} feature traces")
    else:
        print("  [i] No feature traces selected for plotting")

    if args.layer_split_mode == "explicit":
        if args.layer_split_depth is None:
            raise ValueError("--layer-split-depth is required when --layer-split-mode explicit.")
        layer_split_depth = float(args.layer_split_depth)
        layer_split_method = "explicit"
        print(f"\n[2/9] Using explicit layer split depth: {layer_split_depth:.2f} m")
    else:
        if args.layer_split_depth is not None:
            raise ValueError("Use --layer-split-mode explicit when providing --layer-split-depth.")
        if args.layer_split_mode == "global_pycnocline":
            print("\n[2/9] Computing global pycnocline split depth...")
            layer_split_depth = _compute_global_pycnocline_depth(
                df=df,
                profile_cols=profile_cols,
                sal_col=args.salinity_col,
                temp_col=args.temperature_col,
                depth_col=args.depth_col,
                pressure_col=args.pressure_col,
                latitude_col=args.latitude_col,
                longitude_col=args.longitude_col,
                latitude=args.latitude,
                longitude=args.longitude,
            )
            if np.isfinite(layer_split_depth):
                layer_split_method = "global_pycnocline_max"
                print(f"  [i] Using global pycnocline max depth for layer split: {layer_split_depth:.2f} m")
            else:
                layer_split_method = "none"
                print("  [i] No valid pycnocline depths; layer metrics will be NaN.")
        elif args.layer_split_mode == "mld125":
            print("\n[2/9] Computing global mld125 split depth...")
            layer_split_depth = _compute_global_mld_depth(
                df=df,
                profile_cols=profile_cols,
                sal_col=args.salinity_col,
                temp_col=args.temperature_col,
                depth_col=args.depth_col,
                pressure_col=args.pressure_col,
                latitude_col=args.latitude_col,
                longitude_col=args.longitude_col,
                latitude=args.latitude,
                longitude=args.longitude,
                ref_depth=args.mld_reference_depth,
                delta_rho=0.125,
                stat=args.layer_split_stat,
            )
            if np.isfinite(layer_split_depth):
                layer_split_method = f"mld125_{args.layer_split_stat}"
                print(
                    f"  [i] Using global mld125 ({args.layer_split_stat}) depth for layer split: "
                    f"{layer_split_depth:.2f} m"
                )
            else:
                layer_split_method = "none"
                print("  [i] No valid mld125 depths; layer metrics will be NaN.")
        else:
            raise ValueError(f"Unknown --layer-split-mode '{args.layer_split_mode}'.")

    density_rows = []
    n2_rows = []
    summary_rows = []
    mld_rows = []
    warned_lat_zero = False
    fixed_thresholds = _split_float_csv(args.mld_delta_rho)
    fixed_thresholds = [t for t in fixed_thresholds if t > 0]
    if not fixed_thresholds:
        raise ValueError("--mld-delta-rho must include at least one positive threshold.")

    print("\n[3/9] Computing profiles...")
    for keys, group in df.groupby(profile_cols):
        if not isinstance(keys, tuple):
            keys = (keys,)
        key_map = dict(zip(profile_cols, keys))
        profile_label = _profile_label(key_map)

        g = group.copy()
        g["__sal"] = _numeric(g[args.salinity_col])
        g["__temp"] = _numeric(g[args.temperature_col])
        g["__depth"] = _numeric(g[args.depth_col])
        g["__press"] = _numeric(g[args.pressure_col]) if args.pressure_col else np.nan
        lat_vec = _resolve_vector(g, args.latitude_col, args.latitude, "latitude", len(g))
        lon_vec = _resolve_vector(g, args.longitude_col, args.longitude, "longitude", len(g))
        g["__lat"] = lat_vec if lat_vec is not None else np.nan
        g["__lon"] = lon_vec if lon_vec is not None else np.nan

        req_mask = np.isfinite(g["__sal"]) & np.isfinite(g["__temp"]) & np.isfinite(g["__depth"])
        g = g.loc[req_mask]
        if g.empty:
            continue

        depth_df = _aggregate_by_depth(
            g,
            depth_col=args.depth_col,
            cols=["__sal", "__temp", "__press", "__lat", "__lon"],
        )

        depth = _numeric(depth_df[args.depth_col])
        sal = _numeric(depth_df["__sal"])
        temp = _numeric(depth_df["__temp"])
        lat = _numeric(depth_df["__lat"]) if "__lat" in depth_df else None
        lon = _numeric(depth_df["__lon"]) if "__lon" in depth_df else None
        if lat is not None and np.isnan(lat).all():
            lat = None
        if lon is not None and np.isnan(lon).all():
            lon = None
        pressure = _numeric(depth_df["__press"]) if args.pressure_col else None

        if lat is None:
            if args.pressure_col:
                lat = np.zeros(len(depth))
                if not warned_lat_zero:
                    print("[i] Latitude not provided; using 0 deg for SA conversion.")
                    warned_lat_zero = True
            else:
                raise ValueError("Latitude is required when pressure is not provided.")
        if lon is None:
            lon = np.zeros(len(depth))

        pressure = _compute_pressure(depth, pressure, lat)

        SA = gsw.SA_from_SP(sal, pressure, lon, lat)
        CT = gsw.CT_from_t(SA, temp, pressure)
        rho = gsw.rho(SA, CT, pressure)
        sigma0 = gsw.sigma0(SA, CT)

        profile_date = None
        if args.date_col:
            date_vals = g[args.date_col].dropna()
            if not date_vals.empty:
                profile_date = date_vals.iloc[0]

        n2 = np.array([])
        p_mid = np.array([])
        depth_mid = np.array([])
        if len(depth) >= 2:
            n2, p_mid = gsw.Nsquared(SA, CT, pressure, lat)
            lat_mid = float(np.nanmean(lat)) if np.isfinite(np.nanmean(lat)) else 0.0
            depth_mid = -gsw.z_from_p(p_mid, lat_mid)

        mld_summary = {}
        mld_rho_ref = np.nan
        mld_entries = []
        for delta_rho in fixed_thresholds:
            mld_depth, rho_ref, delta_below = _compute_mld(
                depth=depth,
                sigma0=sigma0,
                ref_depth=args.mld_reference_depth,
                delta_rho=delta_rho,
            )
            if np.isnan(mld_rho_ref):
                mld_rho_ref = rho_ref
            tag = _format_threshold_tag(delta_rho)
            mld_summary[f"mld_depth_m_dr{tag}"] = mld_depth
            mld_summary[f"mld_delta_rho_below_kg_m3_dr{tag}"] = delta_below
            mld_entries.append(
                {
                    **key_map,
                    "profile_label": profile_label,
                    "profile_date": profile_date,
                    "threshold_type": "fixed",
                    "threshold_label": f"dr{tag}",
                    "delta_rho_kg_m3": float(delta_rho),
                    "mld_depth_m": mld_depth,
                    "mld_rho_ref_kg_m3": rho_ref,
                    "mld_delta_rho_below_kg_m3": delta_below,
                }
            )

        adaptive_delta = _adaptive_delta_rho(
            depth=depth,
            sigma0=sigma0,
            layer_max_depth=args.adaptive_layer_max_depth,
            percentile=args.adaptive_percentile,
        )
        if np.isfinite(adaptive_delta):
            mld_depth, rho_ref, delta_below = _compute_mld(
                depth=depth,
                sigma0=sigma0,
                ref_depth=args.mld_reference_depth,
                delta_rho=adaptive_delta,
            )
        else:
            mld_depth, rho_ref, delta_below = (np.nan, np.nan, np.nan)
        mld_entries.append(
            {
                **key_map,
                "profile_label": profile_label,
                "profile_date": profile_date,
                "threshold_type": "adaptive",
                "threshold_label": "adaptive",
                "delta_rho_kg_m3": adaptive_delta,
                "mld_depth_m": mld_depth,
                "mld_rho_ref_kg_m3": rho_ref,
                "mld_delta_rho_below_kg_m3": delta_below,
            }
        )
        lat_mean = float(np.nanmean(lat)) if np.isfinite(np.nanmean(lat)) else np.nan
        pea = _compute_pea(depth=depth, sigma0=sigma0, lat=lat_mean if np.isfinite(lat_mean) else None)
        n2_mean = float(np.nanmean(n2)) if n2.size else np.nan
        n2_max = float(np.nanmax(n2)) if n2.size else np.nan
        if n2.size:
            idx_max = int(np.nanargmax(n2))
            n2_max_depth = float(depth_mid[idx_max]) if depth_mid.size else np.nan
        else:
            n2_max_depth = np.nan
        pycnocline_depth = n2_max_depth

        if n2.size:
            if np.isfinite(layer_split_depth):
                n2_layer = np.where(depth_mid <= layer_split_depth, "upper", "lower")
            else:
                n2_layer = np.array(["unknown"] * len(n2))
            for i in range(len(n2)):
                row = {
                    **key_map,
                    "profile_label": profile_label,
                    "profile_date": profile_date,
                    "pressure_mid_dbar": p_mid[i],
                    "depth_mid_m": depth_mid[i],
                    "N2_s-2": n2[i],
                    "pycnocline_depth_m": pycnocline_depth,
                    "layer_split_depth_m": layer_split_depth,
                    "layer_split_method": layer_split_method,
                    "layer": n2_layer[i],
                }
                n2_rows.append(row)

        if np.isfinite(layer_split_depth):
            density_layer = np.where(depth <= layer_split_depth, "upper", "lower")
        else:
            density_layer = np.array(["unknown"] * len(depth))

        for i in range(len(depth)):
            row = {
                **key_map,
                "profile_label": profile_label,
                "profile_date": profile_date,
                args.depth_col: depth[i],
                "pycnocline_depth_m": pycnocline_depth,
                "layer_split_depth_m": layer_split_depth,
                "layer_split_method": layer_split_method,
                "layer": density_layer[i],
                "pressure_dbar": pressure[i],
                "salinity_sp": sal[i],
                "temperature_c": temp[i],
                "SA": SA[i],
                "CT": CT[i],
                "rho_kg_m3": rho[i],
                "sigma0_kg_m3": sigma0[i],
                "n_samples_at_depth": int(depth_df["__n_samples"].iloc[i]),
            }
            density_rows.append(row)

        upper_mask_depth, lower_mask_depth = _split_layers(depth, layer_split_depth)
        upper_mask_n2, lower_mask_n2 = _split_layers(depth_mid, layer_split_depth)

        sigma0_upper_mean = float(np.nanmean(sigma0[upper_mask_depth])) if upper_mask_depth.any() else np.nan
        sigma0_lower_mean = float(np.nanmean(sigma0[lower_mask_depth])) if lower_mask_depth.any() else np.nan
        if np.isfinite(sigma0_upper_mean) and np.isfinite(sigma0_lower_mean):
            sigma0_upper_lower_diff = sigma0_lower_mean - sigma0_upper_mean
        else:
            sigma0_upper_lower_diff = np.nan

        pea_upper = (
            _compute_pea(depth[upper_mask_depth], sigma0[upper_mask_depth], lat_mean)
            if upper_mask_depth.any()
            else np.nan
        )
        pea_lower = (
            _compute_pea(depth[lower_mask_depth], sigma0[lower_mask_depth], lat_mean)
            if lower_mask_depth.any()
            else np.nan
        )

        n2_mean_upper = float(np.nanmean(n2[upper_mask_n2])) if upper_mask_n2.any() else np.nan
        n2_mean_lower = float(np.nanmean(n2[lower_mask_n2])) if lower_mask_n2.any() else np.nan
        n2_max_upper = float(np.nanmax(n2[upper_mask_n2])) if upper_mask_n2.any() else np.nan
        n2_max_lower = float(np.nanmax(n2[lower_mask_n2])) if lower_mask_n2.any() else np.nan

        mld_extra = {
            "pycnocline_depth_m": pycnocline_depth,
            "layer_split_depth_m": layer_split_depth,
            "layer_split_method": layer_split_method,
            "sigma0_upper_mean_kg_m3": sigma0_upper_mean,
            "sigma0_lower_mean_kg_m3": sigma0_lower_mean,
            "sigma0_upper_lower_diff_kg_m3": sigma0_upper_lower_diff,
            "pea_J_m3": pea,
            "pea_upper_J_m3": pea_upper,
            "pea_lower_J_m3": pea_lower,
            "n2_mean_s-2": n2_mean,
            "n2_max_s-2": n2_max,
            "n2_mean_upper_s-2": n2_mean_upper,
            "n2_mean_lower_s-2": n2_mean_lower,
            "n2_max_upper_s-2": n2_max_upper,
            "n2_max_lower_s-2": n2_max_lower,
        }
        for entry in mld_entries:
            entry.update(mld_extra)
        mld_rows.extend(mld_entries)

        summary = {
            **key_map,
            "profile_label": profile_label,
            "profile_date": profile_date,
            "n_depths": int(len(depth)),
            "depth_min_m": float(np.nanmin(depth)),
            "depth_max_m": float(np.nanmax(depth)),
            "pycnocline_depth_m": pycnocline_depth,
            "layer_split_depth_m": layer_split_depth,
            "layer_split_method": layer_split_method,
            "sigma0_upper_mean_kg_m3": sigma0_upper_mean,
            "sigma0_lower_mean_kg_m3": sigma0_lower_mean,
            "sigma0_upper_lower_diff_kg_m3": sigma0_upper_lower_diff,
            "mld_reference_depth_m": float(args.mld_reference_depth),
            "mld_delta_rho_fixed_list_kg_m3": ",".join([str(t) for t in fixed_thresholds]),
            "mld_delta_rho_adaptive_kg_m3": adaptive_delta,
            "mld_rho_ref_kg_m3": mld_rho_ref,
            "adaptive_layer_max_depth_m": float(args.adaptive_layer_max_depth),
            "adaptive_percentile": float(args.adaptive_percentile),
            "pea_J_m3": pea,
            "pea_upper_J_m3": pea_upper,
            "pea_lower_J_m3": pea_lower,
            "n2_mean_s-2": n2_mean,
            "n2_max_s-2": n2_max,
            "n2_max_depth_m": n2_max_depth,
            "n2_mean_upper_s-2": n2_mean_upper,
            "n2_mean_lower_s-2": n2_mean_lower,
            "n2_max_upper_s-2": n2_max_upper,
            "n2_max_lower_s-2": n2_max_lower,
        }
        summary.update(mld_summary)
        summary_rows.append(summary)

    density_df = pd.DataFrame(density_rows)
    n2_df = pd.DataFrame(n2_rows)
    summary_df = pd.DataFrame(summary_rows)
    mld_long_df = pd.DataFrame(mld_rows)

    oxygen_lower_col = "deep_intrusion_oxygen_lower_mean"
    if args.deep_intrusion_oxygen_col in df.columns and np.isfinite(layer_split_depth):
        oxygen_work = df[profile_cols + [args.depth_col, args.deep_intrusion_oxygen_col]].copy()
        oxygen_work[args.depth_col] = pd.to_numeric(oxygen_work[args.depth_col], errors="coerce")
        oxygen_work[args.deep_intrusion_oxygen_col] = pd.to_numeric(
            oxygen_work[args.deep_intrusion_oxygen_col], errors="coerce"
        )
        oxygen_lower = (
            oxygen_work.loc[oxygen_work[args.depth_col] > layer_split_depth]
            .groupby(profile_cols, sort=False)[args.deep_intrusion_oxygen_col]
            .mean()
            .rename(oxygen_lower_col)
            .reset_index()
        )
        summary_df = summary_df.merge(oxygen_lower, on=profile_cols, how="left")

    print("\n[4/9] Preparing feature traces...")
    feature_trace_df = None
    if plot_feature_cols:
        if not np.isfinite(layer_split_depth):
            print("  [i] Skipping feature traces (no valid layer split depth).")
        else:
            feature_trace_df = df[profile_cols + [args.depth_col] + plot_feature_cols].copy()
            for col in plot_feature_cols + [args.depth_col]:
                feature_trace_df[col] = pd.to_numeric(feature_trace_df[col], errors="coerce")
            depth_vals = feature_trace_df[args.depth_col]
            upper_df = feature_trace_df[depth_vals <= layer_split_depth]
            lower_df = feature_trace_df[depth_vals > layer_split_depth]
            print(
                "  [i] Feature trace split counts: "
                f"upper_rows={len(upper_df)} lower_rows={len(lower_df)}"
            )
            upper_means = upper_df.groupby(profile_cols, sort=False)[plot_feature_cols].mean()
            lower_means = lower_df.groupby(profile_cols, sort=False)[plot_feature_cols].mean()
            feature_trace_df = upper_means.add_suffix("_upper").join(
                lower_means.add_suffix("_lower"), how="outer"
            ).reset_index()
            feature_trace_df["profile_label"] = (
                feature_trace_df[profile_cols].astype(str).agg("|".join, axis=1)
            )
            if feature_trace_df.empty:
                print("  [i] Feature trace table is empty after aggregation.")
            else:
                nonnull_counts = feature_trace_df.notna().sum().to_dict()
                total_profiles = len(feature_trace_df)
                print(f"  [i] Feature trace table rows: {total_profiles}")
                for col in plot_feature_cols[:5]:
                    upper_col = f"{col}_upper"
                    lower_col = f"{col}_lower"
                    upper_n = nonnull_counts.get(upper_col, 0)
                    lower_n = nonnull_counts.get(lower_col, 0)
                    print(f"  [i] {col}: upper_nonnull={upper_n} lower_nonnull={lower_n}")

    print("\n[5/9] Classifying PEA...")
    if args.pea_class_mode != "none":
        summary_df, pea_cluster_df, pea_bootstrap_df = _classify_pea(
            summary_df=summary_df,
            mode=args.pea_class_mode,
            bootstrap_iterations=args.pea_bootstrap_iterations,
            random_state=args.pea_random_state,
        )
        lower_input = summary_df.copy()
        lower_input["pea_J_m3"] = lower_input["pea_lower_J_m3"]
        lower_classified, pea_lower_cluster_df, pea_lower_bootstrap_df = _classify_pea(
            summary_df=lower_input,
            mode=args.pea_class_mode,
            bootstrap_iterations=args.pea_bootstrap_iterations,
            random_state=args.pea_random_state,
        )
        lower_column_map = {
            "pea_class": "pea_lower_class",
            "pea_class_method": "pea_lower_class_method",
            "pea_threshold_low": "pea_lower_threshold_low",
            "pea_threshold_high": "pea_lower_threshold_high",
            "pea_threshold_low_ci_lower": "pea_lower_threshold_low_ci_lower",
            "pea_threshold_low_ci_upper": "pea_lower_threshold_low_ci_upper",
            "pea_threshold_high_ci_lower": "pea_lower_threshold_high_ci_lower",
            "pea_threshold_high_ci_upper": "pea_lower_threshold_high_ci_upper",
        }
        for source, target in lower_column_map.items():
            summary_df[target] = lower_classified[source].to_numpy()
        if not mld_long_df.empty:
            mld_long_df = mld_long_df.merge(
                summary_df[
                    [
                        "profile_label",
                        "pea_class",
                        "pea_threshold_low",
                        "pea_threshold_high",
                        "pea_lower_class",
                    ]
                ],
                on="profile_label",
                how="left",
            )

    print("\n[5b/9] Classifying multimetric physical regimes...")
    summary_df, physical_selection_df, physical_clusters_df, physical_projection_df = (
        _classify_physical_regimes(
            summary_df=summary_df,
            random_state=args.pea_random_state,
            k_max=args.physical_regime_k_max,
        )
    )
    if not mld_long_df.empty:
        mld_long_df = mld_long_df.merge(
            summary_df[
                [
                    "profile_label",
                    "physical_regime_class",
                    "physical_regime_max_probability",
                    "physical_regime_entropy",
                ]
            ],
            on="profile_label",
            how="left",
        )

    print("\n[6/9] Classifying deep intrusion...")
    if args.deep_intrusion_mode != "none":
        summary_df, intrusion_diagnostics_df = _classify_deep_intrusion(
            summary_df=summary_df,
            mode=args.deep_intrusion_mode,
            quantile=args.deep_intrusion_quantile,
            value_col="sigma0_lower_mean_kg_m3",
            date_col="profile_date",
        )
        summary_df, oxygen_intrusion_events_df = _classify_oxygen_intrusions(
            summary_df=summary_df,
            observations=df,
            profile_cols=profile_cols,
            depth_col=args.depth_col,
            oxygen_col=args.deep_intrusion_oxygen_col,
            split_depth=layer_split_depth,
            low_oxygen_max=args.oxygen_low_compartment_max,
            bottom_n=args.oxygen_intrusion_bottom_n,
            onset_threshold=args.oxygen_intrusion_onset_threshold,
            persistence_threshold=args.oxygen_intrusion_persistence_threshold,
            end_consecutive=args.oxygen_intrusion_end_consecutive,
        )
        summary_df, nitrate_renewal_events_df, oxygen_anomalies_df = (
            _classify_nitrate_renewal(
                summary_df=summary_df,
                observations=df,
                profile_cols=profile_cols,
                depth_col=args.depth_col,
                nitrate_col=args.renewal_nitrate_col,
                split_depth=layer_split_depth,
                bottom_n=args.renewal_nitrate_bottom_n,
                min_depths=args.renewal_nitrate_min_depths,
                detection_limit=args.renewal_nitrate_detection_limit,
                bridge_enabled=args.renewal_bridge_enabled,
                bridge_max_cruises=args.renewal_bridge_max_cruises,
                bridge_max_days=args.renewal_bridge_max_days,
            )
        )
        if not mld_long_df.empty:
            mld_long_df = mld_long_df.merge(
                summary_df[
                    [
                        "profile_label",
                        "deep_intrusion_class",
                        "deep_intrusion_threshold",
                        "deep_intrusion_score",
                        "deep_intrusion_expected_sigma0_kg_m3",
                        "deep_intrusion_density_anomaly_kg_m3",
                        "deep_intrusion_density_change_kg_m3",
                    ]
                ],
                on="profile_label",
                how="left",
            )

    print("\n[7/9] Writing outputs...")
    density_df.to_csv(out_dir / "stratification_density_profiles.tsv", sep="\t", index=False)
    n2_df.to_csv(out_dir / "stratification_n2_profiles.tsv", sep="\t", index=False)
    summary_df.to_csv(out_dir / "stratification_summary.tsv", sep="\t", index=False)
    mld_long_df.to_csv(out_dir / "stratification_mld_timeseries.tsv", sep="\t", index=False)
    if feature_trace_df is not None and not feature_trace_df.empty:
        feature_trace_df.to_csv(out_dir / "stratification_feature_traces.tsv", sep="\t", index=False)
    if args.pea_class_mode != "none":
        pea_cluster_df.to_csv(out_dir / "stratification_pea_clusters.tsv", sep="\t", index=False)
        pea_bootstrap_df.to_csv(out_dir / "stratification_pea_threshold_bootstrap.tsv", sep="\t", index=False)
        pea_lower_cluster_df.to_csv(
            out_dir / "stratification_pea_lower_clusters.tsv", sep="\t", index=False
        )
        pea_lower_bootstrap_df.to_csv(
            out_dir / "stratification_pea_lower_threshold_bootstrap.tsv", sep="\t", index=False
        )
        _pea_class_validation(summary_df).to_csv(
            out_dir / "stratification_pea_class_validation.tsv", sep="\t", index=False
        )
    if args.deep_intrusion_mode != "none":
        intrusion_diagnostics_df.to_csv(
            out_dir / "stratification_deep_intrusion_model.tsv", sep="\t", index=False
        )
        if not oxygen_intrusion_events_df.empty:
            oxygen_intrusion_events_df.to_csv(
                out_dir / "stratification_oxygen_intrusion_events.tsv",
                sep="\t", index=False,
            )
        nitrate_renewal_events_df.to_csv(
            out_dir / "stratification_nitrate_renewal_events.tsv",
            sep="\t", index=False,
        )
        oxygen_anomalies_df.to_csv(
            out_dir / "stratification_oxygen_anomalies.tsv",
            sep="\t", index=False,
        )
    physical_selection_df.to_csv(
        out_dir / "stratification_physical_regime_selection.tsv", sep="\t", index=False
    )
    physical_clusters_df.to_csv(
        out_dir / "stratification_physical_regime_clusters.tsv", sep="\t", index=False
    )
    physical_projection_df.to_csv(
        out_dir / "stratification_physical_regime_projection.tsv", sep="\t", index=False
    )

    if not args.no_plots:
        print("\n[8/9] Plotting review figure...")
        _plot_review(
            summary_df=summary_df,
            mld_long_df=mld_long_df,
            feature_df=feature_trace_df,
            feature_cols=plot_feature_cols,
            out_path=out_dir / "stratification_review.pdf",
            date_col=args.date_col,
            layer_split_depth=layer_split_depth,
        )
        if args.pea_class_mode != "none":
            _plot_pea_diagnostics(summary_df, pea_cluster_df, out_dir)
        if args.deep_intrusion_mode != "none":
            _plot_intrusion_diagnostics(summary_df, out_dir)
            _plot_oxygen_intrusion_diagnostics(summary_df, out_dir)
        _plot_physical_regime_selection(physical_selection_df, out_dir)
        _plot_physical_regime_clusters(
            physical_projection_df, physical_clusters_df, summary_df, out_dir
        )

    print("\n[9/9] Done.")
    print(f"Outputs saved to: {out_dir}\n")


if __name__ == "__main__":
    main()
