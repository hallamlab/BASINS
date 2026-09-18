# Path: stratification_timeseries_anomaly.py
# Full replacement script (single-input biochem/metadata table)
# Outputs (exact same filenames as before):
#   - stratification_timeseries.tsv
#   - annual_extremes.tsv
#   - stratification_monthly_profile.pdf
#   - stratification_vs_pea_timeseries.pdf (when --pea-metrics is supplied)
#   - stratification_physical_biochem_timeseries.{pdf,png,svg} (when --pea-metrics is supplied)
#   - stratification_physical_biochem_timeseries_complete_cases.{pdf,png,svg}
#   - stratification_physical_biochem_timeseries.tsv (when --pea-metrics is supplied)

#!/usr/bin/env python3
"""
stratification_timeseries_anomaly.py

Single-table stratification time-series analysis using a biochem/physical table
that also contains all metadata required for plotting.

Produces ONLY:
  - stratification_timeseries.tsv
  - annual_extremes.tsv
  - stratification_monthly_profile.pdf

Feature selection:
  - If --features is provided: use that comma-separated list.
  - Else: use all columns AFTER --features-after-col (default: Depth_anchored).
  - Then apply coverage filtering with COVERAGE_THRESHOLD (same behavior as before).

Example:
  python stratification_timeseries_anomaly.py \
    --input biochem_table.tsv \
    --sample-id-col cruise_year_month_depth \
    --date-col date \
    --month-col Month \
    --year-col Year \
    --depth-col Depth \
    --output-dir strat_monthly_out
"""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path
from itertools import combinations
from typing import List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import zscore
from scipy.ndimage import gaussian_filter1d
from sklearn.ensemble import IsolationForest
from sklearn.neighbors import LocalOutlierFactor

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shared_plot_style import install_publication_style

warnings.filterwarnings("ignore")

install_publication_style()

COVERAGE_THRESHOLD = 0.51


# ============================================================================
# Helpers
# ============================================================================

def _safe_euclidean(a: np.ndarray, b: np.ndarray) -> float | None:
    mask = np.isfinite(a) & np.isfinite(b)
    if not mask.any():
        return None
    diff = a[mask] - b[mask]
    return np.linalg.norm(diff)


def _split_features_arg(s: str) -> List[str]:
    if s is None:
        return []
    parts = [p.strip() for p in s.split(",")]
    return [p for p in parts if p]


# ============================================================================
# Feature selection
# ============================================================================

def select_feature_columns(
    df: pd.DataFrame,
    features_csv: str | None,
    features_after_col: str,
    coverage_threshold: float,
) -> List[str]:
    """
    Select candidate feature columns, then filter by coverage threshold.
    """
    if features_csv:
        candidates = _split_features_arg(features_csv)
        missing = [c for c in candidates if c not in df.columns]
        if missing:
            raise ValueError(f"--features includes columns not found in input: {missing}")
    else:
        if features_after_col not in df.columns:
            raise ValueError(
                f"--features-after-col '{features_after_col}' not found in input columns."
            )
        idx = list(df.columns).index(features_after_col)
        candidates = list(df.columns)[idx + 1 :]
        if not candidates:
            raise ValueError(
                f"No columns found after --features-after-col '{features_after_col}'."
            )

    # Drop ASV-like cols if present (keeps prior behavior)
    candidates = [c for c in candidates if not str(c).startswith("ASV")]

    # Coverage filter (same concept as before: fraction of rows with non-null)
    cov = df[candidates].notna().sum() / max(len(df), 1)
    keep = [c for c in candidates if cov.get(c, 0.0) >= coverage_threshold]
    drop = [c for c in candidates if c not in keep]

    print(f"  [i] Keeping {len(keep)} feature columns for stratification: {keep}")
    if drop:
        print(f"  [i] Dropping undersampled columns (<{coverage_threshold*100:.0f}% coverage): {drop}")
    if not keep:
        raise ValueError("No feature columns meet the coverage requirement.")
    return keep


# ============================================================================
# Core computations
# ============================================================================

def calculate_stratification_score_timeseries(
    integrated_data: pd.DataFrame,
    metadata: pd.DataFrame,
    cruise_col: str,
    date_col: str,
    depth_col: str,
    feature_cols: List[str],
) -> pd.DataFrame:
    """
    For each unique date:
      - compute mean pairwise distance between depth centroids (features)
    """
    print("  [i] Calculating stratification score time series...")

    results = []
    unique_dates = sorted(metadata[date_col].unique())

    for date in unique_dates:
        date_mask = metadata[date_col] == date
        date_meta = metadata.loc[date_mask]
        date_data = integrated_data.loc[date_mask, feature_cols]

        cruise = date_meta[cruise_col].unique()

        if date_data.empty or date_meta.empty:
            continue

        depths = sorted(date_meta[depth_col].unique())
        depth_centroids = {}

        for depth in depths:
            depth_mask = date_meta[depth_col] == depth
            depth_samples = date_data.loc[depth_mask]
            if depth_samples.empty:
                continue
            centroid = depth_samples.mean(axis=0, skipna=True)
            if centroid.notna().sum() == 0:
                continue
            depth_centroids[depth] = centroid.values

        distances = []
        for depth1, depth2 in combinations(depth_centroids.keys(), 2):
            dist = _safe_euclidean(depth_centroids[depth1], depth_centroids[depth2])
            if dist is not None:
                distances.append(dist)

        if not distances:
            continue

        mean_dist = np.mean(distances)
        total_cells = date_data.size
        non_na = np.isfinite(date_data.values).sum()
        coverage = non_na / total_cells if total_cells > 0 else 0.0

        results.append(
            {
                'Cruise': cruise[0] if len(cruise) == 1 else "multiple",
                "date": date,
                "stratification_score": mean_dist,
                "n_depths": len(depths),
                "depths_present": ",".join(map(str, depths)),
                "n_samples": len(date_data),
                "coverage": coverage,
            }
        )

    timeseries_df = pd.DataFrame(results)
    print(f"      Computed {len(timeseries_df)} time points")
    return timeseries_df


def normalize_to_centered_scale(timeseries_df: pd.DataFrame) -> pd.DataFrame:
    """
    Normalize stratification scores to [-1, 1] centered at 0 (median at 0).
    """
    print("  [i] Normalizing to centered scale...")

    min_score = timeseries_df["stratification_score"].min()
    max_score = timeseries_df["stratification_score"].max()
    median_score = timeseries_df["stratification_score"].median()

    timeseries_df = timeseries_df.copy()
    normalized_scores = []
    for score in timeseries_df["stratification_score"]:
        if score < median_score:
            norm = -1 + (score - min_score) / (median_score - min_score)
        else:
            norm = (score - median_score) / (max_score - median_score)
        normalized_scores.append(norm)

    timeseries_df["normalized_score"] = normalized_scores
    print(
        f"      Normalized: min={min(normalized_scores):.3f}, median=0, max={max(normalized_scores):.3f}"
    )
    return timeseries_df


def detect_anomalies_consensus(
    timeseries_df: pd.DataFrame,
    consensus_threshold: int = 2,
    date_col: str = "date",
    window_months: int = 12,
    min_points_in_window: int = 5,
) -> pd.DataFrame:
    """
    Detect anomalies using a sliding N-month window around each point.

    Requires:
      - timeseries_df['normalized_score']
      - timeseries_df[date_col] parseable as datetime

    For each point i:
      - Build a local window of points with dates in [date_i - window_months, date_i + window_months]
      - If window has < min_points_in_window points (including i): no anomaly (n_votes=0)
      - Run 4 methods on the WINDOW distribution, then evaluate whether point i is an outlier
        relative to that window.
      - Output schema preserved:
          is_anomaly (bool), anomaly_type (str), n_votes (int)

    Notes:
      - Z-score + IQR are computed from the window stats.
      - IsolationForest + LOF are fit on the window, then we use their predictions for point i.
      - This yields a *local* notion of anomaly (good for regime shifts / nonstationarity).
    """
    print(
        f"  [i] Detecting anomalies (±{window_months} month sliding window, consensus ≥{consensus_threshold})..."
    )

    if "normalized_score" not in timeseries_df.columns:
        raise ValueError("timeseries_df missing required column: normalized_score")
    if date_col not in timeseries_df.columns:
        raise ValueError(f"timeseries_df missing required date column: {date_col}")

    out = timeseries_df.copy()

    # Parse dates and keep original order stable
    dates = pd.to_datetime(out[date_col], errors="coerce")
    scores = pd.to_numeric(out["normalized_score"], errors="coerce").to_numpy(dtype=float)

    out["is_anomaly"] = False
    out["anomaly_type"] = "normal"
    out["n_votes"] = 0

    finite = np.isfinite(scores) & dates.notna().to_numpy()

    if finite.sum() < 3:
        return out

    # Work in date order for predictable windows
    order = np.argsort(dates.to_numpy())
    dates_s = dates.to_numpy()[order]
    scores_s = scores[order]
    finite_s = finite[order]

    # Precompute month window as a DateOffset
    halfwin = pd.DateOffset(months=int(window_months))

    # Map back to original row indices
    inv_order = np.empty_like(order)
    inv_order[order] = np.arange(len(order))

    for pos in range(len(scores_s)):
        if not finite_s[pos]:
            continue

        center_date = pd.Timestamp(dates_s[pos])

        left = center_date - halfwin
        right = center_date + halfwin

        in_window = (dates_s >= np.datetime64(left)) & (dates_s <= np.datetime64(right)) & finite_s
        w_idx = np.where(in_window)[0]
        n = len(w_idx)

        if n < int(min_points_in_window):
            continue

        y = scores_s[w_idx]
        # Identify where the center point sits inside the window
        # (pos is in sorted space)
        center_in_window = np.where(w_idx == pos)[0]
        if center_in_window.size != 1:
            # Shouldn't happen, but be defensive
            continue
        j = int(center_in_window[0])
        x0 = float(y[j])

        votes = 0

        # Method 1: Z-score (window)
        mu = float(np.nanmean(y))
        sd = float(np.nanstd(y, ddof=0))
        if sd > 0:
            z = abs((x0 - mu) / sd)
            votes += int(z > 1.5)

        # Method 2: IQR (window)
        Q1 = float(np.nanquantile(y, 0.25))
        Q3 = float(np.nanquantile(y, 0.75))
        IQR = Q3 - Q1
        if IQR > 0:
            lower = Q1 - 1.5 * IQR
            upper = Q3 + 1.5 * IQR
            votes += int((x0 < lower) or (x0 > upper))

        # Method 3: Isolation Forest (fit on window, predict all window points)
        # (contamination is still a fixed fraction; keep as-is for now)
        try:
            iso = IsolationForest(contamination=0.15, random_state=42)
            iso_pred = iso.fit_predict(y.reshape(-1, 1))
            votes += int(iso_pred[j] == -1)
        except Exception:
            pass

        # Method 4: LOF (fit on window, predict all window points)
        try:
            if n >= 3:
                n_neighbors = min(10, n - 1)
                lof = LocalOutlierFactor(n_neighbors=n_neighbors, contamination=0.15)
                lof_pred = lof.fit_predict(y.reshape(-1, 1))
                votes += int(lof_pred[j] == -1)
        except Exception:
            pass

        is_anom = votes >= int(consensus_threshold)

        # Write to *original* row index
        orig_i = int(order[pos])
        out.at[orig_i, "n_votes"] = int(votes)
        out.at[orig_i, "is_anomaly"] = bool(is_anom)

        if is_anom:
            out.at[orig_i, "anomaly_type"] = "high_stratification" if x0 > 0 else "mixing_event"
        else:
            out.at[orig_i, "anomaly_type"] = "normal"

    n_anomalies = int(out["is_anomaly"].sum())
    n_mixing = int((out["anomaly_type"] == "mixing_event").sum())
    n_high_strat = int((out["anomaly_type"] == "high_stratification").sum())

    print(f"      Total anomalies: {n_anomalies}/{len(out)} time points")
    print(f"      Mixing events: {n_mixing}")
    print(f"      High stratification: {n_high_strat}")

    return out


def identify_annual_extremes(
    timeseries_df: pd.DataFrame,
    metadata: pd.DataFrame,
    year_col: str,
) -> pd.DataFrame:
    """
    Identify max stratification and max mixing events per year.
    Kept identical to the earlier (quirky) behavior to preserve outputs.
    """
    print("  [i] Identifying annual extremes...")

    timeseries_df = timeseries_df.copy()

    # This mapping logic is intentionally preserved
    date_to_year = metadata.set_index(metadata.columns[0])[year_col].to_dict()

    years = []
    for date in timeseries_df["date"]:
        year_found = False
        for sample, meta_year in date_to_year.items():
            sample_meta = metadata[metadata.index == sample]
            if len(sample_meta) > 0:
                sample_date = sample_meta[sample_meta.columns[0]].iloc[0]
                if sample_date == date:
                    years.append(meta_year)
                    year_found = True
                    break

        if not year_found:
            try:
                if isinstance(date, str):
                    year = int(date.split("-")[0]) if "-" in date else int(date[:4])
                else:
                    year = date.year
                years.append(year)
            except Exception:
                years.append(None)

    timeseries_df["year"] = years

    extremes = []
    for year in sorted(timeseries_df["year"].dropna().unique()):
        year_data = timeseries_df[timeseries_df["year"] == year]

        max_strat_idx = year_data["normalized_score"].idxmax()
        max_strat_row = year_data.loc[max_strat_idx]
        extremes.append(
            {
                "date": max_strat_row["date"],
                "year": year,
                "stratification_score": max_strat_row["stratification_score"],
                "normalized_score": max_strat_row["normalized_score"],
                "extreme_type": "max_stratification",
            }
        )

        min_mix_idx = year_data["normalized_score"].idxmin()
        min_mix_row = year_data.loc[min_mix_idx]
        extremes.append(
            {
                "date": min_mix_row["date"],
                "year": year,
                "stratification_score": min_mix_row["stratification_score"],
                "normalized_score": min_mix_row["normalized_score"],
                "extreme_type": "max_mixing",
            }
        )

    extremes_df = pd.DataFrame(extremes)
    print(
        f"      Found {len(extremes_df)} extremes across {len(extremes_df['year'].unique())} years"
    )
    return extremes_df


# ============================================================================
# PEA comparison plot
# ============================================================================

def _load_pea_timeseries(pea_path: Path, date_col: str) -> pd.DataFrame:
    pea_df = pd.read_csv(pea_path, sep="\t")
    if date_col not in pea_df.columns:
        raise ValueError(f"PEA metrics file missing date column '{date_col}'")
    pea_df = pea_df.copy()
    pea_df["date"] = pd.to_datetime(pea_df[date_col], errors="coerce")
    pea_df = pea_df.dropna(subset=["date"])

    numeric_cols = [
        "pea_J_m3",
        "pea_upper_J_m3",
        "pea_lower_J_m3",
        "pea_threshold_low",
        "pea_threshold_high",
        "deep_intrusion_threshold",
        "deep_intrusion_score",
        "deep_intrusion_expected_sigma0_kg_m3",
        "deep_intrusion_density_anomaly_kg_m3",
        "deep_intrusion_density_change_kg_m3",
        "deep_intrusion_oxygen_lower_mean",
        "oxygen_intrusion_event_id",
        "oxygen_intrusion_bottom_depths_n",
        "oxygen_intrusion_bottom_median_um",
        "oxygen_intrusion_change_median_um",
        "oxygen_intrusion_baseline_change_median_um",
        "oxygen_intrusion_onset_threshold_um",
        "oxygen_intrusion_persistence_threshold_um",
        "oxygen_intrusion_end_consecutive_cruises",
        "renewal_event_id",
        "renewal_nitrate_bottom_depths_n",
        "renewal_nitrate_bottom_median_um",
        "renewal_nitrate_detection_limit_um",
        "renewal_nitrate_min_depths",
        "n2_mean_s-2",
        "n2_max_s-2",
        "n2_max_depth_m",
        "pycnocline_depth_m",
        "sigma0_upper_lower_diff_kg_m3",
        "mld_depth_m_dr0p03",
        "mld_depth_m_dr0p125",
        "pea_threshold_low_ci_lower",
        "pea_threshold_low_ci_upper",
        "pea_threshold_high_ci_lower",
        "pea_threshold_high_ci_upper",
        "pea_lower_threshold_low",
        "pea_lower_threshold_high",
        "pea_lower_threshold_low_ci_lower",
        "pea_lower_threshold_low_ci_upper",
        "pea_lower_threshold_high_ci_lower",
        "pea_lower_threshold_high_ci_upper",
        "physical_regime_probability_mixed",
        "physical_regime_probability_intermediate",
        "physical_regime_probability_stratified",
        "physical_regime_max_probability",
        "physical_regime_entropy",
    ]
    for col in numeric_cols:
        if col in pea_df.columns:
            pea_df[col] = pd.to_numeric(pea_df[col], errors="coerce")

    keep_cols = [c for c in ["pea_J_m3", "pea_upper_J_m3", "pea_lower_J_m3"] if c in pea_df.columns]
    if not keep_cols:
        raise ValueError("PEA metrics file missing PEA columns (pea_J_m3/pea_upper_J_m3/pea_lower_J_m3)")

    retained_numeric = [col for col in numeric_cols if col in pea_df.columns]
    retained_classes = [
        col for col in [
            "pea_class", "pea_lower_class", "physical_regime_class", "deep_intrusion_class",
            "oxygen_intrusion_class", "renewal_phase", "renewal_nitrate_status",
            "renewal_inference_reason", "renewal_nitrate_bottom_depths_m",
            "oxygen_anomaly_class",
        ]
        if col in pea_df.columns
    ]
    retained_flags = [
        col for col in [
            "renewal_onset", "renewal_last_active", "renewal_end_confirmed",
            "renewal_phase_inferred", "renewal_nitrate_supported",
        ]
        if col in pea_df.columns
    ]
    aggregations = {col: "mean" for col in retained_numeric}
    aggregations.update({col: "first" for col in retained_classes})
    aggregations.update({col: "first" for col in retained_flags})
    if "Cruise" in pea_df.columns:
        aggregations["Cruise"] = "first"
    pea_ts = pea_df.groupby("date", sort=True).agg(aggregations).reset_index()
    return pea_ts


def _combine_stratification_and_physical_metrics(
    timeseries_df: pd.DataFrame,
    pea_df: pd.DataFrame,
) -> pd.DataFrame:
    """Create the auditable source table behind the combined time-series plot."""
    strat_df = timeseries_df.copy()
    strat_df["date"] = pd.to_datetime(strat_df["date"], errors="coerce")
    strat_df = strat_df.dropna(subset=["date"])
    physical_df = pea_df.copy()
    if "Cruise" in physical_df.columns:
        physical_df = physical_df.rename(columns={"Cruise": "physical_cruise"})
    combined = pd.merge(physical_df, strat_df, on="date", how="outer", sort=True)
    if "physical_cruise" in combined.columns:
        if "Cruise" in combined.columns:
            combined["Cruise"] = combined["physical_cruise"].combine_first(combined["Cruise"])
        else:
            combined["Cruise"] = combined["physical_cruise"]
        combined = combined.drop(columns=["physical_cruise"])
    combined = combined.sort_values("date").reset_index(drop=True)
    combined["date"] = combined["date"].dt.strftime("%Y-%m-%d")
    return combined


def _year_panels(dates: pd.Series, max_years: int = 2) -> list[list[int]]:
    years = sorted(pd.to_datetime(dates, errors="coerce").dropna().dt.year.unique().tolist())
    return [years[i:i + max_years] for i in range(0, len(years), max_years)]


def _shared_axis_limits(series: list[pd.Series]) -> tuple[float, float] | None:
    """Return one padded finite range for repeated time-series panels."""
    arrays = [pd.to_numeric(values, errors="coerce").to_numpy(dtype=float) for values in series]
    finite = np.concatenate([values[np.isfinite(values)] for values in arrays if values.size])
    if finite.size == 0:
        return None
    lower = float(finite.min())
    upper = float(finite.max())
    if lower == upper:
        pad = max(abs(lower) * 0.05, 0.5)
    else:
        pad = (upper - lower) * 0.05
    return lower - pad, upper + pad


def _plot_with_missing_connectors(
    ax, x, dates, values, *, color: str, linewidth: float, label: str, sequence=None
) -> None:
    """Plot solid observed runs and dashed bridges across missing observations.

    This follows the older PEA comparison plot's gap convention while also
    detecting explicit NA rows in the combined audit table.  A bridge is used
    when observations have missing table rows between them or their calendar
    separation is more than twice the series' median sampling interval.
    """
    y = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(dtype=float)
    dates = pd.to_datetime(pd.Series(dates), errors="coerce")
    valid = np.flatnonzero(np.isfinite(y) & dates.notna().to_numpy())
    if valid.size == 0:
        return
    valid_dates = dates.iloc[valid].to_numpy(dtype="datetime64[ns]")
    source_sequence = np.asarray(sequence if sequence is not None else np.arange(len(y)))[valid]
    date_deltas = np.diff(valid_dates).astype("timedelta64[D]").astype(float)
    positive_deltas = date_deltas[np.isfinite(date_deltas) & (date_deltas > 0)]
    median_delta = float(np.median(positive_deltas)) if positive_deltas.size else np.nan

    run_start = 0
    first_solid = True
    for j in range(1, valid.size):
        explicit_missing = (
            valid[j] - valid[j - 1] > 1
            or source_sequence[j] - source_sequence[j - 1] > 1
        )
        calendar_gap = (
            np.isfinite(median_delta)
            and median_delta > 0
            and date_deltas[j - 1] > (2.0 * median_delta)
        )
        if explicit_missing or calendar_gap:
            run = valid[run_start:j]
            ax.plot(
                np.asarray(x)[run], y[run], color=color, linewidth=linewidth,
                label=label if first_solid else "_nolegend_",
            )
            first_solid = False
            bridge = valid[j - 1:j + 1]
            ax.plot(
                np.asarray(x)[bridge], y[bridge], color=color, linewidth=linewidth,
                linestyle="--", label="_nolegend_",
            )
            run_start = j
    run = valid[run_start:]
    ax.plot(
        np.asarray(x)[run], y[run], color=color, linewidth=linewidth,
        label=label if first_solid else "_nolegend_",
    )


def _date_cruise_labels(df: pd.DataFrame) -> list[str]:
    """Format combined time-series ticks as ISO date plus an integer cruise ID."""
    cruises = df.get("Cruise", pd.Series(pd.NA, index=df.index))
    labels = []
    for date, cruise in zip(df["date"], cruises):
        date_label = pd.Timestamp(date).strftime("%Y-%m-%d")
        numeric_cruise = pd.to_numeric(pd.Series([cruise]), errors="coerce").iloc[0]
        if pd.notna(numeric_cruise) and float(numeric_cruise).is_integer():
            cruise_label = str(int(numeric_cruise))
        elif pd.notna(cruise) and str(cruise).strip():
            cruise_label = str(cruise).strip()
        else:
            cruise_label = "NA"
        labels.append(f"{date_label} [{cruise_label}]")
    return labels


def plot_physical_biochem_timeseries(
    combined_df: pd.DataFrame,
    output_base: Path,
    title: str = "Physical stratification and biochemical depth separation",
    years_per_panel: int = 2,
    align_months: bool = False,
) -> None:
    """Plot physical and biochemical stratification with comparable split panels."""
    plot_df = combined_df.copy()
    plot_df["date"] = pd.to_datetime(plot_df["date"], errors="coerce")
    plot_df = plot_df.dropna(subset=["date"]).sort_values("date")
    if plot_df.empty:
        print("  [i] Skipping combined physical/biochemical plot (no dated data).")
        return

    panels = _year_panels(plot_df["date"], max_years=years_per_panel)
    if not panels:
        return

    fig, axes = plt.subplots(
        len(panels),
        1,
        figsize=(22, max(4.0, 3.8 * len(panels))),
        squeeze=False,
    )
    axes = axes[:, 0]
    pea_lines = [
        ("pea_J_m3", "PEA total", "tab:green"),
        ("pea_upper_J_m3", "PEA upper", "tab:blue"),
        ("pea_lower_J_m3", "PEA lower", "tab:orange"),
    ]
    left_limits = _shared_axis_limits([plot_df["stratification_score"]]) if "stratification_score" in plot_df else None
    right_limits = _shared_axis_limits([plot_df[col] for col, _, _ in pea_lines if col in plot_df])

    legend_handles = None
    legend_labels = None
    for ax_left, years in zip(axes, panels):
        sub = plot_df[plot_df["date"].dt.year.isin(years)].copy().reset_index(drop=True)
        x = sub["date"].dt.month.to_numpy() if align_months else np.arange(len(sub))
        source_sequence = sub.get("_source_order", pd.Series(x, index=sub.index)).to_numpy()
        ax_right = ax_left.twinx()

        if "stratification_score" in sub.columns:
            dcd = pd.to_numeric(sub["stratification_score"], errors="coerce").to_numpy()
            _plot_with_missing_connectors(
                ax_left, x, sub["date"], dcd, color="black", linewidth=1.3,
                label="Depth-centroid distance (D)", sequence=source_sequence,
            )
            mask = np.isfinite(dcd)
            ax_left.scatter(x[mask], dcd[mask], marker="o", s=34, facecolor="black",
                            edgecolor="black", linewidth=0.7, zorder=5, label="_nolegend_")

        for col, label, color in pea_lines:
            if col not in sub.columns:
                continue
            values = pd.to_numeric(sub[col], errors="coerce").to_numpy()
            _plot_with_missing_connectors(
                ax_right, x, sub["date"], values, color=color, linewidth=1.2, label=label,
                sequence=source_sequence,
            )
            mask = np.isfinite(values)
            ax_right.scatter(x[mask], values[mask], marker="o", s=30, facecolor=color,
                             edgecolor=color, linewidth=0.8, zorder=4, label="_nolegend_")

        if align_months:
            ax_left.set_xlim(0.5, 12.5)
        else:
            ax_left.set_xlim(-0.5, max(len(sub) - 0.5, 0.5))
        ax_left.set_ylabel("Depth-centroid distance (D)")
        ax_right.set_ylabel("PEA (J m$^{-3}$)")
        if left_limits is not None:
            ax_left.set_ylim(*left_limits)
        if right_limits is not None:
            ax_right.set_ylim(*right_limits)
        ax_left.grid(axis="y", linestyle="--", alpha=0.30)
        ax_left.set_title(str(years[0]) if len(years) == 1 else f"{years[0]}–{years[-1]}")
        if align_months:
            month_names = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
            labels_by_month = dict(zip(sub["date"].dt.month, _date_cruise_labels(sub)))
            ax_left.set_xticks(np.arange(1, 13))
            ax_left.set_xticklabels(
                [f"{name}\n{labels_by_month.get(month, '—')}" for month, name in enumerate(month_names, 1)],
                rotation=90, fontsize=7,
            )
        else:
            ax_left.set_xticks(x)
            ax_left.set_xticklabels(_date_cruise_labels(sub), rotation=90, fontsize=8)

        if legend_handles is None:
            left_h, left_l = ax_left.get_legend_handles_labels()
            right_h, right_l = ax_right.get_legend_handles_labels()
            legend_handles = left_h + right_h
            legend_labels = left_l + right_l

    fig.suptitle(title, y=0.998)
    if legend_handles:
        fig.legend(legend_handles, legend_labels, loc="upper left", bbox_to_anchor=(0.84, 0.985))
    fig.tight_layout(rect=[0, 0, 0.83, 0.99])
    for suffix in [".pdf", ".png", ".svg"]:
        fig.savefig(output_base.with_suffix(suffix), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  [✓] Saved combined physical/biochemical time-series plot")


def plot_physical_biochem_monthly_profile(combined_df: pd.DataFrame, output_base: Path) -> None:
    """Legacy-style monthly profiles for each current physical/biochemical metric."""
    plot_df = combined_df.copy()
    plot_df["date"] = pd.to_datetime(plot_df["date"], errors="coerce")
    plot_df = plot_df.dropna(subset=["date"]).sort_values("date").reset_index(drop=True)
    metrics = [
        (
            "stratification_score",
            "Depth-centroid distance (D)",
            "Depth-centroid distance",
        ),
        ("pea_J_m3", "PEA (J m$^{-3}$)", "Total potential energy anomaly"),
        ("pea_upper_J_m3", "PEA (J m$^{-3}$)", "Upper-water-column potential energy anomaly"),
        ("pea_lower_J_m3", "PEA (J m$^{-3}$)", "Lower-water-column potential energy anomaly"),
    ]
    metrics = [item for item in metrics if item[0] in plot_df.columns]
    if not metrics:
        return

    fig, axes = plt.subplots(len(metrics), 1, figsize=(17, 3.4 * len(metrics)), sharex=True, squeeze=False)
    axes = axes[:, 0]
    plot_range = np.linspace(0.8, 12.2, 300)
    for panel_index, (ax, (column, ylabel, panel_title)) in enumerate(zip(axes, metrics)):
        retained_columns = ["date", column]
        if column == "stratification_score" and "coverage" in plot_df.columns:
            retained_columns.append("coverage")
        metric_df = plot_df[retained_columns].copy()
        metric_df[column] = pd.to_numeric(metric_df[column], errors="coerce")
        metric_df = metric_df.dropna(subset=[column])
        metric_df["month"] = metric_df["date"].dt.month.astype(int)
        metric_df["year"] = metric_df["date"].dt.year.astype(int)

        # All observations remain visible, but are deliberately unconnected.
        ax.scatter(
            metric_df["month"], metric_df[column], color="black", s=19,
            alpha=0.38, edgecolor="none", zorder=3,
        )

        monthly = metric_df.groupby("month")[column].agg(
            median="median",
            q25=lambda values: values.quantile(0.25),
            q75=lambda values: values.quantile(0.75),
        ).reindex(range(1, 13))
        smooth = {}
        for statistic in ["median", "q25", "q75"]:
            values = monthly[statistic].interpolate(limit_direction="both")
            smooth[statistic] = np.interp(
                plot_range,
                np.arange(1, 13),
                gaussian_filter1d(values.to_numpy(dtype=float), sigma=1.1),
            )
        ax.fill_between(
            plot_range, smooth["q25"], smooth["q75"],
            color="lightgrey", alpha=0.65, linewidth=0, zorder=1,
        )
        ax.plot(plot_range, smooth["median"], color="black", linewidth=3.0, zorder=2)

        # Match the legacy completeness safeguard: annual extrema are shown
        # only for years with at least seven adequately observed months.
        eligibility_df = metric_df
        if column == "stratification_score" and "coverage" in metric_df.columns:
            eligibility_df = metric_df[
                pd.to_numeric(metric_df["coverage"], errors="coerce") >= COVERAGE_THRESHOLD
            ]
        eligible_counts = eligibility_df.groupby("year")["month"].nunique()
        eligible_years = set(eligible_counts[eligible_counts >= 7].index.tolist())
        for year, annual in metric_df.groupby("year", sort=True):
            if year not in eligible_years:
                continue
            max_row = annual.loc[annual[column].idxmax()]
            min_row = annual.loc[annual[column].idxmin()]
            ax.scatter(
                max_row["month"], max_row[column], marker="^", color="black",
                s=72, linewidth=0.8, zorder=5,
            )
            ax.scatter(
                min_row["month"], min_row[column], marker="v", facecolor="white",
                edgecolor="black", s=72, linewidth=1.2, zorder=5,
            )
        print(
            f"      {panel_title}: annual extrema shown for "
            f"{len(eligible_years)}/{metric_df['year'].nunique()} years (>=7 eligible months)"
        )

        ax.set_title(
            f"{chr(65 + panel_index)}  {panel_title}", loc="left",
            fontsize=12, fontweight="semibold", pad=7,
        )
        ax.set_ylabel(ylabel)
        ax.set_xlim(0.5, 12.5)
        ax.grid(axis="y", linestyle="--", alpha=0.30)

    axes[-1].set_xticks(np.arange(1, 13))
    axes[-1].set_xticklabels(["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
    axes[-1].set_xlabel("Calendar month")
    fig.suptitle("Monthly physical and biochemical profiles", y=0.985, fontsize=13)
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    legend_handles = [
        Line2D([], [], color="black", marker="o", linestyle="None", alpha=0.38, label="Observation"),
        Line2D([], [], color="black", linewidth=3, label="Monthly median"),
        Patch(facecolor="lightgrey", alpha=0.65, label="Monthly IQR"),
        Line2D([], [], color="black", marker="^", linestyle="None", label="Annual maximum"),
        Line2D([], [], marker="v", linestyle="None", markerfacecolor="white",
               markeredgecolor="black", color="black", label="Annual minimum"),
    ]
    fig.legend(
        handles=legend_handles, loc="upper left", bbox_to_anchor=(0.805, 0.92),
        ncol=1, frameon=False,
    )
    fig.tight_layout(rect=[0, 0, 0.79, 0.95])
    for suffix in [".pdf", ".png", ".svg"]:
        fig.savefig(output_base.with_suffix(suffix), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  [✓] Saved all-data monthly physical/biochemical profile")


def plot_stratification_vs_pea_timeseries(
    timeseries_df: pd.DataFrame,
    pea_df: pd.DataFrame,
    output_path: Path,
) -> None:
    print("  [i] Creating stratification vs PEA time-series plot...")
    if timeseries_df.empty or pea_df.empty:
        print("  [i] Skipping PEA comparison plot (no data).")
        return

    strat_df = timeseries_df.copy()
    strat_df["date"] = pd.to_datetime(strat_df["date"], errors="coerce")
    strat_df = strat_df.dropna(subset=["date"]).sort_values("date")
    pea_df = pea_df.sort_values("date")

    def _plot_with_gaps(ax, x, y, color):
        sub = pd.DataFrame({"x": x, "y": y}).dropna().sort_values("x")
        if sub.empty:
            return
        xs = sub["x"].to_numpy()
        ys = sub["y"].to_numpy()
        if len(xs) == 1:
            ax.plot(xs, ys, color=color, label="_nolegend_")
            return

        deltas = np.diff(xs)
        median_delta = np.median(deltas)
        if not np.isfinite(median_delta) or median_delta == 0:
            ax.plot(xs, ys, color=color, label="_nolegend_")
            return

        gap_thresh = median_delta * 2
        seg_start = 0
        for i in range(1, len(xs)):
            if (xs[i] - xs[i - 1]) > gap_thresh:
                ax.plot(xs[seg_start:i], ys[seg_start:i], color=color, label="_nolegend_")
                ax.plot(
                    xs[i - 1:i + 1],
                    ys[i - 1:i + 1],
                    color=color,
                    linestyle=":",
                    label="_nolegend_",
                )
                seg_start = i
        ax.plot(xs[seg_start:], ys[seg_start:], color=color, label="_nolegend_")

    fig, ax_left = plt.subplots(1, 1, figsize=(22, 6))
    ax_right = ax_left.twinx()

    ax_left.scatter(
        strat_df["date"],
        strat_df["stratification_score"],
        color="black",
        s=30,
        label="Depth-centroid distance (D)",
    )
    _plot_with_gaps(
        ax_left,
        strat_df["date"],
        strat_df["stratification_score"],
        color="black",
    )
    ax_left.set_ylabel("Depth-centroid distance (D)")
    ax_left.grid(axis="y", linestyle="--", alpha=0.35)

    if "pea_J_m3" in pea_df.columns:
        ax_right.scatter(pea_df["date"], pea_df["pea_J_m3"], color="tab:green", s=30, label="PEA total")
        _plot_with_gaps(ax_right, pea_df["date"], pea_df["pea_J_m3"], color="tab:green")
    if "pea_upper_J_m3" in pea_df.columns:
        ax_right.scatter(pea_df["date"], pea_df["pea_upper_J_m3"], color="tab:blue", s=30, label="PEA upper")
        _plot_with_gaps(ax_right, pea_df["date"], pea_df["pea_upper_J_m3"], color="tab:blue")
    if "pea_lower_J_m3" in pea_df.columns:
        ax_right.scatter(pea_df["date"], pea_df["pea_lower_J_m3"], color="tab:orange", s=30, label="PEA lower")
        _plot_with_gaps(ax_right, pea_df["date"], pea_df["pea_lower_J_m3"], color="tab:orange")
    ax_right.set_ylabel("PEA (J/m3)")

    ax_left.set_title("Stability index vs PEA")
    ax_left.set_xlabel("Date")

    handles_left, labels_left = ax_left.get_legend_handles_labels()
    handles_right, labels_right = ax_right.get_legend_handles_labels()
    ax_left.legend(
        handles_left + handles_right,
        labels_left + labels_right,
        loc="upper left",
        bbox_to_anchor=(1.02, 1),
        borderaxespad=0,
    )

    fig.tight_layout(rect=[0, 0, 0.85, 1])
    for suffix in [".pdf", ".png", ".svg"]:
        fig.savefig(Path(output_path).with_suffix(suffix), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  [✓] Saved stratification vs PEA time-series plot")


# ============================================================================
# Plot (monthly profile)
# ============================================================================

def plot_stratification_monthly_profile(
    timeseries_df: pd.DataFrame,
    extremes_df: pd.DataFrame,
    metadata: pd.DataFrame,
    date_col: str,
    month_col: str,
    year_col: str,
    output_path: Path,
) -> None:
    """
    Monthly profile plot (kept byte-for-byte logic from prior pruned version).
    """
    print("  [i] Creating enriched monthly stratification profile...")
    plot_df = timeseries_df.copy()

    meta_reset = metadata.reset_index()
    if date_col not in meta_reset.columns:
        raise ValueError(f"Metadata missing date column '{date_col}' for monthly profile")
    if month_col not in meta_reset.columns or year_col not in meta_reset.columns:
        raise ValueError("Metadata must supply month/year columns for monthly profile")

    date_to_month = meta_reset.groupby(date_col)[month_col].first().to_dict()
    date_to_year = meta_reset.groupby(date_col)[year_col].first().to_dict()

    plot_df["month"] = plot_df["date"].map(date_to_month)
    plot_df["year"] = plot_df["date"].map(date_to_year)
    plot_df = plot_df.dropna(subset=["month", "year"]).copy()
    plot_df["coverage"] = plot_df["coverage"].fillna(0.0)
    plot_df["month"] = plot_df["month"].astype(int)
    plot_df["year"] = plot_df["year"].astype(int)
    
    value_col = "stratification_score"
    pivot = (
        plot_df.groupby(["month", "year"])[value_col]
        .mean()
        .unstack(level=1)
        .reindex(range(1, 13))
    )
    pivot = pivot.apply(pd.to_numeric, errors="coerce")
    coverage_pivot = (
        plot_df.groupby(["month", "year"])["coverage"]
        .mean()
        .unstack(level=1)
        .reindex(range(1, 13))
        .apply(pd.to_numeric, errors="coerce")
    )

    plot_range = np.linspace(0.8, 12.2, 300)

    month_means = pivot.mean(axis=1, skipna=True)
    mean_values = (
        month_means.interpolate(limit_direction="both")
        .fillna(method="ffill")
        .fillna(method="bfill")
    )
    mean_smoothed = gaussian_filter1d(mean_values.values, sigma=1.1)
    mean_curve = np.interp(plot_range, np.arange(1, 13), mean_smoothed)

    month_std = pivot.std(axis=1, ddof=0).fillna(0.0)
    std_smoothed = gaussian_filter1d(month_std.values, sigma=1.1)
    lower_curve = np.interp(plot_range, np.arange(1, 13), mean_smoothed - std_smoothed)
    upper_curve = np.interp(plot_range, np.arange(1, 13), mean_smoothed + std_smoothed)

    fig, ax = plt.subplots(figsize=(20, 6))

    ax.fill_between(plot_range, lower_curve, upper_curve, color="lightgrey", alpha=0.6, zorder=1)
    ax.plot(plot_range, mean_curve, color="black", linewidth=3, zorder=2)
    extremes_df = extremes_df.copy()
    extremes_df["month"] = extremes_df["date"].map(date_to_month)
    strat_points = extremes_df[extremes_df["extreme_type"] == "max_stratification"]
    mix_points = extremes_df[extremes_df["extreme_type"] == "max_mixing"]

    coverage_counts = {}
    for year in coverage_pivot.columns:
        coverage_counts[year] = (coverage_pivot[year] >= COVERAGE_THRESHOLD).sum()
    min_months_needed = 7
    eligible_years = {year for year, count in coverage_counts.items() if count >= min_months_needed}

    for _, row in strat_points.iterrows():
        if pd.isna(row["month"]) or row.get("year", None) not in eligible_years:
            continue
        month_val = row["month"]
        ax.scatter(
            month_val,
            row[value_col],
            marker="^",
            color="black",
            s=200,
            edgecolor="black",
            linewidth=2.0,
            zorder=5,
        )

    for _, row in mix_points.iterrows():
        if pd.isna(row["month"]) or row.get("year", None) not in eligible_years:
            continue
        month_val = row["month"]
        ax.scatter(
            month_val,
            row[value_col],
            marker="v",
            color="black",
            s=200,
            edgecolor="black",
            linewidth=2.0,
            zorder=5,
        )

    ax.set_xlim(0.9, 12.1)
    ax.set_xticks(range(1, 13))
    ax.set_xticklabels(
        ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"],
        fontsize=12,
    )

    score_min = timeseries_df[value_col].min()
    score_max = timeseries_df[value_col].max()
    y_pad = (score_max - score_min) * 0.1 if score_max != score_min else 0.5
    ax.set_ylim(score_min - y_pad, score_max + y_pad)

    ax.set_xlabel("Month", fontsize=14, fontweight="bold")
    ax.set_ylabel("Depth-centroid distance (D)", fontsize=14, fontweight="bold")
    ax.text(-0.05, 0.5, "", transform=ax.transAxes, fontsize=12, fontweight="bold", rotation=90, va="center")
    ax.grid(axis="y", linestyle="--", alpha=0.35)
    ax.tick_params(axis="y", which="major", pad=8)
    ax.set_title(
        "Monthly Biochemical Stratification Profile",
        fontsize=16,
        fontweight="bold",
    )

    point_df = pivot.stack().reset_index(name=value_col)
    coverage_df = coverage_pivot.stack().reset_index(name="coverage")
    point_df = point_df.merge(coverage_df, on=["month", "year"], how="left")

    for _, row in point_df.iterrows():
        if pd.isna(row[value_col]) or pd.isna(row["coverage"]):
            continue
        if row["coverage"] < COVERAGE_THRESHOLD:
            continue
        color = "black" #"royalblue" if row["normalized_score"] >= 0 else "darkorange"
        ax.scatter(
            row["month"],
            row[value_col],
            color=color,
            s=48,
            edgecolor="black",
            linewidth=0.6,
            alpha=0.9,
            zorder=5,
        )

    plt.tight_layout()
    for suffix in [".pdf", ".png", ".svg"]:
        plt.savefig(Path(output_path).with_suffix(suffix), dpi=300, bbox_inches="tight")
    plt.close()
    print("  [✓] Saved monthly stratification profile")


# ============================================================================
# Main
# ============================================================================

def parse_args():
    ap = argparse.ArgumentParser(
        description="Stratification monthly profile from a single biochem/metadata TSV",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--input", type=Path, required=True, help="Single input TSV containing metadata + biochem features.")
    ap.add_argument("--sep", default=",", help="Delimiter (default: comma). Use '\t' for TSV.")
    ap.add_argument("--sample-id-col", default="cruise_year_month_depth", help="Unique sample ID column.")
    ap.add_argument("--cruise-col", default="Cruise", help="Cruise column name.")
    ap.add_argument("--date-col", default="date", help="Date column name.")
    ap.add_argument("--month-col", default="Month", help="Month column name.")
    ap.add_argument("--year-col", default="Year", help="Year column name.")
    ap.add_argument("--depth-col", default="Depth_anchored", help="Depth column name.")
    ap.add_argument("--depth-min", type=float, default=None, help="Minimum depth to include.")
    ap.add_argument("--depth-max", type=float, default=None, help="Maximum depth to include.")
    ap.add_argument(
        "--pea-metrics",
        type=Path,
        default=None,
        help="Path to stratification_summary.tsv from env_stratification_metrics.py (for PEA comparison plot).",
    )
    ap.add_argument(
        "--pea-date-col",
        default="profile_date",
        help="Date column in PEA metrics file (default profile_date).",
    )

    ap.add_argument(
        "--features",
        default=None,
        help="Comma-separated list of feature columns to use (overrides --features-after-col).",
    )
    ap.add_argument(
        "--features-after-col",
        default="Depth_anchored",
        help="If --features is not provided, use all columns after this column as features.",
    )

    ap.add_argument("--consensus-threshold", type=int, default=1, help="Minimum methods for anomaly consensus [1].")
    ap.add_argument("--output-dir", type=Path, required=True, help="Output directory.")
    return ap.parse_args()


def main():
    args = parse_args()

    out_dir = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    print("\n" + "=" * 70)
    print("STRATIFICATION MONTHLY PROFILE (SINGLE INPUT)")
    print("=" * 70)

    print("\n[1/5] Loading input table...")
    df = pd.read_csv(args.input, sep=args.sep)

    # Validate required columns
    required_cols = [args.sample_id_col, args.date_col, args.month_col, args.year_col, args.depth_col]
    missing_required = [c for c in required_cols if c not in df.columns]
    if missing_required:
        raise ValueError(f"Input missing required columns: {missing_required}")

    # Select features
    feature_cols = select_feature_columns(
        df=df,
        features_csv=args.features,
        features_after_col=args.features_after_col,
        coverage_threshold=COVERAGE_THRESHOLD,
    )

    # Build "integrated_data" and "metadata" to match existing functions
    print("\n[2/5] Preparing matrices...")
    df_dedup = df.drop_duplicates(subset=[args.sample_id_col]).copy()

    # Metadata indexed by sample id
    metadata = df_dedup.set_index(args.sample_id_col)

    # Feature matrix indexed by sample id
    integrated_data = df_dedup.set_index(args.sample_id_col)[feature_cols].copy()
    integrated_data = integrated_data.apply(pd.to_numeric, errors="coerce")

    # Align indices (defensive)
    common = integrated_data.index.intersection(metadata.index)
    integrated_data = integrated_data.loc[common]
    metadata = metadata.loc[common]

    # Optional depth filtering
    if args.depth_min is not None or args.depth_max is not None:
        depth_vals = pd.to_numeric(metadata[args.depth_col], errors="coerce")
        depth_mask = np.isfinite(depth_vals)
        if args.depth_min is not None:
            depth_mask &= depth_vals >= args.depth_min
        if args.depth_max is not None:
            depth_mask &= depth_vals <= args.depth_max
        integrated_data = integrated_data.loc[depth_mask]
        metadata = metadata.loc[depth_mask]
        print(f"  [i] Depth filter kept {len(metadata)} samples")

    print(f"  Samples: {len(common)}")
    print(f"  Features used: {len(feature_cols)}")
    print(f"  Time points: {len(metadata[args.date_col].unique())}")

    print("\n[3/5] Calculating stratification time series...")
    timeseries_df = calculate_stratification_score_timeseries(
        integrated_data=integrated_data,
        metadata=metadata,
        cruise_col=args.cruise_col,
        date_col=args.date_col,
        depth_col=args.depth_col,
        feature_cols=feature_cols,
    )

    print("\n[4/5] Normalizing + anomaly detection...")
    timeseries_df = normalize_to_centered_scale(timeseries_df)
    timeseries_df = detect_anomalies_consensus(timeseries_df, consensus_threshold=args.consensus_threshold)

    print("\n[5/5] Annual extremes + monthly profile...")
    extremes_df = identify_annual_extremes(timeseries_df, metadata, args.year_col)

    # Write the 2 TSVs (rename columns for output)
    rename_map = {
        "stratification_score": "depth_centroid_distance",
        "normalized_score": "normalized_depth_centroid_distance",
    }
    timeseries_out = timeseries_df.rename(columns=rename_map)
    extremes_out = extremes_df.rename(columns=rename_map)
    timeseries_out.to_csv(out_dir / "stratification_timeseries.tsv", sep="\t", index=False)
    extremes_out.to_csv(out_dir / "annual_extremes.tsv", sep="\t", index=False)

    # Write the PDF
    plot_stratification_monthly_profile(
        timeseries_df=timeseries_df,
        extremes_df=extremes_df,
        metadata=metadata,
        date_col=args.date_col,
        month_col=args.month_col,
        year_col=args.year_col,
        output_path=out_dir / "stratification_monthly_profile.pdf",
    )

    # Optional PEA comparison plot
    if args.pea_metrics is not None:
        pea_ts = _load_pea_timeseries(args.pea_metrics, args.pea_date_col)
        plot_stratification_vs_pea_timeseries(
            timeseries_df=timeseries_df,
            pea_df=pea_ts,
            output_path=out_dir / "stratification_vs_pea_timeseries.pdf",
        )
        combined_df = _combine_stratification_and_physical_metrics(timeseries_df, pea_ts)
        combined_out = combined_df.rename(
            columns={
                "stratification_score": "depth_centroid_distance",
                "normalized_score": "normalized_depth_centroid_distance",
            }
        )
        combined_out.to_csv(
            out_dir / "stratification_physical_biochem_timeseries.tsv",
            sep="\t",
            index=False,
        )
        plot_physical_biochem_timeseries(
            combined_df=combined_df,
            output_base=out_dir / "stratification_physical_biochem_timeseries",
            title="Physical stratification and biochemical depth separation — all available data",
        )
        plot_physical_biochem_timeseries(
            combined_df=combined_df,
            output_base=out_dir / "stratification_physical_biochem_timeseries_yearly_month_aligned",
            title="Physical stratification and biochemical depth separation — yearly month-aligned, all data",
            years_per_panel=1,
            align_months=True,
        )
        plot_physical_biochem_monthly_profile(
            combined_df=combined_df,
            output_base=out_dir / "stratification_physical_biochem_monthly_profile_all_data",
        )
        complete_columns = [
            "stratification_score",
            "pea_J_m3",
            "pea_upper_J_m3",
            "pea_lower_J_m3",
        ]
        if all(col in combined_df.columns for col in complete_columns):
            complete_source = combined_df.copy()
            complete_source["_source_order"] = np.arange(len(complete_source))
            complete_df = complete_source.dropna(subset=complete_columns).copy()
            plot_physical_biochem_timeseries(
                combined_df=complete_df,
                output_base=out_dir / "stratification_physical_biochem_timeseries_complete_cases",
                title="Physical stratification and biochemical depth separation — complete cases",
            )
            plot_physical_biochem_timeseries(
                combined_df=complete_df,
                output_base=out_dir / "stratification_physical_biochem_timeseries_yearly_month_aligned_complete_cases",
                title="Physical stratification and biochemical depth separation — yearly month-aligned, complete cases",
                years_per_panel=1,
                align_months=True,
            )
            print(
                f"  [i] Complete-case time series retained {len(complete_df)}/{len(combined_df)} cruises"
            )

    print("\n" + "=" * 70)
    print("DONE")
    print("=" * 70)
    print(f"Outputs saved to: {out_dir}\n")


if __name__ == "__main__":
    main()
