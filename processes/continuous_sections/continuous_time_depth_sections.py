#!/usr/bin/env python3
"""Create auditable continuous time-depth sections from the final BASINS matrix."""

from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, to_hex
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import cKDTree

PROCESS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(PROCESS_DIR.parent))
from shared_plot_export import save_figure_all_formats  # noqa: E402
from shared_plot_style import install_publication_style  # noqa: E402


DEFAULT_VARIABLES = (
    "Oxygen,Nitrate,Nitrite,Nitrous Oxide,Ammonium,Hydrogen Sulfide,Methane,"
    "Phosphate,Silicate,Temperature,Salinity,Density,Iron,Dimethyl Sulfide"
)
# Main low-to-high shading sequence recovered from the ODV rainbow color bar in
# the reference SI interpolation collection. A single colormap object is used
# for continuous measurements and sampled for every compartment surface; this
# deliberately prevents system-specific alternative palettes from appearing in
# otherwise directly comparable time-depth sections.
ODV_RAINBOW_COLORS = (
    "#070268", "#16096D", "#211273", "#2A1A79", "#32217E", "#3A2985",
    "#41308A", "#483790", "#504096", "#57489C", "#5D4FA2", "#5B55A5",
    "#595BA7", "#5660AA", "#5465AC", "#506CAF", "#4C71B2", "#4677B5",
    "#417CB7", "#3A81BA", "#3386BC", "#398CBB", "#4092B8", "#4697B7",
    "#4C9CB4", "#50A1B3", "#55A8B0", "#5AAEAD", "#5EB4AA", "#61B9A8",
    "#64BFA6", "#69C2A5", "#71C5A5", "#78C8A5", "#7FCAA5", "#87CDA5",
    "#8DD0A5", "#94D3A4", "#9AD5A4", "#A0D8A4", "#A6DAA4", "#ACDDA3",
    "#B3E0A2", "#B9E3A1", "#BFE4A0", "#C4E79F", "#CAE99E", "#CFEB9D",
    "#D5EE9B", "#DBF09A", "#E0F299", "#E5F498", "#E8F69B", "#EBF79F",
    "#EEF8A4", "#EFF9A7", "#F2FAAB", "#F4FBAE", "#F7FCB2", "#F9FDB6",
    "#FBFEBA", "#FEFFBE", "#FFFCBA", "#FFF9B6", "#FFF6B0", "#FFF3AB",
    "#FFF0A6", "#FFEDA1", "#FFEB9D", "#FFE797", "#FEE491", "#FEE08C",
    "#FEDD88", "#FED884", "#FED480", "#FECE7C", "#FECA78", "#FEC574",
    "#FEC06F", "#FEBB6B", "#FEB567", "#FDB163", "#FDAC5F", "#FCA65D",
    "#FBA05A", "#FB9956", "#FA9454", "#F98D51", "#F8874D", "#F7804B",
    "#F67A48", "#F57345", "#F46D43", "#F06844", "#EE6445", "#EB6047",
    "#E75A48", "#E4564A", "#E1524B", "#DE4D4C", "#DC494D", "#D8434E",
    "#D53F4F", "#D13A4E", "#CB364D", "#C5314B", "#C02C4A", "#BA2649",
    "#B52147", "#B01C46", "#AB1645", "#A60F44", "#A00642",
)
ODV_RAINBOW_CMAP = LinearSegmentedColormap.from_list(
    "basin_odv_rainbow", ODV_RAINBOW_COLORS, N=256
)
UNITS = {
    "Oxygen": "µM",
    "Nitrate": "µM",
    "Nitrite": "µM",
    "Nitrous Oxide": "nM",
    "Ammonium": "µM",
    "Hydrogen Sulfide": "µM",
    "Methane": "nM",
    "Phosphate": "µM",
    "Silicate": "µM",
    "Temperature": "°C",
    "Salinity": "PSU",
    "Density": r"kg m$^{-3}$",
    "Iron": "nM",
    "Dimethyl Sulfide": "nM",
}
O2_COMPARTMENT_PALETTE = {
    "oxic": "#FF0000",
    "dysoxic": "#008000",
    "suboxic": "#ADD8E6",
    "anoxic": "#800080",
}
HYBRID_COMPARTMENT_PALETTE = {
    "oxic-GMM0": "#990000", "oxic-GMM1": "#D40000", "oxic-GMM2": "#FF0F0F",
    "oxic-GMM3": "#FF4A4A", "oxic-GMM4": "#FF8585",
    "dysoxic-GMM0": "#009900", "dysoxic-GMM1": "#00D400", "dysoxic-GMM2": "#0FFF0F",
    "dysoxic-GMM3": "#4AFF4A", "dysoxic-GMM4": "#85FF85",
    "suboxic-GMM0": "#206379", "suboxic-GMM1": "#2C89A7", "suboxic-GMM2": "#42ABCD",
    "suboxic-GMM3": "#70BFD9", "suboxic-GMM4": "#9ED4E5",
    "anoxic-GMM0": "#990099", "anoxic-GMM1": "#D400D4", "anoxic-GMM2": "#FF0FFF",
    "anoxic-GMM3": "#FF4AFF", "anoxic-GMM4": "#FF85FF",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--variables", default=DEFAULT_VARIABLES)
    p.add_argument("--date-col", default="date")
    p.add_argument("--depth-col", default="Depth_anchored")
    p.add_argument("--renewal-events")
    p.add_argument("--renewal-date-col", default="start_date")
    p.add_argument("--maximum-depth", type=float, default=210.0)
    p.add_argument("--depth-step", type=float, default=1.0)
    p.add_argument("--time-step-days", type=int, default=7)
    p.add_argument("--max-time-support-days", type=float, default=90.0)
    p.add_argument("--max-depth-support-m", type=float, default=30.0)
    p.add_argument("--minimum-samples", type=int, default=4)
    p.add_argument("--levels", type=int, default=64)
    p.add_argument("--formats", default="pdf,png,svg")
    p.add_argument("--o2-assignments")
    p.add_argument("--gmm-assignments")
    p.add_argument("--hybrid-assignments")
    return p.parse_args()


def slug(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_")


def renewal_dates(path: str | None, column: str) -> pd.DatetimeIndex:
    if not path or not Path(path).is_file():
        return pd.DatetimeIndex([])
    table = pd.read_csv(path, sep=None, engine="python")
    if column not in table:
        return pd.DatetimeIndex([])
    return pd.DatetimeIndex(pd.to_datetime(table[column], errors="coerce").dropna().unique()).sort_values()


def interpolate_variable(obs: pd.DataFrame, dates: pd.DatetimeIndex, depths: np.ndarray,
                         max_days: float, max_depth: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    origin = dates.min()
    obs_days = (obs["date"] - origin).dt.total_seconds().to_numpy() / 86400.0
    grid_days = (dates - origin).total_seconds().to_numpy() / 86400.0
    # Median sampling intervals provide a dimensionless geometry without allowing
    # one coordinate's units to dominate the two-dimensional triangulation.
    unique_days = np.unique(obs_days)
    unique_depths = np.unique(obs["depth"].to_numpy(float))
    time_scale = max(float(np.median(np.diff(unique_days))) if len(unique_days) > 1 else 1.0, 1.0)
    depth_scale = max(float(np.median(np.diff(unique_depths))) if len(unique_depths) > 1 else 1.0, 1.0)
    points = np.column_stack((obs_days / time_scale, obs["depth"].to_numpy(float) / depth_scale))
    xx, yy = np.meshgrid(grid_days, depths)
    query = np.column_stack((xx.ravel() / time_scale, yy.ravel() / depth_scale))
    try:
        values = LinearNDInterpolator(points, obs["value"].to_numpy(float), fill_value=np.nan)(query)
    except Exception:
        values = np.full(len(query), np.nan)
    values = np.asarray(values, dtype=float).reshape(xx.shape)

    # Explicit rectangular support limits prevent long temporal or vertical gaps
    # from being painted merely because they fall inside the global convex hull.
    time_dist = np.min(np.abs(grid_days[:, None] - unique_days[None, :]), axis=1)
    depth_dist = np.min(np.abs(depths[:, None] - unique_depths[None, :]), axis=1)
    support = np.isfinite(values) & (time_dist[None, :] <= max_days) & (depth_dist[:, None] <= max_depth)

    # Nearest two-dimensional observation distance is retained for auditing.
    nearest_scaled, _ = cKDTree(points).query(query, k=1)
    nearest_scaled = nearest_scaled.reshape(xx.shape)
    values[~support] = np.nan
    return (
        values,
        support,
        np.broadcast_to(time_dist, xx.shape),
        np.broadcast_to(depth_dist[:, None], xx.shape),
        nearest_scaled,
    )


def nice_color_scale(lo: float, hi: float, target_intervals: int = 7) -> tuple[float, float, np.ndarray, int]:
    """Return round plotting limits and readable 1/2/5-based colorbar ticks."""
    if not np.isfinite(lo) or not np.isfinite(hi):
        return lo, hi, np.array([], dtype=float), 0
    if math.isclose(lo, hi):
        hi = lo + max(abs(lo) * 0.01, 1e-9)
    span = hi - lo
    raw_step = span / max(target_intervals, 1)
    exponent = math.floor(math.log10(raw_step)) if raw_step > 0 else 0
    magnitude = 10.0 ** exponent
    normalized = raw_step / magnitude
    multiplier = next(value for value in (1.0, 2.0, 5.0, 10.0) if normalized <= value)
    step = multiplier * magnitude
    # For ordinary measurement ranges, fractional labels are harder to scan and
    # provide no useful extra precision. Sub-unit ranges retain a suitable scale.
    if span >= 2.0 and step < 1.0:
        step = 1.0
    scale_lo = math.floor(lo / step + 1e-12) * step
    scale_hi = math.ceil(hi / step - 1e-12) * step
    if lo >= 0 and scale_lo < 0:
        scale_lo = 0.0
    ticks = np.arange(scale_lo, scale_hi + step * 0.5, step)
    decimals = max(0, -int(math.floor(math.log10(step)))) if step < 1 else 0
    return scale_lo, scale_hi, ticks, decimals


def format_colorbar(colorbar, ticks: np.ndarray, decimals: int, unit: str) -> None:
    colorbar.set_ticks(ticks)
    colorbar.set_ticklabels([f"{value:.{decimals}f}" for value in ticks])
    colorbar.set_label(unit)


def add_section(ax: plt.Axes, dates: pd.DatetimeIndex, depths: np.ndarray, values: np.ndarray,
                obs: pd.DataFrame, variable: str, renewals: pd.DatetimeIndex, levels: int,
                show_title: bool = True):
    # All panels retain the same complete study extent. Unsupported cells show
    # through as neutral gray rather than being mistaken for zero-valued water.
    ax.set_facecolor("#E6E6E6")
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        ax.text(0.5, 0.5, "Insufficient supported data", ha="center", va="center", transform=ax.transAxes)
        mappable = None
        ticks = np.array([], dtype=float)
        decimals = 0
    else:
        lo, hi = float(np.min(finite)), float(np.max(finite))
        scale_lo, scale_hi, ticks, decimals = nice_color_scale(lo, hi)
        # Put the first contour boundary infinitesimally below the displayed
        # minimum. Matplotlib otherwise leaves values exactly equal to the first
        # boundary (notably measured non-detect zeros) visually unfilled, making
        # them indistinguishable from unsupported/missing grid cells.
        boundary_epsilon = max((scale_hi - scale_lo) * 1e-7, 1e-12)
        level_lo = scale_lo - boundary_epsilon if lo <= scale_lo else scale_lo
        lvls = np.linspace(level_lo, scale_hi, max(8, levels))
        mappable = ax.contourf(dates, depths, values, levels=lvls, cmap=ODV_RAINBOW_CMAP, extend="neither")
    ax.scatter(obs["date"], obs["depth"], s=5, c="black", linewidths=0, zorder=5, rasterized=True)
    for event_date in renewals:
        ax.axvline(event_date, color="black", linestyle="--", linewidth=0.65, alpha=0.8, zorder=4)
    ax.set_xlim(dates.min(), dates.max())
    ax.set_ylim(float(np.max(depths)), 0.0)
    ax.set_ylabel("Depth (m)")
    ax.xaxis.set_major_locator(mdates.YearLocator(2))
    ax.xaxis.set_minor_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    if show_title:
        ax.set_title(variable)
    return mappable, ticks, decimals


def fixed_collection_layout(n_rows: int) -> tuple[plt.Figure, list[plt.Axes], list[plt.Axes]]:
    """Return identical fixed-size main and annotation facets for both collections."""
    if n_rows < 1:
        raise ValueError("A collection must contain at least one row")
    figure_width = 14.5
    row_pitch = 3.8
    panel_height = 3.35
    vertical_margin = (row_pitch - panel_height) / 2.0
    main_left = 0.75
    main_width = 9.50
    side_left = 10.45
    side_width = 3.80
    figure_height = row_pitch * n_rows
    fig = plt.figure(figsize=(figure_width, figure_height))
    main_axes: list[plt.Axes] = []
    side_axes: list[plt.Axes] = []
    shared_axis = None
    for row in range(n_rows):
        bottom_inches = figure_height - vertical_margin - panel_height - row * row_pitch
        bounds = [
            main_left / figure_width,
            bottom_inches / figure_height,
            main_width / figure_width,
            panel_height / figure_height,
        ]
        ax = fig.add_axes(bounds, sharex=shared_axis, sharey=shared_axis)
        if shared_axis is None:
            shared_axis = ax
        side = fig.add_axes([
            side_left / figure_width,
            bottom_inches / figure_height,
            side_width / figure_width,
            panel_height / figure_height,
        ])
        side.set_axis_off()
        main_axes.append(ax)
        side_axes.append(side)
    return fig, main_axes, side_axes


def assignment_palette(
    system: str,
    labels: list[str],
) -> tuple[dict[str, str], str]:
    """Use the exact categorical palettes established by the other BASINS plots."""
    if system == "O2":
        return ({label: O2_COMPARTMENT_PALETTE[label.lower()] for label in labels}, "BASINS O2 palette")
    if system == "GMM":
        # env_compare_compartments.py uses evenly spaced grayscale values from
        # 0.15 to 0.85 in sorted component order, avoiding pure black and white.
        ordered = sorted(labels, key=lambda label: int(re.search(r"(\d+)$", label).group(1)))
        gray_values = np.linspace(0.15, 0.85, len(ordered))
        palette = {
            label: to_hex((gray, gray, gray), keep_alpha=False).upper()
            for label, gray in zip(ordered, gray_values)
        }
        return palette, "BASINS GMM grayscale palette"
    if system == "Hybrid":
        missing = [label for label in labels if label not in HYBRID_COMPARTMENT_PALETTE]
        if missing:
            raise ValueError(f"Hybrid palette lacks observed compartments: {missing}")
        return (
            {label: HYBRID_COMPARTMENT_PALETTE[label] for label in labels},
            "BASINS complete hybrid parent palette",
        )
    raise ValueError(f"Unsupported compartment system: {system}")


def render_compartment_system(path: str, system: str, dates: pd.DatetimeIndex, depths: np.ndarray,
                              args: argparse.Namespace, renewals: pd.DatetimeIndex,
                              plots: Path, tables: Path, formats: tuple[str, ...]):
    source = pd.read_csv(path)
    required = {args.date_col, args.depth_col, "component"}
    if not required.issubset(source.columns):
        raise ValueError(f"{system} assignments lack columns: {sorted(required - set(source.columns))}")
    source[args.date_col] = pd.to_datetime(source[args.date_col], errors="coerce")
    source[args.depth_col] = pd.to_numeric(source[args.depth_col], errors="coerce")
    resp_cols = sorted(
        [c for c in source if re.fullmatch(r"resp_\d+", c)],
        key=lambda c: int(c.split("_")[1]),
    )
    source = source.dropna(subset=[args.date_col, args.depth_col, "component"])
    observed_components = sorted(pd.to_numeric(source["component"], errors="coerce").dropna().astype(int).unique())
    resp_cols = [f"resp_{component}" for component in observed_components if f"resp_{component}" in resp_cols]
    if len(resp_cols) < 2:
        raise ValueError(f"{system} assignments contain fewer than two observed responsibility columns")
    work = source[[args.date_col, args.depth_col, "component", *resp_cols]].copy()
    for column in resp_cols:
        work[column] = pd.to_numeric(work[column], errors="coerce")
    grouped = work.groupby([args.date_col, args.depth_col], as_index=False)[resp_cols].median()
    probability_surfaces = []
    common_support = None
    time_distance = depth_distance = None
    for column in resp_cols:
        obs = grouped[[args.date_col, args.depth_col, column]].dropna().copy()
        obs.columns = ["date", "depth", "value"]
        surface, support, td, dd, _ = interpolate_variable(
            obs, dates, depths, args.max_time_support_days, args.max_depth_support_m
        )
        probability_surfaces.append(surface)
        common_support = support if common_support is None else (common_support & support)
        time_distance, depth_distance = td, dd
    cube = np.stack(probability_surfaces, axis=0)
    cube = np.clip(cube, 0.0, None)
    denominator = np.nansum(cube, axis=0)
    valid = common_support & np.isfinite(denominator) & (denominator > 0)
    cube[:, valid] /= denominator[valid]
    winner_index = np.argmax(np.where(np.isfinite(cube), cube, -np.inf), axis=0)
    winner_index = np.ma.array(winner_index, mask=~valid)
    maximum_probability = np.where(valid, np.max(cube, axis=0), np.nan)

    component_ids = [int(column.split("_")[1]) for column in resp_cols]
    if system == "O2":
        lookup = source.drop_duplicates("component").set_index("component").get("compartment_name", pd.Series(dtype=str)).to_dict()
        labels = [str(lookup.get(component, f"O2-{component}")) for component in component_ids]
    elif system == "Hybrid":
        oxygen_names = ("oxic", "dysoxic", "suboxic", "anoxic")
        labels = [f"{oxygen_names[component // 5]}-GMM{component % 5}" for component in component_ids]
    else:
        labels = [f"GMM{component}" for component in component_ids]
    palette, palette_source = assignment_palette(system, labels)

    fig, ax = plt.subplots(figsize=(10.5, 4.2))
    ax.set_facecolor("#E6E6E6")
    ax.contourf(
        dates, depths, winner_index,
        levels=np.arange(len(labels) + 1) - 0.5,
        cmap=ListedColormap([palette[label] for label in labels]),
        antialiased=True,
    )
    actual_components = pd.to_numeric(source["component"], errors="coerce").astype("Int64")
    component_to_index = {component: index for index, component in enumerate(component_ids)}
    actual_indices = actual_components.map(component_to_index)
    keep = actual_indices.notna()
    ax.scatter(
        source.loc[keep, args.date_col], source.loc[keep, args.depth_col],
        c=[palette[labels[int(index)]] for index in actual_indices[keep]],
        s=10, edgecolors="black", linewidths=0.25, zorder=5,
    )
    for event_date in renewals:
        ax.axvline(event_date, color="black", linestyle="--", linewidth=0.65, alpha=0.8, zorder=4)
    ax.set_xlim(dates.min(), dates.max())
    ax.set_ylim(float(np.max(depths)), 0.0)
    ax.set_xlabel("Sampling date")
    ax.set_ylabel("Depth (m)")
    ax.set_title(f"{system} compartments")
    ax.xaxis.set_major_locator(mdates.YearLocator(2))
    ax.xaxis.set_minor_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax.legend(
        handles=[Patch(facecolor=palette[label], edgecolor="black", linewidth=0.4, label=label) for label in labels],
        loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False, ncol=1 if len(labels) < 8 else 2,
    )
    fig.tight_layout()
    output_slug = system.lower().replace("o2", "o2")
    save_figure_all_formats(fig, plots / f"continuous_time_depth_compartments_{output_slug}", dpi=300, formats=formats)
    plt.close(fig)

    grid = pd.DataFrame({
        "compartment_system": system,
        "date": np.tile(dates.to_numpy(), len(depths)),
        "depth_m": np.repeat(depths, len(dates)),
        "interpolated_compartment": np.asarray([labels[int(x)] if not masked else pd.NA for x, masked in zip(winner_index.data.ravel(), np.ma.getmaskarray(winner_index).ravel())], dtype=object),
        "maximum_interpolated_probability": maximum_probability.ravel(),
        "supported": valid.ravel(),
        "nearest_observed_date_days": time_distance.ravel(),
        "nearest_observed_depth_m": depth_distance.ravel(),
    })
    for index, label in enumerate(labels):
        grid[f"probability_{label}"] = cube[index].ravel()
    grid.to_csv(tables / f"continuous_compartment_grid_{output_slug}.tsv.gz", sep="\t", index=False, compression="gzip")
    pd.DataFrame({
        "component": component_ids,
        "compartment": labels,
        "color_hex": [palette[x] for x in labels],
        "palette_source": palette_source,
    }).to_csv(
        tables / f"continuous_compartment_palette_{output_slug}.tsv", sep="\t", index=False
    )
    return system, labels, palette, winner_index, source, actual_indices


def main() -> None:
    args = parse_args()
    if args.depth_step <= 0 or args.time_step_days <= 0 or args.maximum_depth <= 0:
        raise ValueError("Grid steps and maximum depth must be positive")
    install_publication_style()
    out = Path(args.outdir)
    plots, tables = out / "plots", out / "tables"
    plots.mkdir(parents=True, exist_ok=True)
    tables.mkdir(parents=True, exist_ok=True)
    # Remove the superseded oxygen-ranked color audit from earlier contour
    # implementations so a rerun cannot leave a misleading stale product.
    (tables / "continuous_gmm_posthoc_oxygen_color_order.tsv").unlink(missing_ok=True)
    formats = tuple(x.strip().lower() for x in args.formats.split(",") if x.strip())
    variables = [x.strip() for x in args.variables.split(",") if x.strip()]
    data = pd.read_csv(args.input, sep=None, engine="python")
    missing = [x for x in (args.date_col, args.depth_col) if x not in data]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    data[args.date_col] = pd.to_datetime(data[args.date_col], errors="coerce")
    data[args.depth_col] = pd.to_numeric(data[args.depth_col], errors="coerce")
    data = data.dropna(subset=[args.date_col, args.depth_col])
    dates = pd.date_range(
        data[args.date_col].min().floor("D"),
        data[args.date_col].max().ceil("D"),
        freq=f"{args.time_step_days}D",
    ).union(pd.DatetimeIndex(data[args.date_col].dropna().unique())).sort_values()
    depths = np.arange(0.0, args.maximum_depth + args.depth_step * 0.5, args.depth_step)
    renewals = renewal_dates(args.renewal_events, args.renewal_date_col)
    audits, grids, panels = [], [], []
    for variable in variables:
        if variable not in data:
            audits.append({"variable": variable, "status": "missing_column"})
            continue
        frame = data[[args.date_col, args.depth_col, variable]].copy()
        frame[variable] = pd.to_numeric(frame[variable], errors="coerce")
        frame = frame.dropna().groupby([args.date_col, args.depth_col], as_index=False)[variable].median()
        frame.columns = ["date", "depth", "value"]
        if len(frame) < args.minimum_samples or frame["date"].nunique() < 2 or frame["depth"].nunique() < 2:
            audits.append({"variable": variable, "status": "insufficient_data", "n_observations": len(frame)})
            continue
        values, support, time_distance, depth_distance, nearest_scaled = interpolate_variable(
            frame, dates, depths, args.max_time_support_days, args.max_depth_support_m
        )
        grid = pd.DataFrame({
            "environmental_variable": variable,
            "date": np.tile(dates.to_numpy(), len(depths)),
            "depth_m": np.repeat(depths, len(dates)),
            "interpolated_value": values.ravel(),
            "supported": support.ravel(),
            "nearest_observed_date_days": time_distance.ravel(),
            "nearest_observed_depth_m": depth_distance.ravel(),
            "nearest_observation_scaled_distance": nearest_scaled.ravel(),
        })
        grids.append(grid)
        frame.assign(environmental_variable=variable).to_csv(tables / f"{slug(variable)}_observations.tsv", sep="\t", index=False)
        fig, ax = plt.subplots(figsize=(10.5, 4.2))
        mappable, ticks, decimals = add_section(
            ax, dates, depths, values, frame, variable, renewals, args.levels
        )
        ax.set_xlabel("Sampling date")
        if mappable is not None:
            colorbar = fig.colorbar(mappable, ax=ax, pad=0.015)
            format_colorbar(colorbar, ticks, decimals, UNITS.get(variable, "units not specified"))
        fig.tight_layout()
        save_figure_all_formats(fig, plots / f"continuous_time_depth_{slug(variable)}", dpi=300, formats=formats)
        plt.close(fig)
        panels.append((variable, frame, values))
        audits.append({
            "variable": variable, "status": "rendered", "n_observations": len(frame),
            "n_dates": frame.date.nunique(), "n_depths": frame.depth.nunique(),
            "first_date": frame.date.min().date(), "last_date": frame.date.max().date(),
            "minimum_value": frame.value.min(), "maximum_value": frame.value.max(),
            "n_zero_observations": int(frame.value.eq(0).sum()),
            "n_positive_observations": int(frame.value.gt(0).sum()),
            "display_unit": UNITS.get(variable, "units not specified"),
            "colorbar_minimum": ticks[0] if len(ticks) else np.nan,
            "colorbar_maximum": ticks[-1] if len(ticks) else np.nan,
            "colorbar_tick_interval": ticks[1] - ticks[0] if len(ticks) > 1 else np.nan,
            "supported_grid_fraction": float(support.mean()),
        })
    pd.DataFrame(audits).to_csv(tables / "continuous_section_variable_audit.tsv", sep="\t", index=False)
    if grids:
        pd.concat(grids, ignore_index=True).to_csv(tables / "continuous_time_depth_grids.tsv.gz", sep="\t", index=False, compression="gzip")
    if panels:
        fig, collection_axes, collection_side_axes = fixed_collection_layout(len(panels))
        for row, (ax, side_ax, (variable, frame, values)) in enumerate(
            zip(collection_axes, collection_side_axes, panels)
        ):
            mappable, ticks, decimals = add_section(
                ax, dates, depths, values, frame, variable, renewals, args.levels
            )
            if mappable is not None:
                side_bounds = side_ax.get_position().bounds
                colorbar_ax = fig.add_axes([
                    side_bounds[0] + 0.025,
                    side_bounds[1] + side_bounds[3] * 0.12,
                    0.018,
                    side_bounds[3] * 0.76,
                ])
                colorbar = fig.colorbar(mappable, cax=colorbar_ax)
                format_colorbar(colorbar, ticks, decimals, UNITS.get(variable, "units not specified"))
            ax.tick_params(axis="x", which="both", labelbottom=row == len(panels) - 1)
        collection_axes[-1].set_xlabel("Sampling date")
        save_figure_all_formats(fig, plots / "continuous_time_depth_collection", dpi=300, formats=formats)
        plt.close(fig)
    compartment_inputs = (
        ("O2", args.o2_assignments),
        ("GMM", args.gmm_assignments),
        ("Hybrid", args.hybrid_assignments),
    )
    compartment_manifest = []
    compartment_panels = []
    for system, path in compartment_inputs:
        if not path:
            compartment_manifest.append({"compartment_system": system, "status": "not_configured"})
            continue
        if not Path(path).is_file():
            compartment_manifest.append({"compartment_system": system, "status": "missing_input", "input": path})
            continue
        rendered = render_compartment_system(
            path, system, dates, depths, args, renewals, plots, tables, formats,
        )
        _, labels, palette, winner_index, source, actual_indices = rendered
        compartment_panels.append(rendered)
        compartment_manifest.append({
            "compartment_system": system, "status": "rendered", "input": path,
            "number_of_observed_compartments": len(labels), "number_of_anchor_samples": len(source),
        })
    pd.DataFrame(compartment_manifest).to_csv(
        tables / "continuous_compartment_manifest.tsv", sep="\t", index=False
    )
    if compartment_panels:
        fig, collection_axes, collection_side_axes = fixed_collection_layout(len(compartment_panels))
        for row, (ax, side_ax, (system, labels, palette, winner_index, source, actual_indices)) in enumerate(
            zip(collection_axes, collection_side_axes, compartment_panels)
        ):
            ax.set_facecolor("#E6E6E6")
            ax.contourf(
                dates, depths, winner_index,
                levels=np.arange(len(labels) + 1) - 0.5,
                cmap=ListedColormap([palette[label] for label in labels]),
                antialiased=True,
            )
            keep = actual_indices.notna()
            ax.scatter(
                source.loc[keep, args.date_col], source.loc[keep, args.depth_col],
                c=[palette[labels[int(index)]] for index in actual_indices[keep]],
                s=8, edgecolors="black", linewidths=0.2, zorder=5,
            )
            for event_date in renewals:
                ax.axvline(event_date, color="black", linestyle="--", linewidth=0.65, alpha=0.8, zorder=4)
            ax.set_xlim(dates.min(), dates.max())
            ax.set_ylim(float(np.max(depths)), 0.0)
            ax.set_ylabel("Depth (m)")
            ax.set_title(f"{system} compartments")
            side_ax.legend(
                handles=[Patch(facecolor=palette[label], edgecolor="black", linewidth=0.4, label=label) for label in labels],
                loc="upper left", bbox_to_anchor=(0.0, 1.0), frameon=False,
                ncol=1 if len(labels) < 8 else 2,
            )
            ax.tick_params(axis="x", which="both", labelbottom=row == len(compartment_panels) - 1)
        collection_axes[-1].set_xlabel("Sampling date")
        collection_axes[-1].xaxis.set_major_locator(mdates.YearLocator(2))
        collection_axes[-1].xaxis.set_minor_locator(mdates.YearLocator())
        collection_axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
        save_figure_all_formats(
            fig, plots / "continuous_time_depth_compartment_collection", dpi=300, formats=formats
        )
        plt.close(fig)
    params = vars(args).copy()
    params["variables"] = args.variables
    pd.DataFrame({"parameter": list(params), "value": [str(v) for v in params.values()]}).to_csv(
        tables / "continuous_section_parameters.tsv", sep="\t", index=False
    )


if __name__ == "__main__":
    main()
