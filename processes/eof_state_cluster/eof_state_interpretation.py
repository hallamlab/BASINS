#!/usr/bin/env python3
"""Interpret neutral cruise groups without using interpretation metrics to fit them."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shared_plot_style import install_publication_style

install_publication_style()


def read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t" if path.suffix.lower() in {".tsv", ".txt"} else ",")


def norm_id(series: pd.Series) -> pd.Series:
    def one(value: object) -> str:
        if pd.isna(value):
            return ""
        try:
            number = float(value)
            return str(int(number)) if number.is_integer() else str(value)
        except (TypeError, ValueError):
            return str(value).strip()
    return series.map(one)


def cruise_group_grayscale(groups) -> dict[str, object]:
    """Match the light-to-dark mapping used by season redundancy plots."""
    ordered = sorted(pd.Series(list(groups)).dropna().astype(str).unique())
    shades = plt.cm.Greys(np.linspace(0.35, 0.85, len(ordered)))
    return dict(zip(ordered, shades))


def cruise_group_marker_edge(
    color: object, phase: str = "baseline"
) -> str:
    """Outline light group markers without changing dark-group styling."""
    rgb = np.asarray(matplotlib.colors.to_rgb(color), dtype=float)
    if float(rgb.mean()) >= 0.5:
        return "0.25"
    return "black" if phase != "baseline" else "white"


def save(fig: plt.Figure, base: Path, formats: list[str]) -> None:
    for fmt in formats:
        fig.savefig(base.with_suffix(f".{fmt}"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def renewal_phase(
    frame: pd.DataFrame,
    renewal_col: str = "oxygen_intrusion_class",
    onset_col: str = "oxygen_intrusion_onset",
    event_col: str = "oxygen_intrusion_event_id",
) -> pd.Series:
    """Return nitrate-qualified phases, with legacy O2 fallback."""
    if "renewal_phase" in frame:
        direct = (
            frame["renewal_phase"].fillna("baseline").astype(str).str.strip().str.lower()
            .str.replace("_", "-", regex=False)
        )
        # Unknown nitrate coverage and oxygen-only anomalies are deliberately
        # not shown as renewal states on the primary interpretation figures.
        return direct.where(
            direct.isin({"baseline", "renewal", "post-renewal"}), "baseline"
        )

    phase = pd.Series("baseline", index=frame.index, dtype=object)
    active = frame.get(
        renewal_col, pd.Series("baseline", index=frame.index)
    ).eq("intrusion")
    phase.loc[active] = "post-renewal"

    if onset_col in frame:
        onset_raw = frame[onset_col]
        if pd.api.types.is_bool_dtype(onset_raw):
            onset = onset_raw.fillna(False)
        else:
            onset = (
                onset_raw.astype(str).str.strip().str.lower()
                .isin({"true", "1", "yes", "y"})
            )
        phase.loc[active & onset] = "renewal"
        return phase

    # Older interpreted tables did not retain the explicit onset flag. The
    # first active cruise in each event is equivalent to the source onset flag.
    if event_col in frame:
        dated = frame.loc[active & frame[event_col].notna()].copy()
        dated["_source_index"] = dated.index
        dated["_event_date"] = pd.to_datetime(dated.get("date"), errors="coerce")
        dated = dated.sort_values(
            [event_col, "_event_date", "_source_index"], kind="mergesort"
        )
        onset_indices = dated.groupby(event_col, sort=False).head(1)["_source_index"]
        phase.loc[onset_indices] = "renewal"
        phase.loc[active & frame[event_col].isna()] = "renewal"
    else:
        phase.loc[active] = "renewal"
    return phase


def monthly_group_profile(
    frame: pd.DataFrame, metrics: list[tuple[str, str, str]], group_col: str,
    colors: dict[str, object], base: Path, formats: list[str],
    renewal_col: str = "oxygen_intrusion_class",
) -> None:
    """Plot monthly PEA and D without encoding cruise-group membership."""
    plot = frame.copy()
    plot["date"] = pd.to_datetime(plot["date"], errors="coerce")
    plot = plot.dropna(subset=["date"])
    plot["month"] = plot["date"].dt.month
    plot["_renewal_phase"] = renewal_phase(plot, renewal_col=renewal_col)
    fig, axes = plt.subplots(len(metrics), 1, figsize=(17, 3.8 * len(metrics)), sharex=True, squeeze=False)
    curve_x = np.linspace(.8, 12.2, 300)
    for panel_index, (ax, (column, ylabel, title)) in enumerate(zip(axes[:, 0], metrics)):
        values = plot.dropna(subset=[column]).copy()
        monthly = values.groupby("month")[column].agg(
            median="median", q25=lambda x: x.quantile(.25), q75=lambda x: x.quantile(.75)
        ).reindex(range(1, 13))
        smooth = {}
        for statistic in ("median", "q25", "q75"):
            series = monthly[statistic].interpolate(limit_direction="both")
            smooth[statistic] = np.interp(
                curve_x, np.arange(1, 13), gaussian_filter1d(series.to_numpy(float), sigma=1.1)
            )
        ax.fill_between(curve_x, smooth["q25"], smooth["q75"], color="lightgrey", alpha=.65, zorder=1)
        ax.plot(curve_x, smooth["median"], color="black", linewidth=3, zorder=2)
        ax.axhline(
            values[column].median(),
            color="0.30",
            linewidth=1.6,
            linestyle="--",
            zorder=2,
        )
        renewal = values["_renewal_phase"].eq("renewal")
        post_renewal = values["_renewal_phase"].eq("post-renewal")
        baseline = ~(renewal | post_renewal)
        ax.scatter(
            values.loc[baseline, "month"], values.loc[baseline, column],
            marker="o", s=50, facecolor="0.35", edgecolor="white",
            linewidth=.35, alpha=.82, zorder=3,
        )
        ax.scatter(
            values.loc[renewal, "month"], values.loc[renewal, column],
            marker="s", s=72, facecolor="white", edgecolor="black",
            linewidth=1.1, alpha=1, zorder=4,
        )
        ax.scatter(
            values.loc[post_renewal, "month"], values.loc[post_renewal, column],
            marker="D", s=66, facecolor="white", edgecolor="black",
            linewidth=1.1, alpha=1, zorder=4,
        )
        ax.set(ylabel=ylabel, xlim=(.5, 12.5))
        ax.set_title(
            f"{chr(65 + panel_index)}  {title}", loc="left", fontsize=12,
            fontweight="semibold", pad=7,
        )
        ax.grid(axis="y", linestyle="--", alpha=.3)
    axes[-1, 0].set_xticks(range(1, 13),
                           labels=["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
    axes[-1, 0].set_xlabel("Calendar month")
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    handles = [
        Line2D([], [], marker="o", linestyle="None", markerfacecolor="0.35",
               markeredgecolor="white", color="0.35", label="Cruise observation"),
        Line2D([], [], marker="s", linestyle="None", markerfacecolor="white",
               markeredgecolor="black", markeredgewidth=1.1, color="white",
               label="Renewal"),
        Line2D([], [], marker="D", linestyle="None", markerfacecolor="white",
               markeredgecolor="black", markeredgewidth=1.1, color="white",
               label="Post-renewal"),
        Line2D([], [], color="black", linewidth=3, label="Monthly median"),
        Line2D([], [], color="0.30", linewidth=1.6, linestyle="--",
               label="Overall median"),
        Patch(facecolor="lightgrey", alpha=.65, label="Monthly IQR"),
    ]
    fig.suptitle(
        "Monthly physical and biochemical profiles by nitrate-qualified renewal phase",
        y=.985, fontsize=13,
    )
    fig.legend(
        handles=handles, loc="upper left", bbox_to_anchor=(.805, .92),
        ncol=1, frameon=False,
    )
    fig.tight_layout(rect=[0, 0, .79, .94])
    save(fig, base, formats)


def _monthly_envelope(ax: plt.Axes, values: pd.DataFrame, column: str) -> None:
    curve_x = np.linspace(.8, 12.2, 300)
    monthly = values.groupby("month")[column].agg(
        median="median", q25=lambda x: x.quantile(.25), q75=lambda x: x.quantile(.75)
    ).reindex(range(1, 13))
    smooth = {}
    for statistic in ("median", "q25", "q75"):
        series = monthly[statistic].interpolate(limit_direction="both")
        smooth[statistic] = np.interp(
            curve_x, np.arange(1, 13), gaussian_filter1d(series.to_numpy(float), sigma=1.1)
        )
    ax.fill_between(curve_x, smooth["q25"], smooth["q75"], color="lightgrey", alpha=.6, zorder=1)
    ax.plot(curve_x, smooth["median"], color="black", linewidth=2.5, zorder=2)


def maintext_monthly_profile(
    interpreted: pd.DataFrame,
    all_stratification: pd.DataFrame,
    centroid_col: str,
    group_col: str,
    colors: dict[str, object],
    base: Path,
    formats: list[str],
    renewal_col: str = "oxygen_intrusion_class",
) -> None:
    """Publication-oriented comparison of the two centroid presentations."""
    grouped = interpreted.copy()
    grouped["date"] = pd.to_datetime(grouped["date"], errors="coerce")
    grouped = grouped.dropna(subset=["date"])
    grouped["month"] = grouped["date"].dt.month
    grouped["year"] = grouped["date"].dt.year
    grouped["_renewal_phase"] = renewal_phase(grouped, renewal_col=renewal_col)

    physical = all_stratification.copy()
    date_col = "date" if "date" in physical else "profile_date"
    physical["date"] = pd.to_datetime(physical[date_col], errors="coerce")
    physical = physical.dropna(subset=["date"])
    physical["month"] = physical["date"].dt.month
    physical["year"] = physical["date"].dt.year

    if centroid_col not in grouped or centroid_col not in physical:
        return

    panels = [
        ("All available cruises", False),
        ("Cruise groups and nitrate-qualified renewal phase", True),
    ]
    fig, axes = plt.subplots(2, 1, figsize=(17, 7.6), sharex=True, sharey=True)
    for panel_index, (ax, (title, grouped_panel)) in enumerate(zip(axes, panels)):
        column = centroid_col
        source = grouped if grouped_panel else physical
        values = source.copy()
        values[column] = pd.to_numeric(values[column], errors="coerce")
        values = values.dropna(subset=[column])
        _monthly_envelope(ax, values, column)
        if grouped_panel:
            for group, sub in values.groupby(group_col):
                for phase, marker, size in [
                    ("baseline", "o", 56),
                    ("renewal", "s", 72),
                    ("post-renewal", "D", 66),
                ]:
                    selected = sub["_renewal_phase"].eq(phase)
                    ax.scatter(
                        sub.loc[selected, "month"], sub.loc[selected, column], marker=marker,
                        s=size, facecolor=colors[group],
                        edgecolor=cruise_group_marker_edge(colors[group], phase),
                        linewidth=.9 if phase != "baseline" else .35, alpha=.85, zorder=4,
                    )
        else:
            ax.scatter(values["month"], values[column], color="black", s=18, alpha=.3,
                       edgecolor="none", zorder=3)
            eligibility_values = values
            if "coverage" in values:
                eligibility_values = values[pd.to_numeric(values["coverage"], errors="coerce") >= .51]
            eligible = eligibility_values.groupby("year")["month"].nunique()
            for year, annual in values.groupby("year", sort=True):
                if eligible.get(year, 0) < 7:
                    continue
                maximum = annual.loc[annual[column].idxmax()]
                minimum = annual.loc[annual[column].idxmin()]
                ax.scatter(maximum["month"], maximum[column], marker="^", color="black", s=58, zorder=5)
                ax.scatter(minimum["month"], minimum[column], marker="v", facecolor="white",
                           edgecolor="black", linewidth=1, s=58, zorder=5)
        ax.set_title(f"{chr(65 + panel_index)}  {title}", loc="left", fontsize=12, fontweight="semibold")
        ax.set_ylabel("Depth-centroid distance (D)")
        ax.set_xlim(.5, 12.5)
        ax.grid(axis="y", linestyle="--", alpha=.25)

    month_labels = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    for ax in axes:
        ax.set_xticks(range(1, 13), labels=month_labels)
    axes[-1].set_xlabel("Calendar month")

    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    group_handles = [
        Line2D(
            [], [], marker="o", linestyle="None",
            markerfacecolor=color,
            markeredgecolor=cruise_group_marker_edge(color),
            markeredgewidth=0.35,
            color=color,
            label=group,
        )
        for group, color in colors.items()
    ]
    group_handles += [
        Line2D([], [], marker="s", linestyle="None", markerfacecolor="white",
               markeredgecolor="black", markeredgewidth=.6, color="white",
               label="Renewal"),
        Line2D([], [], marker="D", linestyle="None", markerfacecolor="white",
               markeredgecolor="black", markeredgewidth=.6, color="white",
               label="Post-renewal"),
        Line2D([], [], color="black", linewidth=2.5, label="Monthly median"),
        Patch(facecolor="lightgrey", alpha=.6, label="Monthly IQR"),
    ]
    all_cruise_handles = [
        Line2D([], [], marker="o", linestyle="None", color="0.45", alpha=.45,
               label="All-cruise observation"),
        Line2D([], [], color="black", linewidth=2.5, label="Monthly median"),
        Patch(facecolor="lightgrey", alpha=.6, label="Monthly IQR"),
        Line2D([], [], marker="^", linestyle="None", color="black", label="Annual maximum"),
        Line2D([], [], marker="v", linestyle="None", markerfacecolor="white",
               markeredgecolor="black", color="black", label="Annual minimum"),
    ]
    axes[0].legend(handles=all_cruise_handles, loc="upper left", bbox_to_anchor=(1.02, 1),
                   ncol=1, frameon=False, borderaxespad=0)
    axes[1].legend(handles=group_handles, loc="upper left", bbox_to_anchor=(1.02, 1),
                   ncol=1, frameon=False, borderaxespad=0)
    fig.tight_layout(rect=[0, 0, .79, .98])
    save(fig, base, formats)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--assignments", type=Path, required=True)
    ap.add_argument("--cruise-profiles", type=Path, required=True)
    ap.add_argument("--stratification", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--cruise-col", default="Cruise")
    ap.add_argument("--pea-col", default="pea_J_m3")
    ap.add_argument("--centroid-col", default="depth_centroid_distance")
    ap.add_argument("--formats", default="pdf,png,svg")
    args = ap.parse_args()

    tables = args.outdir / "tables"
    plots = args.outdir / "plots"
    tables.mkdir(parents=True, exist_ok=True)
    plots.mkdir(parents=True, exist_ok=True)
    formats = [x.strip() for x in args.formats.split(",") if x.strip()]

    assignments = read_table(args.assignments)
    profiles = read_table(args.cruise_profiles)
    strat = read_table(args.stratification)
    if "cruise_group" not in assignments.columns:
        source = "component" if "component" in assignments.columns else "state"
        offset = 1 if source == "component" else 0
        assignments["cruise_group"] = [
            f"cruise_group_{int(value) + offset}" for value in assignments[source]
        ]
    for frame in (assignments, profiles, strat):
        frame[args.cruise_col] = norm_id(frame[args.cruise_col])

    keep_strat = [args.cruise_col] + [
        c for c in strat.columns
        if c != args.cruise_col and c not in assignments.columns
    ]
    interpreted = assignments.merge(strat[keep_strat], on=args.cruise_col, how="left")
    interpreted.to_csv(tables / "cruise_group_assignments_interpreted.tsv", sep="\t", index=False)

    group_col = "cruise_group"
    exclude = {args.cruise_col, "state", "component", "Year", "Month", "Day"}
    metrics = [
        c for c in interpreted.columns
        if c not in exclude and pd.to_numeric(interpreted[c], errors="coerce").notna().sum() >= 3
        and not c.startswith(("PC", "resp_", "p_state"))
        and c not in {"max_prob", "resp_entropy", "resp_entropy_normalized", "assignment_uncertain"}
    ]
    summary_rows = []
    for group, sub in interpreted.groupby(group_col, sort=True):
        for metric in metrics:
            values = pd.to_numeric(sub[metric], errors="coerce").dropna().astype(float)
            if values.empty:
                continue
            summary_rows.append({
                group_col: group, "metric": metric, "n": len(values),
                "mean": values.mean(), "sd": values.std(), "median": values.median(),
                "q25": values.quantile(.25), "q75": values.quantile(.75),
            })
    pd.DataFrame(summary_rows).to_csv(tables / "cruise_group_oceanographic_summary.tsv", sep="\t", index=False)

    profile_cols = [c for c in profiles.columns if "@" in c]
    profile = profiles.merge(assignments[[args.cruise_col, group_col]], on=args.cruise_col, how="inner")
    long = profile.melt(id_vars=[args.cruise_col, group_col], value_vars=profile_cols,
                        var_name="feature_depth", value_name="value")
    split = long["feature_depth"].str.rsplit("@", n=1, expand=True)
    long["feature"] = split[0]
    long["depth_m"] = pd.to_numeric(split[1], errors="coerce")
    profile_summary = long.groupby([group_col, "feature", "depth_m"], as_index=False).agg(
        n=("value", "count"), mean=("value", "mean"), median=("value", "median"),
        q25=("value", lambda x: x.quantile(.25)), q75=("value", lambda x: x.quantile(.75)),
    )
    profile_summary.to_csv(tables / "cruise_group_vertical_profiles.tsv", sep="\t", index=False)

    group_values = sorted(interpreted[group_col].dropna().unique())
    colors = {
        group: plt.get_cmap("tab10")(index)
        for index, group in enumerate(group_values)
    }
    grayscale_colors = cruise_group_grayscale(group_values)
    if {"PC1", "PC2"}.issubset(interpreted.columns):
        fig, ax = plt.subplots(figsize=(8, 6))
        for group, sub in interpreted.groupby(group_col):
            ax.scatter(sub.PC1, sub.PC2, label=group, color=colors[group], s=45, alpha=.8)
        uncertain = interpreted.get("assignment_uncertain", pd.Series(False, index=interpreted.index)).astype(bool)
        ax.scatter(interpreted.loc[uncertain, "PC1"], interpreted.loc[uncertain, "PC2"],
                   facecolors="none", edgecolors="black", s=90, linewidths=1.2, label="uncertain")
        ax.set(xlabel="EOF PC1", ylabel="EOF PC2", title="Neutral cruise groups in fitted EOF space")
        ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
        save(fig, plots / "cruise_groups_eof_space", formats)

    if {args.pea_col, args.centroid_col}.issubset(interpreted.columns):
        fig, ax = plt.subplots(figsize=(8, 6))
        for group, sub in interpreted.groupby(group_col):
            ax.scatter(sub[args.pea_col], sub[args.centroid_col], label=group,
                       color=colors[group], s=48, alpha=.8)
        ax.set(xlabel="PEA (J m$^{-3}$)", ylabel="Biochemical depth-centroid distance (D)",
               title="Post hoc physical and biochemical interpretation")
        ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
        save(fig, plots / "cruise_groups_pea_vs_depth_centroid", formats)

        dated = interpreted.copy()
        date_col = "date" if "date" in dated.columns else None
        if date_col:
            dated[date_col] = pd.to_datetime(dated[date_col], errors="coerce")
            dated = dated.dropna(subset=[date_col]).sort_values(date_col)
            fig, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=True)
            for ax, metric, ylabel in zip(
                axes, [args.pea_col, args.centroid_col],
                ["PEA (J m$^{-3}$)", "Depth-centroid distance (D)"],
            ):
                ax.plot(dated[date_col], dated[metric], color="0.75", linewidth=1, zorder=0)
                for group, sub in dated.groupby(group_col):
                    ax.scatter(sub[date_col], sub[metric], color=colors[group], label=group, s=42)
                ax.set_ylabel(ylabel)
                ax.grid(axis="y", alpha=.2)
            axes[0].set_title("Neutral cruise groups overlaid on interpretation metrics")
            axes[0].legend(frameon=False, ncol=max(1, len(colors)), fontsize=8)
            axes[-1].set_xlabel("Cruise date")
            save(fig, plots / "cruise_groups_physical_biochem_timeseries", formats)

        monthly_group_profile(
            interpreted,
            [
                (
                    args.centroid_col,
                    "Depth-centroid distance (D)",
                    "Biochemical depth-centroid distance",
                ),
                (args.pea_col, "PEA (J m$^{-3}$)", "Potential energy anomaly"),
            ],
            group_col, colors,
            plots / "cruise_groups_physical_biochem_monthly_profile_all_data", formats,
        )
        maintext_monthly_profile(
            interpreted=interpreted,
            all_stratification=strat,
            centroid_col=args.centroid_col,
            group_col=group_col,
            colors=grayscale_colors,
            base=plots / "cruise_groups_physical_biochem_monthly_profile_maintext",
            formats=formats,
        )

    z = profile[profile_cols].apply(pd.to_numeric, errors="coerce")
    z = (z - z.mean()) / z.std(ddof=0).replace(0, np.nan)
    means = z.assign(**{group_col: profile[group_col]}).groupby(group_col).mean()
    fig, ax = plt.subplots(figsize=(18, max(3, .7 * len(means))))
    image = ax.imshow(means.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.set_yticks(range(len(means)), labels=means.index)
    step = max(1, len(profile_cols) // 25)
    ticks = np.arange(0, len(profile_cols), step)
    ax.set_xticks(ticks, labels=[profile_cols[i] for i in ticks], rotation=90, fontsize=7)
    ax.set_title("Standardized mean water-column profile by cruise group")
    fig.colorbar(image, ax=ax, label="Standardized profile value")
    save(fig, plots / "cruise_group_profile_heatmap", formats)


if __name__ == "__main__":
    main()
