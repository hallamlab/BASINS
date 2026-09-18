# File: BASINS/processes/split_o2_by_gmm/env_split_o2_by_gmm.py
#
# What this does
# --------------
# 1) Compute O2 compartments from matrix_cleaned (oxic/dysoxic/suboxic/anoxic).
# 2) Define subcompartments as the intersection: O2 compartment × GMM component.
# 3) Optionally collapse tiny intersections to "<o2>__other" (reporting convenience).
# 4) Optionally reassign "borderline" samples (typically in "<o2>__other") into the nearest
#    *core* subcompartment within the SAME O2 compartment using standardized PC space centroids,
#    with a conservative radius gate derived from within-core dispersion.
#
# Inputs
# ------
# - matrix_cleaned.csv : Oxygen + metadata (for O2 compartment labeling)
# - eigenvectors_scores.csv : PC columns used for centroid distance calculations
# - compartments_assignments_smoothed.csv : GMM component labels + max_prob (optional)
#
# Outputs
# -------
# --outdir/
#   run_config.json
#   tables/
#     merged_o2_split_by_gmm.csv
#     o2_subcompartment_counts_before.csv
#     o2_subcompartment_counts_after.csv
#     o2_by_gmm_confusion_raw.csv
#     o2_by_gmm_confusion_row_norm.csv
#     o2_by_gmm_confusion_col_norm.csv
#     reassignment_qc_summary.csv
#     reassignment_centroids.csv
#     reassignment_cluster_radii.csv
#
# Notes
# -----
# - No re-fitting of clustering models. This is post-hoc hierarchical labeling only.
# - Reassignment is constrained strictly within each O2 compartment.
# - Reassignment is gated by an empirical radius: the q-quantile of within-core distances.
#
# Extensions in this script compared to env_split_o2_by_gmm.py
# Adds to the previous script:
#  1) Summary plots per O2 compartment (PC-space + biochem overlays)
#  2) Depth-profile visualizations colored by subcompartment (y-axis inverted: shallow at top)
#  3) Within-O2 silhouettes (unweighted + weighted by max_prob) computed in PC space
#
# Usage: see example command at bottom of file docstring.

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from dataclasses import dataclass
from typing import List, Optional, Tuple, Dict

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from scipy.ndimage import gaussian_filter

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shared_plot_style import install_publication_style

install_publication_style()

from sklearn.metrics import silhouette_score
from sklearn.metrics import pairwise_distances


# ----------------------------
# Palettes (user-specified)
# ----------------------------

O2_COMPARTMENT_PALETTE = {
    "oxic": "red",
    "dysoxic": "green",
    "suboxic": "lightblue",
    "anoxic": "purple",
}

BIOCHEM_COLOR_MAP = {
    "Oxygen": "black",
    "Nitrogen Oxides": "#E7298A",
    "Nitrate": "#1B9E77",
    "Nitrite": "#66A61E",
    "Nitrous Oxide": "#0C5196",
    "Ammonium": "#7570B3",
    "Hydrogen Sulfide": "#D95F02",
    "Methane": "violet",
}


# ----------------------------
# Config
# ----------------------------

@dataclass
class Config:
    matrix_cleaned: str
    eigenvectors: str
    umap_embedding: str
    assignments: str
    o2_assignments: Optional[str]
    outdir: str
    sep_matrix: str
    sep_eig: str
    sep_assign: str
    sep_o2_assign: str

    # keying
    id_col: str
    key_mode: str
    key_cols: List[str]
    key_sep: str
    derived_key_col: str

    # core cols
    oxygen_col: str
    o2_compartment_col: str
    cruise_col: str
    depth_col: str
    depth_anchored_col: str
    date_col: str

    # O2 thresholds (µM)
    o2_oxic_gt: float
    o2_dysoxic_hi: float
    o2_dysoxic_lo: float
    o2_suboxic_hi: float
    o2_suboxic_lo: float

    # GMM columns
    gmm_component_col: str
    max_prob_col: str

    # PC columns
    pc_cols: List[str]

    # labeling
    sub_label_sep: str
    prefix_gmm: str

    # collapse tiny intersections
    min_subcluster_size: int

    # reassignment
    do_reassign: bool
    borderline_mode: str
    borderline_max_prob: Optional[float]
    core_min_prob: Optional[float]
    reassign_radius_quantile: float
    reassign_min_core_n: int

    # plotting
    do_plots: bool
    plot_formats: List[str]
    png_dpi: int
    point_size: float
    alpha: float

    # time-depth curtain (display only; assignments are not changed)
    curtain_hide_other: bool
    curtain_time_subdivisions_per_month: int
    curtain_time_sigma_months: float
    curtain_depth_step_m: float
    curtain_maximum_depth_m: float
    curtain_contour_visual_depth_sigma_m: float
    curtain_renewal_events: Optional[str]
    curtain_renewal_date_col: str


def parse_args() -> Config:
    ap = argparse.ArgumentParser(
        description="Hierarchical labeling: O2 compartments split by GMM components with optional centroid reassignment, plots, and within-O2 silhouettes."
    )

    ap.add_argument("--matrix-cleaned", required=True)
    ap.add_argument("--eigenvectors", required=True)
    ap.add_argument("--umap-embedding", required=True)
    ap.add_argument("--assignments", required=True)
    ap.add_argument(
        "--o2-assignments",
        default=None,
        help="Optional path to o2_compartments_assignments_{base|smoothed}.csv. If provided, uses these O2 labels instead of thresholding Oxygen.",
    )
    ap.add_argument("--outdir", required=True)

    ap.add_argument("--sep-matrix", default=",")
    ap.add_argument("--sep-eig", default=",")
    ap.add_argument("--sep-assign", default=",")
    ap.add_argument("--sep-o2-assign", default=",")

    ap.add_argument("--id-col", default="cruise_year_month_depth")
    ap.add_argument("--key-mode", choices=["composite", "id"], default="composite")
    ap.add_argument("--key-cols", default="Cruise,Year,Month,Day,Depth")
    ap.add_argument("--key-sep", default="|")

    ap.add_argument("--oxygen-col", default="Oxygen")
    ap.add_argument(
        "--o2-compartment-col",
        default="compartment_name",
        help="Column in --o2-assignments containing O2 compartment labels (default compartment_name).",
    )
    ap.add_argument("--cruise-col", default="Cruise")
    ap.add_argument("--depth-col", default="Depth")
    ap.add_argument("--depth-anchored-col", default="Depth_anchored")
    ap.add_argument("--date-col", default="date")

    ap.add_argument("--o2-oxic-gt", type=float, default=90.0)
    ap.add_argument("--o2-dysoxic-hi", type=float, default=90.0)
    ap.add_argument("--o2-dysoxic-lo", type=float, default=20.0)
    ap.add_argument("--o2-suboxic-hi", type=float, default=20.0)
    ap.add_argument("--o2-suboxic-lo", type=float, default=1.0)

    ap.add_argument("--gmm-component-col", default="component")
    ap.add_argument("--max-prob-col", default="max_prob")

    ap.add_argument("--pc-cols", default="PC1,PC2,PC3")

    ap.add_argument("--sub-label-sep", default="__")
    ap.add_argument("--prefix-gmm", default="gmm")

    ap.add_argument("--min-subcluster-size", type=int, default=10)

    ap.add_argument("--reassign", action="store_true")
    ap.add_argument(
        "--borderline-mode",
        choices=["other_only", "low_conf_only", "other_or_low_conf"],
        default="other_only",
    )
    ap.add_argument("--borderline-max-prob", type=float, default=None)
    ap.add_argument("--core-min-prob", type=float, default=None)
    ap.add_argument("--reassign-radius-quantile", type=float, default=0.90)
    ap.add_argument("--reassign-min-core-n", type=int, default=20)

    ap.add_argument("--plots", action="store_true", help="Write summary plots per O2 compartment.")
    ap.add_argument(
        "--plot-formats",
        default="png,pdf,svg",
        help="Comma-separated list: png,pdf,svg (default: all three).",
    )
    ap.add_argument("--png-dpi", type=int, default=300)
    ap.add_argument("--point-size", type=float, default=18.0)
    ap.add_argument("--alpha", type=float, default=0.65)
    ap.add_argument(
        "--curtain-hide-other", action="store_true",
        help="Exclude '<oxygen>__other' states from the curtain support surface and legend.",
    )
    ap.add_argument("--curtain-time-subdivisions-per-month", type=int, default=4)
    ap.add_argument("--curtain-time-sigma-months", type=float, default=0.75)
    ap.add_argument("--curtain-depth-step-m", type=float, default=1.0)
    ap.add_argument("--curtain-maximum-depth-m", type=float, default=210.0)
    ap.add_argument("--curtain-contour-visual-depth-sigma-m", type=float, default=5.0)
    ap.add_argument(
        "--curtain-renewal-events", default=None,
        help=(
            "Optional nitrate-qualified renewal-event table. Only event-onset "
            "dates are drawn as vertical dashed lines on the curtain."
        ),
    )
    ap.add_argument("--curtain-renewal-date-col", default="start_date")

    ns = ap.parse_args()
    key_cols = [c.strip() for c in ns.key_cols.split(",") if c.strip()]
    pc_cols = [c.strip() for c in ns.pc_cols.split(",") if c.strip()]
    plot_formats = [f.strip().lower() for f in ns.plot_formats.split(",") if f.strip()]

    return Config(
        matrix_cleaned=ns.matrix_cleaned,
        eigenvectors=ns.eigenvectors,
        assignments=ns.assignments,
        o2_assignments=ns.o2_assignments,
        umap_embedding=ns.umap_embedding,
        outdir=ns.outdir,
        sep_matrix=ns.sep_matrix,
        sep_eig=ns.sep_eig,
        sep_assign=ns.sep_assign,
        sep_o2_assign=ns.sep_o2_assign,
        id_col=ns.id_col,
        key_mode=ns.key_mode,
        key_cols=key_cols,
        key_sep=ns.key_sep,
        derived_key_col="__merge_key__",
        oxygen_col=ns.oxygen_col,
        o2_compartment_col=ns.o2_compartment_col,
        cruise_col=ns.cruise_col,
        depth_col=ns.depth_col,
        depth_anchored_col=ns.depth_anchored_col,
        date_col=ns.date_col,
        o2_oxic_gt=ns.o2_oxic_gt,
        o2_dysoxic_hi=ns.o2_dysoxic_hi,
        o2_dysoxic_lo=ns.o2_dysoxic_lo,
        o2_suboxic_hi=ns.o2_suboxic_hi,
        o2_suboxic_lo=ns.o2_suboxic_lo,
        gmm_component_col=ns.gmm_component_col,
        max_prob_col=ns.max_prob_col,
        pc_cols=pc_cols,
        sub_label_sep=ns.sub_label_sep,
        prefix_gmm=ns.prefix_gmm,
        min_subcluster_size=int(ns.min_subcluster_size),
        do_reassign=bool(ns.reassign),
        borderline_mode=str(ns.borderline_mode),
        borderline_max_prob=ns.borderline_max_prob,
        core_min_prob=ns.core_min_prob,
        reassign_radius_quantile=float(ns.reassign_radius_quantile),
        reassign_min_core_n=int(ns.reassign_min_core_n),
        do_plots=bool(ns.plots),
        plot_formats=plot_formats,
        png_dpi=int(ns.png_dpi),
        point_size=float(ns.point_size),
        alpha=float(ns.alpha),
        curtain_hide_other=bool(ns.curtain_hide_other),
        curtain_time_subdivisions_per_month=max(1, int(ns.curtain_time_subdivisions_per_month)),
        curtain_time_sigma_months=max(0.0, float(ns.curtain_time_sigma_months)),
        curtain_depth_step_m=max(0.1, float(ns.curtain_depth_step_m)),
        curtain_maximum_depth_m=max(0.1, float(ns.curtain_maximum_depth_m)),
        curtain_contour_visual_depth_sigma_m=max(
            0.0, float(ns.curtain_contour_visual_depth_sigma_m)
        ),
        curtain_renewal_events=ns.curtain_renewal_events,
        curtain_renewal_date_col=str(ns.curtain_renewal_date_col),
    )


# ----------------------------
# IO / helpers
# ----------------------------

def ensure_dirs(outdir: str) -> Tuple[str, str]:
    tables = os.path.join(outdir, "tables")
    plots = os.path.join(outdir, "plots")
    os.makedirs(tables, exist_ok=True)
    os.makedirs(plots, exist_ok=True)
    return tables, plots


def read_table_dedup_cols(path: str, sep: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep=sep)
    if df.columns.duplicated().any():
        df = df.loc[:, ~df.columns.duplicated()].copy()
    return df

def read_umap_embedding(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    required = {"__merge_key__", "UMAP1", "UMAP2"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"UMAP embedding missing columns: {missing}")
    return df

def _normalize_str_series(s: pd.Series) -> pd.Series:
    return s.astype("object").fillna("NA").astype(str).str.strip()

def _series_equal_tolerant(a: pd.Series, b: pd.Series, rtol=1e-6, atol=1e-9) -> bool:
    # Try numeric compare first
    an = pd.to_numeric(a, errors="coerce")
    bn = pd.to_numeric(b, errors="coerce")

    both_numeric = an.notna() & bn.notna()
    both_nan = an.isna() & bn.isna()

    # For numeric-overlap rows, compare with tolerance
    if both_numeric.any():
        if not np.allclose(an[both_numeric].to_numpy(), bn[both_numeric].to_numpy(), rtol=rtol, atol=atol):
            return False

    # For rows not both numeric, compare normalized strings
    other = ~(both_numeric | both_nan)
    if other.any():
        if not _normalize_str_series(a[other]).equals(_normalize_str_series(b[other])):
            return False

    return True

def coalesce_merge_suffix_columns(
    df: pd.DataFrame,
    suffixes=("_x", "_y"),
    prefer="x",
    rtol=1e-6,
    atol=1e-9,
    fail_on_mismatch=True,
) -> pd.DataFrame:
    """
    Coalesce columns produced by pandas merge suffixes.
    For each base col where base_x and base_y exist:
      - if identical (tolerant), keep single base col (prefer x or y)
      - else: raise (default) or keep both.
    """
    out = df.copy()
    sx, sy = suffixes

    # find candidate bases
    bases = []
    for c in out.columns:
        if c.endswith(sx):
            base = c[: -len(sx)]
            if (base + sy) in out.columns:
                bases.append(base)

    for base in bases:
        cx = base + sx
        cy = base + sy

        equal = _series_equal_tolerant(out[cx], out[cy], rtol=rtol, atol=atol)
        if not equal:
            if fail_on_mismatch:
                # show a few example mismatches to debug quickly
                ax = out[cx]
                ay = out[cy]
                # build a mismatch mask (string compare fallback)
                nx = pd.to_numeric(ax, errors="coerce")
                ny = pd.to_numeric(ay, errors="coerce")
                both_num = nx.notna() & ny.notna()
                mm = pd.Series(False, index=out.index)
                mm[both_num] = ~np.isclose(nx[both_num], ny[both_num], rtol=rtol, atol=atol)
                other = ~both_num & ~(nx.isna() & ny.isna())
                mm[other] = _normalize_str_series(ax[other]) != _normalize_str_series(ay[other])
                examples = out.loc[mm, [cx, cy]].head(10)
                raise ValueError(
                    f"Merge produced non-identical duplicate columns for base='{base}': {cx} vs {cy}\n"
                    f"Examples (first 10 mismatches):\n{examples.to_string(index=False)}"
                )
            else:
                # keep both; continue
                continue

        keep_col = cx if prefer == "x" else cy
        out[base] = out[keep_col]
        out = out.drop(columns=[cx, cy])

    return out

def build_merge_key(df: pd.DataFrame, cfg: Config) -> pd.DataFrame:
    out = df.copy()
    if cfg.key_mode == "id":
        if cfg.id_col not in out.columns:
            raise ValueError(f"key-mode=id but id-col not found: {cfg.id_col}")
        out[cfg.derived_key_col] = out[cfg.id_col].astype(str)
        return out

    missing = [c for c in cfg.key_cols if c not in out.columns]
    if missing:
        raise ValueError(f"key-mode=composite but missing columns: {missing}")

    parts = [out[c].astype(str).fillna("NA") for c in cfg.key_cols]
    key = parts[0]
    for p in parts[1:]:
        key = key + cfg.key_sep + p
    out[cfg.derived_key_col] = key
    return out


def label_o2_compartment(o2_uM: pd.Series, cfg: Config) -> pd.Series:
    x = pd.to_numeric(o2_uM, errors="coerce")
    out = pd.Series(["NA"] * len(x), index=x.index, dtype="object")
    out[x > cfg.o2_oxic_gt] = "oxic"
    out[(x <= cfg.o2_dysoxic_hi) & (x >= cfg.o2_dysoxic_lo)] = "dysoxic"
    out[(x < cfg.o2_suboxic_hi) & (x >= cfg.o2_suboxic_lo)] = "suboxic"
    out[x < cfg.o2_suboxic_lo] = "anoxic"
    return out


def normalize_o2_labels(raw: pd.Series) -> pd.Series:
    out = pd.Series(["NA"] * len(raw), index=raw.index, dtype="object")

    # numeric map
    num = pd.to_numeric(raw, errors="coerce")
    num_map = {0: "oxic", 1: "dysoxic", 2: "suboxic", 3: "anoxic"}
    for k, v in num_map.items():
        out[num == float(k)] = v

    # text map
    txt = raw.astype("object").fillna("NA").astype(str).str.strip().str.lower()
    txt_map = {
        "oxic": "oxic",
        "dysoxic": "dysoxic",
        "suboxic": "suboxic",
        "anoxic": "anoxic",
    }
    mapped_txt = txt.map(txt_map)
    out[mapped_txt.notna()] = mapped_txt[mapped_txt.notna()]
    return out


def o2_labels_from_assignments(df_o2: pd.DataFrame, cfg: Config) -> pd.Series:
    candidates = [cfg.o2_compartment_col, "compartment_name", "o2_compartment", "component"]
    col = next((c for c in candidates if c in df_o2.columns), None)
    if col is None:
        raise ValueError(
            "Could not find an O2 compartment column in --o2-assignments. "
            f"Tried: {candidates}"
        )
    return normalize_o2_labels(df_o2[col])


def confusion_tables(y_true: pd.Series, y_pred: pd.Series) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    yt = y_true.astype("object").fillna("NA")
    yp = y_pred.astype("object").fillna("NA")
    raw = pd.crosstab(yt, yp, rownames=["O2_compartment"], colnames=["GMM_component"], dropna=False)
    row_norm = raw.div(raw.sum(axis=1).replace(0, np.nan), axis=0)
    col_norm = raw.div(raw.sum(axis=0).replace(0, np.nan), axis=1)
    return raw, row_norm, col_norm


def standardize_pc_space(df: pd.DataFrame, pc_cols: List[str]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    X = df[pc_cols].to_numpy(dtype=float)
    mu = np.nanmean(X, axis=0)
    sd = np.nanstd(X, axis=0, ddof=0)
    sd = np.where(sd == 0, 1.0, sd)
    Xz = (X - mu) / sd
    return Xz, mu, sd


def _make_sub_label(o2: str, gmm: str, cfg: Config) -> str:
    return f"{o2}{cfg.sub_label_sep}{cfg.prefix_gmm}{gmm}"


# ----------------------------
# Color utilities
# ----------------------------

def _hex_from_rgb01(rgb: Tuple[float, float, float]) -> str:
    r, g, b = rgb
    r = int(np.clip(round(r * 255), 0, 255))
    g = int(np.clip(round(g * 255), 0, 255))
    b = int(np.clip(round(b * 255), 0, 255))
    return f"#{r:02X}{g:02X}{b:02X}"


def _rgb01_from_any(color: str) -> Tuple[float, float, float]:
    # Let matplotlib do the parsing (named colors + hex)
    import matplotlib.colors as mcolors
    rgb = mcolors.to_rgb(color)
    return float(rgb[0]), float(rgb[1]), float(rgb[2])


def _mix(rgb_a: Tuple[float, float, float], rgb_b: Tuple[float, float, float], t: float) -> Tuple[float, float, float]:
    # linear interpolation
    t = float(np.clip(t, 0.0, 1.0))
    return (rgb_a[0] * (1 - t) + rgb_b[0] * t,
            rgb_a[1] * (1 - t) + rgb_b[1] * t,
            rgb_a[2] * (1 - t) + rgb_b[2] * t)


import colorsys

def _parse_gmm_index_from_label(label: str, cfg: Config) -> Optional[int]:
    """
    Expected labels like: '<o2>__gmm0', '<o2>__gmm12', '<o2>__other'
    Returns integer index if parseable, else None.
    """
    s = str(label)
    if cfg.sub_label_sep not in s:
        return None
    tok = s.split(cfg.sub_label_sep, 1)[1]  # e.g. 'gmm3' or 'other'
    if tok.startswith(cfg.prefix_gmm):
        tail = tok[len(cfg.prefix_gmm):]
        try:
            return int(tail)
        except Exception:
            return None
    return None

def _rgb01(color: str) -> Tuple[float, float, float]:
    import matplotlib.colors as mcolors
    return mcolors.to_rgb(color)

def _rgb01_to_hex(rgb: Tuple[float, float, float]) -> str:
    r, g, b = rgb
    r = int(np.clip(round(r * 255), 0, 255))
    g = int(np.clip(round(g * 255), 0, 255))
    b = int(np.clip(round(b * 255), 0, 255))
    return f"#{r:02X}{g:02X}{b:02X}"

def _variant_hls(base_color: str, t: float, sat_floor: float = 0.50) -> str:
    """
    Hue-locked variants (stays within O2 color family):
      - Hue fixed
      - Saturation fixed but floored (dysoxic can be darker and not washed out)
      - Lightness ramps from dark -> light across t in [0,1]
    """
    r, g, b = _rgb01(base_color)
    h, l, s = colorsys.rgb_to_hls(r, g, b)

    s2 = max(s, sat_floor)

    # darker + lighter bounds around base lightness
    l_dark = max(0.08, l * 0.40)  # darker (dysoxic can be darker)
    l_light = min(0.92, 1 - (1 - l) * 0.30)

    t = float(np.clip(t, 0.0, 1.0))
    l2 = (1 - t) * l_dark + t * l_light

    r2, g2, b2 = colorsys.hls_to_rgb(h, l2, s2)
    return _rgb01_to_hex((r2, g2, b2))

def build_full_subcompartment_palette(
    m: pd.DataFrame,
    cfg: Config,
    sat_floor: float = 0.50,
) -> pd.DataFrame:
    """
    Build a COMPLETE, deterministic palette for every:
      '<o2>__gmm0' ... '<o2>__gmmMax' plus '<o2>__other'

    gmm0 is always darkest; gmmMax always lightest.
    Returns a dataframe with columns: o2_compartment, label, gmm_index, color_hex
    """
    rows = []

    if "o2_compartment" not in m.columns:
        return pd.DataFrame(columns=["o2_compartment", "label", "gmm_index", "color_hex"])

    # Determine max gmm index per O2 from labels if present; fallback to gmm_component column if needed.
    for o2 in sorted(m["o2_compartment"].astype(str).unique(), key=str):
        base = O2_COMPARTMENT_PALETTE.get(o2, "gray")

        # collect gmm indices observed in final labels
        observed = (
            m.loc[m["o2_compartment"].astype(str) == o2, "o2_subcompartment_final"]
            if "o2_subcompartment_final" in m.columns
            else pd.Series([], dtype="object")
        )
        idxs = []
        for lab in observed.astype(str).unique().tolist():
            gi = _parse_gmm_index_from_label(lab, cfg)
            if gi is not None:
                idxs.append(gi)

        # fallback: if none found, try gmm_component numeric parse (optional)
        if len(idxs) == 0 and "gmm_component" in m.columns:
            tmp = m.loc[m["o2_compartment"].astype(str) == o2, "gmm_component"].astype(str).unique().tolist()
            for t in tmp:
                try:
                    idxs.append(int(t))
                except Exception:
                    pass

        max_g = int(max(idxs)) if len(idxs) > 0 else 0

        # build labels gmm0..gmmMax
        if max_g < 0:
            max_g = 0
        gmm_labels = [f"{o2}{cfg.sub_label_sep}{cfg.prefix_gmm}{k}" for k in range(0, max_g + 1)]

        # color ramp across gmm labels
        n = max(len(gmm_labels), 1)
        ts = np.linspace(0.0, 1.0, num=n).tolist()  # gmm0 darkest -> gmmMax lightest
        for lab, t in zip(gmm_labels, ts):
            rows.append(
                {
                    "o2_compartment": o2,
                    "label": lab,
                    "gmm_index": int(_parse_gmm_index_from_label(lab, cfg)),
                    "color_hex": _variant_hls(base, float(t), sat_floor=sat_floor),
                }
            )

        # add 'other' as lightest variant (t=1.0)
        other_lab = f"{o2}{cfg.sub_label_sep}other"
        rows.append(
            {
                "o2_compartment": o2,
                "label": other_lab,
                "gmm_index": np.nan,
                "color_hex": "#FFFFFF",
            }
        )

    return pd.DataFrame(rows)

def palette_df_to_dict(pal_df: pd.DataFrame) -> Dict[str, str]:
    return {str(r["label"]): str(r["color_hex"]) for _, r in pal_df.iterrows()}


# ----------------------------
# Intersection collapse + reassignment selection
# ----------------------------

def collapse_small_intersections(m: pd.DataFrame, cfg: Config) -> pd.DataFrame:
    counts = m["o2_subcompartment"].value_counts(dropna=False)
    small = set(counts[counts < cfg.min_subcluster_size].index.astype(str).tolist())
    if len(small) == 0:
        return m

    def _collapse(lbl: str) -> str:
        if lbl not in small:
            return lbl
        o2 = str(lbl).split(cfg.sub_label_sep)[0]
        return f"{o2}{cfg.sub_label_sep}other"

    out = m.copy()
    out["o2_subcompartment"] = out["o2_subcompartment"].astype(str).map(_collapse)
    return out


def determine_borderline_mask(m: pd.DataFrame, cfg: Config) -> pd.Series:
    is_other = m["o2_subcompartment"].astype(str).str.endswith(f"{cfg.sub_label_sep}other")

    if cfg.borderline_mode == "other_only":
        return is_other

    if cfg.borderline_max_prob is None:
        raise ValueError("--borderline-mode uses low_conf but --borderline-max-prob was not provided.")
    if cfg.max_prob_col not in m.columns:
        raise ValueError("--borderline-mode uses low_conf but max_prob column not present in merged table.")
    mp = pd.to_numeric(m[cfg.max_prob_col], errors="coerce")
    is_low = mp < float(cfg.borderline_max_prob)

    if cfg.borderline_mode == "low_conf_only":
        return is_low

    return is_other | is_low


def compute_core_centroids_and_radii(
    m_pc: pd.DataFrame,
    Xz: np.ndarray,
    cfg: Config
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, np.ndarray], Dict[str, float]]:
    sub = m_pc["o2_subcompartment"].astype(str)
    not_other = ~sub.str.endswith(f"{cfg.sub_label_sep}other")

    core_mask = not_other.copy()
    if cfg.core_min_prob is not None:
        if cfg.max_prob_col not in m_pc.columns:
            raise ValueError("--core-min-prob set but max_prob column not present.")
        mp = pd.to_numeric(m_pc[cfg.max_prob_col], errors="coerce")
        core_mask = core_mask & (mp >= float(cfg.core_min_prob))

    centroid_map: Dict[str, np.ndarray] = {}
    radius_map: Dict[str, float] = {}

    rows_centroids = []
    rows_radii = []

    for lab in sorted(sub[not_other].unique(), key=str):
        idx_all = np.where((sub == lab).to_numpy())[0]
        idx_core = np.where((sub == lab).to_numpy() & core_mask.to_numpy())[0]

        n_all = int(len(idx_all))
        n_core = int(len(idx_core))

        if n_core < cfg.reassign_min_core_n:
            continue

        c = np.nanmean(Xz[idx_core, :], axis=0)
        d = np.sqrt(np.sum((Xz[idx_core, :] - c) ** 2, axis=1))
        r = float(np.nanquantile(d, cfg.reassign_radius_quantile))

        centroid_map[str(lab)] = c
        radius_map[str(lab)] = r

        o2 = str(lab).split(cfg.sub_label_sep)[0]
        rows_radii.append(
            {
                "o2_compartment": o2,
                "o2_subcompartment": str(lab),
                "n_all": n_all,
                "n_core": n_core,
                "radius_quantile": cfg.reassign_radius_quantile,
                "radius": r,
            }
        )
        for j, pc in enumerate(cfg.pc_cols):
            rows_centroids.append(
                {
                    "o2_compartment": o2,
                    "o2_subcompartment": str(lab),
                    "pc": pc,
                    "centroid_z": float(c[j]),
                    "n_core": n_core,
                }
            )

    return pd.DataFrame(rows_centroids), pd.DataFrame(rows_radii), centroid_map, radius_map


def reassign_borderline(
    m_pc: pd.DataFrame,
    Xz: np.ndarray,
    cfg: Config,
    centroid_map: Dict[str, np.ndarray],
    radius_map: Dict[str, float],
    borderline_mask_pc: pd.Series
) -> pd.DataFrame:
    out = m_pc.copy()
    out["o2_subcompartment_before_reassign"] = out["o2_subcompartment"].astype(str)
    out["o2_subcompartment_after_reassign"] = out["o2_subcompartment"].astype(str)
    out["o2_subcompartment_final"] = out["o2_subcompartment"].astype(str)

    out["reassigned"] = False
    out["reassign_target"] = ""
    out["reassign_dist"] = np.nan
    out["reassign_radius"] = np.nan
    out["reassign_accept"] = False

    by_o2: Dict[str, List[str]] = {}
    for lab in centroid_map.keys():
        o2 = str(lab).split(cfg.sub_label_sep)[0]
        by_o2.setdefault(o2, []).append(lab)

    b_idx = np.where(borderline_mask_pc.to_numpy())[0]
    if len(b_idx) == 0:
        return out

    o2_series = out["o2_compartment"].astype(str).to_numpy()

    for i in b_idx:
        o2 = str(o2_series[i])
        candidates = by_o2.get(o2, [])
        if len(candidates) == 0:
            continue

        x = Xz[i, :]
        best_lab = None
        best_dist = None

        for lab in candidates:
            c = centroid_map[lab]
            d = float(np.sqrt(np.sum((x - c) ** 2)))
            if best_dist is None or d < best_dist:
                best_dist = d
                best_lab = lab

        if best_lab is None or best_dist is None:
            continue

        r = float(radius_map.get(best_lab, np.nan))
        accept = bool(np.isfinite(r) and (best_dist <= r))

        out.iloc[i, out.columns.get_loc("reassign_target")] = str(best_lab)
        out.iloc[i, out.columns.get_loc("reassign_dist")] = float(best_dist)
        out.iloc[i, out.columns.get_loc("reassign_radius")] = float(r)
        out.iloc[i, out.columns.get_loc("reassign_accept")] = bool(accept)

        if accept:
            out.iloc[i, out.columns.get_loc("o2_subcompartment_after_reassign")] = str(best_lab)
            out.iloc[i, out.columns.get_loc("o2_subcompartment_final")] = str(best_lab)
            out.iloc[i, out.columns.get_loc("reassigned")] = True

    return out


# ----------------------------
# Silhouette utilities
# ----------------------------

def weighted_silhouette_precomputed(X: np.ndarray, labels: np.ndarray, weights: np.ndarray) -> float:
    """
    Robust weighted silhouette:
      - compute Euclidean distance matrix
      - silhouette_score(metric="precomputed", sample_weight=weights)

    Returns np.nan if not computable.
    """
    lab = pd.Series(labels).astype(str).fillna("NA").to_numpy()
    if len(set(lab)) < 2:
        return np.nan

    w = pd.to_numeric(pd.Series(weights), errors="coerce").to_numpy(dtype=float)
    if not np.isfinite(w).any():
        return np.nan
    w = np.where(np.isfinite(w), w, 0.0)
    if np.all(w <= 0):
        return np.nan

    try:
        D = pairwise_distances(X, metric="euclidean")
        return float(silhouette_score(D, lab, metric="precomputed", sample_weight=w))
    except Exception:
        return np.nan


def unweighted_silhouette(X: np.ndarray, labels: np.ndarray) -> float:
    lab = pd.Series(labels).astype(str).fillna("NA").to_numpy()
    if len(set(lab)) < 2:
        return np.nan
    try:
        return float(silhouette_score(X, lab, metric="euclidean"))
    except Exception:
        return np.nan


def compute_silhouette_bundle_pcspace(
    m_pc: pd.DataFrame,
    Xz: np.ndarray,
    label_col: str,
    weights: Optional[np.ndarray] = None,
) -> Dict[str, object]:
    """
    Computes unweighted + (optional) weighted silhouette for a given label column in PC space.
    Returns a dict with metrics + a short note when not computable.
    """
    if label_col not in m_pc.columns:
        return {
            "label_col": label_col,
            "n_rows_used": int(len(m_pc)),
            "n_labels": np.nan,
            "silhouette_unweighted": np.nan,
            "silhouette_weighted_max_prob": np.nan,
            "note": "label_col_missing",
        }

    lab = m_pc[label_col].astype("object").fillna("NA").astype(str).to_numpy()
    n_labels = len(set(lab))
    if len(m_pc) < 5:
        return {
            "label_col": label_col,
            "n_rows_used": int(len(m_pc)),
            "n_labels": int(n_labels),
            "silhouette_unweighted": np.nan,
            "silhouette_weighted_max_prob": np.nan,
            "note": "too_few_rows",
        }
    if n_labels < 2:
        return {
            "label_col": label_col,
            "n_rows_used": int(len(m_pc)),
            "n_labels": int(n_labels),
            "silhouette_unweighted": np.nan,
            "silhouette_weighted_max_prob": np.nan,
            "note": "only_one_label",
        }

    sil_u = unweighted_silhouette(Xz, lab)

    sil_w = np.nan
    note = ""
    if weights is not None:
        sil_w = weighted_silhouette_precomputed(Xz, lab, weights)
    else:
        note = "weights_missing"

    return {
        "label_col": label_col,
        "n_rows_used": int(len(m_pc)),
        "n_labels": int(n_labels),
        "silhouette_unweighted": float(sil_u) if np.isfinite(sil_u) else np.nan,
        "silhouette_weighted_max_prob": float(sil_w) if np.isfinite(sil_w) else np.nan,
        "note": note,
    }


# ----------------------------
# Plotting utilities
# ----------------------------

def _savefig_all(fig: plt.Figure, outbase: str, cfg: Config) -> None:
    for fmt in cfg.plot_formats:
        path = f"{outbase}.{fmt}"
        if fmt == "png":
            fig.savefig(path, dpi=cfg.png_dpi, bbox_inches="tight")
        else:
            fig.savefig(path, bbox_inches="tight")


def _hybrid_display_label(label: object, cfg: Config) -> str:
    raw = str(label)
    if cfg.sub_label_sep not in raw:
        return raw
    oxygen, suffix = raw.split(cfg.sub_label_sep, 1)
    if suffix.startswith(cfg.prefix_gmm):
        return f"{oxygen}-GMM{suffix[len(cfg.prefix_gmm):]}"
    return f"{oxygen}-{suffix}"


def _hybrid_label_sort_key(label: object, cfg: Config) -> Tuple[int, int, str]:
    raw = str(label)
    oxygen = raw.split(cfg.sub_label_sep, 1)[0]
    oxygen_order = {"oxic": 0, "dysoxic": 1, "suboxic": 2, "anoxic": 3}
    gmm_index = _parse_gmm_index_from_label(raw, cfg)
    return (oxygen_order.get(oxygen, 99), 999 if gmm_index is None else gmm_index, raw)


def _draw_categorical_contours(
    ax: plt.Axes,
    x: np.ndarray,
    depth: np.ndarray,
    codes: np.ndarray,
    colors: List[str],
) -> None:
    """Render each categorical region separately without ordinal interpolation."""
    for code, color in enumerate(colors):
        mask = (codes == code).astype(float).T
        if not mask.any():
            continue
        ax.contourf(
            x, depth, mask,
            levels=[0.5, 1.5], colors=[color],
            antialiased=True, zorder=1,
        )
        if mask.min() < 0.5 < mask.max():
            ax.contour(
                x, depth, mask,
                levels=[0.5], colors=["#4D4D4D"],
                linewidths=0.25, alpha=0.55, zorder=2,
            )


def build_smoothed_categorical_curtain(
    profile: pd.DataFrame,
    cfg: Config,
    maximum_depth: float,
) -> Tuple[pd.DataFrame, List[str], np.ndarray, np.ndarray, np.ndarray, np.ndarray, pd.Timestamp, pd.Timestamp]:
    """Build a display-only categorical surface from observed hybrid labels.

    Each cruise is first represented by midpoint-bounded nearest-depth states.
    State-specific one-hot support is interpolated between actual cruise dates
    on a grid containing every calendar month in every represented year, then
    Gaussian-smoothed along time only. No smoothing is applied across depth.
    The displayed state is the state with maximum longitudinally smoothed
    support. Labels ending in ``__other`` can be excluded from support without
    altering the underlying sample assignments.
    """
    source = profile.copy()
    source["_label"] = source["o2_subcompartment_final"].astype(str)
    source["_is_other"] = source["_label"].str.endswith(
        f"{cfg.sub_label_sep}other"
    ) | source["_label"].str.lower().eq("outlier")
    eligible = source.loc[~source["_is_other"]].copy() if cfg.curtain_hide_other else source
    if eligible.empty:
        raise ValueError("No non-'other' hybrid states were available for curtain smoothing")

    labels = sorted(
        eligible["_label"].unique(),
        key=lambda label: _hybrid_label_sort_key(label, cfg),
    )
    label_to_code = {label: index for index, label in enumerate(labels)}
    cruise_dates = source[["cruise_index", "_date"]].drop_duplicates("cruise_index")
    cruise_dates = cruise_dates.sort_values("cruise_index").set_index("cruise_index")["_date"]
    cruise_dates = pd.to_datetime(cruise_dates, errors="coerce")
    if cruise_dates.isna().any():
        raise ValueError("All cruises require valid dates for calendar-time curtain smoothing")
    n_cruises = len(cruise_dates)
    calendar_start = pd.Timestamp(year=int(cruise_dates.dt.year.min()), month=1, day=1)
    calendar_end = pd.Timestamp(year=int(cruise_dates.dt.year.max()) + 1, month=1, day=1)
    month_starts = pd.date_range(calendar_start, calendar_end, freq="MS")
    subdivisions = cfg.curtain_time_subdivisions_per_month
    fine_dates_list: List[pd.Timestamp] = []
    for left, right in zip(month_starts[:-1], month_starts[1:]):
        width = right - left
        fine_dates_list.extend(
            left + width * ((index + 0.5) / subdivisions)
            for index in range(subdivisions)
        )
    fine_dates = pd.DatetimeIndex(fine_dates_list)
    fine_x = (fine_dates - calendar_start).total_seconds().to_numpy() / 86400.0
    coarse_x = (cruise_dates - calendar_start).dt.total_seconds().to_numpy() / 86400.0
    depth_centers = np.arange(
        0.0, maximum_depth + cfg.curtain_depth_step_m * 0.5,
        cfg.curtain_depth_step_m,
    )
    depth_centers[-1] = min(depth_centers[-1], maximum_depth)
    coarse = np.full((n_cruises, len(depth_centers)), np.nan, dtype=float)

    for cruise_index, frame in eligible.groupby("cruise_index", sort=True):
        frame = frame.sort_values("_depth").drop_duplicates("_depth", keep="first")
        depths = frame["_depth"].to_numpy(float)
        codes = frame["_label"].map(label_to_code).to_numpy(int)
        if len(depths) == 1:
            selected = np.zeros(len(depth_centers), dtype=int)
        else:
            midpoints = (depths[:-1] + depths[1:]) / 2.0
            selected = np.searchsorted(midpoints, depth_centers, side="right")
        coarse[int(cruise_index), :] = codes[selected]

    support = np.zeros((len(labels), len(fine_x), len(depth_centers)), dtype=float)
    for label_code in range(len(labels)):
        for depth_index in range(len(depth_centers)):
            available = np.isfinite(coarse[:, depth_index])
            if not available.any():
                continue
            values = (coarse[available, depth_index] == label_code).astype(float)
            dated = pd.DataFrame({"x": coarse_x[available], "value": values})
            dated = dated.groupby("x", as_index=False)["value"].mean().sort_values("x")
            support[label_code, :, depth_index] = np.interp(
                fine_x,
                dated["x"],
                dated["value"],
            )
        support[label_code] = gaussian_filter(
            support[label_code],
            sigma=(
                cfg.curtain_time_sigma_months * subdivisions,
                0.0,
            ),
            mode="nearest",
        )

    total_support = support.sum(axis=0)
    display_codes = support.argmax(axis=0)
    maximum_support = support.max(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        support_fraction = np.divide(
            maximum_support,
            total_support,
            out=np.zeros_like(maximum_support),
            where=total_support > 0,
        )
    grid = pd.DataFrame({
        "calendar_date": np.repeat(fine_dates.to_numpy(), len(depth_centers)),
        "calendar_month": np.repeat(fine_dates.to_period("M").astype(str), len(depth_centers)),
        "display_time_coordinate_days": np.repeat(fine_x, len(depth_centers)),
        "depth_m": np.tile(depth_centers, len(fine_x)),
        "display_state_code": display_codes.ravel(),
        "o2_subcompartment_display": np.asarray(labels, dtype=object)[display_codes.ravel()],
        "maximum_smoothed_support_fraction": support_fraction.ravel(),
    })
    return grid, labels, fine_x, depth_centers, display_codes, support, calendar_start, calendar_end


def enforce_observed_curtain_anchors(
    render_codes: np.ndarray,
    fine_x: np.ndarray,
    depth_centers: np.ndarray,
    profile: pd.DataFrame,
    labels: List[str],
    cfg: Config,
    calendar_start: pd.Timestamp,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Force the rendered surface to agree with observed cruise-depth states.

    Each observed cruise profile owns the fine-grid columns immediately
    bracketing its sampling date (and an adjacent column when the date lies on
    a grid node). Competing anchors are resolved by distance to the sampling
    date. The entire midpoint-bounded observed depth profile is imposed on
    those columns after all visual smoothing.
    """
    anchored = np.asarray(render_codes, dtype=int).copy()
    anchor_mask = np.zeros_like(anchored, dtype=bool)
    anchor_owner = np.full_like(anchored, -1, dtype=int)
    owner_distance = np.full(len(fine_x), np.inf, dtype=float)
    label_to_code = {label: index for index, label in enumerate(labels)}

    source = profile.copy()
    source["_label"] = source["o2_subcompartment_final"].astype(str)
    source["_is_other"] = source["_label"].str.endswith(
        f"{cfg.sub_label_sep}other"
    ) | source["_label"].str.lower().eq("outlier")
    if cfg.curtain_hide_other:
        source = source.loc[~source["_is_other"]].copy()

    for cruise_index, frame in source.groupby("cruise_index", sort=True):
        frame = frame.sort_values("_depth").drop_duplicates("_depth", keep="first")
        if frame.empty:
            continue
        depths = frame["_depth"].to_numpy(float)
        observed_codes = frame["_label"].map(label_to_code).to_numpy(int)
        if len(depths) == 1:
            selected = np.zeros(len(depth_centers), dtype=int)
        else:
            selected = np.searchsorted(
                (depths[:-1] + depths[1:]) / 2.0,
                depth_centers,
                side="right",
            )
        depth_profile_codes = observed_codes[selected]
        sample_date = pd.to_datetime(frame["_date"].iloc[0], errors="coerce")
        if pd.isna(sample_date):
            continue
        sample_x = float((sample_date - calendar_start).total_seconds() / 86400.0)
        insertion = int(np.searchsorted(fine_x, sample_x, side="left"))
        candidate_columns = {
            max(0, min(len(fine_x) - 1, insertion - 1)),
            max(0, min(len(fine_x) - 1, insertion)),
        }
        if insertion < len(fine_x) and np.isclose(fine_x[insertion], sample_x):
            candidate_columns.add(min(len(fine_x) - 1, insertion + 1))
        for column in sorted(candidate_columns):
            distance = abs(float(fine_x[column]) - sample_x)
            if distance > owner_distance[column]:
                continue
            anchored[column, :] = depth_profile_codes
            anchor_mask[column, :] = True
            anchor_owner[column, :] = int(cruise_index)
            owner_distance[column] = distance
    return anchored, anchor_mask, anchor_owner


def plot_hybrid_time_depth_curtain(
    m: pd.DataFrame,
    depth_col: str,
    sub_palette: Dict[str, str],
    plots_dir: str,
    tables_dir: str,
    cfg: Config,
    maximum_depth: float = 210.0,
) -> pd.DataFrame:
    """Draw a smoothed categorical time-depth surface with sampled-depth audit points."""
    required = {cfg.cruise_col, depth_col, "o2_subcompartment_final"}
    missing = sorted(required.difference(m.columns))
    if missing:
        raise ValueError(
            "Hybrid time-depth curtain is missing required columns: "
            + ", ".join(missing)
        )
    data = m.copy()
    data["_depth"] = pd.to_numeric(data[depth_col], errors="coerce")
    if cfg.date_col in data:
        data["_date"] = pd.to_datetime(data[cfg.date_col], errors="coerce")
    elif {"Year", "Month", "Day"}.issubset(data.columns):
        data["_date"] = pd.to_datetime(
            data[["Year", "Month", "Day"]].rename(
                columns={"Year": "year", "Month": "month", "Day": "day"}
            ),
            errors="coerce",
        )
    else:
        data["_date"] = pd.NaT
    data = data.loc[
        data["_depth"].between(0.0, maximum_depth, inclusive="both")
        & data["o2_subcompartment_final"].notna()
        & data["_date"].notna()
    ].copy()
    if data.empty:
        raise ValueError("No classified observations were available for the hybrid time-depth curtain")

    cruise_table = data[[cfg.cruise_col, "_date"]].drop_duplicates().copy()
    cruise_table["_cruise_numeric"] = pd.to_numeric(
        cruise_table[cfg.cruise_col], errors="coerce"
    )
    cruise_table = cruise_table.sort_values(
        ["_date", "_cruise_numeric", cfg.cruise_col], na_position="last"
    ).reset_index(drop=True)
    cruise_table["cruise_index"] = np.arange(len(cruise_table), dtype=int)
    calendar_start = pd.Timestamp(
        year=int(cruise_table["_date"].dt.year.min()), month=1, day=1
    )
    cruise_table["time_coordinate_days"] = (
        cruise_table["_date"] - calendar_start
    ).dt.total_seconds() / 86400.0
    data = data.merge(
        cruise_table[[cfg.cruise_col, "_date", "cruise_index", "time_coordinate_days"]],
        on=[cfg.cruise_col, "_date"], how="inner", validate="many_to_one",
    )

    confidence = pd.to_numeric(data.get(cfg.max_prob_col), errors="coerce")
    data["_confidence"] = confidence if confidence is not None else np.nan
    data = data.sort_values(
        ["cruise_index", "_depth", "_confidence", "o2_subcompartment_final"],
        ascending=[True, True, False, True], na_position="last",
    )
    # Multiple bottles at the same cruise and anchored depth represent one
    # curtain cell; the most confident final assignment is used for its fill.
    profile = data.drop_duplicates(["cruise_index", "_depth"], keep="first")

    cells = []
    for cruise_index, frame in profile.groupby("cruise_index", sort=True):
        frame = frame.sort_values("_depth").reset_index(drop=True)
        depths = frame["_depth"].to_numpy(float)
        if len(depths) == 1:
            boundaries = np.array([0.0, maximum_depth], dtype=float)
        else:
            boundaries = np.concatenate((
                [0.0],
                (depths[:-1] + depths[1:]) / 2.0,
                [maximum_depth],
            ))
        for row_index, row in frame.iterrows():
            label = str(row["o2_subcompartment_final"])
            cells.append({
                "cruise_index": int(cruise_index),
                cfg.cruise_col: row[cfg.cruise_col],
                "date": row["_date"],
                "time_coordinate_days": float(row["time_coordinate_days"]),
                "measurement_depth_m": float(row["_depth"]),
                "cell_top_depth_m": float(boundaries[row_index]),
                "cell_bottom_depth_m": float(boundaries[row_index + 1]),
                "o2_subcompartment_final": label,
                "display_label": _hybrid_display_label(label, cfg),
                "color_hex": sub_palette.get(label, "#BDBDBD"),
            })
    cell_table = pd.DataFrame(cells)
    cell_table["is_other_assignment"] = cell_table[
        "o2_subcompartment_final"
    ].astype(str).str.endswith(f"{cfg.sub_label_sep}other")
    cell_table["included_in_smoothed_display_support"] = (
        ~cell_table["is_other_assignment"]
        if cfg.curtain_hide_other
        else True
    )
    grid, observed_labels, fine_x, depth_centers, display_codes, support, calendar_start, calendar_end = (
        build_smoothed_categorical_curtain(profile, cfg, maximum_depth)
    )
    contour_support = np.stack([
        gaussian_filter(
            state_support,
            sigma=(
                0.0,
                cfg.curtain_contour_visual_depth_sigma_m / cfg.curtain_depth_step_m,
            ),
            mode="nearest",
        )
        for state_support in support
    ])
    contour_codes = contour_support.argmax(axis=0)
    contour_codes, anchor_mask, anchor_owner = enforce_observed_curtain_anchors(
        contour_codes,
        fine_x,
        depth_centers,
        profile,
        observed_labels,
        cfg,
        calendar_start,
    )
    grid["render_state_code"] = contour_codes.ravel()
    grid["o2_subcompartment_rendered"] = np.asarray(
        observed_labels, dtype=object
    )[contour_codes.ravel()]
    grid["observed_anchor_enforced"] = anchor_mask.ravel()
    grid["anchor_cruise_index"] = np.where(
        anchor_mask, anchor_owner, np.nan
    ).ravel()

    rendered_at_measurement = []
    for row in cell_table.itertuples(index=False):
        time_index = int(np.argmin(np.abs(fine_x - float(row.time_coordinate_days))))
        depth_index = int(np.argmin(np.abs(depth_centers - float(row.measurement_depth_m))))
        rendered_at_measurement.append(observed_labels[int(contour_codes[time_index, depth_index])])
    cell_table["rendered_state_at_measurement"] = rendered_at_measurement
    cell_table["rendered_state_matches_observation"] = (
        cell_table["rendered_state_at_measurement"].astype(str)
        == cell_table["o2_subcompartment_final"].astype(str)
    )
    included_mismatch = (
        cell_table["included_in_smoothed_display_support"].astype(bool)
        & ~cell_table["rendered_state_matches_observation"]
    )
    if included_mismatch.any():
        raise RuntimeError(
            "Observed-anchor enforcement failed for "
            f"{int(included_mismatch.sum())} curtain cells"
        )
    cell_table.to_csv(
        os.path.join(tables_dir, "hybrid_compartment_time_depth_curtain_cells.csv"),
        index=False,
    )
    grid.to_csv(
        os.path.join(tables_dir, "hybrid_compartment_time_depth_curtain_grid.csv"),
        index=False,
    )

    renewal_rows = pd.DataFrame(columns=["renewal_date", "time_coordinate_days"])
    renewal_path = getattr(cfg, "curtain_renewal_events", None)
    renewal_date_col = getattr(cfg, "curtain_renewal_date_col", "start_date")
    if renewal_path:
        renewal_sep = "\t" if str(renewal_path).lower().endswith((".tsv", ".txt")) else ","
        renewal_source = pd.read_csv(renewal_path, sep=renewal_sep)
        if renewal_date_col not in renewal_source.columns:
            raise ValueError(
                f"Renewal-event table lacks configured onset column {renewal_date_col!r}"
            )
        renewal_dates = pd.to_datetime(
            renewal_source[renewal_date_col], errors="coerce"
        ).dropna().drop_duplicates().sort_values()
        renewal_rows = pd.DataFrame({"renewal_date": renewal_dates})
        renewal_rows["time_coordinate_days"] = (
            renewal_rows["renewal_date"] - calendar_start
        ).dt.total_seconds() / 86400.0
        renewal_rows = renewal_rows.loc[
            renewal_rows["time_coordinate_days"].between(
                0.0, float((calendar_end - calendar_start).days), inclusive="both"
            )
        ].reset_index(drop=True)
    renewal_rows.to_csv(
        os.path.join(tables_dir, "hybrid_compartment_time_depth_curtain_renewals.csv"),
        index=False,
    )
    pd.DataFrame([{
        "hide_other": cfg.curtain_hide_other,
        "calendar_start": calendar_start,
        "calendar_end_exclusive": calendar_end,
        "calendar_months_n": (calendar_end.year - calendar_start.year) * 12,
        "time_subdivisions_per_calendar_month": cfg.curtain_time_subdivisions_per_month,
        "time_gaussian_sigma_months": cfg.curtain_time_sigma_months,
        "depth_grid_step_m": cfg.curtain_depth_step_m,
        "maximum_display_depth_m": maximum_depth,
        "depth_gaussian_sigma_m": 0.0,
        "smoothing_axes": "calendar_time_only",
        "rendering": "categorical_filled_contours",
        "contour_visual_depth_sigma_m": cfg.curtain_contour_visual_depth_sigma_m,
        "contour_visual_smoothing_changes_assignments": False,
        "observed_anchor_constraint": True,
        "anchor_grid_cells_n": int(anchor_mask.sum()),
        "included_observation_mismatches_n": int(included_mismatch.sum()),
        "source_samples_n": len(profile),
        "other_source_samples_n": int(
            profile["o2_subcompartment_final"].astype(str).str.endswith(
                f"{cfg.sub_label_sep}other"
            ).sum()
        ),
        "display_states_n": len(observed_labels),
        "renewal_onsets_plotted_n": len(renewal_rows),
        "renewal_onset_source": renewal_path or "none",
        "renewal_onset_date_column": renewal_date_col,
    }]).to_csv(
        os.path.join(tables_dir, "hybrid_compartment_time_depth_curtain_smoothing.csv"),
        index=False,
    )

    fig, ax = plt.subplots(figsize=(18.0, 7.5))
    plot_x = np.concatenate((
        [0.0], fine_x, [float((calendar_end - calendar_start).days)]
    ))
    plot_codes = np.concatenate(
        (contour_codes[:1], contour_codes, contour_codes[-1:]), axis=0
    )
    _draw_categorical_contours(
        ax,
        plot_x,
        depth_centers,
        plot_codes,
        [sub_palette.get(label, "#BDBDBD") for label in observed_labels],
    )
    ax.scatter(
        profile["time_coordinate_days"], profile["_depth"],
        s=7.5, facecolor="black", edgecolor="none", alpha=1.0, zorder=4,
    )
    for renewal_x in renewal_rows["time_coordinate_days"]:
        ax.axvline(
            renewal_x, color="black", linestyle="--", linewidth=0.9,
            alpha=0.9, zorder=3,
        )
    year_starts = pd.date_range(calendar_start, calendar_end, freq="YS", inclusive="left")
    year_x = (year_starts - calendar_start).total_seconds().to_numpy() / 86400.0
    month_starts = pd.date_range(calendar_start, calendar_end, freq="MS", inclusive="left")
    month_x = (month_starts - calendar_start).total_seconds().to_numpy() / 86400.0
    ax.set_xticks(year_x)
    ax.set_xticklabels([str(date.year) for date in year_starts], ha="left", fontsize=8)
    ax.set_xticks(month_x, minor=True)
    ax.tick_params(axis="x", which="minor", length=2.5)
    ax.set_xlim(0.0, float((calendar_end - calendar_start).days))
    ax.set_ylim(maximum_depth, 0.0)
    ax.set_yticks(np.arange(0, maximum_depth + 1, 25))
    ax.set_ylabel("Depth (m)")
    ax.set_xlabel("Calendar time (monthly intervals)")
    ax.set_title(r"Hybrid O$_2$-GMM compartments through time and depth")

    handles = [
        Patch(
            facecolor=sub_palette.get(label, "#BDBDBD"), edgecolor="none",
            label=_hybrid_display_label(label, cfg),
        )
        for label in observed_labels
    ]
    handles.append(Line2D(
        [], [], marker="o", linestyle="", markersize=3.5,
        markerfacecolor="black", markeredgecolor="none", label="Sampled depth",
    ))
    if not renewal_rows.empty:
        handles.append(Line2D(
            [], [], color="black", linestyle="--", linewidth=0.9,
            label="Predicted renewal onset",
        ))
    ax.legend(
        handles=handles, title=r"Hybrid O$_2$-GMM compartment",
        loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False,
        fontsize=8, title_fontsize=9,
    )
    fig.subplots_adjust(left=0.07, right=0.82, bottom=0.20, top=0.91)
    _savefig_all(
        fig,
        os.path.join(plots_dir, "hybrid_compartment_time_depth_curtain"),
        cfg,
    )
    plt.close(fig)
    return cell_table


def _smooth_profile_by_depth_bins(
    depth: np.ndarray,
    x: np.ndarray,
    n_bins: int = 40,
    smooth_window: int = 7,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Smooth mean profile x(depth), per subcompartment:
      - quantile-bin depths
      - compute mean depth + mean x per bin
      - rolling mean smooth across bins

    Returns (x_smooth, depth_bin_centers)
    """
    d = pd.to_numeric(pd.Series(depth), errors="coerce").to_numpy(dtype=float)
    v = pd.to_numeric(pd.Series(x), errors="coerce").to_numpy(dtype=float)
    ok = np.isfinite(d) & np.isfinite(v)
    d = d[ok]
    v = v[ok]
    if len(d) < 10:
        return np.array([]), np.array([])

    df = pd.DataFrame({"depth": d, "x": v}).sort_values("depth")

    # robust bins even with uneven depth sampling
    try:
        q = min(n_bins, max(5, len(df) // 5))
        df["bin"] = pd.qcut(df["depth"], q=q, duplicates="drop")
    except Exception:
        return np.array([]), np.array([])

    g = df.groupby("bin", observed=True).agg(
        depth_mean=("depth", "mean"),
        x_mean=("x", "mean"),
        n=("x", "size"),
    ).reset_index(drop=True)

    if len(g) < 5:
        return np.array([]), np.array([])

    w = int(smooth_window)
    if w < 3:
        w = 3
    if w % 2 == 0:
        w += 1
    if w > len(g):
        w = max(3, (len(g) // 2) * 2 + 1)

    x_smooth = g["x_mean"].rolling(window=w, center=True, min_periods=max(3, w // 3)).mean().to_numpy()
    d_smooth = g["depth_mean"].to_numpy()
    return x_smooth, d_smooth


def plot_pc_scatter_subcompartments(
    df: pd.DataFrame,
    o2: str,
    pc_x: str,
    pc_y: str,
    sub_palette: Dict[str, str],
    plots_dir: str,
    cfg: Config,
    xlim: Tuple[float, float],
    ylim: Tuple[float, float],
) -> None:
    d = df.loc[df["o2_compartment"].astype(str) == o2].copy()
    if len(d) < 3:
        return
    if pc_x not in d.columns or pc_y not in d.columns:
        return

    x = pd.to_numeric(d[pc_x], errors="coerce")
    y = pd.to_numeric(d[pc_y], errors="coerce")
    ok = x.notna() & y.notna()
    d = d.loc[ok].copy()
    if len(d) < 3:
        return

    labs = d["o2_subcompartment_final"].astype(str).to_numpy()

    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111)

    # draw each subcompartment separately for a clean legend
    unique_labs = sorted(pd.unique(labs), key=str)
    for lab in unique_labs:
        msk = (labs == lab)
        col = sub_palette.get(str(lab), "gray")
        ax.scatter(
            d.loc[msk, pc_x],
            d.loc[msk, pc_y],
            s=cfg.point_size,
            alpha=cfg.alpha,
            label=str(lab),
            color=col,            # <-- REQUIRED: force palette color for points
            edgecolors="gray",
            linewidths=0.5,
        )

    # overlay reassigned points, colored by FINAL subcompartment
    if "reassigned" in d.columns:
        rr = d["reassigned"].astype(bool)
        if rr.any():
            d_rr = d.loc[rr].copy()
            rr_labs = d_rr["o2_subcompartment_final"].astype(str).to_numpy()
            for lab in sorted(pd.unique(rr_labs), key=str):
                msk = (rr_labs == lab)
                col_rr = sub_palette.get(str(lab), "gray")
                ax.scatter(
                    d_rr.loc[msk, pc_x],
                    d_rr.loc[msk, pc_y],
                    s=cfg.point_size * 1.2,
                    alpha=1.0,
                    marker="D",
                    linewidths=0.5,
                    color=col_rr,
                    edgecolors="gray",
                    label=f"{lab} (reassigned)",
                )

    ax.set_autoscale_on(False)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel(pc_x)
    ax.set_ylabel(pc_y)
    ax.set_title(f"{o2}: PC scatter by subcompartment")
    present = set(d["o2_subcompartment_final"].astype(str).unique())

    handles = [
        Line2D(
            [],
            [],
            marker="o",
            linestyle="None",
            markersize=6,
            markerfacecolor=color,
            markeredgecolor="gray",
            markeredgewidth=0.5,
            label=str(key),
        )
        for key, color in sub_palette.items()
        if str(key) in present
    ]

    ax.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        fontsize=7,
        frameon=False,
    )
    fig.subplots_adjust(right=0.75)

    outbase = os.path.join(plots_dir, f"o2_{o2}_pc_scatter_subcompartments_{pc_x}_vs_{pc_y}")
    _savefig_all(fig, outbase, cfg)
    plt.close(fig)


def plot_pc_scatter_all_subcompartments(
    df: pd.DataFrame,
    pc_x: str,
    pc_y: str,
    sub_palette: Dict[str, str],
    plots_dir: str,
    cfg: Config,
) -> None:
    if pc_x not in df.columns or pc_y not in df.columns:
        return

    d = df.copy()
    x = pd.to_numeric(d[pc_x], errors="coerce")
    y = pd.to_numeric(d[pc_y], errors="coerce")
    ok = x.notna() & y.notna()
    d = d.loc[ok].copy()
    if len(d) < 5:
        return

    labs = d["o2_subcompartment_final"].astype(str).to_numpy()
    unique_labs = sorted(pd.unique(labs), key=str)

    fig = plt.figure(figsize=(8.5, 6))
    ax = fig.add_subplot(111)

    for lab in unique_labs:
        msk = (labs == lab)
        col = sub_palette.get(str(lab), "gray")
        ax.scatter(
            d.loc[msk, pc_x],
            d.loc[msk, pc_y],
            s=cfg.point_size,
            alpha=cfg.alpha,
            color=col,
            label=str(lab),
            edgecolors="gray",
            linewidths=0.5,
        )

    ax.set_xlabel(pc_x)
    ax.set_ylabel(pc_y)
    ax.set_title("All subcompartments: PC scatter (EDA)")
    present = set(d["o2_subcompartment_final"].astype(str).unique())

    handles = [
        Line2D(
            [],
            [],
            marker="o",
            linestyle="None",
            markersize=6,
            markerfacecolor=color,
            markeredgecolor="gray",
            markeredgewidth=0.5,
            label=str(key),
            alpha=cfg.alpha,
        )
        for key, color in sub_palette.items()
        if str(key) in present
    ]

    ax.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        fontsize=7,
        frameon=False,
    )
    fig.subplots_adjust(right=0.70)

    outbase = os.path.join(plots_dir, f"EDA_all_subcompartments_pc_scatter_{pc_x}_vs_{pc_y}")
    _savefig_all(fig, outbase, cfg)
    plt.close(fig)


def plot_depth_profile_all_subcompartments(
    df: pd.DataFrame,
    value_col: str,
    depth_col: str,
    sub_palette: Dict[str, str],
    plots_dir: str,
    cfg: Config,
) -> None:
    if value_col not in df.columns or depth_col not in df.columns:
        return

    d = df.copy()
    x = pd.to_numeric(d[value_col], errors="coerce")
    dep = pd.to_numeric(d[depth_col], errors="coerce")
    ok = x.notna() & dep.notna()
    d = d.loc[ok].copy()
    if len(d) < 5:
        return

    labs = d["o2_subcompartment_final"].astype(str).to_numpy()
    unique_labs = sorted(pd.unique(labs), key=str)

    fig = plt.figure(figsize=(8.5, 6))
    ax = fig.add_subplot(111)

    for lab in unique_labs:
        msk = (labs == lab)
        col = sub_palette.get(str(lab), "gray")
        ax.scatter(
            d.loc[msk, value_col],
            d.loc[msk, depth_col],
            s=cfg.point_size,
            alpha=cfg.alpha,
            color=col,
            label=str(lab),
            edgecolors="gray",
            linewidths=0.5,
        )

        # Smooth mean "fitted" profile per subcompartment (same color as points)
        xs, ds = _smooth_profile_by_depth_bins(
            d.loc[msk, depth_col].to_numpy(),
            d.loc[msk, value_col].to_numpy(),
            n_bins=40,
            smooth_window=7,
        )
        if len(xs) > 0:
            ax.plot(xs, ds, linewidth=3.0, alpha=0.95, color=col)

    ax.invert_yaxis()
    ax.set_xlabel(value_col)
    ax.set_ylabel(depth_col)
    ax.set_title(f"All subcompartments: depth profile of {value_col} (EDA)")
    present = set(d["o2_subcompartment_final"].astype(str).unique())

    handles = [
        Line2D(
            [],
            [],
            marker="o",
            linestyle="None",
            markersize=6,
            markerfacecolor=color,
            markeredgecolor="gray",
            markeredgewidth=0.5,
            label=str(key),
        )
        for key, color in sub_palette.items()
        if str(key) in present
    ]

    ax.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        fontsize=7,
        frameon=False,
    )
    fig.subplots_adjust(right=0.70)

    outbase = os.path.join(plots_dir, f"EDA_all_subcompartments_depth_profile_{value_col}")
    _savefig_all(fig, outbase, cfg)
    plt.close(fig)


def plot_pc_scatter_biochem_overlay(
    df: pd.DataFrame,
    o2: str,
    pc_x: str,
    pc_y: str,
    feature: str,
    plots_dir: str,
    cfg: Config,
) -> None:
    d = df.loc[df["o2_compartment"].astype(str) == o2].copy()
    if len(d) < 3:
        return
    if pc_x not in d.columns or pc_y not in d.columns or feature not in d.columns:
        return

    x = pd.to_numeric(d[pc_x], errors="coerce")
    y = pd.to_numeric(d[pc_y], errors="coerce")
    v = pd.to_numeric(d[feature], errors="coerce")
    ok = x.notna() & y.notna() & v.notna()
    d = d.loc[ok].copy()
    if len(d) < 3:
        return

    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111)

    sc = ax.scatter(
        d[pc_x],
        d[pc_y],
        c=d[feature],
        s=cfg.point_size,
        alpha=cfg.alpha,
    )

    ax.set_xlabel(pc_x)
    ax.set_ylabel(pc_y)
    ax.set_title(f"{o2}: {feature} over PC space")

    cb = fig.colorbar(sc, ax=ax)
    cb.set_label(feature)

    outbase = os.path.join(plots_dir, f"o2_{o2}_pc_scatter_{feature}_{pc_x}_vs_{pc_y}")
    _savefig_all(fig, outbase, cfg)
    plt.close(fig)


def plot_depth_profile_subcompartments(
    df: pd.DataFrame,
    o2: str,
    value_col: str,
    depth_col: str,
    sub_palette: Dict[str, str],
    plots_dir: str,
    cfg: Config,
) -> None:
    """
    Depth profile scatter: x=value_col, y=depth (inverted so shallow at top).
    Colored by o2_subcompartment_final.
    """
    d = df.loc[df["o2_compartment"].astype(str) == o2].copy()
    if len(d) < 3:
        return
    if value_col not in d.columns or depth_col not in d.columns:
        return

    x = pd.to_numeric(d[value_col], errors="coerce")
    dep = pd.to_numeric(d[depth_col], errors="coerce")
    ok = x.notna() & dep.notna()
    d = d.loc[ok].copy()
    if len(d) < 3:
        return

    labs = d["o2_subcompartment_final"].astype(str).to_numpy()
    unique_labs = sorted(pd.unique(labs), key=str)

    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111)

    for lab in unique_labs:
        msk = (labs == lab)
        col = sub_palette.get(str(lab), "gray")
        ax.scatter(
            d.loc[msk, value_col],
            d.loc[msk, depth_col],
            s=cfg.point_size,
            alpha=cfg.alpha,
            label=str(lab),
            color=col,           # <-- REQUIRED
            edgecolors="gray",
            linewidths=0.5,
        )

        # Smooth mean "fitted" profile per subcompartment (same color as points)
        xs, ds = _smooth_profile_by_depth_bins(
            d.loc[msk, depth_col].to_numpy(),
            d.loc[msk, value_col].to_numpy(),
            n_bins=40,
            smooth_window=7,
        )
        if len(xs) > 0:
            ax.plot(xs, ds, linewidth=3.0, alpha=0.95, color=col)   # <-- same palette color

    # shallow at top
    ax.invert_yaxis()

    ax.set_xlabel(value_col)
    ax.set_ylabel(depth_col)
    ax.set_title(f"{o2}: depth profile colored by subcompartment")
    present = set(d["o2_subcompartment_final"].astype(str).unique())

    handles = [
        Line2D(
            [],
            [],
            marker="o",
            linestyle="None",
            markersize=6,
            markerfacecolor=color,
            markeredgecolor="gray",
            markeredgewidth=0.5,
            label=str(key),
        )
        for key, color in sub_palette.items()
        if str(key) in present
    ]

    ax.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        fontsize=8,
        frameon=False,
    )
    fig.subplots_adjust(right=0.75)

    outbase = os.path.join(plots_dir, f"o2_{o2}_depth_profile_{value_col}_by_subcompartment")
    _savefig_all(fig, outbase, cfg)
    plt.close(fig)

def plot_umap_subcompartments(
    m: pd.DataFrame,
    umap_df: pd.DataFrame,
    sub_palette: Dict[str, str],
    plots_dir: str,
    cfg: Config,
) -> None:
    # merge UMAP with final labels
    d = umap_df.merge(
        m[["__merge_key__", "o2_subcompartment_final"]],
        on="__merge_key__",
        how="left",
    )

    d["UMAP1"] = pd.to_numeric(d["UMAP1"], errors="coerce")
    d["UMAP2"] = pd.to_numeric(d["UMAP2"], errors="coerce")
    d = d.dropna(subset=["UMAP1", "UMAP2", "o2_subcompartment_final"])

    if len(d) < 5:
        return

    labs = d["o2_subcompartment_final"].astype(str).to_numpy()
    unique_labs = sorted(pd.unique(labs), key=str)

    fig = plt.figure(figsize=(8.5, 6.5))
    ax = fig.add_subplot(111)

    for lab in unique_labs:
        msk = (labs == lab)
        col = sub_palette.get(str(lab), "gray")
        ax.scatter(
            d.loc[msk, "UMAP1"],
            d.loc[msk, "UMAP2"],
            s=cfg.point_size,
            alpha=cfg.alpha,
            color=col,
            label=str(lab),
            edgecolors="gray",
            linewidths=0.4,
        )

    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title("UMAP embedding colored by final O₂ × GMM subcompartments")
    present = set(d["o2_subcompartment_final"].astype(str).unique())

    handles = [
        Line2D(
            [],
            [],
            marker="o",
            linestyle="None",
            markersize=6,
            markerfacecolor=color,
            markeredgecolor="gray",
            markeredgewidth=0.4,
            label=str(key),
        )
        for key, color in sub_palette.items()
        if str(key) in present
    ]

    ax.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        fontsize=8,
        frameon=False,
    )

    fig.subplots_adjust(right=0.72)

    outbase = os.path.join(plots_dir, "umap_final_subcompartments")
    _savefig_all(fig, outbase, cfg)
    plt.close(fig)


# ----------------------------
# Main
# ----------------------------

def main() -> None:
    cfg = parse_args()
    tables_dir, plots_dir = ensure_dirs(cfg.outdir)

    df_matrix = read_table_dedup_cols(cfg.matrix_cleaned, cfg.sep_matrix)
    df_eig = read_table_dedup_cols(cfg.eigenvectors, cfg.sep_eig)
    df_assign = read_table_dedup_cols(cfg.assignments, cfg.sep_assign)
    df_o2_assign = (
        read_table_dedup_cols(cfg.o2_assignments, cfg.sep_o2_assign)
        if cfg.o2_assignments
        else None
    )

    df_matrix = build_merge_key(df_matrix, cfg)
    df_eig = build_merge_key(df_eig, cfg)
    df_assign = build_merge_key(df_assign, cfg)
    if df_o2_assign is not None:
        df_o2_assign = build_merge_key(df_o2_assign, cfg)

    if (df_o2_assign is None) and (cfg.oxygen_col not in df_matrix.columns):
        raise ValueError(f"matrix_cleaned missing oxygen col: {cfg.oxygen_col}")
    if cfg.gmm_component_col not in df_assign.columns:
        raise ValueError(f"assignments missing GMM component col: {cfg.gmm_component_col}")

    missing_pc = [c for c in cfg.pc_cols if c not in df_eig.columns]
    if missing_pc:
        raise ValueError(f"eigenvectors missing requested pc-cols: {missing_pc}")

    # merge: assignment rows authoritative
    keep_matrix_cols = [cfg.derived_key_col]
    for c in [cfg.oxygen_col, cfg.cruise_col, cfg.depth_col, cfg.depth_anchored_col, cfg.date_col]:
        if c in df_matrix.columns and c not in keep_matrix_cols:
            keep_matrix_cols.append(c)

    keep_eig_cols = [cfg.derived_key_col] + cfg.pc_cols

    m = df_assign.merge(df_matrix[keep_matrix_cols], on=cfg.derived_key_col, how="left")
    m = coalesce_merge_suffix_columns(m, prefer="x")  # keeps assignment-side cols if duplicates

    m = m.merge(df_eig[keep_eig_cols], on=cfg.derived_key_col, how="left")
    m = coalesce_merge_suffix_columns(m, prefer="x")

    # O2 compartments (prefer external soft/smoothed assignments when provided)
    if df_o2_assign is not None:
        o2_tbl = df_o2_assign[[cfg.derived_key_col]].copy()
        o2_tbl["o2_compartment"] = o2_labels_from_assignments(df_o2_assign, cfg)
        o2_tbl = o2_tbl.dropna(subset=[cfg.derived_key_col]).drop_duplicates(subset=[cfg.derived_key_col], keep="first")
        o2_lookup = o2_tbl.set_index(cfg.derived_key_col)["o2_compartment"]
        m["o2_compartment"] = m[cfg.derived_key_col].map(o2_lookup).fillna("NA").astype("object")
    else:
        m["o2_compartment"] = label_o2_compartment(m[cfg.oxygen_col], cfg)
    m["gmm_component"] = m[cfg.gmm_component_col].astype("object").fillna("NA").astype(str)

    # intersection label
    m["o2_subcompartment"] = [
        _make_sub_label(str(o2), str(g), cfg) for o2, g in zip(m["o2_compartment"].astype(str), m["gmm_component"].astype(str))
    ]

    # collapse tiny intersections
    if cfg.min_subcluster_size > 1:
        m = collapse_small_intersections(m, cfg)

    # baseline counts
    counts_before = (
        m.groupby(["o2_compartment", "o2_subcompartment"], dropna=False)
         .size()
         .reset_index(name="n")
         .sort_values(["o2_compartment", "n"], ascending=[True, False])
    )
    counts_before.to_csv(os.path.join(tables_dir, "o2_subcompartment_counts_before.csv"), index=False)

    # confusion O2 vs GMM
    raw, row_norm, col_norm = confusion_tables(m["o2_compartment"], m["gmm_component"])
    raw.to_csv(os.path.join(tables_dir, "o2_by_gmm_confusion_raw.csv"))
    row_norm.to_csv(os.path.join(tables_dir, "o2_by_gmm_confusion_row_norm.csv"))
    col_norm.to_csv(os.path.join(tables_dir, "o2_by_gmm_confusion_col_norm.csv"))

    # coerce PCs numeric
    for c in cfg.pc_cols:
        m[c] = pd.to_numeric(m[c], errors="coerce")

    pc_complete_mask = m[cfg.pc_cols].notna().all(axis=1)

    # default final labels (no reassignment)
    m["o2_subcompartment_before_reassign"] = m["o2_subcompartment"].astype(str)
    m["o2_subcompartment_after_reassign"] = m["o2_subcompartment"].astype(str)
    m["o2_subcompartment_final"] = m["o2_subcompartment"].astype(str)
    m["reassigned"] = False
    m["reassign_target"] = ""
    m["reassign_dist"] = np.nan
    m["reassign_radius"] = np.nan
    m["reassign_accept"] = False

    qc = {
        "pc_cols_used": cfg.pc_cols,
        "n_rows_total": int(len(m)),
        "n_rows_pc_complete": int(pc_complete_mask.sum()),
        "do_reassign": bool(cfg.do_reassign),
        "borderline_mode": cfg.borderline_mode,
        "borderline_max_prob": cfg.borderline_max_prob,
        "core_min_prob": cfg.core_min_prob,
        "reassign_radius_quantile": cfg.reassign_radius_quantile,
        "reassign_min_core_n": cfg.reassign_min_core_n,
        "n_borderline_candidates": 0,
        "n_reassignable_pc_complete": 0,
        "n_reassigned": 0,
    }

    centroids_df = pd.DataFrame([])
    radii_df = pd.DataFrame([])
    scaler_mu = None
    scaler_sd = None

    if cfg.do_reassign:
        borderline_mask = determine_borderline_mask(m, cfg)
        qc["n_borderline_candidates"] = int(borderline_mask.sum())
        qc["n_reassignable_pc_complete"] = int((borderline_mask & pc_complete_mask).sum())

        m_pc = m.loc[pc_complete_mask].copy()
        Xz, mu, sd = standardize_pc_space(m_pc, cfg.pc_cols)
        scaler_mu = mu
        scaler_sd = sd

        centroids_df, radii_df, centroid_map, radius_map = compute_core_centroids_and_radii(m_pc, Xz, cfg)

        borderline_pc = borderline_mask.loc[m_pc.index]
        m_pc_reassigned = reassign_borderline(m_pc, Xz, cfg, centroid_map, radius_map, borderline_pc)

        cols_pull = [
            "o2_subcompartment_before_reassign",
            "o2_subcompartment_after_reassign",
            "o2_subcompartment_final",
            "reassigned",
            "reassign_target",
            "reassign_dist",
            "reassign_radius",
            "reassign_accept",
        ]
        m.loc[m_pc_reassigned.index, cols_pull] = m_pc_reassigned[cols_pull]
        qc["n_reassigned"] = int(m["reassigned"].sum())

        qc["pc_scaler_mu"] = {pc: float(mu[i]) for i, pc in enumerate(cfg.pc_cols)}
        qc["pc_scaler_sd"] = {pc: float(sd[i]) for i, pc in enumerate(cfg.pc_cols)}

    # counts after reassignment
    counts_after = (
        m.groupby(["o2_compartment", "o2_subcompartment_final"], dropna=False)
         .size()
         .reset_index(name="n")
         .sort_values(["o2_compartment", "n"], ascending=[True, False])
    )
    counts_after.to_csv(os.path.join(tables_dir, "o2_subcompartment_counts_after.csv"), index=False)

    # within-O2 silhouettes in PC space (PC-complete only)
    rows_sil = []
    if pc_complete_mask.sum() >= 5:
        m_pc_all = m.loc[pc_complete_mask].copy()
        Xz_all, _, _ = standardize_pc_space(m_pc_all, cfg.pc_cols)

        weights_all = None
        if cfg.max_prob_col in m_pc_all.columns:
            weights_all = pd.to_numeric(m_pc_all[cfg.max_prob_col], errors="coerce").to_numpy(dtype=float)

        for o2 in sorted(m_pc_all["o2_compartment"].astype(str).unique(), key=str):
            idx = np.where(m_pc_all["o2_compartment"].astype(str).to_numpy() == o2)[0]
            if len(idx) < 5:
                rows_sil.append({
                    "o2_compartment": o2,
                    "n_rows_used": int(len(idx)),
                    "n_subcompartments": np.nan,
                    "silhouette_unweighted": np.nan,
                    "silhouette_weighted_max_prob": np.nan,
                    "note": "too_few_rows",
                })
                continue

            labs = m_pc_all["o2_subcompartment_final"].astype(str).to_numpy()[idx]
            n_labels = len(set(labs))
            if n_labels < 2:
                rows_sil.append({
                    "o2_compartment": o2,
                    "n_rows_used": int(len(idx)),
                    "n_subcompartments": int(n_labels),
                    "silhouette_unweighted": np.nan,
                    "silhouette_weighted_max_prob": np.nan,
                    "note": "only_one_subcompartment",
                })
                continue

            Xu = Xz_all[idx, :]
            sil_u = unweighted_silhouette(Xu, labs)

            sil_w = np.nan
            note = ""
            if weights_all is not None:
                sil_w = weighted_silhouette_precomputed(Xu, labs, weights_all[idx])
            else:
                note = "max_prob_missing"

            rows_sil.append({
                "o2_compartment": o2,
                "n_rows_used": int(len(idx)),
                "n_subcompartments": int(n_labels),
                "silhouette_unweighted": float(sil_u) if np.isfinite(sil_u) else np.nan,
                "silhouette_weighted_max_prob": float(sil_w) if np.isfinite(sil_w) else np.nan,
                "note": note,
            })

    pd.DataFrame(rows_sil).to_csv(os.path.join(tables_dir, "within_o2_silhouette_pcspace.csv"), index=False)

    # ----------------------------------------
    # GLOBAL silhouettes in PC space (PC-complete only)
    #   1) O2 compartments
    #   2) GMM components
    #   3) Combined hierarchical labels (o2_subcompartment_final)
    # ----------------------------------------
    sil_rows_global = []
    counts_rows_global = []

    if pc_complete_mask.sum() >= 5:
        m_pc_all = m.loc[pc_complete_mask].copy()
        Xz_all, _, _ = standardize_pc_space(m_pc_all, cfg.pc_cols)

        weights_all = None
        if cfg.max_prob_col in m_pc_all.columns:
            weights_all = pd.to_numeric(m_pc_all[cfg.max_prob_col], errors="coerce").to_numpy(dtype=float)

        # (A) O2 compartments
        sil_rows_global.append(
            compute_silhouette_bundle_pcspace(
                m_pc_all, Xz_all, label_col="o2_compartment", weights=weights_all
            )
        )
        counts_rows_global.append(
            m_pc_all["o2_compartment"].astype("object").fillna("NA").astype(str).value_counts().rename_axis("label").reset_index(name="n").assign(label_col="o2_compartment")
        )

        # (B) GMM components
        sil_rows_global.append(
            compute_silhouette_bundle_pcspace(
                m_pc_all, Xz_all, label_col="gmm_component", weights=weights_all
            )
        )
        counts_rows_global.append(
            m_pc_all["gmm_component"].astype("object").fillna("NA").astype(str).value_counts().rename_axis("label").reset_index(name="n").assign(label_col="gmm_component")
        )

        # (C) Combined hierarchical labels (FINAL)
        sil_rows_global.append(
            compute_silhouette_bundle_pcspace(
                m_pc_all, Xz_all, label_col="o2_subcompartment_final", weights=weights_all
            )
        )
        counts_rows_global.append(
            m_pc_all["o2_subcompartment_final"].astype("object").fillna("NA").astype(str).value_counts().rename_axis("label").reset_index(name="n").assign(label_col="o2_subcompartment_final")
        )

    # write summary + counts
    pd.DataFrame(sil_rows_global).to_csv(
        os.path.join(tables_dir, "silhouette_pcspace_global_summary.csv"),
        index=False,
    )

    if len(counts_rows_global) > 0:
        pd.concat(counts_rows_global, ignore_index=True).to_csv(
            os.path.join(tables_dir, "silhouette_pcspace_global_label_counts.csv"),
            index=False,
        )
    else:
        pd.DataFrame(columns=["label_col", "label", "n"]).to_csv(
            os.path.join(tables_dir, "silhouette_pcspace_global_label_counts.csv"),
            index=False,
        )

    # Human-friendly name (based on FINAL labels, within each O2 by abundance)
    name_rows = []
    for o2 in sorted(m["o2_compartment"].astype(str).unique(), key=str):
        sub = m.loc[m["o2_compartment"].astype(str) == o2, "o2_subcompartment_final"].astype(str)
        vc = sub.value_counts()
        for i, lab in enumerate(list(vc.index), start=1):
            name_rows.append(
                {
                    "o2_compartment": o2,
                    "o2_subcompartment_final": lab,
                    "o2_subcompartment_name": f"{o2}: subcompartment {i}",
                }
            )
    name_map = pd.DataFrame(name_rows).drop_duplicates(subset=["o2_compartment", "o2_subcompartment_final"])
    m = m.merge(name_map, on=["o2_compartment", "o2_subcompartment_final"], how="left")

    # write merged output
    m.to_csv(os.path.join(tables_dir, "merged_o2_split_by_gmm.csv"), index=False)

    # write qc + centroid/radius tables
    pd.DataFrame([qc]).to_csv(os.path.join(tables_dir, "reassignment_qc_summary.csv"), index=False)

    if len(centroids_df) > 0:
        centroids_df.to_csv(os.path.join(tables_dir, "reassignment_centroids.csv"), index=False)
    else:
        pd.DataFrame(columns=["o2_compartment", "o2_subcompartment", "pc", "centroid_z", "n_core"]).to_csv(
            os.path.join(tables_dir, "reassignment_centroids.csv"), index=False
        )

    if len(radii_df) > 0:
        radii_df.to_csv(os.path.join(tables_dir, "reassignment_cluster_radii.csv"), index=False)
    else:
        pd.DataFrame(columns=["o2_compartment", "o2_subcompartment", "n_all", "n_core", "radius_quantile", "radius"]).to_csv(
            os.path.join(tables_dir, "reassignment_cluster_radii.csv"), index=False
        )

    # ----------------------------
    # Plots (per O2 compartment)
    # ----------------------------
    if cfg.do_plots:

        # Build COMPLETE deterministic palette (covers every '<o2>__gmmK' and '<o2>__other')
        pal_df = build_full_subcompartment_palette(m, cfg, sat_floor=0.50)
        order_map = {v: i for i, v in enumerate(O2_COMPARTMENT_PALETTE.keys())}
        pal_df = pal_df.sort_values(
            by=["o2_compartment", "label"],
            key=lambda s: s.map(order_map) if s.name == "o2_compartment" else s
        )

        pal_df.to_csv(os.path.join(tables_dir, "subcompartment_palette.csv"), index=False)  # audit/stability
        sub_palette = palette_df_to_dict(pal_df)

        depth_for_curtain = (
            cfg.depth_anchored_col
            if cfg.depth_anchored_col in m.columns
            else cfg.depth_col
        )
        plot_hybrid_time_depth_curtain(
            m=m,
            depth_col=depth_for_curtain,
            sub_palette=sub_palette,
            plots_dir=plots_dir,
            tables_dir=tables_dir,
            cfg=cfg,
            maximum_depth=cfg.curtain_maximum_depth_m,
        )

        # Choose PC axes
        pc_x = cfg.pc_cols[0] if len(cfg.pc_cols) >= 1 else None
        pc_y = cfg.pc_cols[1] if len(cfg.pc_cols) >= 2 else None

        # Depth column to use for profiles (prefer anchored if present)
        depth_for_profile = cfg.depth_anchored_col if cfg.depth_anchored_col in m.columns else cfg.depth_col

        # EDA: all subcompartments together
        if pc_x is not None and pc_y is not None:
            plot_pc_scatter_all_subcompartments(m, pc_x, pc_y, sub_palette, plots_dir, cfg)

        depth_for_profile = cfg.depth_anchored_col if cfg.depth_anchored_col in m.columns else cfg.depth_col

        # EDA: depth profiles for each PC
        for pc in cfg.pc_cols:
            plot_depth_profile_all_subcompartments(m, pc, depth_for_profile, sub_palette, plots_dir, cfg)

        # EDA: depth profiles for each biochem feature present
        for feat in BIOCHEM_COLOR_MAP.keys():
            if feat in m.columns:
                plot_depth_profile_all_subcompartments(m, feat, depth_for_profile, sub_palette, plots_dir, cfg)

        # compute global PC limits across ALL data (or all data you want comparable)
        gx = pd.to_numeric(m[pc_x], errors="coerce")
        gy = pd.to_numeric(m[pc_y], errors="coerce")
        ok = gx.notna() & gy.notna()

        xmin, xmax = gx[ok].min(), gx[ok].max()
        ymin, ymax = gy[ok].min(), gy[ok].max()

        # optional padding
        dx = xmax - xmin
        dy = ymax - ymin
        xlim = (xmin - 0.05 * dx, xmax + 0.05 * dx)
        ylim = (ymin - 0.05 * dy, ymax + 0.05 * dy)

        for o2 in sorted(m["o2_compartment"].astype(str).unique(), key=str):
            # PC scatter by subcompartment
            if pc_x is not None and pc_y is not None:
                plot_pc_scatter_subcompartments(m, o2, pc_x, pc_y, sub_palette, plots_dir, cfg, xlim, ylim)

                # biochem overlays in PC space (only if feature exists)
                for feat in BIOCHEM_COLOR_MAP.keys():
                    if feat in m.columns:
                        plot_pc_scatter_biochem_overlay(m, o2, pc_x, pc_y, feat, plots_dir, cfg)

            # Depth profiles colored by subcompartment:
            # - for each PC used
            for pc in cfg.pc_cols:
                plot_depth_profile_subcompartments(m, o2, pc, depth_for_profile, sub_palette, plots_dir, cfg)

            # - for each biochem feature (if exists)
            for feat in BIOCHEM_COLOR_MAP.keys():
                if feat in m.columns:
                    plot_depth_profile_subcompartments(m, o2, feat, depth_for_profile, sub_palette, plots_dir, cfg)

        # ----------------------------------------
        # UMAP overlay of final subcompartments
        # ----------------------------------------
        if cfg.umap_embedding is not None and os.path.exists(cfg.umap_embedding):
            umap_df = read_umap_embedding(cfg.umap_embedding)
            plot_umap_subcompartments(
                m=m,
                umap_df=umap_df,
                sub_palette=sub_palette,
                plots_dir=plots_dir,
                cfg=cfg,
            )

    # ----------------------------
    # Save config
    # ----------------------------
    with open(os.path.join(cfg.outdir, "run_config.json"), "w") as f:
        json.dump({"config": cfg.__dict__}, f, indent=2)

    # ----------------------------
    # Console summary
    # ----------------------------
    print(f"[OK] Wrote outputs to: {cfg.outdir}")
    print(f"     Tables: {tables_dir}")
    print(f"     Plots : {plots_dir}" if cfg.do_plots else "     Plots : (disabled)")
    print(f"     key_mode={cfg.key_mode}  (derived key col: {cfg.derived_key_col})")
    print(f"     O2 subcompartments: intersection (O2 × GMM), collapse < {cfg.min_subcluster_size} → '<o2>__other'")
    if cfg.do_reassign:
        print(f"     Reassign: enabled  (borderline_mode={cfg.borderline_mode})")
        print(f"       borderline_max_prob={cfg.borderline_max_prob}  core_min_prob={cfg.core_min_prob}")
        print(f"       radius_quantile={cfg.reassign_radius_quantile}  min_core_n={cfg.reassign_min_core_n}")
        print(f"       reassigned={qc['n_reassigned']} of borderline={qc['n_borderline_candidates']} (PC-complete eligible={qc['n_reassignable_pc_complete']})")
    else:
        print("     Reassign: disabled")


if __name__ == "__main__":
    main()


"""
Minimal example command (with reassignment + plots + depth profiles):

python processes/split_o2_by_gmm/env_split_o2_by_gmm.py \
  --matrix-cleaned ../V4_ncbi_output/env_pca/tables/matrix_cleaned.csv \
  --eigenvectors ../V4_ncbi_output/env_pca/tables/eigenvectors_scores.csv \
  --assignments ../V4_ncbi_output/env_compartments_gmm/tables/compartments_assignments_smoothed.csv \
  --outdir ../V4_ncbi_output/env_o2_split_by_gmm \
  --sep-matrix ',' \
  --sep-eig ',' \
  --sep-assign ',' \
  --key-mode composite \
  --key-cols "Cruise,Year,Month,Day,Depth" \
  --pc-cols "PC1,PC2,PC3" \
  --min-subcluster-size 20 \
  --reassign \
  --core-min-prob 0.8 \
  --reassign-radius-quantile 0.95 \
  --reassign-min-core-n 20 \
  --plots \
  --plot-formats "png,pdf,svg"
"""
