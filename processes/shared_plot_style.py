#!/usr/bin/env python3
"""Repository-wide publication figure contract for BASINS.

The typography and vector-export settings follow ASPIRE's publication style.
The save hook is intentional: it reapplies the required font immediately
before export so later seaborn/theme calls or artist-level overrides cannot
silently reintroduce Matplotlib's default DejaVu Sans.

BASINS's scientific color mappings are left unchanged.  In particular, this
contract does not recolor heatmaps or categorical compartment figures.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import matplotlib as mpl
import numpy as np
from matplotlib import font_manager
from matplotlib.figure import Figure


FONT_FAMILY = "Times New Roman"
_INSTALLED = False
_ORIGINAL_SAVEFIG = Figure.savefig


def _register_system_font() -> list[Path]:
    """Register Times files even when a task-local Matplotlib cache omits them."""
    candidates: set[Path] = set()
    try:
        result = subprocess.run(
            ["fc-list", ":family=Times New Roman", "-f", "%{file}\n"],
            check=True,
            capture_output=True,
            text=True,
        )
        candidates.update(
            Path(line.strip())
            for line in result.stdout.splitlines()
            if line.strip()
        )
    except (FileNotFoundError, subprocess.CalledProcessError):
        pass

    for root in (Path("/usr/share/fonts"), Path("/usr/local/share/fonts")):
        if not root.is_dir():
            continue
        candidates.update(root.rglob("*Times*Roman*.ttf"))
        candidates.update(root.rglob("times*.ttf"))

    registered: list[Path] = []
    for path in sorted(candidates):
        if not path.is_file():
            continue
        try:
            font_manager.fontManager.addfont(str(path))
            registered.append(path)
        except (OSError, RuntimeError):
            continue
    return registered


def require_font() -> Path:
    """Resolve the exact required font, registering system files if necessary."""
    try:
        resolved = Path(
            font_manager.findfont(FONT_FAMILY, fallback_to_default=False)
        )
    except ValueError:
        _register_system_font()
        try:
            resolved = Path(
                font_manager.findfont(FONT_FAMILY, fallback_to_default=False)
            )
        except ValueError as exc:
            raise RuntimeError(
                f"Required publication font '{FONT_FAMILY}' is unavailable. "
                "Install the Microsoft core fonts and rebuild the font cache."
            ) from exc
    if not resolved.is_file():
        raise RuntimeError(f"Resolved publication font does not exist: {resolved}")
    return resolved


def apply_publication_contract(fig: Figure) -> None:
    """Enforce Times New Roman on every visible text artist before export."""
    for text in fig.findobj(match=mpl.text.Text):
        text.set_fontfamily(FONT_FAMILY)


def _spread_pixel_positions(
    desired: list[float],
    lower: float,
    upper: float,
    minimum_gap: float,
) -> list[float]:
    """Return ordered positions constrained to a collision-free pixel interval."""
    if not desired:
        return []
    if len(desired) == 1:
        return [min(max(desired[0], lower), upper)]

    available = max(upper - lower, 1.0)
    gap = min(float(minimum_gap), available / (len(desired) - 1))
    placed = [max(float(desired[0]), lower)]
    for value in desired[1:]:
        placed.append(max(float(value), placed[-1] + gap))

    if placed[-1] > upper:
        shift = placed[-1] - upper
        placed = [value - shift for value in placed]
    if placed[0] < lower:
        shift = lower - placed[0]
        placed = [value + shift for value in placed]

    # A backward pass protects the upper bound after the lower-bound shift.
    placed[-1] = min(placed[-1], upper)
    for index in range(len(placed) - 2, -1, -1):
        placed[index] = min(placed[index], placed[index + 1] - gap)
    return placed


def layout_biplot_labels(
    ax,
    text_artists: list,
    anchor_xy: list[tuple[float, float]],
    *,
    lane_padding_fraction: float = 0.04,
    vertical_padding_points: float = 3.0,
) -> None:
    """Place biplot labels in non-overlapping lanes outside all vector arrows.

    Labels are divided by the sign of their vector tip, distributed vertically
    in display coordinates, and connected back to their own vector.  Right-lane
    text extends only to the right and left-lane text only to the left, so the
    opaque label boxes cannot cover any arrow, whose x extent ends before its
    corresponding lane.
    """
    if not text_artists:
        return
    if len(text_artists) != len(anchor_xy):
        raise ValueError("Biplot text and anchor collections must have equal length.")

    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    x_limits = ax.get_xlim()
    x_span = max(abs(x_limits[1] - x_limits[0]), 1.0)
    lane_padding = lane_padding_fraction * x_span

    groups = {
        "left": [i for i, (x, _) in enumerate(anchor_xy) if x < 0],
        "right": [i for i, (x, _) in enumerate(anchor_xy) if x >= 0],
    }
    lane_x = {
        "left": min((anchor_xy[i][0] for i in groups["left"]), default=0.0)
        - lane_padding,
        "right": max((anchor_xy[i][0] for i in groups["right"]), default=0.0)
        + lane_padding,
    }

    axes_box = ax.get_window_extent(renderer=renderer)
    point_to_pixel = fig.dpi / 72.0
    margin_px = 5.0 * point_to_pixel
    lower = axes_box.y0 + margin_px
    upper = axes_box.y1 - margin_px

    for side, indices in groups.items():
        if not indices:
            continue
        indices = sorted(indices, key=lambda i: anchor_xy[i][1])
        desired = [
            float(ax.transData.transform(anchor_xy[index])[1])
            for index in indices
        ]
        max_height = max(
            text_artists[index].get_window_extent(renderer=renderer).height
            for index in indices
        )
        minimum_gap = max_height + vertical_padding_points * point_to_pixel
        placed = _spread_pixel_positions(desired, lower, upper, minimum_gap)

        for index, display_y in zip(indices, placed):
            _, data_y = ax.transData.inverted().transform(
                (ax.transData.transform((lane_x[side], 0.0))[0], display_y)
            )
            text = text_artists[index]
            text.set_position((lane_x[side], float(data_y)))
            text.set_ha("right" if side == "left" else "left")
            text.set_va("center")
            text.set_clip_on(False)
            patch = text.get_bbox_patch()
            if patch is not None:
                patch.set_facecolor("white")
                patch.set_edgecolor("none")
                patch.set_alpha(1.0)

            tip_x, tip_y = anchor_xy[index]
            ax.plot(
                [tip_x, lane_x[side]],
                [tip_y, float(data_y)],
                color=text.get_color(),
                linewidth=0.7,
                alpha=0.75,
                zorder=2.5,
            )


def calculate_biplot_vector_scale(
    x_values,
    y_values,
    vector_norms,
) -> tuple[float, float]:
    """Return the canonical BASINS biplot cloud and vector scaling factors."""
    x = np.asarray(x_values, dtype=float)
    y = np.asarray(y_values, dtype=float)
    xr = np.nanpercentile(x, 99) - np.nanpercentile(x, 1)
    yr = np.nanpercentile(y, 99) - np.nanpercentile(y, 1)
    cloud_scale = (
        0.35 * float(min(xr, yr))
        if np.isfinite(xr) and np.isfinite(yr)
        else 1.0
    )
    if not np.isfinite(cloud_scale) or cloud_scale <= 0:
        cloud_scale = 1.0

    norms = np.asarray(vector_norms, dtype=float)
    finite = norms[np.isfinite(norms)]
    denominator = float(np.max(finite)) if finite.size else 1.0
    if denominator <= 0:
        denominator = 1.0
    return cloud_scale, cloud_scale / denominator


def draw_biplot_vector(
    ax,
    tip_x: float,
    tip_y: float,
    color,
    vector_type: str,
    cloud_scale: float,
) -> None:
    """Draw the canonical outlined PCA-loading or sparse-correlation vector."""
    if vector_type == "Core loading":
        # Light under-stroke keeps the arrow distinct over dense point clouds.
        ax.arrow(
            0,
            0,
            tip_x,
            tip_y,
            length_includes_head=True,
            head_width=0.032 * cloud_scale,
            linewidth=3.2,
            color="white",
            zorder=3,
        )
        ax.arrow(
            0,
            0,
            tip_x,
            tip_y,
            length_includes_head=True,
            head_width=0.03 * cloud_scale,
            linewidth=2.2,
            color=color,
            zorder=4,
        )
        return

    # Sparse features use the PCA biplot's outlined dashed shaft and round tip.
    ax.plot(
        [0, tip_x],
        [0, tip_y],
        linestyle="--",
        linewidth=3.2,
        color="white",
        zorder=3,
    )
    ax.plot(
        [0, tip_x],
        [0, tip_y],
        linestyle="--",
        linewidth=2.2,
        color=color,
        zorder=4,
    )
    ax.scatter([tip_x], [tip_y], s=22, color=color, zorder=4)


def _contract_savefig(self: Figure, *args, **kwargs):
    apply_publication_contract(self)
    kwargs.setdefault("facecolor", "white")
    return _ORIGINAL_SAVEFIG(self, *args, **kwargs)


def install_publication_style() -> None:
    """Install the ASPIRE-compatible BASINS publication style once per process."""
    global _INSTALLED
    if _INSTALLED:
        return

    require_font()
    mpl.rcParams.update(
        {
            "font.family": FONT_FAMILY,
            "font.serif": [FONT_FAMILY],
            "font.size": 11,
            "axes.titlesize": 13,
            "axes.labelsize": 12,
            "axes.linewidth": 0.8,
            "axes.facecolor": "white",
            "axes.grid": False,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 10,
            "legend.title_fontsize": 11,
            "legend.frameon": False,
            "figure.titlesize": 14,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.dpi": 300,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "mathtext.fontset": "stix",
        }
    )
    Figure.savefig = _contract_savefig
    _INSTALLED = True
