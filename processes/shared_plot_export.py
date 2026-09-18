#!/usr/bin/env python3
"""Shared BASINS helpers for consistent static-figure exports."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable


STATIC_PLOT_FORMATS = ("png", "pdf", "svg")


def save_figure_all_formats(
    figure: Any,
    path: str | Path,
    *,
    dpi: int = 300,
    formats: Iterable[str] = STATIC_PLOT_FORMATS,
    **kwargs: Any,
) -> tuple[Path, ...]:
    """Save a Matplotlib figure beside ``path`` in every requested format."""
    requested = tuple(dict.fromkeys(str(fmt).lower().lstrip(".") for fmt in formats))
    unsupported = sorted(set(requested) - set(STATIC_PLOT_FORMATS))
    if unsupported:
        raise ValueError(f"Unsupported static plot formats: {unsupported}")

    base = Path(path)
    if base.suffix.lower().lstrip(".") in STATIC_PLOT_FORMATS:
        base = base.with_suffix("")
    outputs = tuple(base.with_suffix(f".{fmt}") for fmt in requested)
    for output in outputs:
        save_kwargs = dict(kwargs)
        if output.suffix.lower() == ".png":
            save_kwargs.setdefault("dpi", dpi)
        figure.savefig(output, **save_kwargs)
    return outputs
