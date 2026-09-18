#!/usr/bin/env python3
"""Build cruise EOF scores from departures around a local multi-year baseline."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--profiles", type=Path, required=True)
    ap.add_argument("--metadata", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--cruise-col", default="Cruise")
    ap.add_argument("--date-col", default="date")
    ap.add_argument("--baseline-months", type=int, default=24)
    ap.add_argument("--baseline-min-cruises", type=int, default=5)
    ap.add_argument("--min-column-coverage", type=float, default=0.50)
    ap.add_argument("--n-components", type=int, default=10)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    profiles = pd.read_csv(args.profiles, sep="\t")
    metadata = pd.read_csv(args.metadata)
    metadata = metadata.drop_duplicates(args.cruise_col)
    meta_cols = [c for c in [args.cruise_col, args.date_col, "Season", "Year", "Month", "Day"] if c in metadata]
    frame = profiles.merge(metadata[meta_cols], on=args.cruise_col, how="left")
    frame[args.date_col] = pd.to_datetime(frame[args.date_col], errors="coerce")
    frame = frame.dropna(subset=[args.date_col]).sort_values(args.date_col).reset_index(drop=True)

    feature_cols = [c for c in profiles.columns if c != args.cruise_col]
    coverage = frame[feature_cols].notna().mean()
    retained = coverage[coverage >= args.min_column_coverage].index.tolist()
    if len(retained) < 2:
        raise ValueError("Fewer than two cruise-profile columns passed coverage filtering")

    dated = frame.set_index(args.date_col)
    window_days = max(1, int(round(args.baseline_months * 365.25 / 12)))
    baseline = dated[retained].rolling(
        f"{window_days}D", center=True, min_periods=args.baseline_min_cruises
    ).median()
    residual = dated[retained] - baseline
    valid_baseline = baseline.notna().sum(axis=1)
    keep = valid_baseline >= max(2, int(np.ceil(0.50 * len(retained))))
    dated = dated.loc[keep].copy()
    baseline = baseline.loc[keep].copy()
    residual = residual.loc[keep].copy()

    impute = residual.median(axis=0)
    residual_imputed = residual.fillna(impute)
    scaler = StandardScaler()
    scaled = scaler.fit_transform(residual_imputed)
    n_components = min(args.n_components, scaled.shape[0], scaled.shape[1])
    pca = PCA(n_components=n_components, random_state=args.seed)
    scores = pca.fit_transform(scaled)

    metadata_out = dated[[c for c in meta_cols if c != args.date_col]].copy()
    metadata_out.insert(1 if args.cruise_col in metadata_out else 0, args.date_col, dated.index)
    metadata_out = metadata_out.reset_index(drop=True)
    for idx in range(n_components):
        metadata_out[f"PC{idx + 1}"] = scores[:, idx]
    metadata_out.to_csv(args.outdir / "local_eof_scores_by_cruise.csv", index=False)

    keys = dated[[args.cruise_col]].reset_index()
    residual_out = keys.merge(residual.reset_index(), on=args.date_col, how="left")
    residual_out.to_csv(args.outdir / "local_cruise_profile_residuals.tsv", sep="\t", index=False)
    keys.merge(baseline.reset_index(), on=args.date_col, how="left").to_csv(
        args.outdir / "local_cruise_profile_baselines.tsv", sep="\t", index=False
    )
    pd.DataFrame(pca.components_.T, index=retained,
                 columns=[f"PC{i + 1}" for i in range(n_components)]).to_csv(
        args.outdir / "local_eof_loadings.tsv", sep="\t"
    )
    pd.DataFrame({
        "PC": [f"PC{i + 1}" for i in range(n_components)],
        "explained_variance_ratio": pca.explained_variance_ratio_,
        "cumulative_ratio": np.cumsum(pca.explained_variance_ratio_),
    }).to_csv(args.outdir / "local_eof_explained_variance.tsv", sep="\t", index=False)
    pd.DataFrame({"feature_depth": coverage.index, "coverage": coverage.values,
                  "retained": coverage.index.isin(retained)}).to_csv(
        args.outdir / "local_profile_column_coverage.tsv", sep="\t", index=False
    )
    pd.DataFrame([{
        "input_cruises": len(frame), "retained_cruises": len(metadata_out),
        "profile_columns_input": len(feature_cols), "profile_columns_retained": len(retained),
        "baseline_months": args.baseline_months,
        "baseline_min_cruises": args.baseline_min_cruises,
    }]).to_csv(args.outdir / "local_eof_input_audit.tsv", sep="\t", index=False)
    (args.outdir / "local_eof_run_config.json").write_text(json.dumps(vars(args), default=str, indent=2) + "\n")


if __name__ == "__main__":
    main()
