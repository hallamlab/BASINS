#!/usr/bin/env python3
"""Benchmark BASINS cruise groups against season using independent chemistry."""

from __future__ import annotations

import argparse
from itertools import combinations
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import normalized_mutual_info_score

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shared_plot_export import save_figure_all_formats
from shared_plot_style import install_publication_style

install_publication_style()


SEASON_ORDER = {
    "djf": 0,
    "winter": 0,
    "mam": 1,
    "spring": 1,
    "jja": 2,
    "summer": 2,
    "son": 3,
    "fall": 3,
    "autumn": 3,
}


def parse_csv(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def ordered_season_labels(values) -> list[object]:
    """Return conventional winter-to-fall order with stable unknown labels."""
    return sorted(
        list(values),
        key=lambda value: (
            SEASON_ORDER.get(str(value).strip().lower(), len(SEASON_ORDER)),
            str(value),
        ),
    )


def normalize_id(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    try:
        number = float(text)
        if np.isfinite(number) and number.is_integer():
            return str(int(number))
    except ValueError:
        pass
    return text


def truthy(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def bias_corrected_cramers_v(left: pd.Series, right: pd.Series) -> float:
    table = pd.crosstab(left.astype(str), right.astype(str), dropna=False).to_numpy(dtype=float)
    n = float(table.sum())
    if n <= 1 or min(table.shape) < 2:
        return float("nan")
    expected = np.outer(table.sum(axis=1), table.sum(axis=0)) / n
    valid = expected > 0
    chi2 = float(np.sum(np.square(table[valid] - expected[valid]) / expected[valid]))
    rows, cols = table.shape
    phi2 = chi2 / n
    phi2_corrected = max(0.0, phi2 - ((cols - 1) * (rows - 1)) / (n - 1))
    rows_corrected = rows - np.square(rows - 1) / (n - 1)
    cols_corrected = cols - np.square(cols - 1) / (n - 1)
    denominator = min(rows_corrected - 1, cols_corrected - 1)
    return float(np.sqrt(phi2_corrected / denominator)) if denominator > 0 else float("nan")


def season_group_redundancy(
    frame: pd.DataFrame,
    season_col: str,
    group_col: str,
    year_col: str,
    permutations: int,
    seed: int,
) -> pd.DataFrame:
    matched = frame[[season_col, group_col, year_col]].dropna().reset_index(drop=True)
    season = matched[season_col].astype(str)
    group = matched[group_col].astype(str)
    observed_nmi = float(normalized_mutual_info_score(season, group))
    observed_v = bias_corrected_cramers_v(season, group)

    rng = np.random.default_rng(seed)
    years = matched[year_col].astype(str).to_numpy()
    group_values = group.to_numpy()
    null_values = []
    for _ in range(permutations):
        permuted = group_values.copy()
        for year in pd.unique(years):
            indices = np.flatnonzero(years == year)
            if len(indices) > 1:
                permuted[indices] = rng.permutation(permuted[indices])
        null_values.append(bias_corrected_cramers_v(season, pd.Series(permuted)))
    null = np.asarray(null_values, dtype=float)
    valid = np.isfinite(null)
    p_value = (
        float((1 + np.sum(null[valid] >= observed_v)) / (1 + valid.sum()))
        if np.isfinite(observed_v) and valid.any()
        else np.nan
    )
    return pd.DataFrame([{
        "n_cruises": int(len(matched)),
        "n_seasons": int(season.nunique()),
        "n_cruise_groups": int(group.nunique()),
        "normalized_mutual_information": observed_nmi,
        "cramers_v_bias_corrected": observed_v,
        "within_year_permutation_p_value": p_value,
        "permutations_completed": int(valid.sum()),
    }])


def aggregate_sparse_profiles(
    matrix: pd.DataFrame,
    features: list[str],
    cruise_col: str,
    depth_col: str,
) -> pd.DataFrame:
    """Give each observed anchored depth equal weight within each cruise."""
    parts = []
    for feature in features:
        if feature not in matrix:
            continue
        current = matrix[[cruise_col, depth_col, feature]].copy()
        current[feature] = pd.to_numeric(current[feature], errors="coerce")
        current[depth_col] = pd.to_numeric(current[depth_col], errors="coerce")
        current = current.replace([np.inf, -np.inf], np.nan).dropna()
        if current.empty:
            continue
        by_depth = (
            current.groupby([cruise_col, depth_col], as_index=False)[feature]
            .median()
        )
        summary = (
            by_depth.groupby(cruise_col)
            .agg(
                **{
                    feature: (feature, "median"),
                    f"{feature}__n_depths": (depth_col, "nunique"),
                    f"{feature}__min_depth_m": (depth_col, "min"),
                    f"{feature}__max_depth_m": (depth_col, "max"),
                }
            )
            .reset_index()
        )
        parts.append(summary)
    if not parts:
        return pd.DataFrame(columns=[cruise_col])
    merged = parts[0]
    for part in parts[1:]:
        merged = merged.merge(part, on=cruise_col, how="outer")
    return merged


def train_test_design(
    train: pd.DataFrame,
    test: pd.DataFrame,
    categorical: list[str],
) -> tuple[np.ndarray, np.ndarray]:
    train_parts = [np.ones((len(train), 1), dtype=float)]
    test_parts = [np.ones((len(test), 1), dtype=float)]
    for column in categorical:
        train_values = train[column].astype(str)
        test_values = test[column].astype(str)
        levels = sorted(train_values.unique().tolist())
        for level in levels[1:]:
            train_parts.append(train_values.eq(level).to_numpy(dtype=float)[:, None])
            test_parts.append(test_values.eq(level).to_numpy(dtype=float)[:, None])
    return np.column_stack(train_parts), np.column_stack(test_parts)


def leave_year_out_cv(
    frame: pd.DataFrame,
    feature: str,
    season_col: str,
    group_col: str,
    year_col: str,
    analysis_set: str,
) -> pd.DataFrame:
    models = {
        "Season": [season_col],
        "Cruise group": [group_col],
        "Season + cruise group": [season_col, group_col],
    }
    y_raw = pd.to_numeric(frame[feature], errors="coerce").to_numpy(dtype=float)
    y = np.sign(y_raw) * np.log1p(np.abs(y_raw))
    years = frame[year_col].astype(str)
    rows = []
    for fold_year in sorted(years.unique()):
        test = years.eq(fold_year).to_numpy()
        train = ~test
        if test.sum() < 1 or train.sum() < 5:
            continue
        train_mean = float(np.mean(y[train]))
        baseline_sse = float(np.square(y[test] - train_mean).sum())
        for model, categorical in models.items():
            x_train, x_test = train_test_design(
                frame.loc[train],
                frame.loc[test],
                categorical,
            )
            coefficients = np.linalg.pinv(x_train) @ y[train]
            predicted = x_test @ coefficients
            sse = float(np.square(y[test] - predicted).sum())
            rows.append({
                "analysis_set": analysis_set,
                "feature": feature,
                "fold_year": fold_year,
                "model": model,
                "n_train": int(train.sum()),
                "n_test": int(test.sum()),
                "sse": sse,
                "intercept_baseline_sse": baseline_sse,
                "cv_r2": 1.0 - sse / baseline_sse if baseline_sse > 0 else np.nan,
                "rmse_log1p": float(np.sqrt(sse / test.sum())),
            })
    return pd.DataFrame(rows)


def summarize_cv(folds: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for keys, group in folds.groupby(["analysis_set", "feature", "model"], sort=False):
        analysis_set, feature, model = keys
        total_sse = float(group["sse"].sum())
        total_baseline = float(group["intercept_baseline_sse"].sum())
        total_n = int(group["n_test"].sum())
        rows.append({
            "analysis_set": analysis_set,
            "feature": feature,
            "model": model,
            "folds": int(group["fold_year"].nunique()),
            "heldout_cruises": total_n,
            "pooled_cv_r2": 1.0 - total_sse / total_baseline if total_baseline > 0 else np.nan,
            "fold_cv_r2_mean": float(group["cv_r2"].mean()),
            "fold_cv_r2_median": float(group["cv_r2"].median()),
            "pooled_rmse_log1p": float(np.sqrt(total_sse / total_n)) if total_n else np.nan,
        })
    return pd.DataFrame(rows)


def paired_comparisons(
    folds: pd.DataFrame,
    bootstrap: int,
    permutations: int,
    seed: int,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    for keys, current in folds.groupby(["analysis_set", "feature"]):
        analysis_set, feature = keys
        pivot = current.pivot(index="fold_year", columns="model", values="cv_r2")
        for model_a, model_b in combinations(pivot.columns.tolist(), 2):
            paired = pivot[[model_a, model_b]].dropna()
            if paired.empty:
                continue
            differences = (paired[model_a] - paired[model_b]).to_numpy(dtype=float)
            boot = np.array([
                rng.choice(differences, size=len(differences), replace=True).mean()
                for _ in range(bootstrap)
            ])
            null = np.array([
                np.mean(differences * rng.choice([-1.0, 1.0], size=len(differences)))
                for _ in range(permutations)
            ])
            observed = float(differences.mean())
            rows.append({
                "analysis_set": analysis_set,
                "feature": feature,
                "model_a": model_a,
                "model_b": model_b,
                "metric": "fold_cv_r2",
                "n_paired_years": int(len(differences)),
                "mean_difference_a_minus_b": observed,
                "ci_lower": float(np.quantile(boot, 0.025)),
                "ci_upper": float(np.quantile(boot, 0.975)),
                "sign_flip_p_value": float(
                    (1 + np.sum(np.abs(null) >= abs(observed))) / (permutations + 1)
                ),
            })
    return pd.DataFrame(rows)


def bh_adjust(values: pd.Series) -> pd.Series:
    result = pd.Series(np.nan, index=values.index, dtype=float)
    valid = values.dropna().astype(float)
    if valid.empty:
        return result
    ordered = valid.sort_values()
    adjusted = (
        ordered.to_numpy() * len(ordered) / np.arange(1, len(ordered) + 1)
    )[::-1]
    adjusted = np.minimum.accumulate(adjusted)[::-1]
    result.loc[ordered.index] = np.minimum(adjusted, 1.0)
    return result


def plot_redundancy(
    counts: pd.DataFrame,
    out_base: Path,
) -> None:
    counts = counts.reindex(ordered_season_labels(counts.index))
    proportions = counts.div(counts.sum(axis=1).replace(0, np.nan), axis=0)
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.2))
    image = axes[0].imshow(
        proportions.to_numpy(), vmin=0, vmax=1, cmap="Greys"
    )
    axes[0].set_xticks(range(len(proportions.columns)))
    axes[0].set_xticklabels(proportions.columns, rotation=35, ha="right")
    axes[0].set_yticks(range(len(proportions.index)))
    axes[0].set_yticklabels(proportions.index)
    axes[0].set_xlabel("Cruise group")
    axes[0].set_ylabel("Season")
    axes[0].set_title("Within-season cruise-group proportions")
    for i in range(len(proportions.index)):
        for j in range(len(proportions.columns)):
            value = proportions.iloc[i, j]
            axes[0].text(
                j, i, f"{int(counts.iloc[i, j])}", ha="center", va="center",
                color="white" if value >= 0.55 else "black",
            )
    fig.colorbar(image, ax=axes[0], fraction=0.046, pad=0.04)

    bottom = np.zeros(len(proportions))
    group_colors = dict(zip(
        proportions.columns,
        plt.cm.Greys(np.linspace(0.35, 0.85, len(proportions.columns))),
    ))
    for group in reversed(proportions.columns.tolist()):
        values = proportions[group].to_numpy(dtype=float)
        axes[1].bar(
            proportions.index,
            values,
            bottom=bottom,
            label=group,
            color=group_colors[group],
            edgecolor="black",
            linewidth=0.5,
        )
        bottom += values
    axes[1].set_ylim(0, 1)
    axes[1].set_ylabel("Fraction of cruises")
    axes[1].set_title("Cruise-group composition by season")
    handles, labels = axes[1].get_legend_handles_labels()
    handle_by_label = dict(zip(labels, handles))
    axes[1].legend(
        [handle_by_label[group] for group in proportions.columns],
        proportions.columns.tolist(),
        frameon=False,
        title="Cruise group",
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        borderaxespad=0.0,
    )
    fig.tight_layout(rect=(0.0, 0.0, 0.86, 1.0))
    save_figure_all_formats(fig, out_base, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_cv(summary: pd.DataFrame, out_base: Path) -> None:
    if summary.empty:
        return
    analysis_sets = summary["analysis_set"].drop_duplicates().tolist()
    features = summary["feature"].drop_duplicates().tolist()
    models = ["Season", "Cruise group", "Season + cruise group"]
    fig, axes = plt.subplots(
        len(analysis_sets),
        len(features),
        figsize=(5.1 * len(features), 4.7 * len(analysis_sets)),
        squeeze=False,
        sharey=True,
    )
    colors = ["#BDBDBD", "#737373", "#252525"]
    for row, analysis_set in enumerate(analysis_sets):
        for col, feature in enumerate(features):
            ax = axes[row, col]
            current = summary[
                summary["analysis_set"].eq(analysis_set)
                & summary["feature"].eq(feature)
            ].set_index("model")
            values = [
                current.loc[model, "pooled_cv_r2"] if model in current.index else np.nan
                for model in models
            ]
            ax.bar(range(len(models)), values, color=colors, edgecolor="black", linewidth=0.5)
            ax.axhline(0, color="black", linewidth=0.8)
            ax.set_xticks(range(len(models)))
            ax.set_xticklabels(models, rotation=30, ha="right")
            ax.set_title(f"{feature}\n{analysis_set.replace('_', ' ')}")
            ax.set_ylabel("Leave-one-year-out pooled R²")
    fig.suptitle("Season and cruise-group prediction of PCA-excluded chemistry")
    fig.tight_layout()
    save_figure_all_formats(fig, out_base, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--assignments", type=Path, required=True)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--assignments-sep", default="\t")
    parser.add_argument("--matrix-sep", default=",")
    parser.add_argument("--cruise-col", default="Cruise")
    parser.add_argument("--season-col", default="Season")
    parser.add_argument("--year-col", default="Year")
    parser.add_argument("--group-col", default="cruise_group")
    parser.add_argument("--uncertain-col", default="assignment_uncertain")
    parser.add_argument("--depth-col", default="Depth_anchored")
    parser.add_argument(
        "--sparse-features",
        default="Nitrous Oxide,Hydrogen Sulfide,Methane,Iron,Dimethyl Sulfide",
    )
    parser.add_argument("--min-depths", type=int, default=3)
    parser.add_argument("--bootstrap", type=int, default=5000)
    parser.add_argument("--permutations", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    if args.min_depths < 1 or args.bootstrap < 1 or args.permutations < 1:
        raise ValueError("Minimum depths, bootstrap, and permutation counts must be positive.")

    tables = args.outdir / "tables"
    plots = args.outdir / "plots"
    tables.mkdir(parents=True, exist_ok=True)
    plots.mkdir(parents=True, exist_ok=True)

    assignments = pd.read_csv(args.assignments, sep=args.assignments_sep)
    matrix = pd.read_csv(args.matrix, sep=args.matrix_sep)
    required_assignments = [
        args.cruise_col, args.season_col, args.year_col, args.group_col
    ]
    missing = [column for column in required_assignments if column not in assignments]
    if missing:
        raise ValueError(f"Cruise assignments missing required columns: {missing}")
    if args.cruise_col not in matrix or args.depth_col not in matrix:
        raise ValueError("Sparse matrix is missing the cruise or anchored-depth column.")

    assignments = assignments.drop_duplicates(args.cruise_col).copy()
    assignments["__cruise_id__"] = assignments[args.cruise_col].map(normalize_id)
    matrix = matrix.copy()
    matrix["__cruise_id__"] = matrix[args.cruise_col].map(normalize_id)

    redundancy = season_group_redundancy(
        assignments,
        args.season_col,
        args.group_col,
        args.year_col,
        args.permutations,
        args.seed,
    )
    counts = pd.crosstab(
        assignments[args.season_col].astype(str),
        assignments[args.group_col].astype(str),
    )
    proportions = counts.div(counts.sum(axis=1).replace(0, np.nan), axis=0)
    redundancy.to_csv(tables / "cruise_group_season_redundancy.tsv", sep="\t", index=False)
    counts.to_csv(tables / "cruise_group_by_season_counts.tsv", sep="\t")
    proportions.to_csv(tables / "cruise_group_by_season_row_proportions.tsv", sep="\t")
    plot_redundancy(counts, plots / "cruise_group_season_redundancy")

    features = parse_csv(args.sparse_features)
    profiles = aggregate_sparse_profiles(
        matrix,
        features,
        "__cruise_id__",
        args.depth_col,
    )
    profiles.to_csv(tables / "cruise_sparse_chemistry_profiles.tsv", sep="\t", index=False)

    metadata_cols = [
        "__cruise_id__", args.cruise_col, args.season_col, args.year_col,
        args.group_col,
    ]
    if args.uncertain_col in assignments:
        metadata_cols.append(args.uncertain_col)
    analysis = assignments[metadata_cols].merge(
        profiles, on="__cruise_id__", how="left"
    )

    fold_parts = []
    audit_rows = []
    for feature in features:
        depth_count_col = f"{feature}__n_depths"
        if feature not in analysis or depth_count_col not in analysis:
            audit_rows.append({
                "analysis_set": "all_assignments",
                "feature": feature,
                "status": "skipped_missing_column",
                "n_cruises": 0,
                "n_years": 0,
            })
            continue
        base = analysis[
            pd.to_numeric(analysis[depth_count_col], errors="coerce") >= args.min_depths
        ].dropna(
            subset=[feature, args.season_col, args.year_col, args.group_col]
        ).copy()
        cohorts = {"all_assignments": base}
        if args.uncertain_col in base:
            cohorts["confident_assignments"] = base.loc[~truthy(base[args.uncertain_col])].copy()
        for analysis_set, cohort in cohorts.items():
            status = "ok" if cohort[args.year_col].nunique() >= 3 else "insufficient_years"
            audit_rows.append({
                "analysis_set": analysis_set,
                "feature": feature,
                "status": status,
                "n_cruises": int(len(cohort)),
                "n_years": int(cohort[args.year_col].nunique()),
                "n_seasons": int(cohort[args.season_col].nunique()),
                "n_cruise_groups": int(cohort[args.group_col].nunique()),
                "minimum_observed_depths": args.min_depths,
                "response_transform": "signed_log1p_of_equal_depth_weighted_median",
            })
            if status == "ok":
                folds = leave_year_out_cv(
                    cohort,
                    feature,
                    args.season_col,
                    args.group_col,
                    args.year_col,
                    analysis_set,
                )
                if not folds.empty:
                    fold_parts.append(folds)

    cv_folds = pd.concat(fold_parts, ignore_index=True) if fold_parts else pd.DataFrame()
    cv_summary = summarize_cv(cv_folds) if not cv_folds.empty else pd.DataFrame()
    comparisons = (
        paired_comparisons(
            cv_folds,
            args.bootstrap,
            args.permutations,
            args.seed,
        )
        if not cv_folds.empty
        else pd.DataFrame()
    )
    if not comparisons.empty:
        comparisons["sign_flip_q_value"] = comparisons.groupby(
            "analysis_set", group_keys=False
        )["sign_flip_p_value"].transform(bh_adjust)

    cv_folds.to_csv(tables / "cruise_group_season_sparse_cv_folds.tsv", sep="\t", index=False)
    cv_summary.to_csv(tables / "cruise_group_season_sparse_cv_summary.tsv", sep="\t", index=False)
    comparisons.to_csv(
        tables / "cruise_group_season_sparse_paired_comparisons.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame(audit_rows).to_csv(
        tables / "cruise_group_season_sparse_cohort_audit.tsv",
        sep="\t",
        index=False,
    )
    plot_cv(cv_summary, plots / "cruise_group_season_sparse_cv_performance")
    (tables / "cruise_group_season_benchmark_run_config.json").write_text(
        json.dumps(vars(args), default=str, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
