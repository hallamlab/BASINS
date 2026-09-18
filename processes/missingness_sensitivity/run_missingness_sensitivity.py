#!/usr/bin/env python3
"""Run and summarize a BASINS feature-missingness threshold sensitivity analysis."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shared_plot_export import save_figure_all_formats
from shared_plot_style import install_publication_style

install_publication_style()


DEFAULT_SELECTION_COVARIANCE_TYPE = "tied"
DEFAULT_FINAL_COVARIANCE_TYPE = "tied"


def parse_thresholds(raw: str) -> list[float]:
    values = sorted({float(piece.strip()) for piece in raw.split(",") if piece.strip()})
    if not values:
        raise ValueError("At least one missingness threshold is required.")
    if any(not np.isfinite(value) or value <= 0 or value >= 1 for value in values):
        raise ValueError(
            "Missingness thresholds must be finite and strictly between 0 and 1."
        )
    return values


def threshold_label(value: float) -> str:
    return f"cutoff_{value:.3f}".rstrip("0").rstrip(".").replace(".", "p")


def bool_series(values: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(values):
        return values.fillna(False)
    return values.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def run_command(command: list[str], stdout_path: Path, stderr_path: Path) -> int:
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    with stdout_path.open("w") as stdout, stderr_path.open("w") as stderr:
        result = subprocess.run(command, stdout=stdout, stderr=stderr, check=False)
    return int(result.returncode)


def required_outputs_exist(paths: Iterable[Path]) -> bool:
    return all(path.is_file() and path.stat().st_size > 0 for path in paths)


def read_json(path: Path) -> dict:
    with path.open() as handle:
        return json.load(handle)


def selected_row(path: Path) -> pd.Series | None:
    if not path.is_file():
        return None
    table = pd.read_csv(path)
    selected = table[bool_series(table["SELECTED"])]
    return selected.iloc[0] if len(selected) == 1 else None


def assignment_table(path: Path, id_col: str) -> pd.DataFrame | None:
    if not path.is_file():
        return None
    table = pd.read_csv(path)
    if id_col not in table or "component" not in table:
        return None
    out = table[[id_col, "component"]].dropna().copy()
    out[id_col] = out[id_col].astype(str)
    out = out[~out[id_col].duplicated(keep=False)]
    return out


def retained_loading_basis(run_dir: Path) -> tuple[list[str], np.ndarray] | None:
    loadings_path = run_dir / "pca" / "tables" / "pca_loadings.csv"
    keep_path = run_dir / "pca" / "tables" / "pc_keep_decision.csv"
    if not loadings_path.is_file() or not keep_path.is_file():
        return None
    loadings = pd.read_csv(loadings_path)
    keep = pd.read_csv(keep_path)
    if "feature" not in loadings or "PC" not in keep or "KEEP" not in keep:
        return None
    pcs = keep.loc[bool_series(keep["KEEP"]), "PC"].astype(str).tolist()
    pcs = [pc for pc in pcs if pc in loadings]
    if not pcs:
        return None
    return loadings["feature"].astype(str).tolist(), loadings[pcs].to_numpy(float)


def principal_angles(
    left: tuple[list[str], np.ndarray] | None,
    right: tuple[list[str], np.ndarray] | None,
) -> tuple[float, float, int]:
    if left is None or right is None:
        return np.nan, np.nan, 0
    left_features, left_matrix = left
    right_features, right_matrix = right
    common = sorted(set(left_features).intersection(right_features))
    if len(common) < 2:
        return np.nan, np.nan, len(common)
    left_index = {feature: index for index, feature in enumerate(left_features)}
    right_index = {feature: index for index, feature in enumerate(right_features)}
    a = left_matrix[[left_index[feature] for feature in common], :]
    b = right_matrix[[right_index[feature] for feature in common], :]
    rank = min(a.shape[1], b.shape[1], len(common))
    if rank == 0:
        return np.nan, np.nan, len(common)
    qa, _ = np.linalg.qr(a)
    qb, _ = np.linalg.qr(b)
    singular = np.linalg.svd(qa[:, :rank].T @ qb[:, :rank], compute_uv=False)
    angles = np.degrees(np.arccos(np.clip(singular, -1, 1)))
    return float(np.mean(angles)), float(np.max(angles)), len(common)


def summarize_runs(
    outdir: Path,
    thresholds: Iterable[float],
    statuses: dict[float, dict],
    id_col: str,
    *,
    min_cluster_frac: float,
    min_gmm_stability_ari: float,
) -> float:
    tables_dir = outdir / "tables"
    plots_dir = outdir / "plots"
    tables_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)

    summary_rows: list[dict] = []
    feature_frames: list[pd.DataFrame] = []
    assignments: dict[float, pd.DataFrame] = {}
    bases: dict[float, tuple[list[str], np.ndarray]] = {}

    for threshold in thresholds:
        run_dir = outdir / "runs" / threshold_label(threshold)
        pca_dir = run_dir / "pca"
        select_dir = run_dir / "selectk"
        gmm_dir = run_dir / "gmm"
        state = statuses[threshold]
        row = {
            "missingness_cutoff": threshold,
            "run_label": threshold_label(threshold),
            "status": state["status"],
            "failed_stage": state.get("failed_stage", ""),
            "pca_exit_code": state.get("pca_exit_code", np.nan),
            "selectk_exit_code": state.get("selectk_exit_code", np.nan),
            "gmm_exit_code": state.get("gmm_exit_code", np.nan),
        }

        feature_path = pca_dir / "tables" / "core_vs_sparse_features.csv"
        if feature_path.is_file():
            feature = pd.read_csv(feature_path)
            feature.insert(0, "missingness_cutoff", threshold)
            feature_frames.append(feature)
            row["n_core_features"] = int((feature["status"] == "core").sum())
            row["n_sparse_features"] = int((feature["status"] == "sparse").sum())
            row["core_features"] = "|".join(
                feature.loc[feature["status"] == "core", "feature"].astype(str)
            )
            row["sparse_features"] = "|".join(
                feature.loc[feature["status"] == "sparse", "feature"].astype(str)
            )

        qc_path = pca_dir / "qc_summary.json"
        if qc_path.is_file():
            qc = read_json(qc_path)
            for source, destination in (
                ("n_rows_input", "n_rows_input"),
                ("n_rows_kept", "n_rows_pca"),
                ("n_rows_dropped_pre", "n_rows_dropped_pre"),
                ("n_rows_dropped_post", "n_rows_dropped_post"),
            ):
                row[destination] = qc.get(source, np.nan)

        keep_path = pca_dir / "tables" / "pc_keep_decision.csv"
        if keep_path.is_file():
            keep = pd.read_csv(keep_path)
            kept_pcs = keep.loc[bool_series(keep["KEEP"]), "PC"].astype(str).tolist()
            row["n_pcs_kept"] = len(kept_pcs)
            row["pcs_kept"] = "|".join(kept_pcs)

        variance_path = pca_dir / "tables" / "pca_explained_variance.csv"
        if variance_path.is_file():
            variance = pd.read_csv(variance_path)
            row["pc1_variance_ratio"] = (
                float(variance.iloc[0]["explained_variance_ratio"]) if len(variance) else np.nan
            )
            row["pc1_pc2_cumulative_ratio"] = (
                float(variance.iloc[min(1, len(variance) - 1)]["cumulative_ratio"])
                if len(variance)
                else np.nan
            )

        choice = selected_row(select_dir / "tables" / "gmm_k_selection_decision.csv")
        if choice is not None:
            for source, destination in (
                ("K", "selected_k"),
                ("BIC", "selected_bic"),
                ("ICL", "selected_icl"),
                ("mean_resp_entropy", "selected_mean_resp_entropy"),
                ("min_cluster_frac", "selected_min_cluster_frac"),
                ("stability_median_ARI", "selected_stability_median_ari"),
                ("stability_n_reps", "selected_stability_n_reps"),
            ):
                row[destination] = choice.get(source, np.nan)

        if state["status"] == "complete":
            assignment = assignment_table(
                gmm_dir / "tables" / "compartments_assignments_smoothed.csv", id_col
            )
            if assignment is not None:
                assignments[threshold] = assignment
                row["n_unique_assignment_ids"] = len(assignment)
            basis = retained_loading_basis(run_dir)
            if basis is not None:
                bases[threshold] = basis
        summary_rows.append(row)

    summary = pd.DataFrame(summary_rows).sort_values("missingness_cutoff")
    summary.to_csv(tables_dir / "threshold_run_summary.tsv", sep="\t", index=False)

    feature_long = (
        pd.concat(feature_frames, ignore_index=True)
        if feature_frames
        else pd.DataFrame(columns=["missingness_cutoff", "feature", "status"])
    )
    feature_long.to_csv(tables_dir / "threshold_feature_status.tsv", sep="\t", index=False)

    pairwise_columns = [
        "cutoff_left",
        "cutoff_right",
        "adjacent_cutoffs",
        "n_shared_assignments",
        "adjusted_rand_index",
        "n_common_loading_features",
        "mean_principal_angle_degrees",
        "max_principal_angle_degrees",
    ]
    pairwise_rows = []
    successful = sorted(assignments)
    for left_index, left_threshold in enumerate(successful):
        for right_threshold in successful[left_index + 1 :]:
            shared = assignments[left_threshold].merge(
                assignments[right_threshold],
                on=id_col,
                suffixes=("_left", "_right"),
            )
            mean_angle, max_angle, n_common_features = principal_angles(
                bases.get(left_threshold), bases.get(right_threshold)
            )
            pairwise_rows.append(
                {
                    "cutoff_left": left_threshold,
                    "cutoff_right": right_threshold,
                    "adjacent_cutoffs": bool(
                        successful.index(right_threshold)
                        == successful.index(left_threshold) + 1
                    ),
                    "n_shared_assignments": len(shared),
                    "adjusted_rand_index": (
                        adjusted_rand_score(
                            shared["component_left"], shared["component_right"]
                        )
                        if len(shared) >= 2
                        else np.nan
                    ),
                    "n_common_loading_features": n_common_features,
                    "mean_principal_angle_degrees": mean_angle,
                    "max_principal_angle_degrees": max_angle,
                }
            )
    pairwise = pd.DataFrame(pairwise_rows, columns=pairwise_columns)
    pairwise.to_csv(tables_dir / "threshold_pairwise_stability.tsv", sep="\t", index=False)

    failures = summary[summary["status"] != "complete"].copy()
    failures.to_csv(tables_dir / "threshold_failures.tsv", sep="\t", index=False)
    selected_cutoff = select_cutoff(
        summary,
        pairwise,
        tables_dir,
        min_cluster_frac=min_cluster_frac,
        min_gmm_stability_ari=min_gmm_stability_ari,
    )
    plot_summary(summary, feature_long, pairwise, plots_dir)
    return selected_cutoff


def fail_reasons(row: pd.Series) -> str:
    checks = (
        ("passes_run", "run_failed"),
        ("passes_pc_selection", "no_retained_pcs"),
        ("passes_gmm_stability", "gmm_stability"),
        ("passes_min_cluster", "min_cluster_fraction"),
    )
    return ",".join(reason for column, reason in checks if not bool(row[column]))


def select_cutoff(
    summary: pd.DataFrame,
    pairwise: pd.DataFrame,
    tables_dir: Path,
    *,
    min_cluster_frac: float,
    min_gmm_stability_ari: float,
) -> float:
    """Choose the most sample-rich cutoff yielding the largest supported feature set."""
    decision = summary.copy()
    numeric_columns = (
        "n_core_features",
        "n_rows_pca",
        "n_pcs_kept",
        "selected_stability_median_ari",
        "selected_min_cluster_frac",
    )
    for column in numeric_columns:
        if column not in decision:
            decision[column] = np.nan
        decision[column] = pd.to_numeric(decision[column], errors="coerce")

    maximum_rows = float(decision["n_rows_pca"].max())
    decision["relative_rows_retained"] = (
        decision["n_rows_pca"] / maximum_rows if maximum_rows > 0 else np.nan
    )
    decision["lower_successful_cutoff"] = np.nan
    decision["ari_to_lower_successful_cutoff"] = np.nan
    completed = sorted(
        decision.loc[decision["status"] == "complete", "missingness_cutoff"].astype(float)
    )
    pair_lookup = {}
    if not pairwise.empty:
        pair_lookup = {
            (float(row.cutoff_left), float(row.cutoff_right)): float(
                row.adjusted_rand_index
            )
            for row in pairwise.itertuples()
        }
    for index, cutoff in enumerate(completed):
        mask = decision["missingness_cutoff"].eq(cutoff)
        if index == 0:
            continue
        lower = completed[index - 1]
        decision.loc[mask, "lower_successful_cutoff"] = lower
        decision.loc[mask, "ari_to_lower_successful_cutoff"] = pair_lookup.get(
            (lower, cutoff), np.nan
        )
    decision["passes_run"] = decision["status"].eq("complete")
    decision["passes_pc_selection"] = decision["n_pcs_kept"].ge(1)
    decision["passes_gmm_stability"] = decision[
        "selected_stability_median_ari"
    ].ge(min_gmm_stability_ari)
    decision["passes_min_cluster"] = decision["selected_min_cluster_frac"].ge(
        min_cluster_frac
    )
    gate_columns = [
        "passes_run",
        "passes_pc_selection",
        "passes_gmm_stability",
        "passes_min_cluster",
    ]
    decision["feasible"] = decision[gate_columns].all(axis=1)
    decision["FAIL_REASONS"] = decision.apply(fail_reasons, axis=1)
    decision["SELECTED"] = False

    feasible = decision[decision["feasible"]].copy()
    if feasible.empty:
        decision.to_csv(
            tables_dir / "missingness_cutoff_selection_decision.tsv",
            sep="\t",
            index=False,
        )
        raise RuntimeError(
            "No missingness cutoff passed the PC, GMM-stability, and "
            "cluster-size gates."
        )

    maximum_features = int(feasible["n_core_features"].max())
    eligible = feasible[feasible["n_core_features"].eq(maximum_features)]
    maximum_eligible_rows = int(eligible["n_rows_pca"].max())
    eligible = eligible[eligible["n_rows_pca"].eq(maximum_eligible_rows)]
    selected = float(eligible["missingness_cutoff"].min())
    decision.loc[decision["missingness_cutoff"].eq(selected), "SELECTED"] = True
    decision.to_csv(
        tables_dir / "missingness_cutoff_selection_decision.tsv",
        sep="\t",
        index=False,
    )
    (tables_dir.parent / "SELECTED_MISSINGNESS_CUTOFF.txt").write_text(
        f"{selected:.6g}\n"
    )
    rationale = {
        "selected_missingness_cutoff": selected,
        "selection_rule": (
            "Among cutoffs passing all feasibility gates, maximize the number "
            "of core features, then maximize retained PCA samples, and finally "
            "choose the smallest cutoff."
        ),
        "gates": {
            "minimum_selected_model_median_oob_ari": min_gmm_stability_ari,
            "minimum_cluster_fraction": min_cluster_frac,
            "minimum_retained_pcs": 1,
        },
        "maximum_supported_core_features": maximum_features,
        "maximum_pca_samples_among_maximum_feature_candidates": maximum_eligible_rows,
        "sample_retention_use": "reported_and_used_only_as_a_tiebreaker",
        "cross_cutoff_ari_use": "review_only_not_used_for_feasibility_or_selection",
    }
    (tables_dir / "missingness_cutoff_selection_rationale.json").write_text(
        json.dumps(rationale, indent=2) + "\n"
    )
    return selected


def plot_summary(
    summary: pd.DataFrame,
    feature_long: pd.DataFrame,
    pairwise: pd.DataFrame,
    plots_dir: Path,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    valid = summary[summary["status"] == "complete"].copy()

    axes[0, 0].plot(
        summary["missingness_cutoff"] * 100,
        summary.get("n_core_features", pd.Series(index=summary.index, dtype=float)),
        marker="o",
        color="#0072B2",
    )
    axes[0, 0].set(xlabel="Allowed not-measured values (%)", ylabel="Core features")
    axes[0, 0].grid(alpha=0.25)

    axes[0, 1].plot(
        summary["missingness_cutoff"] * 100,
        summary.get("n_rows_pca", pd.Series(index=summary.index, dtype=float)),
        marker="o",
        color="#009E73",
    )
    axes[0, 1].set(xlabel="Allowed not-measured values (%)", ylabel="Rows retained for PCA")
    axes[0, 1].grid(alpha=0.25)

    if not valid.empty and "selected_k" in valid:
        axes[1, 0].plot(
            valid["missingness_cutoff"] * 100,
            valid["selected_k"],
            marker="o",
            color="#D55E00",
            label="Selected K",
        )
        ari = valid.get("selected_stability_median_ari")
        if ari is not None:
            secondary = axes[1, 0].twinx()
            secondary.plot(
                valid["missingness_cutoff"] * 100,
                ari,
                marker="s",
                color="#CC79A7",
                label="Median OOB ARI",
            )
            secondary.axhline(0.70, color="#CC79A7", linestyle="--", alpha=0.5)
            secondary.set_ylabel("Selected-model median OOB ARI")
    axes[1, 0].set(xlabel="Allowed not-measured values (%)", ylabel="Selected K")
    axes[1, 0].grid(alpha=0.25)

    if not pairwise.empty:
        adjacent = pairwise[pairwise["adjacent_cutoffs"]].copy()
        x = (adjacent["cutoff_left"] + adjacent["cutoff_right"]) * 50
        axes[1, 1].plot(
            x,
            adjacent["adjusted_rand_index"],
            marker="o",
            color="#56B4E9",
        )
    axes[1, 1].set(
        xlabel="Midpoint between adjacent cutoffs (%)",
        ylabel="Assignment ARI on shared samples (review only)",
        ylim=(-0.05, 1.05),
    )
    axes[1, 1].grid(alpha=0.25)
    fig.suptitle("Missingness-threshold sensitivity")
    fig.tight_layout()
    save_figure_all_formats(fig, plots_dir / "missingness_threshold_summary.png", dpi=200)
    plt.close(fig)

    if feature_long.empty:
        return
    matrix = feature_long.pivot(
        index="feature", columns="missingness_cutoff", values="status"
    )
    numeric = matrix.apply(
        lambda column: column.map({"sparse": 0.0, "core": 1.0})
    ).astype(float)
    fig_height = max(4.5, 0.38 * len(numeric) + 1.8)
    fig, ax = plt.subplots(figsize=(9, fig_height))
    image = ax.imshow(numeric.to_numpy(), aspect="auto", cmap="YlGnBu", vmin=0, vmax=1)
    ax.set_yticks(np.arange(len(numeric)), labels=numeric.index)
    ax.set_xticks(
        np.arange(len(numeric.columns)),
        labels=[f"{value * 100:g}%" for value in numeric.columns],
    )
    ax.set_xlabel("Allowed not-measured values")
    ax.set_title("Feature eligibility across missingness thresholds")
    colorbar = fig.colorbar(image, ax=ax, ticks=[0, 1], fraction=0.03, pad=0.02)
    colorbar.ax.set_yticklabels(["Sparse", "Core"])
    fig.tight_layout()
    save_figure_all_formats(fig, plots_dir / "feature_eligibility_heatmap.png", dpi=200)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--feature-cols", required=True)
    parser.add_argument("--thresholds", default="0.10,0.15,0.20,0.25,0.30,0.35,0.40")
    parser.add_argument("--eigenvectors-script", required=True, type=Path)
    parser.add_argument("--selectk-script", required=True, type=Path)
    parser.add_argument("--gmm-script", required=True, type=Path)
    parser.add_argument("--id-col", default="cruise_year_month_depth")
    parser.add_argument("--pc-parallel-replicates", type=int, default=100)
    parser.add_argument("--pc-stability-replicates", type=int, default=100)
    parser.add_argument("--gmm-stability-replicates", type=int, default=100)
    parser.add_argument("--min-cluster-frac", type=float, default=0.02)
    parser.add_argument("--min-gmm-stability-ari", type=float, default=0.70)
    parser.add_argument(
        "--selection-covariance-type",
        choices=["full", "tied", "diag", "spherical"],
        default=DEFAULT_SELECTION_COVARIANCE_TYPE,
    )
    parser.add_argument(
        "--final-covariance-type",
        choices=["full", "tied", "diag", "spherical"],
        default=DEFAULT_FINAL_COVARIANCE_TYPE,
    )
    args = parser.parse_args()

    thresholds = parse_thresholds(args.thresholds)
    labels = [threshold_label(value) for value in thresholds]
    if len(labels) != len(set(labels)):
        raise ValueError(
            "Threshold labels collide at three-decimal precision; use more widely "
            "separated candidate values."
        )
    for value in (
        args.pc_parallel_replicates,
        args.pc_stability_replicates,
        args.gmm_stability_replicates,
    ):
        if value < 1:
            raise ValueError("Sensitivity replicate counts must be positive.")
    if not 0 <= args.min_gmm_stability_ari <= 1:
        raise ValueError("--min-gmm-stability-ari must be between 0 and 1.")
    if not 0 < args.min_cluster_frac < 1:
        raise ValueError("--min-cluster-frac must be strictly between 0 and 1.")
    required_inputs = (
        args.input,
        args.eigenvectors_script,
        args.selectk_script,
        args.gmm_script,
    )
    missing_inputs = [str(path) for path in required_inputs if not path.is_file()]
    if missing_inputs:
        raise FileNotFoundError(f"Missing required input/script files: {missing_inputs}")

    args.outdir.mkdir(parents=True, exist_ok=True)
    config = {
        "input": str(args.input),
        "feature_cols": args.feature_cols,
        "thresholds": thresholds,
        "measurement_semantics": (
            "null_blank_non_numeric_nonfinite_or_negative_is_not_measured;"
            "zero_is_measured_non_detect;positive_finite_is_measured"
        ),
        "pc_parallel_replicates": args.pc_parallel_replicates,
        "pc_stability_replicates": args.pc_stability_replicates,
        "gmm_stability_replicates": args.gmm_stability_replicates,
        "min_cluster_frac": args.min_cluster_frac,
        "min_gmm_stability_ari": args.min_gmm_stability_ari,
        "selection_covariance_type": args.selection_covariance_type,
        "final_covariance_type": args.final_covariance_type,
        "selection_note": (
            "The selected cutoff is passed to the production PCA when the "
            "Nextflow sensitivity stage is enabled."
        ),
    }
    (args.outdir / "run_config.json").write_text(json.dumps(config, indent=2) + "\n")

    statuses: dict[float, dict] = {}
    for threshold in thresholds:
        run_dir = args.outdir / "runs" / threshold_label(threshold)
        if run_dir.is_dir():
            shutil.rmtree(run_dir)
        pca_dir = run_dir / "pca"
        select_dir = run_dir / "selectk"
        gmm_dir = run_dir / "gmm"
        logs_dir = run_dir / "logs"

        pca_command = [
            sys.executable,
            str(args.eigenvectors_script),
            "--input",
            str(args.input),
            "--outdir",
            str(pca_dir),
            "--feature-cols",
            args.feature_cols,
            "--dropna-col-thresh",
            str(threshold),
            "--pc-selection",
            "--anchor-depths",
            "--pcsel-parallel-B",
            str(args.pc_parallel_replicates),
            "--pcsel-stability-R",
            str(args.pc_stability_replicates),
        ]
        pca_code = run_command(
            pca_command, logs_dir / "pca.stdout.log", logs_dir / "pca.stderr.log"
        )
        pca_outputs = (
            pca_dir / "tables" / "eigenvectors_scores.csv",
            pca_dir / "tables" / "pc_keep_decision.csv",
            pca_dir / "tables" / "core_vs_sparse_features.csv",
            pca_dir / "tables" / "matrix_cleaned_with_sparse.csv",
            pca_dir / "qc_summary.json",
        )
        if pca_code == 0 and not required_outputs_exist(pca_outputs):
            pca_code = 90
        state = {
            "status": "pca_failed" if pca_code else "pca_complete",
            "failed_stage": "pca" if pca_code else "",
            "pca_exit_code": pca_code,
        }
        statuses[threshold] = state
        if pca_code:
            continue

        select_command = [
            sys.executable,
            str(args.selectk_script),
            "--eigenvectors",
            str(pca_dir / "tables" / "eigenvectors_scores.csv"),
            "--pc-keep",
            str(pca_dir / "tables" / "pc_keep_decision.csv"),
            "--outdir",
            str(select_dir),
            "--sep",
            ",",
            "--stability-block-col",
            "Cruise",
            "--min-cluster-frac",
            str(args.min_cluster_frac),
            "--stability-R",
            str(args.gmm_stability_replicates),
            "--stability-min-ari",
            str(args.min_gmm_stability_ari),
            "--covariance-type",
            args.selection_covariance_type,
        ]
        select_code = run_command(
            select_command,
            logs_dir / "selectk.stdout.log",
            logs_dir / "selectk.stderr.log",
        )
        select_outputs = (
            select_dir / "SELECTED_K.txt",
            select_dir / "tables" / "gmm_k_selection_decision.csv",
        )
        if select_code == 0 and not required_outputs_exist(select_outputs):
            select_code = 91
        state.update(
            {
                "status": "selectk_failed" if select_code else "selectk_complete",
                "failed_stage": "selectk" if select_code else "",
                "selectk_exit_code": select_code,
            }
        )
        if select_code:
            continue

        selected_k = (select_dir / "SELECTED_K.txt").read_text().strip()
        if not selected_k.isdigit():
            state.update({"status": "selectk_failed", "failed_stage": "selected_k_parse"})
            continue
        gmm_command = [
            sys.executable,
            str(args.gmm_script),
            "--eigenvectors",
            str(pca_dir / "tables" / "eigenvectors_scores.csv"),
            "--pc-keep",
            str(pca_dir / "tables" / "pc_keep_decision.csv"),
            "--outdir",
            str(gmm_dir),
            "--sep",
            ",",
            "--pc-use-mode",
            "keep",
            "--standardize-pc-space",
            "--covariance-type",
            args.final_covariance_type,
            "--episodic-smoothing",
            "--random-state",
            "42",
            "--K",
            selected_k,
            "--matrix-cleaned",
            str(pca_dir / "tables" / "matrix_cleaned_with_sparse.csv"),
        ]
        gmm_code = run_command(
            gmm_command, logs_dir / "gmm.stdout.log", logs_dir / "gmm.stderr.log"
        )
        gmm_output = gmm_dir / "tables" / "compartments_assignments_smoothed.csv"
        if gmm_code == 0 and not required_outputs_exist((gmm_output,)):
            gmm_code = 92
        state.update(
            {
                "status": "gmm_failed" if gmm_code else "complete",
                "failed_stage": "gmm" if gmm_code else "",
                "gmm_exit_code": gmm_code,
            }
        )

    selected = summarize_runs(
        args.outdir,
        thresholds,
        statuses,
        args.id_col,
        min_cluster_frac=args.min_cluster_frac,
        min_gmm_stability_ari=args.min_gmm_stability_ari,
    )
    print(f"[OK] Missingness sensitivity outputs: {args.outdir}")
    print(f"     Selected missingness cutoff = {selected:g}")


if __name__ == "__main__":
    main()
