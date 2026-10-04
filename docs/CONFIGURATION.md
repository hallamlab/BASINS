# BASINS Configuration Reference

BASINS is the **Biochemical Analysis Suite for Investigating Niche Spaces**.
Run configurations with `./run_basins_pipeline.sh <config.yml>`. Existing
`basin_pipeline_nextflow.yml` settings, `BASIN_*` environment variables, and
`.basin` runtime paths remain supported; no migration of saved runs is required.

This reference covers every public key in `basin_pipeline_nextflow.yml`. Copy
that template for each run; do not edit the repository template in place. Empty
YAML values mean “not set.” Unknown keys are not a supported extension
interface and may be ignored.

## From raw tables to a first run

1. Normalize both CSV inputs as described by the README input contract.
2. Copy `basin_pipeline_nextflow.yml` to a run-specific YAML.
3. Set all `paths` values and `biochem.table_a`/`table_b`.
4. Adapt `clean_rename_map`, `clean_keep_cols`, and `feature_cols` together.
5. Run through `BIOCHEM_EIGENVECTORS` on representative profiles.
6. Inspect merged, density, cleaned, and PCA matrix tables.
7. Correct columns, units, and missing values before clustering the full data.

## `paths`

| Key | Type/default | Meaning |
|---|---|---|
| `output_dir` | path; required | Top-level public output containing modules, summary, report, and logs. |
| `runtime_dir` | path; `<output_dir>/.basin` | Persistent work, cache, and publication-staging root. |
| `keep_runtime_dir` | boolean, `true` | Retain runtime state for resume. Set false only when successful-run caches should be discarded. |
| `work_dir` | path; `<runtime_dir>/nf_work` | Optional Nextflow task-cache override. Retain it for resume and use ample storage. |
| `conda_cache_dir` | path; `<runtime_dir>/conda_cache` | Optional package/environment-cache override. Avoid unsafe concurrent writers. |

## `resources`

| Key | Type/default | Meaning |
|---|---|---|
| `threads` | positive integer; host CPU count if absent | CPU allocation passed to stages supporting parallel work. It does not limit disk use. |
| `math_threads` | positive integer; `4` if absent | Workflow-enforced per-task cap for native BLAS, OpenMP, NumExpr, Numba, BLIS, and Accelerate thread pools. |
| `max_concurrent_tasks` | positive integer; `floor(threads / math_threads)`, minimum `1` | Workflow-enforced total Nextflow executor queue cap, preventing bursts of simultaneously launched tasks. |

## `biochem`

| Key | Type/default | Meaning |
|---|---|---|
| `enabled` | boolean, `true` | Enables the BASINS branch. `false` produces no scientific analysis. |
| `table_a` | CSV path; required | Geochemistry observations with compatible profile keys and numeric depth. |
| `table_b` | CSV path; required | CTD/profile observations used for nearest-depth matching, oxygen, and density. |
| `ctd_max_depth_diff_m` | nonnegative number; `10.0` | Maximum permitted absolute depth difference for geochemistry-to-CTD matching. More distant CTD fields remain missing. |
| `output_root` | relative path or path below `output_dir`; empty | Optional staging namespace. Leave empty for the canonical ASPIRE-style layout. External paths are rejected so publication remains atomic. |
| `cleaned_density_filename` | filename | Name of the cleaned density table; this is not a directory. |
| `clean_keep_cols` | list of exact column names | Ordered allow-list retained by custom cleaning. Needed identifiers and features must survive it. |
| `clean_drop_cols` | list, empty | Columns removed after selection. |
| `clean_rename_map` | map `source: canonical` | Renames source measurements. Targets must agree with `feature_cols`; avoid duplicate targets. |
| `feature_cols` | comma-separated names | Numeric features supplied to PCA/EOF and compartment modeling. Exclude IDs, dates, labels, and text. |
| `missingness_cutoff` | number in `(0,1)`; `0.20` | Fixed maximum not-measured fraction for a core feature when sensitivity selection is disabled. Zero is measured; null-like, non-finite, and negative assay values are not measured. |
| `missingness_sensitivity_enabled` | boolean; `false` | Run threshold-specific PCA, PC selection, select-K, and final GMM analyses, select a supported cutoff, and pass it into the production PCA. |
| `missingness_sensitivity_thresholds` | list; `0.10` through `0.40` by `0.05` | Candidate maximum not-measured fractions. Each run is isolated under the missingness-sensitivity module. |
| `missingness_sensitivity_pc_parallel_replicates` | positive integer; `100` | Parallel-analysis permutations in each cutoff run. |
| `missingness_sensitivity_pc_stability_replicates` | positive integer; `100` | Block-bootstrap PCA stability replicates in each cutoff run. |
| `missingness_sensitivity_gmm_stability_replicates` | positive integer; `100` | Cruise-block GMM stability replicates in each cutoff run. |
| `missingness_sensitivity_min_gmm_ari` | 0–1; `0.70` | Minimum selected-model median out-of-bag ARI for a feasible cutoff. |
| `gmm_k` | `auto` or positive integer | GMM component count. `auto` uses select-K diagnostics; an integer forces the count. |
| `eof_pcs` | comma-separated positive integers | Supported cruise-profile EOF modes used for neutral grouping; the SI configuration uses the independently retained stable modes `1,2`. |
| `eof_k_min`, `eof_k_max` | integers, `2 <= min <= max` | Candidate number of neutral cruise groups. |
| `eof_covariance_type` | `full`, `tied`, `diag`, or `spherical` | GMM covariance model in retained EOF space; SI uses `tied` to avoid unstable component-specific covariance estimates. |
| `eof_stability_min_ari` | 0–1 | Required median out-of-bag bootstrap adjusted Rand index; SI uses `0.50` as a moderate-stability discovery threshold and reports the achieved value. |
| `eof_assignment_prob_threshold` | 0–1 | Maximum posterior membership probability below which a cruise assignment is flagged uncertain. Default `0.80`; the best-fitting group is retained for plotting and sensitivity analyses. |
| `eof_min_cluster_frac` | 0–1 | Required minimum fraction of cruises in every group. |
| `eof_min_cluster_n` | positive integer | Required minimum absolute number of cruises in every group. |
| `eof_baseline_months` | positive integer | Width of the centered rolling-median baseline removed from each feature-at-depth profile before cruise EOF analysis; SI uses `24` months to retain seasonal-to-annual departures. |
| `eof_baseline_min_cruises` | positive integer | Minimum neighboring cruises required to estimate each local baseline value. |
| `pea_bootstrap_iterations` | nonnegative integer; `500` | Profile-bootstrap fits used for 95% uncertainty intervals around the ordered three-class k-means boundaries. Higher values improve interval stability and increase runtime. |
| `pea_random_state` | integer; `42` | Reproducible PEA k-means and bootstrap seed. |
| `physical_regime_k_max` | integer at least 3; `6` | Largest candidate K evaluated by BIC, AIC, ICL, silhouette, entropy, and cluster-size diagnostics. The comparison classification retains three ordered physical regimes. |
| `deep_intrusion_quantile` | number in `(0,1)`; `0.90` | Upper threshold for lower-layer density residuals relative to the centered three-calendar-month median. Consecutive cruises may all be intrusions. |
| `oxygen_low_compartment_max` | number; `90.0` | Existing oxic/dysoxic boundary in µM. An O₂ intrusion is only retained while at least one sampled depth remains at or below this value. |
| `oxygen_intrusion_bottom_n` | positive integer; `3` | Number of deepest matched samples whose median O₂ change defines onset and persistence. |
| `oxygen_intrusion_onset_threshold` | nonnegative number; `4.5` | Required bottom-layer median increase from the preceding cruise, in µM, to start an event. |
| `oxygen_intrusion_persistence_threshold` | nonnegative number; `4.5` | Required bottom-layer median elevation above the pre-event profile, in µM, to remain active. |
| `oxygen_intrusion_end_consecutive` | positive integer; `2` | Consecutive cruises below the persistence threshold required to confirm event termination. |
| `renewal_nitrate_col` | column name; `Nitrate` | Chemistry column used to qualify candidate O₂ events as renewals and track post-renewal persistence. |
| `renewal_nitrate_bottom_n` | positive integer; `3` | Number of deepest sampled depths considered for nitrate qualification. SI normally evaluates 165, 185, and 200 m, excluding 150 m. |
| `renewal_nitrate_min_depths` | integer from 1 through `renewal_nitrate_bottom_n`; `2` | Minimum valid nitrate measurements required among the nitrate qualification depths. SI therefore permits one missing value among 165, 185, and 200 m. |
| `renewal_nitrate_detection_limit` | nonnegative number; `0.0` | Deep median nitrate must be greater than this value to support renewal or post-renewal. With the default, zero remains a measured non-detect. |
| `renewal_bridge_enabled` | boolean; `true` | Infer post-renewal for short insufficient-coverage blocks directly bracketed by nitrate-supported cruises from the same event. Measured non-detects and O₂ anomalies are not bridged. |
| `renewal_bridge_max_cruises` | positive integer; `2` | Maximum consecutive insufficient-coverage cruises eligible for bracket inference. |
| `renewal_bridge_max_days` | positive number; `150` | Maximum elapsed days between the two nitrate-supported flanking cruises. |

The template's `clean_rename_map` entries (`PO4`, `NO2`, `NO3`, `NH4`, `H2S`,
`N2O`, `CH4`, and `density_kg_m3`) are SI source-column examples, not reserved
BASINS parameters. Replace or remove them when the new source tables use
different names.

`biochem_pre_asv` is an internal legacy alias for the complete `biochem` map.
New configurations should use `biochem`.

When missingness sensitivity is enabled, BASINS writes
`SELECTED_MISSINGNESS_CUTOFF.txt` and
`missingness_cutoff_selection_decision.tsv`. A cutoff is feasible only when
the run completes with at least one retained PC and passes the configured GMM
stability and minimum-component-size gates. Among feasible cutoffs, BASINS
maximizes the number of eligible core features, then the retained PCA sample
count, and finally chooses the smallest cutoff. Sample retention is reported
but is not a feasibility gate. The selected value replaces
`missingness_cutoff` for the production PCA and all downstream analyses.
Adjusted Rand indices among every pair of successful cutoffs, including
neighboring cutoffs, are reported as review-only sensitivity diagnostics and
do not affect feasibility or selection.

## `environments`

Each value is a Conda/Mamba YAML path. Supported keys are `biochem`,
`biochem_merge`, `biochem_density`, `biochem_strat_metrics`,
`biochem_custom_clean`, `biochem_eigenvectors`, `biochem_selectk`,
`biochem_missingness_sensitivity`, `biochem_gmm`, `biochem_o2_soft`, `biochem_hybrid`, `biochem_compare`,
`biochem_split_o2_by_gmm`, `biochem_strat_anomaly`,
`biochem_state_transitions`, `biochem_succession`, `biochem_feature_assoc`,
`biochem_eof_pipeline`, `biochem_eof_state_cluster`,
`biochem_eof_mode_plots`, and `biochem_within_gmm`.
The `master_summary` key selects the summary builder environment.

Stage-specific entries fall back to `environments.biochem`. Keep the committed
environments for normal runs. Overrides are for dependency development and
change the reproducibility environment.

## Dependencies and interpretation

Stages execute in the order printed by `--list-stages`. Cleaning and feature
choices affect every PCA, GMM, oxygen/hybrid, anomaly, EOF, and within-GMM
result downstream. `gmm_k: auto` depends on SELECTK output. EOF stages require
valid eigenvectors and every requested `eof_pcs` mode. A successful exit proves
technical execution, not that units, scaling, cluster count, or ecological
interpretation are appropriate.

## `master_summary`

`master_summary.enabled` is a boolean defaulting to `true`. It adds integrated
run overview, key-output accounting, and module inventory tables before output
publication. The wrapper independently creates the checksum manifests, module
output visualization, execution links, and `BASIN_run_report.html`.

## Preflight checklist

- Both input paths exist and are readable CSV files.
- Profile identifiers, dates, depth units, and missing-value conventions agree.
- Density inputs are numeric and coordinates are valid.
- Every `feature_cols` value exists after renaming and is numeric.
- Output, work, and Conda-cache locations have sufficient space.
- A representative test produces sensible merged and cleaned tables before the
  full clustering run.
