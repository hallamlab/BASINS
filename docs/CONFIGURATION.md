# BASIN Configuration Reference

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
| `output_dir` | path; required | Top-level public output containing modules, summary, report, logs, and compatibility links. |
| `runtime_dir` | path; `<output_dir>/.basin` | Persistent work, cache, and publication-staging root. |
| `keep_runtime_dir` | boolean, `true` | Retain runtime state for resume. Set false only when successful-run caches should be discarded. |
| `work_dir` | path; `<runtime_dir>/nf_work` | Optional Nextflow task-cache override. Retain it for resume and use ample storage. |
| `conda_cache_dir` | path; `<runtime_dir>/conda_cache` | Optional package/environment-cache override. Avoid unsafe concurrent writers. |

## `resources`

| Key | Type/default | Meaning |
|---|---|---|
| `threads` | positive integer; host CPU count if absent | CPU allocation passed to stages supporting parallel work. It does not limit disk use. |

## `biochem`

| Key | Type/default | Meaning |
|---|---|---|
| `enabled` | boolean, `true` | Enables the BASIN branch. `false` produces no scientific analysis. |
| `table_a` | CSV path; required | Geochemistry observations with compatible profile keys and numeric depth. |
| `table_b` | CSV path; required | CTD/profile observations used for nearest-depth matching, oxygen, and density. |
| `output_root` | relative path or path below `output_dir` | Staging namespace for legacy scientific directories. Use `biochem_pipeline` for compatibility with ASPIRE. External paths are rejected so publication remains atomic. |
| `cleaned_density_filename` | filename | Name of the cleaned density table; this is not a directory. |
| `clean_keep_cols` | list of exact column names | Ordered allow-list retained by custom cleaning. Needed identifiers and features must survive it. |
| `clean_drop_cols` | list, empty | Columns removed after selection. |
| `clean_rename_map` | map `source: canonical` | Renames source measurements. Targets must agree with `feature_cols`; avoid duplicate targets. |
| `feature_cols` | comma-separated names | Numeric features supplied to PCA/EOF and compartment modeling. Exclude IDs, dates, labels, and text. |
| `gmm_k` | `auto` or positive integer | GMM component count. `auto` uses select-K diagnostics; an integer forces the count. |
| `eof_pcs` | comma-separated positive integers | One-based PCA/EOF modes used downstream, for example `1,2,4`; requested modes must exist. |

The template's `clean_rename_map` entries (`PO4`, `NO2`, `NO3`, `NH4`, `H2S`,
`N2O`, `CH4`, and `density_kg_m3`) are SI source-column examples, not reserved
BASIN parameters. Replace or remove them when the new source tables use
different names.

`biochem_pre_asv` is an internal legacy alias for the complete `biochem` map.
New configurations should use `biochem`.

## `environments`

Each value is a Conda/Mamba YAML path. Supported keys are `biochem`,
`biochem_merge`, `biochem_density`, `biochem_strat_metrics`,
`biochem_custom_clean`, `biochem_eigenvectors`, `biochem_selectk`,
`biochem_gmm`, `biochem_o2_soft`, `biochem_hybrid`, `biochem_compare`,
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
