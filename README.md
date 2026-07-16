# BASIN

BASIN stands for **Biogeochemical Analysis Suite for Investigating Niches**.

BASIN is a Nextflow DSL2 workflow for the biogeochemistry branch extracted from SPARK and organized in the same project style as ASPIRE. It builds cleaned environmental feature matrices, computes density and stratification summaries, derives PCA/EOF structure, assigns oxygen/GMM/hybrid compartments, and writes downstream compartment diagnostics.

The supported entrypoint is `run_basin_pipeline.sh`. It bootstraps a small controller environment, launches `basin_pipeline.nf`, manages resume behavior, and supports stage-aware reruns with `--rerun-from`.

## Key Files

- `run_basin_pipeline.sh`: main wrapper for routine runs.
- `basin_pipeline.nf`: current Nextflow workflow.
- `basin_pipeline_nextflow.yml`: full config template.
- `examples/si.local.yml`: Salton Sea config mirroring the SPARK run that used `SI_asv_pipeline_nextflow.yml`.
- `processes/`: scripts and conda environment YAMLs used by individual stages.

## Quick Start

Install runtime prerequisites:

```bash
command -v mamba
```

Create a run config from the full template, then edit all paths for your environment:

```bash
cp basin_pipeline_nextflow.yml my_basin_run.yml
```

At minimum, review:

- `paths.output_dir`
- `paths.work_dir`
- `paths.conda_cache_dir`
- `biochem.table_a`
- `biochem.table_b`
- `biochem.output_root`
- any `/abs/path/...` placeholder

Run the pipeline:

```bash
./run_basin_pipeline.sh my_basin_run.yml
```

List valid stage names for targeted reruns:

```bash
./run_basin_pipeline.sh --list-stages
```

Force a rerun from one stage onward:

```bash
./run_basin_pipeline.sh my_basin_run.yml --rerun-from BIOCHEM_EIGENVECTORS
```

Pass extra Nextflow options after `--`:

```bash
./run_basin_pipeline.sh my_basin_run.yml -- -with-report report.html -with-trace trace.tsv
```

Direct Nextflow invocation is supported, but the wrapper is preferred:

```bash
nextflow run basin_pipeline.nf --params-file my_basin_run.yml --pipeline_config my_basin_run.yml
```

## Inputs

BASIN currently expects two source tables:

- `biochem.table_a`: geochemistry table.
- `biochem.table_b`: CTD table.

The first workflow stage merges these by nearest depth within station/date groups and creates `Oxygen_best_available`. Later stages expect the Salton Sea column conventions used by the current SPARK biochemistry branch. The default `biochem.clean_rename_map`, `biochem.clean_keep_cols`, and `biochem.feature_cols` in the config template match the last successful SPARK configuration.

## Workflow Stages

The wrapper's current stage order is:

1. `BIOCHEM_MERGE`
2. `BIOCHEM_DENSITY`
3. `BIOCHEM_STRAT_METRICS`
4. `BIOCHEM_CUSTOM_CLEAN`
5. `BIOCHEM_EIGENVECTORS`
6. `BIOCHEM_SELECTK`
7. `BIOCHEM_GMM`
8. `BIOCHEM_O2_SOFT`
9. `BIOCHEM_HYBRID`
10. `BIOCHEM_COMPARE`
11. `BIOCHEM_SPLIT_O2_BY_GMM`
12. `BIOCHEM_STRAT_ANOMALY`
13. `BIOCHEM_STATE_TRANSITIONS`
14. `BIOCHEM_SUCCESSION_GRAPH`
15. `BIOCHEM_FEATURE_ASSOC`
16. `BIOCHEM_EOF_PIPELINE`
17. `BIOCHEM_EOF_STATE_CLUSTER`
18. `BIOCHEM_EOF_MODE_PLOTS`
19. `BIOCHEM_WITHIN_GMM_HDBSCAN`

This extraction intentionally stops before ASV-facing modules such as `BIOCHEM_NETWORK_OVERLAY`, metadata merges, sample agreement dendrograms, and master summaries.

## Scripts Used By Stage

| Workflow stage | Script used | Purpose |
|---|---|---|
| Pipeline launch | `run_basin_pipeline.sh` | Initializes the controller environment, resolves work/cache directories, and launches `basin_pipeline.nf` with the selected YAML config. |
| Workflow orchestration | `basin_pipeline.nf` | Defines process order, config parsing, inputs/outputs, conda environments, and stage commands. |
| `BIOCHEM_MERGE` | `processes/merge_tables/merge_tables_ctd_nearest_depth.py` | Merges geochemistry and CTD tables by nearest depth and creates best-available oxygen. |
| `BIOCHEM_DENSITY` | `processes/density/env_calc_density.py` | Computes in-situ density and sigma0 from salinity, temperature, depth, latitude, and longitude. |
| `BIOCHEM_STRAT_METRICS` | `processes/stratification_metrics/env_stratification_metrics.py` | Computes stratification metrics, density profiles, N2 profiles, MLD, and PEA summaries. |
| `BIOCHEM_CUSTOM_CLEAN` | `processes/custom_clean/custom_density_cleaner.py` | Keeps, drops, and renames columns to produce the cleaned density table used by PCA. |
| `BIOCHEM_EIGENVECTORS` | `processes/eigenvectors/env_eigenvectors.py` | Builds the cleaned feature matrix and computes PCA/EOF-oriented environmental eigenvectors. |
| `BIOCHEM_SELECTK` | `processes/selectk/env_compartments_selectk.py` | Selects the GMM compartment count using BIC/ICL/stability diagnostics. |
| `BIOCHEM_GMM` | `processes/gmm/env_compartments_gmm.py` | Assigns GMM environmental compartments and smoothed responsibilities. |
| `BIOCHEM_O2_SOFT` | `processes/o2_soft/env_compartments_o2_soft.py` | Assigns soft oxygen compartments with episodic smoothing. |
| `BIOCHEM_HYBRID` | `processes/hybrid/env_hybrid_compartment_builder.py` | Combines GMM and oxygen compartments into hybrid compartment outputs. |
| `BIOCHEM_COMPARE` | `processes/compare_compartments/env_compare_compartments.py` | Compares oxygen and GMM compartments and generates UMAP/diagnostic outputs. |
| `BIOCHEM_SPLIT_O2_BY_GMM` | `processes/split_o2_by_gmm/env_split_o2_by_gmm.py` | Splits oxygen states by GMM structure and writes final merged compartment assignments. |
| `BIOCHEM_STRAT_ANOMALY` | `processes/stratification_anomaly/env_stratification_anomaly_detection.py` | Builds stratification time series and anomaly labels. |
| `BIOCHEM_STATE_TRANSITIONS` | `processes/state_transitions/env_state_transition_analysis.py` | Summarizes state transitions, changepoints, and coupling metrics. |
| `BIOCHEM_SUCCESSION_GRAPH` | `processes/succession_graph/env_succession_graph.py` | Builds compartment succession graphs. |
| `BIOCHEM_FEATURE_ASSOC` | `processes/feature_assoc/env_compartment_feature_assoc.py` | Quantifies feature associations with compartment assignments. |
| `BIOCHEM_EOF_PIPELINE` | `processes/eof_pipeline/env_eof_pipeline.py` | Runs EOF-mode PCA summaries from the cleaned matrix and core loadings. |
| `BIOCHEM_EOF_STATE_CLUSTER` | `processes/eof_state_cluster/eof_state_clustering.py` | Clusters EOF scores into environmental states. |
| `BIOCHEM_EOF_MODE_PLOTS` | `processes/eof_mode_plots/eof_mode_plots.py` | Plots EOF modes and explained variance. |
| `BIOCHEM_WITHIN_GMM_HDBSCAN` | `processes/within_gmm_hdbscan/env_within_gmm_hdbscan.py` | Detects high-confidence within-GMM HDBSCAN subclusters. |

## Outputs

By default, outputs are written under `paths.output_dir`. If `biochem.output_root` is set, all workflow result directories are written there instead.

Important outputs include:

- `biochem_processing/02_oxygen_best_available.tsv`
- `biochem_processing/02_oxygen_best_available_density.tsv`
- `biochem_processing/02_oxygen_best_available_density_RJM.tsv`
- `biochem_processing/stratification_metrics/stratification_summary.tsv`
- `env_pca/tables/eigenvectors_scores.csv`
- `env_pca/tables/matrix_cleaned.csv`
- `env_compartments_selectk/SELECTED_K.txt`
- `env_compartments_gmm/tables/compartments_assignments_smoothed.csv`
- `env_o2_soft_compartments/tables/o2_compartments_assignments_smoothed.csv`
- `env_hybrid_soft_compartments/tables/compartments_assignments_hybrid.csv`
- `env_o2_split_by_gmm/tables/merged_o2_split_by_gmm.csv`
- `env_stratification_index/stratification_timeseries.tsv`

## Lineage

The stage code was extracted from `SPARK/biochem_modeling` and organized into the ASPIRE-style `processes/` structure. SPARK is treated as read-only source material during this migration.
