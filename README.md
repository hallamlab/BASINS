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
- `paths.runtime_dir`
- `paths.keep_runtime_dir`
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

### Input contract for a new dataset

`table_a` and `table_b` must be CSV files with one row per observation. Column
names are case-sensitive. Before launching a full run, confirm that the two
tables contain compatible cruise/station, date, and numeric depth fields and
that the CTD table contains the temperature, salinity, oxygen, latitude, and
longitude fields needed by the merge and density stages. Environmental
features selected in `biochem.feature_cols` must either already exist or be
created by `biochem.clean_rename_map`.

The implementation was extracted with Salton Sea column conventions; it is not
a schema-free table importer. For a new study, first make a small test pair
containing several profiles and run through `BIOCHEM_EIGENVECTORS`. Inspect the
merged, density, and cleaned tables before running clustering. BASIN does not
convert arbitrary user units, so normalize dates, profile identifiers, depth
units, missing values, oxygen units, and all other measurement units first.

Preparation checklist:

1. Make cruise/station and sampling-date identifiers agree between both files.
2. Make depth numeric and use the same depth unit in both files.
3. Verify fields required for density are numeric and geographically valid.
4. Map source chemistry names to canonical names used in `feature_cols`.
5. Remove identifiers and categorical text from `feature_cols`.
6. Start with representative profiles and inspect intermediate tables.
7. Only then choose `gmm_k`, EOF modes, and the production output path.

## Complete Configuration Reference

The field-by-field reference is in
[`docs/CONFIGURATION.md`](docs/CONFIGURATION.md). It documents every key in
`basin_pipeline_nextflow.yml`, including types, defaults, required inputs,
study-specific fields, and affected stages. Read it before adapting BASIN to a
dataset that does not already follow the SI conventions.

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
20. `MASTER_SUMMARY`

This extraction intentionally stops before ASV-facing analyses such as
`BIOCHEM_NETWORK_OVERLAY`, metadata merges, and sample-agreement dendrograms.
Its BASIN-specific master summary inventories environmental results only.

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
| `MASTER_SUMMARY` | `processes/master_summary/build_basin_summary.py` | Builds integrated run, key-output, and module inventory tables. |
| Output finalization | `processes/output_layout/organize_outputs.py` | Atomically publishes module tables/plots, summaries, reports, logs, and compatibility links. |

## Outputs

Successful wrapper runs publish an ASPIRE-style output tree:

```text
<output_dir>/
├── .basin/                 # persistent Nextflow/runtime staging
├── modules/
│   └── <module>/
│       ├── tables/
│       └── plots/
├── intermediates/
├── references/
├── summary/
│   ├── tables/
│   ├── plots/
│   └── report/BASIN_run_report.html
├── logs/
└── biochem_pipeline/       # compatibility links for legacy consumers
```

Scientific work is written below `.basin/publication_staging` and is published
atomically only after Nextflow succeeds. The `biochem_pipeline` compatibility
tree preserves historical paths used by ASPIRE while the canonical organized
files live below `modules/`.

Important outputs include:

- `modules/biochemical_processing/tables/02_oxygen_best_available.tsv`
- `modules/biochemical_processing/tables/02_oxygen_best_available_density.tsv`
- `modules/biochemical_processing/tables/02_oxygen_best_available_density_RJM.tsv`
- `modules/biochemical_processing/tables/stratification_metrics/stratification_summary.tsv`
- `modules/environmental_pca/tables/eigenvectors_scores.csv`
- `modules/environmental_pca/tables/matrix_cleaned.csv`
- `modules/compartment_selection/tables/SELECTED_K.txt`
- `modules/gmm_compartments/tables/compartments_assignments_smoothed.csv`
- `modules/oxygen_compartments/tables/o2_compartments_assignments_smoothed.csv`
- `modules/hybrid_compartments/tables/compartments_assignments_hybrid.csv`
- `modules/oxygen_gmm_subcompartments/tables/merged_o2_split_by_gmm.csv`
- `modules/stratification/tables/stratification_timeseries.tsv`

Equivalent historical paths are available below `biochem_pipeline/` as
compatibility links. Summary products include
`summary/tables/basin_run_overview.tsv`, `basin_key_outputs.tsv`,
`basin_module_inventory.tsv`, module/checksum manifests, and the HTML report.

## Lineage

The stage code was extracted from `SPARK/biochem_modeling` and organized into the ASPIRE-style `processes/` structure. SPARK is treated as read-only source material during this migration.
