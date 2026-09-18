# BASINS

BASINS stands for **Biochemical Analysis Suite for Investigating Niche Spaces**.

Source repository: https://github.com/hallamlab/BASINS

BASINS is a Nextflow DSL2 workflow for integrating biogeochemical observations and characterizing environmental niche space. It builds cleaned environmental feature matrices, computes density and stratification summaries, derives PCA/EOF structure, assigns oxygen/GMM/hybrid compartments, and writes downstream compartment diagnostics. It can also train gapseq-compatible media from the final core+sparse biochemical matrix, reconstruct fixed genome-scale models with gapseq 2.1.0, and compare those models across supported media.

The supported entrypoint is `run_basins_pipeline.sh`. The existing `run_basin_pipeline.sh` launcher remains available; both use the same workflow, configuration files, and resume/cache paths. It bootstraps a small controller environment, launches `basin_pipeline.nf`, manages resume behavior, and supports stage-aware reruns with `--rerun-from`.

The launcher automatically caps native math and OpenMP thread pools inside every task, preventing libraries such as OpenBLAS and MKL from oversubscribing shared servers. It also limits total concurrent Nextflow tasks. The defaults are four native threads per task and at most `floor(resources.threads / resources.math_threads)` concurrent tasks; both controls are enforced by the workflow without requiring shell setup.

## Key Files

- `run_basins_pipeline.sh`: main wrapper for routine runs.
- `basin_pipeline.nf`: current Nextflow workflow.
- `basin_pipeline_nextflow.yml`: full config template.
- `examples/si.local.yml`: Saanich Inlet example with site-specific paths and settings; adapt these before running elsewhere.
- `examples/genomes_manifest.tsv.example`: labeled multi-genome input template.
- `processes/`: scripts and conda environment YAMLs used by individual stages.

The local genome manifest referenced by the SI example is not distributed. When enabling genome modeling, create a manifest from `examples/genomes_manifest.tsv.example` and set `genome_modeling.genomes_manifest` to its path.

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
./run_basins_pipeline.sh my_basin_run.yml
```

List valid stage names for targeted reruns:

```bash
./run_basins_pipeline.sh --list-stages
```

Force a rerun from one stage onward:

```bash
./run_basins_pipeline.sh my_basin_run.yml --rerun-from BIOCHEM_EIGENVECTORS
```

Pass extra Nextflow options after `--`:

```bash
./run_basins_pipeline.sh my_basin_run.yml -- -with-report report.html -with-trace trace.tsv
```

Direct Nextflow invocation is supported, but the wrapper is preferred:

```bash
NXF_SYNTAX_PARSER=v1 nextflow run basin_pipeline.nf --pipeline_config "$PWD/my_basin_run.yml"
```

## Inputs

BASINS currently expects two source tables:

- `biochem.table_a`: geochemistry table.
- `biochem.table_b`: CTD table.

The first workflow stage merges these by nearest depth within station/date groups and creates `Oxygen_best_available`. Later stages use the Saanich Inlet column conventions represented in the template. Adapt `biochem.clean_rename_map`, `biochem.clean_keep_cols`, and `biochem.feature_cols` together when using other source tables.

### Input contract for a new dataset

`table_a` and `table_b` must be CSV files with one row per observation. Column
names are case-sensitive. Before launching a full run, confirm that the two
tables contain compatible cruise/station, date, and numeric depth fields and
that the CTD table contains the temperature, salinity, oxygen, latitude, and
longitude fields needed by the merge and density stages. Environmental
features selected in `biochem.feature_cols` must either already exist or be
created by `biochem.clean_rename_map`.

The input defaults follow Saanich Inlet column conventions; this is not
a schema-free table importer. For a new study, first make a small test pair
containing several profiles and run through `BIOCHEM_EIGENVECTORS`. Inspect the
merged, density, and cleaned tables before running clustering. BASINS does not
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
study-specific fields, and affected stages. Read it before adapting BASINS to a
dataset that does not already follow the SI conventions.

## Workflow Stages

The wrapper's current stage order is:

1. `BIOCHEM_MERGE`
2. `BIOCHEM_DENSITY`
3. `BIOCHEM_STRAT_METRICS`
4. `BIOCHEM_CUSTOM_CLEAN`
5. `BIOCHEM_MISSINGNESS_SENSITIVITY` (optional cutoff selection)
6. `BIOCHEM_EIGENVECTORS`
7. `BIOCHEM_SELECTK`
8. `BIOCHEM_GMM`
9. `BIOCHEM_O2_SOFT`
10. `BIOCHEM_HYBRID`
11. `BIOCHEM_COMPARE`
12. `BIOCHEM_SPLIT_O2_BY_GMM`
13. `BIOCHEM_STRAT_ANOMALY`
14. `BIOCHEM_STATE_TRANSITIONS`
15. `BIOCHEM_SUCCESSION_GRAPH`
16. `BIOCHEM_FEATURE_ASSOC`
17. `BIOCHEM_EOF_PIPELINE`
18. `BIOCHEM_EOF_STATE_CLUSTER`
19. `BIOCHEM_EOF_MODE_PLOTS`
20. `BIOCHEM_GAPSEQ_MEDIA` (optional)
21. `BIOCHEM_WITHIN_GMM_HDBSCAN`
22. `GAPSEQ_RECONSTRUCT` (optional; scattered by genome)
23. `GAPSEQ_COMPARE_MEDIA` (optional; scattered by genome)
24. `GAPSEQ_COMBINE_RESULTS` (optional)
25. `MASTER_SUMMARY`

This extraction intentionally stops before ASV-facing analyses such as
`BIOCHEM_NETWORK_OVERLAY`, metadata merges, and sample-agreement dendrograms.
Its BASINS-specific master summary inventories environmental results only.

## Scripts Used By Stage

| Workflow stage | Script used | Purpose |
|---|---|---|
| Pipeline launch | `run_basins_pipeline.sh` | Initializes the controller environment, resolves work/cache directories, and launches `basin_pipeline.nf` with the selected YAML config. |
| Workflow orchestration | `basin_pipeline.nf` | Defines process order, config parsing, inputs/outputs, conda environments, and stage commands. |
| `BIOCHEM_MERGE` | `processes/merge_tables/merge_tables_ctd_nearest_depth.py` | Merges geochemistry and CTD tables by nearest depth within the configured 10 m default limit and creates best-available oxygen. |
| `BIOCHEM_DENSITY` | `processes/density/env_calc_density.py` | Computes in-situ density and sigma0 from salinity, temperature, depth, latitude, and longitude. |
| `BIOCHEM_STRAT_METRICS` | `processes/stratification_metrics/env_stratification_metrics.py` | Computes stratification metrics, density profiles, N2 profiles, MLD, and PEA summaries. |
| `BIOCHEM_CUSTOM_CLEAN` | `processes/custom_clean/custom_density_cleaner.py` | Keeps, drops, and renames columns to produce the cleaned density table used by PCA. |
| `BIOCHEM_MISSINGNESS_SENSITIVITY` | `processes/missingness_sensitivity/run_missingness_sensitivity.py` | Repeats PCA, PC selection, select-K, and final GMM fitting across candidate feature-missingness cutoffs, writes the selection decision, and passes the preferred cutoff downstream. |
| `BIOCHEM_EIGENVECTORS` | `processes/eigenvectors/env_eigenvectors.py` | Builds the cleaned feature matrix and computes PCA/EOF-oriented environmental eigenvectors. |
| `BIOCHEM_SELECTK` | `processes/selectk/env_compartments_selectk.py` | Selects the GMM compartment count using BIC/ICL/stability diagnostics. |
| `BIOCHEM_GMM` | `processes/gmm/env_compartments_gmm.py` | Assigns GMM environmental compartments and smoothed responsibilities. |
| `BIOCHEM_O2_SOFT` | `processes/o2_soft/env_compartments_o2_soft.py` | Assigns soft oxygen compartments with episodic smoothing. |
| `BIOCHEM_HYBRID` | `processes/hybrid/env_hybrid_compartment_builder.py` | Combines GMM and oxygen compartments into hybrid compartment outputs. |
| `BIOCHEM_COMPARE` | `processes/compare_compartments/env_compare_compartments.py` | Compares oxygen, GMM, and hybrid compartments; benchmarks redundancy with depth/season and leave-one-cruise-out prediction of PCA-excluded chemistry; and generates UMAP/diagnostic outputs. |
| `BIOCHEM_SPLIT_O2_BY_GMM` | `processes/split_o2_by_gmm/env_split_o2_by_gmm.py` | Splits oxygen states by GMM structure and writes final merged compartment assignments. |
| `BIOCHEM_STRAT_ANOMALY` | `processes/stratification_anomaly/env_stratification_anomaly_detection.py` | Builds stratification time series and anomaly labels. |
| `BIOCHEM_STATE_TRANSITIONS` | `processes/state_transitions/env_state_transition_analysis.py` | Summarizes state transitions, changepoints, and coupling metrics. |
| `BIOCHEM_SUCCESSION_GRAPH` | `processes/succession_graph/env_succession_graph.py` | Builds compartment succession graphs. |
| `BIOCHEM_FEATURE_ASSOC` | `processes/feature_assoc/env_compartment_feature_assoc.py` | Quantifies feature associations with compartment assignments. |
| `BIOCHEM_EOF_PIPELINE` | `processes/eof_pipeline/env_eof_pipeline.py` | Runs EOF-mode PCA summaries from the cleaned matrix and core loadings. |
| `BIOCHEM_EOF_STATE_CLUSTER` | `processes/eof_state_cluster/local_cruise_eof.py`, `eof_state_clustering.py`, `eof_state_interpretation.py`, `cruise_group_season_benchmark.py` | Removes the local multi-year baseline from depth-resolved cruise profiles, selects stable and materially sized neutral groups, interprets them using original oceanographic measurements, and tests their redundancy with season and incremental prediction of PCA-excluded chemistry. |
| `BIOCHEM_EOF_MODE_PLOTS` | `processes/eof_mode_plots/eof_mode_plots.py` | Plots EOF modes and explained variance. |
| `BIOCHEM_GAPSEQ_MEDIA` | `processes/gapseq_media/build_gapseq_media.py` | Layers group-specific chemistry onto one shared built-in or user-supplied basal medium for cruise groups, all legacy O₂/GMM/hybrid compartments, and their cruise-group crossings. |
| `BIOCHEM_WITHIN_GMM_HDBSCAN` | `processes/within_gmm_hdbscan/env_within_gmm_hdbscan.py` | Detects high-confidence within-GMM HDBSCAN subclusters. |
| `GAPSEQ_RECONSTRUCT` | `processes/gapseq_modeling/run_gapseq_reconstruction.sh` | Reconstructs and gap-fills each genome once with gapseq 2.1.0, using threaded DIAMOND by default. |
| `GAPSEQ_COMPARE_MEDIA` | `processes/gapseq_modeling/compare_gapseq_media.py` | Applies every supported BASINS medium to a clean copy of the fixed model and runs FBA, exchange FVA, pairwise comparisons, and counterfactuals. |
| `GAPSEQ_COMBINE_RESULTS` | `processes/gapseq_modeling/combine_gapseq_results.py` | Produces cross-genome growth, exchange, and counterfactual tables. |
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
```

Scientific work is written below `.basin/publication_staging` and is published
atomically only after Nextflow succeeds. Canonical files live below `modules/`.
Every scientific figure is exported with the same basename in PNG, PDF, and
SVG formats so raster review, publication assembly, and vector editing use the
same plotted data. Plot renderers apply the shared ASPIRE-compatible
publication typography: Times New Roman throughout, editable TrueType text in
PDF, and live text in SVG. BASINS's scientific color mappings are retained.

Important outputs include:

- `modules/biochemical_processing/tables/02_oxygen_best_available.tsv`
- `modules/biochemical_processing/tables/02_oxygen_best_available_density.tsv`
- `modules/biochemical_processing/tables/02_oxygen_best_available_density_RJM.tsv`
- `modules/biochemical_processing/tables/stratification_metrics/stratification_summary.tsv`
- `modules/biochemical_processing/tables/stratification_metrics/stratification_pea_clusters.tsv`
- `modules/biochemical_processing/tables/stratification_metrics/stratification_pea_threshold_bootstrap.tsv`
- `modules/biochemical_processing/tables/stratification_metrics/stratification_pea_lower_clusters.tsv`
- `modules/biochemical_processing/tables/stratification_metrics/stratification_pea_lower_threshold_bootstrap.tsv`
- `modules/biochemical_processing/tables/stratification_metrics/stratification_pea_class_validation.tsv`
- `modules/biochemical_processing/tables/stratification_metrics/stratification_deep_intrusion_model.tsv`
- `modules/biochemical_processing/tables/stratification_metrics/stratification_oxygen_intrusion_events.tsv`
- `modules/biochemical_processing/tables/stratification_metrics/stratification_nitrate_renewal_events.tsv`
- `modules/biochemical_processing/tables/stratification_metrics/stratification_oxygen_anomalies.tsv`
- `modules/biochemical_processing/tables/stratification_metrics/stratification_physical_regime_selection.tsv`
- `modules/biochemical_processing/tables/stratification_metrics/stratification_physical_regime_clusters.tsv`
- `modules/biochemical_processing/tables/stratification_metrics/stratification_physical_regime_projection.tsv`
- `modules/biochemical_processing/plots/stratification_metrics/stratification_pea_classification.{pdf,png,svg}`
- `modules/biochemical_processing/plots/stratification_metrics/stratification_deep_intrusion.{pdf,png,svg}`
- `modules/biochemical_processing/plots/stratification_metrics/stratification_oxygen_intrusion_events.{pdf,png,svg}`
- `modules/biochemical_processing/plots/stratification_metrics/stratification_physical_regime_selection.{pdf,png,svg}`
- `modules/biochemical_processing/plots/stratification_metrics/stratification_physical_regime_clusters.{pdf,png,svg}`
- `modules/environmental_pca/tables/eigenvectors_scores.csv`
- `modules/environmental_pca/tables/matrix_cleaned.csv`
- `modules/compartment_selection/tables/SELECTED_K.txt`
- `modules/gmm_compartments/tables/compartments_assignments_smoothed.csv`
- `modules/oxygen_compartments/tables/o2_compartments_assignments_smoothed.csv`
- `modules/hybrid_compartments/tables/compartments_assignments_hybrid.csv`
- `modules/oxygen_gmm_subcompartments/tables/merged_o2_split_by_gmm.csv`
- `modules/stratification/tables/stratification_timeseries.tsv`
- `modules/stratification/tables/stratification_physical_biochem_timeseries.tsv`
- `modules/stratification/plots/stratification_physical_biochem_timeseries.pdf`
- `modules/stratification/plots/stratification_physical_biochem_timeseries.png`
- `modules/stratification/plots/stratification_physical_biochem_timeseries.svg`
- `modules/stratification/plots/stratification_physical_biochem_timeseries_complete_cases.{pdf,png,svg}`
- `modules/stratification/plots/stratification_physical_biochem_timeseries_yearly_month_aligned.{pdf,png,svg}`
- `modules/stratification/plots/stratification_physical_biochem_timeseries_yearly_month_aligned_complete_cases.{pdf,png,svg}`
- `modules/stratification/plots/stratification_physical_biochem_monthly_profile_all_data.{pdf,png,svg}`
- `modules/eof_states/plots/cruise_groups_physical_biochem_monthly_profile_maintext.{pdf,png,svg}`
- `modules/eof_states/tables/cruise_group_season_redundancy.tsv`
- `modules/eof_states/tables/cruise_group_season_sparse_cv_summary.tsv`
- `modules/eof_states/tables/cruise_group_season_sparse_paired_comparisons.tsv`
- `modules/eof_states/plots/cruise_group_season_redundancy.{pdf,png,svg}`
- `modules/eof_states/plots/cruise_group_season_sparse_cv_performance.{pdf,png,svg}`
- `modules/gapseq_media/tables/gapseq_media_manifest.tsv`
- `modules/gapseq_media/tables/gapseq_media_provenance.tsv`
- `modules/gapseq_media/tables/gapseq_media_exclusions.tsv`
- `modules/gapseq_media/tables/gapseq_media_summary.json`
- `modules/gapseq_media/tables/gapseq_media_recipe_audit.tsv`
- `modules/gapseq_media/tables/basal_medium_used.csv`
- `modules/gapseq_media/tables/template_audit.tsv`
- `modules/gapseq_media/tables/media/{depth_baseline,cruise_groups,compartments,cruise_group_x_compartment}/**/*.csv`

The compartment families are named `legacy_o2`, `gmm`, and `hybrid`. Only the
historical oxygen-threshold classification is described as legacy; GMM and
hybrid are BASINS-derived compartment schemes.
- `modules/genome_modeling/tables/reconstructions/<genome_id>/<genome_id>.xml`
- `modules/genome_modeling/tables/comparisons/<genome_id>/`
- `modules/genome_modeling/tables/combined/combined_reconstruction_qc.tsv`
- `modules/genome_modeling/tables/combined/combined_reconstruction_summary.tsv`
- `modules/genome_modeling/tables/combined/combined_growth_comparison.tsv`
- `modules/genome_modeling/tables/combined/combined_genome_performance_summary.tsv`
- `modules/genome_modeling/tables/combined/combined_nutrient_importance_summary.tsv`
- `modules/genome_modeling/tables/combined/combined_counterfactual_taxonomic_summary.tsv`
- `modules/genome_modeling/tables/combined/combined_counterfactual_lineage_summary.tsv`
- `modules/genome_modeling/tables/combined/combined_genome_source_metadata.tsv`
- `modules/genome_modeling/tables/combined/combined_compartment_system_comparison.tsv`
- `modules/genome_modeling/tables/combined/combined_hybrid_incremental_resolution.tsv`

The per-genome performance summary includes one `nutrient_importance_*` column
for each configured counterfactual nutrient. The score is the median
proportional loss of supplemented growth when that nutrient is removed across
media where it was supplied: 0 means no modeled growth effect and 1 means
complete loss of growth. Missing values mean that nutrient was not supplied and
therefore was not tested; they must not be interpreted as zero importance.
`nutrient_importance_summary.tsv` also reports the maximum effect, the fraction
of tested media in which growth was eliminated, and the number of media tested.

The physical-metrics stage partitions whole-profile PEA with ordered
three-class one-dimensional k-means. Midpoints between adjacent cluster centers
define the study-relative mixed, intermediate, and stratified boundaries. It
reports cluster centers, silhouette score, bootstrap boundary intervals, and
physical-metric validation by class. Deep intrusion requires a lower-layer
density residual at or above the 90th percentile relative to a centered
three-calendar-month median baseline. It does not require a positive change, so
consecutive cruises can remain classified as intrusions.
The separate `oxygen_intrusion_class` uses the median O₂ change across the
three deepest matched samples while dysoxic, suboxic, or anoxic water remains
in the profile. The SI configuration starts an event at a 4.5 µM increase and
retains it while bottom O₂ remains at least 4.5 µM above the pre-event profile.
Two consecutive cruises below that persistence threshold confirm termination.
Density support is reported independently; the original density-only
classification is retained unchanged. Cruise-level `oxygen_intrusion_last_active`
and `oxygen_intrusion_end_confirmed` flags distinguish the final active cruise
from the later cruise that supplies the second below-threshold confirmation.

O₂ onsets are candidate oxygenation events rather than final renewal labels.
`renewal_phase` qualifies an onset as `renewal` only when nitrate is valid at
the configured minimum number of nitrate qualification depths and their median
is above the configured detection limit. The nitrate depth window is configured
independently from O₂ and uses the three deepest sampled depths by default
(normally 165, 185, and 200 m in SI). Detected nitrate then maintains
`post-renewal` even after O₂ is no longer elevated; an adequately measured
nitrate non-detect terminates the phase. Missing nitrate produces `unknown`
instead of demonstrating termination. When enabled, short blocks of no more
than two insufficient-coverage cruises are inferred as post-renewal only when
they are directly bracketed within 150 days by nitrate-supported cruises from
the same event. Inferred phases are explicitly flagged; measured non-detects
and O₂ anomalies are never bridged. Unqualified O₂ onsets are retained in
`stratification_oxygen_anomalies.tsv` as `O2-only` or `nitrate-unresolved` and
are excluded from renewal overlays.

`stratification_summary.tsv` also contains `pea_lower_class`, an independent
three-class k-means classification of `pea_lower_J_m3` below the configured
upper/lower layer split. Its midpoint thresholds and bootstrap intervals use
the `pea_lower_*` column prefix.

The separate `physical_regime_class` comparison fits a three-component GMM to
standardized cruise-level log PEA, log maximum N2, upper/lower sigma0 contrast,
relative 0.03-density-threshold MLD, and relative pycnocline depth. Components
are ordered using PEA, N2, density contrast, and inverse relative MLD and named
`mixed`, `intermediate`, and `stratified`. Candidate K diagnostics do not alter
the requested three-class comparison; they disclose whether three components
are well supported. PCA coordinates are generated for display only and are not
used to fit the GMM.

The combined stratification time series overlays total, upper, and lower PEA
with the multivariate biochemical depth-centroid distance. Every split panel
uses the same PEA and centroid-distance axis ranges. Solid lines connect
observed sequences, while dashed lines bridge missing observations or long
sampling gaps. PEA classes and deep-intrusion status remain available in the
source tables and their dedicated plots but are not encoded in this combined
time-series figure. The base figure retains every cruise with any plotted
metric. The `_complete_cases` version retains only cruises having depth-centroid
distance and all three PEA measurements; dashed lines bridge cruises removed
from that complete-case sequence. The two `yearly_month_aligned` variants use
one year per panel and fixed January–December positions so the same months line
up vertically. The all-data monthly-profile variant splits the four metrics
into aligned monochrome panels, shows every available cruise as an unconnected
point, overlays a smoothed monthly median with its interquartile band, and marks
each eligible year's observed maximum and minimum with upward and downward
triangles. As in the legacy profile, extrema require at least seven sampled
months; depth-centroid extrema also require those months to meet the biochemical
coverage threshold.

Summary products include `summary/tables/basin_run_overview.tsv`,
`basin_key_outputs.tsv`,
`basin_module_inventory.tsv`, module/checksum manifests, and the HTML report.
