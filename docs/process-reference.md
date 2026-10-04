# Process reference


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
20. `BIOCHEM_WITHIN_GMM_HDBSCAN`
21. `BIOCHEM_CONTINUOUS_SECTIONS` (optional plotting)
22. `MASTER_SUMMARY`

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
| `BIOCHEM_WITHIN_GMM_HDBSCAN` | `processes/within_gmm_hdbscan/env_within_gmm_hdbscan.py` | Detects high-confidence within-GMM HDBSCAN subclusters. |
| `MASTER_SUMMARY` | `processes/master_summary/build_basin_summary.py` | Builds integrated run, key-output, and module inventory tables. |
| Output finalization | `processes/output_layout/organize_outputs.py` | Atomically publishes module tables/plots, summaries, reports, logs, and compatibility links. |


`BIOCHEM_CONTINUOUS_SECTIONS` runs `processes/continuous_sections/continuous_time_depth_sections.py` when configured, creating continuous environmental section visualizations. The list above is the controller rerun order; the workflow graph connects the enabled environmental stages.
