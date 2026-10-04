# Results and reports


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

## Reading results remotely

Copy the report and its associated output tree to your computer, or serve the output directory over an SSH tunnel. On the server, from the output directory:

```bash
python -m http.server 8765 --bind 127.0.0.1
```

On your own computer, keep a second terminal open:

```bash
ssh -N -L 8765:127.0.0.1:8765 USER@SERVER
```

Open `http://localhost:8765/summary/report/BASIN_run_report.html`. The server command must run on the host reached by that tunnel. BASINS provides an HTML run report and file inventories; do not assume it includes MetaPathways' relational EDA explorer.
