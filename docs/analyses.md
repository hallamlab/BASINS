# Interpreting the environmental analyses

## Physical structure and measurement preparation

Nearest-depth merging aligns chemistry with CTD observations. Density and stratification summarize the physical water column, including buoyancy-frequency structure, mixed-layer depth, and potential energy anomaly (PEA). These calculations depend on consistent location, depth, temperature, salinity, and measurement units. The GSW implementation supports seawater calculations; see [citations](citations.md).

Whole-profile and lower-layer PEA classifications are study-relative. Their labels are not universal numerical thresholds for every site. Renewal interpretation combines oxygen onsets, nitrate qualification, persistence, and missing-data rules. Review the event tables and settings before treating an oxygen increase as a renewal event. [Outputs](outputs.md) explains these classifications in detail.

## Environmental matrix, PCA, and compartments

Cleaning selects and renames numeric features. Missingness handling distinguishes measured zero from unmeasured or invalid assay values. An optional missingness sensitivity stage compares candidate cutoffs before downstream PCA and GMM fitting. Inspect matrix coverage and loadings alongside component scores; a low-dimensional plot alone does not establish robust environmental structure.

The workflow selects a GMM compartment count, assigns GMM responsibilities, calculates soft oxygen compartments, combines the two into hybrid assignments, and compares schemes. Oxygen-threshold compartments are the historical/legacy scheme; GMM and hybrid are BASINS-derived schemes. HDBSCAN explores supported structure within GMM groups. Numerical component labels describe this fitted run and should not be assumed identical across independently fitted studies.

## Time, succession, and EOF structure

Stratification anomalies, state transitions, succession graphs, and feature associations describe temporal/environmental relationships. Cruise-level EOF analyses account for a local multiyear baseline and assess group stability, size, and assignment confidence. Seasonal-redundancy and held-out-chemistry checks help determine whether groups add information beyond season or depth. These are observational associations; group membership does not establish a causal mechanism.

BASINS' public workflow focuses on environmental analyses. ASV ecological overlays and manuscript-specific investigations found in local developer files are not additional public workflow stages.
