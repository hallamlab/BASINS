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

## `gapseq_media`

This optional stage converts the final cleaned core+sparse biochemical matrix
into gapseq medium CSVs. It creates cruise-group profiles, profiles for every
legacy O₂/GMM/hybrid compartment, and the complete cruise-group × compartment
grid. Optional anchored-depth profiles provide a conventional physical
baseline for the compartment systems.

| Key | Type/default | Meaning |
|---|---|---|
| `enabled` | boolean; `false` | Run the gapseq media stage. |
| `nutrients_tsv` | path; required when enabled | Authoritative local gapseq compound vocabulary. |
| `template` | optional CSV path; bundled default | Shared basal medium with `compounds,name,maxFlux`. If omitted, BASINS uses its bundled minimal glucose basal medium. |
| `compound_mapping` | path; bundled mapping | Chemistry-column to gapseq compound mapping. |
| `transform_method` | `ordinal` or `binary`; `ordinal` | Convert observations to permissive flux scenarios. |
| `detection_limit` | nonnegative number; `0` | Group estimates at or below this value receive zero uptake; measured zeroes remain valid non-detections during aggregation. |
| `min_effective_n` | nonnegative number; `1` | Soft-membership support required to emit a medium; unsupported compartment combinations are excluded and audited. |
| `depth_baseline_m` | list; empty | Anchored depths for conventional depth-baseline media, such as `[100, 150, 200]`. |

The ordinal method maps absence to `0`, observations through the compound's
study-wide Q25 to `0.1`, values through Q75 to `1`, and values above Q75 to
`10`. Concentrations are never treated as fluxes directly. The same basal
medium is used for every environmental group, and mapped observations add or
override its corresponding uptake bounds. BASINS writes the exact resolved
template to `basal_medium_used.csv` and records its source path and SHA-256
digest in `template_audit.tsv`.
Each emitted medium receives a content-derived recipe identifier.
`gapseq_media_recipe_audit.tsv` identifies nominal media that resolve to the
same compound bounds. Compound-level provenance reports the measured
observation count, summed measured membership weight, and measured fraction of
the medium's total membership support for every mapped chemical. These fields
are audit information and do not preferentially exclude sparse compounds.
Structurally possible but unobserved cruise-group × compartment cells are not
filled with global chemistry: no CSV is created, and the cell is written to
`gapseq_media_exclusions.tsv`.

## `genome_modeling`

This optional module reconstructs each supplied genome once with gapseq 2.1.0,
then simulates the resulting fixed SBML model over every medium in the current
`gapseq_media_manifest.tsv`. It requires `gapseq_media.enabled: true`.

| Key | Type/default | Meaning |
|---|---|---|
| `enabled` | boolean; `false` | Enable reconstruction and fixed-model simulation. |
| `genomes_manifest` | CSV/TSV path; required | Two columns, `label, filepath`, with one genome per row. Labels become stable genome IDs; supported inputs are `.fna`, `.fa`, `.fasta`, `.faa`, and gzip-compressed equivalents. |
| `genome_metadata_tsv` | optional CSV/TSV path | Genome taxonomy table used to annotate counterfactual responses by lineage. Accepts `genome_id`, `Genome_Id`, or ASPIRE MAG identifier columns and ranks from Domain through Species. |
| `abundance_tsv` | optional long-format table | Genome recruitment counts used to position modeled genomes at their observed abundance-weighted BASINS PC coordinates and compare those positions with model-predicted metabolic positions. |
| `abundance_sample_col` | string; `sample` | Sample identifier column in `abundance_tsv`; Saanich identifiers are parsed as cruise and depth and matched exactly to PCA coordinates. |
| `abundance_genome_col` | string; `genome` | Genome identifier column in `abundance_tsv`. |
| `abundance_value_col` | string; `read_count` | Nonnegative recruitment-count column in `abundance_tsv`. |
| `abundance_normalization` | enum; `auto` | `auto`, `input_fragment_fpm`, `provided_fpkm`, `provided_tpm`, `median_ratio`, `raw`, or `relative`. The `provided_*` modes validate and retain an already normalized value column without another library-size correction. `auto` requires `abundance_seqkit_tsv` and selects FPM. |
| `abundance_seqkit_tsv` | optional SeqKit TSV | Paired-end FASTQ summary containing `file` and `num_seqs`; supplies the metagenome input-fragment denominator. |
| `abundance_seqkit_file_col`, `abundance_seqkit_count_col` | strings; `file`, `num_seqs` | SeqKit filename and sequence-count fields. |
| `transcript_abundance_tsv` | optional long-format table | Metatranscriptome recruitment counts used to generate the parallel observed-expression overlay. |
| `transcript_abundance_*` | analogous to abundance keys | Sample/genome/value fields, normalization, and paired-end SeqKit input for metatranscriptome libraries. |
| `expected_compartment_top_fraction` | fraction; `0.25` | Defines the top-growth recipe set used to test whether a genome's abundance-weighted expected hybrid compartment ranks among its best predicted media. Ranking is across unique supported hybrid recipes and uses average ranks for ties. |
| `taxonomy` | `auto`; enum | `auto`, `Bacteria`, or `Archaea`. |
| `aligner` | `diamond`; enum | `blast`, `diamond`, or `mmseqs2`. |
| `mmseqs2_bin` | optional executable path | Required for `aligner: mmseqs2`; kept external because current MMseqs2 Conda builds conflict with gapseq 2.1.0's R/libSBML stack. |
| `cpus_per_genome` | integer; `8` | Threads passed to gapseq with `-K`. |
| `max_parallel_genomes` | integer; derived | Maximum concurrent reconstruction and simulation tasks. |
| `solver` | `glpk` | COBRApy solver. |
| `fraction_of_optimum` | `1.0` | FVA objective fraction. |
| `flux_threshold` | `1e-9` | Activity and growth threshold. |
| `supplement_max_flux` | positive number; `10.0` | Fixed non-limiting uptake bound for genome-specific gapseq-predicted nutrients absent from every BASINS environmental medium. |
| `run_fva` | boolean; `true` | Run exchange-reaction FVA. |
| `run_counterfactuals` | boolean; `true` | Remove supplied diagnostic compounds one at a time. |
| `counterfactual_compounds` | `all` or comma-separated IDs; `all` | Compounds eligible for removal experiments. `all` tests every supplied compound that maps to a model exchange. |

Genome reconstruction normally used gapseq's default minimum biomass target of
0.01. If, and only if, `gapseq doall` terminated with its explicit numerical
`Final model cannot grow (status=...)` error after producing one draft model and
one predicted-medium file, BASINS retried the existing gap-filling problem once
with GLPK and a minimum biomass target of 0.001. The default and retry logs,
effective target, solver, exit statuses, and reconstruction mode were written to
`reconstruction_retry_audit.tsv` and `reconstruction_provenance.tsv`. Input,
database, annotation, and all other reconstruction failures remained fatal.

When `genome_metadata_tsv` and counterfactual testing are enabled, BASINS writes
both compound-level and long-form lineage summaries. A model is counted as
having a greater-than-1% counterfactual response when its median proportional
growth loss across tested media (`importance_score`) exceeds 0.01. Complete
growth-loss counts are reported separately.

When a paired-end SeqKit table is supplied, BASINS validates exactly one R1 and
R2 row per sample, requires their `num_seqs` values to agree, and uses either
mate's count as the input-fragment denominator. Recruitment is expressed as
fragments per million input fragments (FPM); R1 and R2 are not summed. `auto`
stops when the SeqKit table is absent rather than silently changing the
normalization; `median_ratio` remains available only by explicit selection.
Each modeled genome is positioned at the
abundance-weighted barycenter of exact cruise-depth PCA matches, after which
genome coordinates are summarized by the median within species. BASINS retains
the original predicted-niche figure and additionally writes a paired figure
linking each species' model-predicted and observed-abundance positions. Audit
tables report sample matching, genome- and species-level coordinates, and the
PC-space displacement between paired positions. Metatranscriptome inputs use
the same calculation and produce separately named observed-expression tables
and plots.

When metagenome recruitment is available, BASINS also validates the predicted
media directly against abundance-defined hybrid niches. For each genome,
normalized recruitment at every exactly matched cruise-depth sample is
multiplied by the sample's soft hybrid-compartment responsibilities. The
largest accumulated fraction defines the observed expected compartment among
the hybrid compartments having supported media. Supplemented biomass flux is
ranked across unique hybrid recipes, identical recipes are counted once, and
ties receive average ranks. Outputs report the expected recipe's rank and
percentile, whether it falls within `expected_compartment_top_fraction`, and
its growth difference and ratio relative to nonexpected recipes. Parallel
species summaries first collapse clonal genome records. Lineage enrichment is
then tested with species as the replicate unit using an exact Poisson-binomial
upper-tail null and Benjamini-Hochberg correction.

The principal validation outputs are
`genome_expected_compartment_performance.tsv`,
`species_expected_compartment_performance.tsv`,
`expected_compartment_performance_summary.tsv`, and
`lineage_expected_compartment_enrichment.tsv`. The corresponding
`*_hybrid_recipe_growth_ranks.tsv` files retain every unique recipe rank, and
the `*_observed_hybrid_affinity.tsv` files retain the soft abundance-weighted
compartment masses used to define the expected state.

The genome manifest may be comma- or tab-delimited, but the header must be exactly
`label,filepath` in that order. Blank lines and lines beginning with `#` are
ignored. Relative filepaths are resolved from the directory containing the
manifest, not from the BASINS repository or run directory. Labels must be unique
and may contain letters, numbers, periods, underscores, and hyphens.
See `examples/genomes_manifest.tsv.example`.

The `filepath` field may be one explicit file or a glob containing `*`, `?`,
character classes, or braces. If a glob matches one file, BASINS preserves the
row's label unchanged. If it matches multiple files, BASINS assigns
`<label>_<filename-without-FASTA-extension>` to each genome. Use the reserved
label `auto` to derive IDs directly from filenames in either case. BASINS fails
before reconstruction if a glob matches nothing or produces duplicate genome
IDs.

Nextflow builds the pinned gapseq 2.1.0 environment, and `doall` uses the
reference data packaged at its environment-default location. No separate
`update-sequences`, database path, or reconstruction medium is required.
Reconstruction invokes `gapseq doall` without a user medium, causing gapseq to
run its own medium-prediction step before gap filling. The resulting fixed model
is subsequently tested against every BASINS-generated medium in
`gapseq_media_manifest.tsv`.

Each genome × medium comparison reports two biomass scenarios. The `strict_*`
columns use only the BASINS environmental medium and retain biosynthetic
completeness differences among genomes. The `supplemented_*` columns add
audited genome-specific gapseq-predicted nutrients that are absent from every
BASINS medium. `growth_rescued_by_supplements`, `absolute_rescue_effect`,
`supplement_count`, and `supplement_compounds` disclose the extent of rescue.
Use strict growth for absolute cross-genome comparisons and supplemented growth
for within-genome environmental responses. Reaction flux, exchange flux, FVA,
and counterfactual tables describe the supplemented scenario.

When source paths contain a depth segment such as `100m`, the aggregation
stage records that value as genome-source metadata. Compact comparison tables
summarize depth, legacy O₂, GMM, and hybrid media using chemically unique
recipes and within-genome growth responses, both across the complete genome
panel and separately by source depth. A companion table reports whether GMM
subdivision adds distinct predicted growth states within each legacy O₂ class.

Nextflow parallelizes independent genomes while `-K` threads the sequence
searches within a genome. Set `cpus_per_genome × max_parallel_genomes` no
higher than the CPUs available to the run.

## `environments`

Each value is a Conda/Mamba YAML path. Supported keys are `biochem`,
`biochem_merge`, `biochem_density`, `biochem_strat_metrics`,
`biochem_custom_clean`, `biochem_eigenvectors`, `biochem_selectk`,
`biochem_missingness_sensitivity`, `biochem_gmm`, `biochem_o2_soft`, `biochem_hybrid`, `biochem_compare`,
`biochem_split_o2_by_gmm`, `biochem_strat_anomaly`,
`biochem_state_transitions`, `biochem_succession`, `biochem_feature_assoc`,
`biochem_eof_pipeline`, `biochem_eof_state_cluster`,
`biochem_eof_mode_plots`, `biochem_within_gmm`, and
`biochem_gapseq_media`. Genome reconstruction and simulation use the dedicated
`genome_modeling` and `genome_simulation` environments; these do not fall back
to the biochemical-analysis environment.
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
