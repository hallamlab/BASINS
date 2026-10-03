# Your first BASINS run

BASINS combines water-column chemistry with CTD observations. Prepare two tables before starting: a geochemistry table and a CTD table containing the matching profile identifiers and depths. Read the [input contract](inputs.md); the workflow does not infer arbitrary column names or convert arbitrary units.

## 1. Copy the configuration

From the repository directory:

```bash
cp basin_pipeline_nextflow.yml my_run.yml
```

Use a text editor to set these sections in that copy. Keep the remaining analysis defaults for the first review:

```yaml
paths:
  output_dir: /absolute/path/to/my_basins_results
  keep_runtime_dir: true

resources:
  threads: 12
  math_threads: 4
  max_concurrent_tasks: 3

biochem:
  enabled: true
  table_a: /absolute/path/to/geochemistry.csv
  table_b: /absolute/path/to/ctd.csv
```

This is an excerpt, not a replacement for the full configuration. Delete the template's `runtime_dir`, `work_dir`, and `conda_cache_dir` placeholder entries to use the defaults below your output directory, or replace each with a real path. Leave `biochem.output_root` empty. Adapt `clean_rename_map`, `clean_keep_cols`, and `feature_cols` together to match your source measurements. For the first environmental-only run leave `gapseq_media.enabled` and `genome_modeling.enabled` false.

## 2. Run

```bash
./run_basins_pipeline.sh my_run.yml
```

Keep the controller session alive until the command finishes. Scientific stages report progress through Nextflow. First-run environment creation needs network access; subsequent compatible runs reuse environments.

## 3. Review

Open `summary/report/BASIN_run_report.html` under your output directory. Inspect, in this order: merged CTD matches and oxygen, density/stratification metrics, the cleaned environmental matrix, PCA scores/loadings, model-selection diagnostics, compartment assignments, and final comparisons. Output paths are listed in [outputs](outputs.md).

Several profiles are needed for meaningful clustering and time-series analysis. A tiny input can test formatting without supporting all scientific analyses. A successful exit is not evidence that the dataset supports every compartment or temporal conclusion.

## 4. Repeat or refine

```bash
./run_basins_pipeline.sh my_run.yml
./run_basins_pipeline.sh --list-stages
./run_basins_pipeline.sh my_run.yml --rerun-from BIOCHEM_EIGENVECTORS
```

The first command requests automatic resume. The last forces the named stage and all later stages in controller order to rerun. Read [resume behavior](resources.md) before deleting caches or mixing outputs from different configurations.
