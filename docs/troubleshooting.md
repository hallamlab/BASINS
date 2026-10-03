# Troubleshooting

| Symptom | Check and next step |
|---|---|
| Mamba or Conda not found | Activate the installation providing both commands before invoking the launcher. |
| Environment creation is slow or fails | Read the solver/download error; first-time setup needs internet and enough cache space. Do not delete a cache another run uses. |
| Inputs still point at `/abs/path` | Edit the copied YAML, including optional runtime/cache overrides. |
| CTD fields are missing | Compare exact profile keys and the maximum depth difference; check coordinates and dates. |
| A feature is missing after cleaning | Make `clean_keep_cols`, `clean_rename_map`, and `feature_cols` agree. |
| Clustering fails or groups are unstable | Inspect sample/profile count, variance, missingness, K diagnostics, and coverage before tuning thresholds. |
| Expected public files are absent mid-run | Scientific results are in publication staging until successful finalization. Inspect Nextflow and controller logs. |
| A rerun appears cached | Use `--rerun-from` for deliberate stage recomputation; preserve work for ordinary resume. |
| Modeling cannot find genomes | Supply your own manifest; the private manifest in the SI example is not included. |

Report issues at [GitHub](https://github.com/hallamlab/BASINS/issues). Include the commit, redacted configuration, failing stage, and the actual error before the final Nextflow summary. Include a small input example when the problem depends on table layout.
