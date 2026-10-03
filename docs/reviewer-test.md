# Reviewer and installation checks

## What is included

The repository provides `tests/fixtures/disabled_pipeline.yml`, focused Python test fixtures, the complete configuration template, and a Saanich Inlet example. It does not contain the study's full environmental tables or a packaged end-to-end reviewer dataset. Do not treat the example's private absolute paths as downloadable inputs.

## Lightweight launcher check

After installing the prerequisites, run this from the checkout:

```bash
./run_basins_pipeline.sh --help
./run_basins_pipeline.sh --list-stages
./run_basins_pipeline.sh tests/fixtures/disabled_pipeline.yml --no-resume
```

The disabled fixture writes to `/tmp/basins_disabled_pipeline_test` and disables biochemical analysis, media, and summary generation. Copy it and change the output path if that location is already in use. This checks controller bootstrap and workflow launch only; it does not validate scientific results or exercise all analysis environments.

## A representative scientific review

Use a pair of chemistry/CTD tables you are allowed to distribute, spanning several cruises and depths. Prepare the input contract and complete configuration as described in [quickstart](quickstart.md). Keep optional modeling disabled initially. Inspect merged matches, missingness, density, matrix dimensions, PCA loadings, clustering diagnostics, and assignment coverage. Keep the configuration, software revision, logs, and source-data provenance with the review.

A distributable scientific reviewer bundle needs data provenance, a stable download, and validated expected output checks. Until that bundle exists, do not describe a disabled run as an end-to-end scientific test.
