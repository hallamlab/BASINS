# BASINS

BASINS stands for **Biochemical Analysis Suite for Investigating Niche Spaces**.

Source repository: https://github.com/hallamlab/BASINS

BASINS is a Nextflow DSL2 workflow for integrating biogeochemical observations and characterizing environmental niche space. It builds cleaned environmental feature matrices, computes density and stratification summaries, derives PCA/EOF structure, assigns oxygen/GMM/hybrid compartments, and writes downstream compartment diagnostics. It can also train gapseq-compatible media from the final core+sparse biochemical matrix, reconstruct fixed genome-scale models with gapseq 2.1.0, and compare those models across supported media.

[User guide](docs/index.md) · [Configuration reference](docs/CONFIGURATION.md) · [Issues and feature requests](https://github.com/hallamlab/BASINS/issues)

## Quick start

On Linux, install Mamba and Git, then:

```bash
git clone https://github.com/hallamlab/BASINS.git
cd BASINS
cp basin_pipeline_nextflow.yml my_run.yml
```

Edit `my_run.yml`: set your output directory, geochemistry CSV, and CTD CSV; remove or replace every `/abs/path/` placeholder. See the [first-run walkthrough](docs/quickstart.md) for the small set of fields to configure first.

```bash
./run_basins_pipeline.sh my_run.yml
```

The launcher creates its controller environment and Nextflow creates the analysis environments with Mamba. First installation needs internet access. No manually exported runtime variables are required.

Open `YOUR_OUTPUT/summary/report/BASIN_run_report.html` after successful completion. Start with the summary and module tables before interpreting the figures.

## Workflow

[![BASINS workflow](docs/assets/workflow-main.svg)](docs/assets/workflow-main.svg)

[Detailed workflow, architecture, and data flow](docs/workflow.md) · [SVG](docs/assets/workflow-main.svg) · [PDF](docs/assets/workflow-main.pdf)

## Reviewer and example runs

The repository includes configuration examples and focused test fixtures. It does **not** bundle an end-to-end environmental reviewer dataset. The Saanich Inlet configuration contains study-specific paths and is not a downloadable demo. Follow the [reviewer guide](docs/reviewer-test.md) to validate installation and prepare a representative input pair.

## Full documentation

The [user guide](docs/index.md) covers [input preparation](docs/inputs.md), [configuration](docs/CONFIGURATION.md), [scientific interpretation](docs/analyses.md), [optional genome modeling](docs/genome-modeling.md), and [resources and resuming runs](docs/resources.md).

Documentation source lives in `docs/`. [Report issues or request features](https://github.com/hallamlab/BASINS/issues).
