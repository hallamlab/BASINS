# BASINS: Biochemical Analysis Suite for Investigating Niche Spaces

BASINS is a Nextflow workflow for integrating biogeochemical observations, characterizing environmental niche space.

**[Full user guide](https://hallamlab-basins.readthedocs.io/en/latest/index.html)** · [Test test](https://hallamlab-basins.readthedocs.io/en/latest/test.html) · [Issues and feature requests](https://github.com/hallamlab/BASINS/issues)

## Quick start

On Linux, install Mamba and Git, then:

```bash
git clone https://github.com/hallamlab/BASINS.git
cd BASINS
cp basin_pipeline_nextflow.yml my_run.yml
```

Follow the user guide to set your output directory and geochemistry/CTD input files in `my_run.yml`, replacing the example paths. Then run:

```bash
./run_basins_pipeline.sh my_run.yml
```

The launcher and Nextflow create the required environments with Mamba. Open `YOUR_OUTPUT/summary/report/BASIN_run_report.html` after completion.

## Try the bundled dataset

From the checkout, run the 47-cruise Saanich Inlet example without editing a configuration:

```bash
./run_basins_pipeline.sh examples/test/test.yml
python3 scripts/validate_test.py examples/test/output
```

Open `examples/test/output/summary/report/BASIN_run_report.html`. See the [user guide](https://hallamlab-basins.readthedocs.io/) for the test walkthrough and data provenance.

## Workflow

[![BASINS workflow](docs/assets/workflow-brief.svg?v=mp-tools-20261008)](https://hallamlab-basins.readthedocs.io/)

For installation, input formats, test checks, configuration, scientific interpretation, and HPC execution, see the **[user guide](https://hallamlab-basins.readthedocs.io/)**. A small Saanich Inlet test dataset is bundled (Torres-Beltrán et al., 2017; [DOI](https://doi.org/10.1038/sdata.2017.159)).

[Report issues or request features](https://github.com/hallamlab/BASINS/issues). Documentation source is maintained in `docs/`.
