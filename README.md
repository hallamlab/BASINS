# BASINS

**Biochemical Analysis Suite for Investigating Niche Spaces**

BASINS is a Nextflow workflow for integrating biogeochemical observations, characterizing environmental niche space.

**[Read the full user guide on Read the Docs](https://hallamlab-basins.readthedocs.io/)**

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

## Workflow

[![BASINS workflow](docs/assets/workflow-brief.svg)](https://hallamlab-basins.readthedocs.io/)

For installation, input formats, reviewer checks, configuration, scientific interpretation, and HPC execution, see the **[user guide](https://hallamlab-basins.readthedocs.io/)**. The repository includes focused test fixtures but does not bundle a complete environmental reviewer dataset.

[Report issues or request features](https://github.com/hallamlab/BASINS/issues). Documentation source is maintained in `docs/`.
