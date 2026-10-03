# Installation

BASINS is a Linux workflow launched from its source repository. Mamba supplies the software environments. You need Git, Mamba with its accompanying Conda executable, write access to the checkout and output locations, and internet access for the initial environment creation.

```bash
git clone https://github.com/hallamlab/BASINS.git
cd BASINS
command -v mamba
command -v conda
```

The supported launcher is `./run_basins_pipeline.sh`; `run_basin_pipeline.sh` is a compatible historical name. On first invocation the launcher builds `.controller_env` from `processes/shared_envs/controller.yml`, including Nextflow, Java, Python, jq, and yq. Nextflow builds each scientific environment from the YAML files under `processes/shared_envs/`. You do not need to install each analysis library by hand.

Proceed to the [quickstart](quickstart.md). Environment creation can take longer than the first scientific stage. The shipped environment specifications target Linux and contain platform-specific pins; do not assume they will solve on macOS or Windows.

## Updating

Keep a record of `git rev-parse HEAD` with your run configuration. Stop active runs before updating their checkout, then use `git pull --ff-only`. Use a new output directory when changing scientific settings or implementation if you need a clean comparison. Retain the existing runtime directory for unchanged-run resume.

The repository workflow is the installation route documented here. A BASINS package or container release is not required to follow this guide.
