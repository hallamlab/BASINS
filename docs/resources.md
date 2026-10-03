# Resources, HPC, and resume

## Local resources

Set `resources.threads` for the run, `resources.math_threads` for native numerical-library pools, and optionally `resources.max_concurrent_tasks`. The launcher defaults the task cap to `max(1, floor(threads / math_threads))`. It exports native thread limits for OpenBLAS, MKL, OpenMP, NumExpr, Numba, BLIS, and Accelerate. These controls reduce oversubscription; they do not cap disk usage.

The environmental workflow has many explicit sequential dependencies. Raising the task cap does not make dependent stages simultaneous. Optional genome reconstruction and comparisons scatter by genome; `genome_modeling.cpus_per_genome` and `max_parallel_genomes` control that branch. Resource requests and solver memory still need to fit the available hardware.

## Resume and targeted reruns

The launcher defaults to `--resume-policy last-with-tasks`, selecting a prior Nextflow run with usable task history rather than blindly choosing an empty launch. Alternatives are `--resume-policy latest`, `--resume-run RUN_NAME`, and `--no-resume`.

```bash
./run_basins_pipeline.sh my_run.yml --resume-run RUN_NAME
./run_basins_pipeline.sh my_run.yml --rerun-from BIOCHEM_GMM
```

`--rerun-from` disables caching for that stage and later stages in the controller's ordered list. Use the wrapper's resume switches rather than passing Nextflow `-resume` after `--`: the wrapper sanitizes those arguments.

Keep `.basin`, work directories, and Nextflow history while you need reuse. `paths.keep_runtime_dir: false` removes the configured runtime directory after successful publication, sacrificing that runtime's cached work. The published report is not a replacement for execution checkpoints. Never manually remove active work or caches. The wrapper guards against simultaneous use of the same Conda-cache location.

## Cluster execution

BASINS uses Nextflow, but its wrapper does not expose MetaPathways' Slurm convenience flags. Cluster operation needs a site-specific Nextflow configuration plus shared input/output/work locations and usable software environments on compute nodes. Pass that configuration through the wrapper:

```bash
./run_basins_pipeline.sh my_run.yml -- -c /absolute/path/to/site.config
```

Configure the executor, account, partition, CPU, memory, walltime, and queue policy with your administrator's guidance. See the [Nextflow executor reference](https://www.nextflow.io/docs/latest/executor.html). Do not assume the local YAML CPU budget alone is a complete Slurm resource configuration. Provision scientific environments before using compute nodes without internet access. Cluster deployment is a site integration step, not a separately validated turnkey profile in this repository.
