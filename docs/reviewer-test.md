# Run the reviewer example

The repository includes a small Saanich Inlet dataset drawn from the observations described by [Torres-Beltrán et al. (2017)](https://doi.org/10.1038/sdata.2017.159): **47 cruises in 2009–2012**, **754 chemistry observations**, and **10,997 matching CTD records**. The two CSVs total about **1.06 MiB**. Complete depth profiles and four seasonal cycles support a meaningful environmental workflow test without the full study dataset.

## Run it

After following [installation](installation.md), run from the BASINS checkout:

```bash
./run_basins_pipeline.sh examples/reviewer/reviewer.yml
python3 scripts/validate_reviewer.py examples/reviewer/output
```

No separate dataset download or configuration editing is needed. Mamba creates the required environments on the first run, which requires internet access. The example uses four total CPU threads and two math threads per task. Repeating the command uses the standard resume mechanism.

Open:

```text
examples/reviewer/output/summary/report/BASIN_run_report.html
```

For remote viewing, see [reading results remotely](outputs.md#reading-results-remotely).

## What to inspect

1. **Merge:** all 754 chemistry rows survive; CTD matches respect the configured depth tolerance.
2. **Physical structure:** inspect density, stratification, and renewal diagnostics across the depth profiles.
3. **Feature space:** review missingness, retained features, PCA scores and loadings.
4. **Compartments:** compare GMM, oxygen and hybrid assignments and their membership strengths.
5. **Time structure:** inspect transitions, succession, EOF cruise groups, within-group structure, and continuous sections.
6. **Report and logs:** verify the run completed and follow links to supporting tables and plots.

The validator checks input checksums, chemistry row preservation, finite PCA coordinates, nonempty principal analysis tables and the report. It does not assert a fixed cluster count or cluster numbering. Statistical labels can change across numerical-library versions; they are not permanent biological identifiers.

The configuration enables the environmental workflow and report, with optional missingness sensitivity and experimental stages disabled. To showcase the complete workflow on this short subset, `eof_stability_min_ari` is 0.40 instead of the template’s 0.50. This relaxed demonstration threshold is not a recommendation for research inference. This reviewer subset tests installation and scientific data flow; its shorter temporal baseline cannot reproduce the full study's estimates or support the same inferential claims. Keep the full dataset for research conclusions.

## Data provenance

The bundle selects complete 2009–2012 cruises with matching CTD profiles from the BASINS authors' cleaned chemistry and CTD tables. Source fields and missing-value markers are retained; there is no synthetic filling or selection based on clustering outcomes. `examples/reviewer/data/provenance.json` records the exact source filenames, hashes, selection rule and bundled file hashes. These cleaned tables are not described as identical copies of the original archive.

The [bundle README](https://github.com/hallamlab/BASINS/blob/docs/user-guide/examples/reviewer/README.md) explains how to regenerate the subset using `scripts/build_reviewer_subset.py`. Research manuscript figures, derived study tables and the full SI analysis outputs are not included.

## References

Torres-Beltrán, M. et al. (2017). **A compendium of geochemical information from the Saanich Inlet water column.** *Scientific Data* **4**, 170159. [DOI: 10.1038/sdata.2017.159](https://doi.org/10.1038/sdata.2017.159); [article and methods](https://www.nature.com/articles/sdata2017159).

Torres-Beltrán, M., Hawley, A. K., Capelle, D. et al. (2018). **Data from: A compendium of geochemical information from the Saanich Inlet water column.** Dryad. [DOI: 10.5061/dryad.nh035](https://doi.org/10.5061/dryad.nh035); [dataset and variable descriptions](https://datadryad.org/dataset/doi:10.5061/dryad.nh035).

Three chemistry cruise/date profiles lack matching CTD profiles in the cleaned source tables and are excluded as whole profiles. Their keys are listed in `provenance.json`; remaining measurements are unchanged.

## Validation and resources

The bundled configuration completed all 21 enabled stages. All 754 chemistry observations were retained and matched to CTD; 616 observations entered the PCA and matching GMM/oxygen/hybrid tables, and 37 cruises entered the EOF tables after the normal filtering. These counts describe the validation run, not required cluster labels.

Allow at least **8 GB RAM**. The validation run used about **12 minutes of summed task time** on the four-thread configuration, excluding environment installation, with about **5.9 GB peak task memory**. Actual wall time depends on your system; the bootstrap cluster-selection stage takes most of the compute time. Generated reports, figures and tables occupied about **158 MiB**, in addition to environments and resume caches.
