# Saanich Inlet test dataset

This small, real-data example retains complete depth profiles from 47 cruises in 2009–2012: 754 chemistry observations at 24 recorded depths (10–200 m), plus 10,997 matching CTD records. The two CSVs total 1,112,748 bytes (about 1.06 MiB). No extra download is required after cloning BASINS.

From the repository root:

```bash
./run_basins_pipeline.sh examples/test/test.yml
python3 scripts/validate_test.py examples/test/output
```

Open `examples/test/output/summary/report/BASIN_run_report.html`. The configuration requests four total CPU threads, with two native math threads per task. The launcher creates its controller environment and Nextflow creates the scientific environment with Mamba; first installation needs internet access. Results and caches stay under `examples/test/output/` and are ignored by Git. Repeating the same command uses the normal resume mechanism.

The example exercises chemistry/CTD matching, physical metrics, matrix preparation, PCA, GMM/oxygen/hybrid compartments, temporal analyses, EOF cruise grouping, within-group structure, continuous sections, and the report. Optional missingness sensitivity and experimental stages are disabled. The cruise-group stability threshold is 0.40 for this demonstration (the normal template uses 0.50); other statistical settings follow the normal template. This is an installation and interpretation exercise, not a reproduction of the full Saanich Inlet study: short time coverage changes baseline estimates, cluster selection and inference. Cluster numbers are not stable biological identifiers.

## Provenance and citation

These are row subsets of the BASINS authors' cleaned SI chemistry and CTD tables. All source columns, field values and missing-value markers are retained; only rows outside the selected cruises are removed and CSV line endings normalized. The source filenames and SHA-256 checksums, selection rule, and distributed file checksums are recorded in `data/provenance.json`. The cleaned source files are not claimed to be identical to the original archival release.

Cite the underlying observations:

Torres-Beltrán, M. et al. (2017). **A compendium of geochemical information from the Saanich Inlet water column.** *Scientific Data* **4**, 170159. DOI: [10.1038/sdata.2017.159](https://doi.org/10.1038/sdata.2017.159). [Article and methods](https://www.nature.com/articles/sdata2017159).

Torres-Beltrán, M., Hawley, A. K., Capelle, D. et al. (2018). **Data from: A compendium of geochemical information from the Saanich Inlet water column.** Dryad dataset. DOI: [10.5061/dryad.nh035](https://doi.org/10.5061/dryad.nh035). [Dataset and variable descriptions](https://datadryad.org/dataset/doi:10.5061/dryad.nh035).

To reproduce the subset from the named cleaned source files:

```bash
python3 scripts/build_test_subset.py \
  --chemistry /path/to/SI_JA_Compiled_Geochem_Dec_09_Outlier_RM.csv \
  --ctd /path/to/SI_JA_Compiled_CTD_Data_Dec_18_2025_Outlier_RM.csv \
  --output /path/to/rebuilt-test-data
```

The builder selects all chemistry rows from 2009–2012 cruises with matching CTD profiles and CTD rows matching their cruise/year/month/day. It does not select on cluster labels, significance, or expected outcomes. The full SI study data and manuscript products are not bundled.

Three chemistry cruise/date profiles lack matching CTD profiles in the cleaned source tables and are excluded as whole profiles. Their keys are listed in `provenance.json`; remaining measurements are unchanged.

## Validation and resources

The bundled configuration completed all 21 enabled stages. All 754 chemistry observations were retained and matched to CTD; 616 observations entered the PCA and matching GMM/oxygen/hybrid tables, and 37 cruises entered the EOF tables after the normal filtering. These counts describe the validation run, not required cluster labels.

Allow at least **8 GB RAM**. The validation run used about **12 minutes of summed task time** on the four-thread configuration, excluding environment installation, with about **5.9 GB peak task memory**. Actual wall time depends on your system; the bootstrap cluster-selection stage takes most of the compute time. Generated reports, figures and tables occupied about **158 MiB**, in addition to environments and resume caches.
