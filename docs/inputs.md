# Preparing input tables


BASINS currently expects two source tables:

- `biochem.table_a`: geochemistry table.
- `biochem.table_b`: CTD table.

The first workflow stage merges these by nearest depth within station/date groups and creates `Oxygen_best_available`. Later stages use the Saanich Inlet column conventions represented in the template. Adapt `biochem.clean_rename_map`, `biochem.clean_keep_cols`, and `biochem.feature_cols` together when using other source tables.

## Input contract for a new dataset

`table_a` and `table_b` must be CSV files with one row per observation. Column
names are case-sensitive. Before launching a full run, confirm that the two
tables contain compatible cruise/station, date, and numeric depth fields and
that the CTD table contains the temperature, salinity, oxygen, latitude, and
longitude fields needed by the merge and density stages. Environmental
features selected in `biochem.feature_cols` must either already exist or be
created by `biochem.clean_rename_map`.

The input defaults follow Saanich Inlet column conventions; this is not
a schema-free table importer. For a new study, first make a small test pair
containing several profiles and run through `BIOCHEM_EIGENVECTORS`. Inspect the
merged, density, and cleaned tables before running clustering. BASINS does not
convert arbitrary user units, so normalize dates, profile identifiers, depth
units, missing values, oxygen units, and all other measurement units first.

Preparation checklist:

1. Make cruise/station and sampling-date identifiers agree between both files.
2. Make depth numeric and use the same depth unit in both files.
3. Verify fields required for density are numeric and geographically valid.
4. Map source chemistry names to canonical names used in `feature_cols`.
5. Remove identifiers and categorical text from `feature_cols`.
6. Start with representative profiles and inspect intermediate tables.
7. Only then choose `gmm_k`, EOF modes, and the production output path.


## Exact profile keys

The merge script requires `Latitude`, `Longitude`, `Cruise`, `Year`, `Month`, `Day`, and `Depth` in both tables. Matching is within the first six keys; depth chooses the nearest CTD observation. The default maximum difference is 10 m. Equidistant CTD matches are averaged. CTD `Oxygen` is preferred where present, with table-A `O2` as fallback. This is not an oxygen calibration or batch-correction step.

Use consistent coordinate values and dates: nearly identical but unequal profile keys can prevent a match. Review `CTD_Depth_Used` and missing CTD fields in the merged table. Do not replace missing measurements with zero; zero has a distinct measured/non-detect meaning in downstream chemistry handling.
