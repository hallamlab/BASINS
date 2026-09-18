import subprocess
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "processes/merge_tables/merge_tables_ctd_nearest_depth.py"


def test_default_depth_limit_keeps_distant_ctd_values_missing(tmp_path):
    shared = {
        "Latitude": [48.5, 48.5],
        "Longitude": [-123.5, -123.5],
        "Cruise": [1, 1],
        "Year": [2020, 2020],
        "Month": [1, 1],
        "Day": [1, 1],
    }
    table_a = pd.DataFrame({
        **shared,
        "Depth": [100, 200],
        "O2": [5.0, 2.0],
    })
    table_b = pd.DataFrame({
        **shared,
        "Depth": [110, 185],
        "Oxygen": [50.0, 15.0],
        "Temperature": [9.0, 8.0],
    })
    a_path = tmp_path / "a.csv"
    b_path = tmp_path / "b.csv"
    outdir = tmp_path / "out"
    table_a.to_csv(a_path, index=False)
    table_b.to_csv(b_path, index=False)

    subprocess.run([
        sys.executable,
        str(SCRIPT),
        "--table-a",
        str(a_path),
        "--table-b",
        str(b_path),
        "--outdir",
        str(outdir),
    ], check=True)

    merged = pd.read_csv(outdir / "01_merged_nearest_depth.tsv", sep="\t")
    oxygen = pd.read_csv(outdir / "02_oxygen_best_available.tsv", sep="\t")

    assert merged.loc[0, "CTD_Depth_Used"] == 110
    assert pd.isna(merged.loc[1, "CTD_Depth_Used"])
    assert pd.isna(merged.loc[1, "Temperature"])
    assert oxygen.loc[0, "Oxygen_best_available"] == 50
    assert oxygen.loc[1, "Oxygen_best_available"] == 2
