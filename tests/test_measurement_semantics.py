from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import numpy as np
import pandas as pd


MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "processes"
    / "eigenvectors"
    / "env_eigenvectors.py"
)
SPEC = importlib.util.spec_from_file_location("env_eigenvectors", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_assay_measurement_semantics_preserve_zero_and_reject_negative_values():
    raw = pd.DataFrame(
        {
            "feature": ["", "NA", "NaN", "null", None, "0", "0.0", "2", "3.5", "-0.1", "inf"]
        }
    )

    numeric = MODULE.coerce_numeric(raw, ["feature"])
    cleaned, n_negative, n_nonfinite = MODULE.invalid_measurements_to_nan(
        numeric[["feature"]]
    )

    assert cleaned["feature"].isna().tolist() == [
        True,
        True,
        True,
        True,
        True,
        False,
        False,
        False,
        False,
        True,
        True,
    ]
    assert cleaned.loc[5, "feature"] == 0
    assert cleaned.loc[6, "feature"] == 0
    assert n_negative == 1
    assert n_nonfinite == 1


def test_missingness_audit_counts_zero_as_measured_non_detect():
    frame = pd.DataFrame({"feature": [np.nan, 0.0, 0.0, 2.0]})

    audit = MODULE.basic_missingness_stats(frame, ["feature"]).iloc[0]

    assert audit["n_rows"] == 4
    assert audit["n_measured"] == 3
    assert audit["n_missing"] == 1
    assert audit["n_zero_measured"] == 2
    assert audit["frac_measured"] == 0.75
    assert audit["frac_missing"] == 0.25
    assert audit["frac_zero_of_measured"] == 2 / 3
