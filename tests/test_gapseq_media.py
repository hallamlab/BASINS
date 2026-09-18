import json
import subprocess
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "processes/gapseq_media/build_gapseq_media.py"
DEFAULT_TEMPLATE = ROOT / "processes/gapseq_media/default_basal_medium.csv"


def write_csv(path, rows):
    pd.DataFrame(rows).to_csv(path, index=False)


def test_builds_complete_group_and_legacy_grid(tmp_path):
    matrix = tmp_path / "matrix.csv"
    write_csv(matrix, {
        "cruise_year_month_depth": ["a", "b", "c", "d"],
        "Cruise": ["C1", "C1", "C2", "C2"],
        "Depth": [100, 200, 100, 200],
        "Oxygen": [100, 50, 1, 0],
        "Nitrate": [20, 10, 1, 0],
    })
    cruise = tmp_path / "cruise.tsv"
    pd.DataFrame({
        "Cruise": ["C1", "C2"], "resp_0": [0.9, 0.1], "resp_1": [0.1, 0.9]
    }).to_csv(cruise, sep="\t", index=False)

    def assignments(path, n, hybrid=False):
        rows = {"cruise_year_month_depth": ["a", "b", "c", "d"]}
        for i in range(n):
            rows[f"resp_{i}"] = [1.0 if j % n == i else 0.0 for j in range(4)]
        if hybrid:
            rows["component"] = [j % n for j in range(4)]
            rows["compartment_name"] = [f"hyb_C{j % 2}_G{j % 2}" for j in range(4)]
        write_csv(path, rows)

    o2, gmm, hybrid = (tmp_path / name for name in ("o2.csv", "gmm.csv", "hybrid.csv"))
    assignments(o2, 4); assignments(gmm, 2); assignments(hybrid, 4, True)
    nutrients = tmp_path / "nutrients.tsv"
    basal = pd.read_csv(DEFAULT_TEMPLATE)
    vocabulary = pd.concat([
        basal.rename(columns={"compounds": "id"})[["id", "name"]],
        pd.DataFrame({"id": ["cpd00007", "cpd00209"], "name": ["O2", "Nitrate"]}),
    ]).drop_duplicates("id")
    vocabulary.to_csv(nutrients, sep="\t", index=False)
    mapping = tmp_path / "mapping.tsv"
    pd.DataFrame({
        "chemistry_column": ["Oxygen", "Nitrate"],
        "compound_id": ["cpd00007", "cpd00209"],
        "compound_name": ["O2", "Nitrate"],
        "measurement_unit": ["uM", "uM"],
    }).to_csv(mapping, sep="\t", index=False)
    out = tmp_path / "out"
    (out / "media").mkdir(parents=True)
    (out / "media" / "stale.csv").write_text("stale\n")
    subprocess.run([
        sys.executable, str(SCRIPT), "--matrix", str(matrix), "--cruise-groups", str(cruise),
        "--o2-assignments", str(o2), "--gmm-assignments", str(gmm),
        "--hybrid-assignments", str(hybrid), "--nutrients", str(nutrients),
        "--mapping", str(mapping), "--outdir", str(out),
        "--min-effective-n", "0.5",
        "--depth-baseline-m", "100,200",
    ], check=True)
    summary = json.loads((out / "tables/gapseq_media_summary.json").read_text())
    # All 10 legacy cells are supported. Within each cruise group, only the
    # compartments represented by its two samples are emitted.
    assert summary["n_media"] == 26
    assert summary["depth_baseline_m"] == [100.0, 200.0]
    assert summary["n_excluded_unsupported_media"] == 8
    manifest = pd.read_csv(out / "tables/gapseq_media_manifest.tsv", sep="\t")
    assert len(manifest) == 26
    assert not (out / "media" / "stale.csv").exists()
    assert set(manifest.loc[manifest["family"].eq("depth"), "compartment"]) == {
        "depth_100m", "depth_200m"
    }
    assert manifest["recipe_id"].str.fullmatch(r"[0-9a-f]{12}").all()
    recipe_audit = pd.read_csv(
        out / "tables/gapseq_media_recipe_audit.tsv", sep="\t"
    )
    assert len(recipe_audit) == len(manifest)
    provenance = pd.read_csv(
        out / "tables/gapseq_media_provenance.tsv", sep="\t"
    )
    assert {
        "measured_observations_n",
        "measured_membership_weight",
        "measured_membership_fraction",
    }.issubset(provenance.columns)
    assert provenance["measured_membership_fraction"].between(0, 1).all()
    exclusions = pd.read_csv(out / "tables/gapseq_media_exclusions.tsv", sep="\t")
    assert len(exclusions) == 8
    assert set(exclusions["reason"]) == {"insufficient_soft_membership_support"}
    assert all((out / path).is_file() for path in manifest["medium_file"])
    assert set(manifest["template_source"]) == {"basin_builtin_default"}
    pd.testing.assert_frame_equal(
        pd.read_csv(out / "tables/basal_medium_used.csv"), basal
    )
    audit = pd.read_csv(out / "tables/template_audit.tsv", sep="\t")
    assert audit.loc[0, "source"] == "basin_builtin_default"
    assert len(audit.loc[0, "sha256"]) == 64
    expected_compounds = set(basal["compounds"]) | {"cpd00007", "cpd00209"}
    for path in manifest["medium_file"]:
        medium = pd.read_csv(out / path)
        assert list(medium) == ["compounds", "name", "maxFlux"]
        assert set(medium["compounds"]) == expected_compounds


def test_rejects_unknown_mapping_compound(tmp_path):
    # Unit-level validation is intentionally strict and happens before grouping.
    sys.path.insert(0, str(SCRIPT.parent))
    from build_gapseq_media import load_mapping

    mapping = tmp_path / "mapping.tsv"
    pd.DataFrame({
        "chemistry_column": ["Oxygen"], "compound_id": ["bad"],
        "compound_name": ["O2"],
    }).to_csv(mapping, sep="\t", index=False)
    nutrients = pd.DataFrame({"id": ["cpd00007"], "name": ["O2"]})
    try:
        load_mapping(mapping, pd.DataFrame({"Oxygen": [1]}), nutrients)
    except ValueError as exc:
        assert "absent from nutrients.tsv" in str(exc)
    else:
        raise AssertionError("unknown compound id was accepted")


def test_custom_shared_template_is_loaded_and_identified(tmp_path):
    sys.path.insert(0, str(SCRIPT.parent))
    from build_gapseq_media import load_template

    template = tmp_path / "custom.csv"
    pd.DataFrame({
        "compounds": ["cpd00001"],
        "name": ["H2O"],
        "maxFlux": [50],
    }).to_csv(template, index=False)
    nutrients = pd.DataFrame({"id": ["cpd00001"], "name": ["H2O"]})

    frame, audit = load_template(template, nutrients, "custom")

    assert frame.loc[0, "maxFlux"] == 50
    assert audit["template"] == "custom"
    assert audit["source"] == "custom"
    assert len(audit["sha256"]) == 64
