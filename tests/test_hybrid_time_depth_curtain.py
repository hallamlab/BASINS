from __future__ import annotations

import importlib.util
import sys
import tempfile
import types
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "processes" / "split_o2_by_gmm" / "env_split_o2_by_gmm.py"
SPEC = importlib.util.spec_from_file_location("basin_split_o2_by_gmm", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class HybridTimeDepthCurtainTests(unittest.TestCase):
    def test_current_hybrid_palette_contract(self) -> None:
        data = pd.DataFrame({
            "o2_compartment": ["oxic", "dysoxic", "suboxic", "anoxic"],
            "o2_subcompartment_final": [
                "oxic__gmm3", "dysoxic__gmm1", "suboxic__gmm4", "anoxic__gmm4",
            ],
        })
        cfg = types.SimpleNamespace(sub_label_sep="__", prefix_gmm="gmm")
        observed = MODULE.build_full_subcompartment_palette(data, cfg).set_index("label")["color_hex"].to_dict()
        expected = {
            "oxic__gmm0": "#660000", "oxic__gmm1": "#D40000",
            "oxic__gmm2": "#FF4444", "oxic__gmm3": "#FFB2B2",
            "dysoxic__gmm0": "#003300", "dysoxic__gmm1": "#8CFF8C",
            "suboxic__gmm0": "#26667C", "suboxic__gmm1": "#3897B7",
            "suboxic__gmm2": "#6AB8D1", "suboxic__gmm3": "#A5D4E3",
            "suboxic__gmm4": "#E0F0F5", "anoxic__gmm0": "#330033",
            "anoxic__gmm1": "#890089", "anoxic__gmm2": "#DF00DF",
            "anoxic__gmm3": "#FF36FF", "anoxic__gmm4": "#FF8CFF",
        }
        for label, color in expected.items():
            self.assertEqual(observed[label], color)

    def test_pipeline_uses_the_activated_conda_python(self) -> None:
        pipeline = (ROOT / "basin_pipeline.nf").read_text()
        self.assertIn(
            '"\\${CONDA_PREFIX}/bin/python" "${biochemSplitScriptPath}"',
            pipeline,
        )
        self.assertNotIn('\npython "${', pipeline)
        self.assertIn("--curtain-maximum-depth-m 210.0", pipeline)

    def test_cells_span_zero_to_200_and_retain_sampled_depths(self) -> None:
        data = pd.DataFrame(
            {
                "Cruise": [1, 1, 1, 2, 2, 2],
                "date": ["2020-01-01"] * 3 + ["2020-02-01"] * 3,
                "Depth_anchored": [10, 50, 200, 10, 100, 200],
                "max_prob": [0.9, 0.8, 0.95, 0.85, 0.75, 0.9],
                "o2_subcompartment_final": [
                    "oxic__gmm0", "dysoxic__gmm1", "suboxic__gmm4",
                    "oxic__gmm0", "suboxic__gmm1", "suboxic__gmm4",
                ],
            }
        )
        palette = {
            "oxic__gmm0": "#990000",
            "dysoxic__gmm1": "#00D400",
            "suboxic__gmm1": "#2C89A7",
            "suboxic__gmm4": "#9ED4E5",
        }
        cfg = types.SimpleNamespace(
            cruise_col="Cruise",
            date_col="date",
            max_prob_col="max_prob",
            sub_label_sep="__",
            prefix_gmm="gmm",
            plot_formats=["png"],
            png_dpi=100,
            curtain_hide_other=True,
            curtain_time_subdivisions_per_month=2,
            curtain_time_sigma_months=0.5,
            curtain_depth_step_m=5.0,
            curtain_contour_visual_depth_sigma_m=5.0,
        )
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "plots").mkdir()
            (root / "tables").mkdir()
            cells = MODULE.plot_hybrid_time_depth_curtain(
                data,
                "Depth_anchored",
                palette,
                str(root / "plots"),
                str(root / "tables"),
                cfg,
                maximum_depth=200.0,
            )
            self.assertTrue((root / "plots" / "hybrid_compartment_time_depth_curtain.png").is_file())
            self.assertTrue((root / "tables" / "hybrid_compartment_time_depth_curtain_cells.csv").is_file())
            grid = pd.read_csv(
                root / "tables" / "hybrid_compartment_time_depth_curtain_grid.csv"
            )

        self.assertEqual(len(cells), 6)
        self.assertEqual(cells["cruise_index"].nunique(), 2)
        for _, frame in cells.groupby("cruise_index"):
            self.assertEqual(frame["cell_top_depth_m"].min(), 0.0)
            self.assertEqual(frame["cell_bottom_depth_m"].max(), 200.0)
            self.assertEqual(
                sorted(frame["measurement_depth_m"].tolist()),
                sorted(data.loc[data["Cruise"].eq(frame["Cruise"].iloc[0]), "Depth_anchored"].tolist()),
            )
        self.assertNotIn("other", " ".join(grid["o2_subcompartment_display"].astype(str)))
        self.assertIn("o2_subcompartment_rendered", grid.columns)
        self.assertIn("render_state_code", grid.columns)
        self.assertIn("observed_anchor_enforced", grid.columns)
        self.assertTrue(cells["rendered_state_matches_observation"].all())
        self.assertGreater(int(grid["observed_anchor_enforced"].sum()), 0)

    def test_other_assignments_are_audited_but_not_displayed(self) -> None:
        profile = pd.DataFrame({
            "cruise_index": [0, 0, 0, 1, 1],
            "_date": pd.to_datetime([
                "2020-02-01", "2020-02-01", "2020-02-01",
                "2020-11-01", "2020-11-01",
            ]),
            "_depth": [10.0, 100.0, 200.0, 10.0, 200.0],
            "o2_subcompartment_final": [
                "oxic__gmm0", "oxic__other", "suboxic__gmm4",
                "oxic__gmm0", "suboxic__gmm4",
            ],
        })
        cfg = types.SimpleNamespace(
            sub_label_sep="__",
            prefix_gmm="gmm",
            curtain_hide_other=True,
            curtain_time_subdivisions_per_month=2,
            curtain_time_sigma_months=0.5,
            curtain_depth_step_m=5.0,
            curtain_contour_visual_depth_sigma_m=5.0,
        )
        grid, labels, *_ = MODULE.build_smoothed_categorical_curtain(
            profile, cfg, maximum_depth=200.0
        )

        self.assertEqual(labels, ["oxic__gmm0", "suboxic__gmm4"])
        self.assertFalse(grid["o2_subcompartment_display"].str.contains("other").any())
        self.assertEqual(grid["calendar_month"].nunique(), 12)


if __name__ == "__main__":
    unittest.main()
