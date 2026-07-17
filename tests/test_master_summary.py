import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "processes" / "master_summary" / "build_basin_summary.py"


class BasinMasterSummaryTests(unittest.TestCase):
    def test_builds_inventory_key_outputs_and_overview(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source = root / "source"
            output = root / "summary"
            processing = source / "biochem_processing"
            selection = source / "env_compartments_selectk"
            processing.mkdir(parents=True)
            selection.mkdir(parents=True)
            (processing / "02_oxygen_best_available_density_RJM.tsv").write_text(
                "Sample\tOxygen\na\t1\nb\t2\n"
            )
            (selection / "SELECTED_K.txt").write_text("4\n")

            subprocess.run(
                [sys.executable, str(SCRIPT), "--input-root", str(source), "--output-dir", str(output)],
                check=True,
            )

            with (output / "basin_run_overview.tsv").open() as handle:
                overview = {row["metric"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}
            self.assertEqual("4", overview["selected_gmm_k"])
            self.assertEqual("2", overview["modules_with_outputs"])
            with (output / "basin_key_outputs.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            density = next(row for row in rows if row["description"].startswith("Best available"))
            self.assertEqual("True", density["present"])
            self.assertEqual("2", density["data_rows"])


if __name__ == "__main__":
    unittest.main()
