import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "processes" / "output_layout" / "organize_outputs.py"
SPEC = importlib.util.spec_from_file_location("basin_output_layout", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class BasinOutputLayoutTests(unittest.TestCase):
    def test_publish_organizes_modules_and_report(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            staging = root / "staging"
            output = root / "published"
            scientific = staging
            gmm = scientific / "env_compartments_gmm" / "tables"
            processing = scientific / "biochem_processing"
            master = scientific / "master_summary"
            logs = staging / "logs"
            for directory in (gmm, processing, master, logs):
                directory.mkdir(parents=True, exist_ok=True)
            (gmm / "compartments_assignments_smoothed.csv").write_text("sample,cluster\na,0\n")
            (processing / "02_oxygen_best_available_density_RJM.tsv").write_text("Sample\tOxygen\na\t1\n")
            (master / "basin_run_overview.tsv").write_text("metric\tvalue\nmodules\t2\n")
            (logs / "nextflow_trace.tsv").write_text("task_id\tstatus\n1\tCOMPLETED\n")
            config = root / "run.yml"
            config.write_text("paths:\n  output_dir: published\n")

            MODULE.publish(staging, output, config)

            self.assertTrue((output / "modules" / "gmm_compartments" / "tables" / "compartments_assignments_smoothed.csv").is_file())
            self.assertTrue((output / "summary" / "tables" / "basin_run_overview.tsv").is_file())
            self.assertTrue((output / "summary" / "report" / "BASIN_run_report.html").is_file())
            with (output / "summary" / "tables" / "module_output_manifest.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(any(row["module"] == "gmm_compartments" for row in rows))

if __name__ == "__main__":
    unittest.main()
