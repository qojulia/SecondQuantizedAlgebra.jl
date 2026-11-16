"""Publication integrity and hard-deadline checks, without scientific dependencies."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

from comparison_support import CONTRACT
from run import execute, report, validate_dataset


class WorkflowTests(unittest.TestCase):
    def archive(self, directory):
        source = directory / "source"
        source.mkdir()
        (source / "julia_bench.jl").write_text("fixture")
        checksum = hashlib.sha256(b"fixture").hexdigest()
        manifest = dict(run_id="example", contract=CONTRACT, contract_sha256="contract",
                        packages=["sqa", "quantumalgebra"], jobs={"julia": "passed"},
                        keys=["jc_heisenberg"])
        (directory / "run.json").write_text(json.dumps(manifest))
        for name in manifest["packages"]:
            data = dict(meta=dict(run_id="example", contract=CONTRACT,
                                  contract_sha256="contract", benchmark_script_sha256=checksum),
                        results={"jc_heisenberg": {"time_ns": 100}}, capped={})
            (directory / f"{name}.json").write_text(json.dumps(data))
        return manifest

    def test_mixed_run_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            self.archive(directory)
            path = directory / "quantumalgebra.json"
            data = json.loads(path.read_text())
            data["meta"]["run_id"] = "different"
            path.write_text(json.dumps(data))
            with self.assertRaisesRegex(ValueError, "Mixed run"):
                validate_dataset(directory)

    def test_missing_scenario_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            self.archive(directory)
            path = directory / "sqa.json"
            data = json.loads(path.read_text())
            data["results"] = {}
            path.write_text(json.dumps(data))
            with self.assertRaisesRegex(ValueError, "Missing or unexpected"):
                validate_dataset(directory)

    def test_smoke_run_cannot_replace_publication(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            self.archive(directory)
            with self.assertRaisesRegex(ValueError, "complete core suite"):
                report(argparse.Namespace(run=directory, publish=True, julia="unused"))

    def test_hung_job_has_hard_deadline(self):
        with tempfile.TemporaryDirectory() as temporary:
            status = execute([sys.executable, "-c", "import time; time.sleep(10)"],
                             {}, Path(temporary) / "job.log", 0.1)
            self.assertEqual(status, "timed_out")

    def test_unrecorded_package_job_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            manifest = self.archive(directory)
            manifest["jobs"] = {}
            (directory / "run.json").write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError, "failed or timed-out"):
                validate_dataset(directory)

    def test_nonfinite_timing_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            self.archive(directory)
            path = directory / "sqa.json"
            data = json.loads(path.read_text())
            data["results"]["jc_heisenberg"]["time_ns"] = float("inf")
            path.write_text(json.dumps(data))
            with self.assertRaisesRegex(ValueError, "Invalid timing"):
                validate_dataset(directory)

    def test_changed_script_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            self.archive(directory)
            (directory / "source/julia_bench.jl").write_text("changed")
            with self.assertRaisesRegex(ValueError, "Script hash mismatch"):
                validate_dataset(directory)


if __name__ == "__main__":
    unittest.main()
