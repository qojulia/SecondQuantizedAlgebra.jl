import copy
import json
from pathlib import Path
import tempfile
import unittest

from migrate_history import RENAMES, PREFIX, migrate, migrate_file


def dashboard(names):
    return {"lastUpdate": 123, "repoUrl": "repo", "entries": {"Benchmark Results": [
        {"commit": {"id": "abc"}, "date": 42, "tool": "julia", "benches": [
            {"name": name, "value": i + 0.125, "unit": "ns", "range": "± 1", "extra": "allocs"}
            for i, name in enumerate(names)
        ]}
    ]}}


def group(data, tags=None):
    return ["BenchmarkGroup", {"tags": tags or [], "data": data}]


class HistoryTests(unittest.TestCase):
    def test_every_rename_preserves_measurements_and_metadata(self):
        source = dashboard([*RENAMES, "Simplify/Jaynes-Cummings/H"])
        original = copy.deepcopy(source)
        result = migrate(source)
        self.assertEqual(source, original)
        self.assertEqual(result, migrate(result))
        for before, after in zip(source["entries"]["Benchmark Results"][0]["benches"],
                                 result["entries"]["Benchmark Results"][0]["benches"]):
            expected = copy.deepcopy(before)
            expected["name"] = "/".join(RENAMES.get(before["name"], before["name"].split("/")))
            self.assertEqual(expected, after)
            after["name"] = before["name"]
        self.assertEqual(result, source)

    def test_rejects_unknown_depth_and_collisions(self):
        for names in [["Unknown/case"], ["Too/many/name/parts"], ["Empty//legend"],
                      [next(iter(RENAMES)), "/".join(next(iter(RENAMES.values())))]]:
            with self.assertRaises(ValueError):
                migrate(dashboard(names))

    def test_benchmarktools_metadata_and_leaves(self):
        leaf = ["TrialEstimate", {"time": 12.5, "memory": 8, "allocs": 1}]
        source = [{"Julia": "1.11", "BenchmarkTools": "1.8"}, [group({
            "Formal operator functions": group({"Mixed polynomial/formal sum 32": leaf}, ["tag"]),
            "Simplify": group({"Jaynes-Cummings": group({"H": leaf})})
        }, ["root"])]]
        result = migrate(source)
        topic = result[1][0][1]["data"]["Formal operator functions"]
        self.assertEqual(topic[1]["tags"], ["tag"])
        self.assertEqual(topic[1]["data"]["Arithmetic"][1]["data"]["Mixed polynomial and formal sum 32"], leaf)
        self.assertEqual(result[0], source[0])
        self.assertEqual(result[1][0][1]["tags"], ["root"])
        self.assertEqual(result, migrate(result))
        duplicate = copy.deepcopy(source)
        duplicate[1][0][1]["data"]["Formal operator functions"][1]["data"]["Arithmetic"] = group({"Mixed polynomial and formal sum 32": leaf})
        with self.assertRaises(ValueError):
            migrate(duplicate)

    def test_files_backups_idempotence_and_failure(self):
        for prefix in ["", PREFIX]:
            with tempfile.TemporaryDirectory() as directory:
                source = Path(directory) / "data.js"
                original = prefix + json.dumps(dashboard(list(RENAMES)))
                source.write_text(original)
                self.assertTrue(migrate_file(source, source))
                self.assertEqual(source.with_name("data.js.bak").read_text(), original)
                canonical = source.read_bytes()
                self.assertFalse(migrate_file(source, source))
                self.assertEqual(source.read_bytes(), canonical)
                source.write_text(prefix + json.dumps(dashboard(["Unknown/case"])))
                previous = source.read_bytes()
                with self.assertRaises(ValueError):
                    migrate_file(source, source)
                self.assertEqual(source.read_bytes(), previous)


if __name__ == "__main__":
    unittest.main()
