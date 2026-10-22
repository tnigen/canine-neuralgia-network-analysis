import gzip
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import openpyxl

import run_proximity


class ScoreCutoffTests(unittest.TestCase):
    def test_none_is_distinct_from_zero(self):
        self.assertIsNone(run_proximity.parse_score_cutoff("none"))
        self.assertEqual(run_proximity.parse_score_cutoff("0"), 0)

    def test_numeric_cutoff(self):
        self.assertEqual(run_proximity.parse_score_cutoff("700"), 700)
        self.assertEqual(run_proximity.parse_score_cutoff("900"), 900)

    def test_public_command_removes_absolute_paths(self):
        command = run_proximity.public_command(
            [
                "/tmp/project/run_proximity.py",
                "--network",
                "/tmp/project/network.txt.gz",
                "--output=/tmp/project/results/run",
            ]
        )
        self.assertEqual(
            command,
            "run_proximity.py --network network.txt.gz --output=run",
        )


class PipelineIntegrationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.network = self.root / "toy.links.txt.gz"
        with gzip.open(self.network, "wt", encoding="utf-8") as handle:
            handle.write("protein1 protein2 combined_score\n")
            edges = [
                ("9615.ENSCAFP00000000001", "9615.ENSCAFP00000000002", 950),
                ("9615.ENSCAFP00000000002", "9615.ENSCAFP00000000001", 950),
                ("9615.ENSCAFP00000000002", "9615.ENSCAFP00000000003", 800),
                ("9615.ENSCAFP00000000003", "9615.ENSCAFP00000000002", 800),
                ("9615.ENSCAFP00000000003", "9615.ENSCAFP00000000004", 650),
                ("9615.ENSCAFP00000000004", "9615.ENSCAFP00000000003", 650),
                ("9615.ENSCAFP00000000004", "9615.ENSCAFP00000000005", 900),
                ("9615.ENSCAFP00000000005", "9615.ENSCAFP00000000004", 900),
            ]
            for a, b, score in edges:
                handle.write(f"{a} {b} {score}\n")

        self.disease = self.root / "disease.txt"
        self.disease.write_text(
            "9615.ENSCAFP00000000001\n"
            "9615.ENSCAFP00000000002\n"
            "9615.ENSCAFP00000000002\n\n",
            encoding="utf-8",
        )
        self.targets = self.root / "targets.xlsx"
        workbook = openpyxl.Workbook()
        sheet = workbook.active
        sheet.append(["drug", "dog_target_protein", "string_protein_id"])
        sheet.append(["Drug A", "T3", "9615.ENSCAFP00000000003"])
        sheet.append(["Drug B", "T4", "9615.ENSCAFP00000000004"])
        sheet.append(["Drug B", "UNKNOWN", None])
        workbook.save(self.targets)

    def tearDown(self):
        self.temp.cleanup()

    def run_cli(self, output, cutoff, *extra):
        command = [
            sys.executable,
            str(Path(run_proximity.__file__).resolve()),
            "--network",
            str(self.network),
            "--disease",
            str(self.disease),
            "--drug-targets",
            str(self.targets),
            "--score-cutoff",
            cutoff,
            "--output",
            str(output),
            "--seed",
            "123",
            *extra,
        ]
        env = os.environ.copy()
        return subprocess.run(command, env=env, text=True, capture_output=True)

    def test_dry_run_qc_and_none_filter(self):
        out_700 = self.root / "dry700"
        result = self.run_cli(out_700, "700", "--dry-run")
        self.assertEqual(result.returncode, 0, result.stderr)
        network_qc = (out_700 / "network_qc.txt").read_text()
        self.assertIn("score_cutoff=700", network_qc)
        self.assertIn("retained_unique_edges=3", network_qc)
        self.assertIn("duplicate_edge_records=3", network_qc)
        excluded_targets = (out_700 / "excluded_drug_targets.tsv").read_text()
        self.assertIn("missing_string_protein_id", excluded_targets)
        disease_qc = (out_700 / "disease_set_qc.tsv").read_text()
        self.assertIn("duplicate_occurrences_beyond_first\t1", disease_qc)
        parameters = (out_700 / "analysis_parameters.txt").read_text()
        self.assertNotIn(str(self.root), parameters)

        out_none = self.root / "drynone"
        result = self.run_cli(out_none, "none", "--dry-run")
        self.assertEqual(result.returncode, 0, result.stderr)
        network_qc = (out_none / "network_qc.txt").read_text()
        self.assertIn("score_cutoff=none", network_qc)
        self.assertIn("retained_unique_edges=4", network_qc)

    def test_small_netmedpy_analysis(self):
        output = self.root / "analysis"
        result = self.run_cli(
            output,
            "none",
            "--n-randomizations",
            "20",
            "--n-processors",
            "1",
            "--no-save-distance-matrix",
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        results = (output / "proximity_results.tsv").read_text()
        self.assertIn("observed_dc", results)
        self.assertIn("Drug A", results)
        self.assertIn("analyzed", results)
        parameters = (output / "analysis_parameters.txt").read_text()
        self.assertIn("netmedpy_version=0.2.0.1", parameters)
        self.assertIn("analysis_status=complete", parameters)
        self.assertTrue((output / "netmedpy_null_distributions.tsv.gz").is_file())

        repeat_output = self.root / "analysis_repeat"
        repeat = self.run_cli(
            repeat_output,
            "none",
            "--n-randomizations",
            "20",
            "--n-processors",
            "1",
            "--no-save-distance-matrix",
        )
        self.assertEqual(repeat.returncode, 0, repeat.stderr)
        self.assertEqual(
            (output / "proximity_results.tsv").read_text(),
            (repeat_output / "proximity_results.tsv").read_text(),
        )
        with gzip.open(
            output / "netmedpy_null_distributions.tsv.gz", "rt"
        ) as first, gzip.open(
            repeat_output / "netmedpy_null_distributions.tsv.gz", "rt"
        ) as second:
            self.assertEqual(first.read(), second.read())


if __name__ == "__main__":
    unittest.main()
