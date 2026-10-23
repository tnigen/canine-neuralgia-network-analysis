"""Check command-line safeguards and lossless reference subsetting."""

import contextlib
import gzip
import io
from pathlib import Path
import runpy
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]


class EntrypointTests(unittest.TestCase):
    def test_conflicting_modes_stop_before_writing(self):
        runner = runpy.run_path(str(ROOT / "scripts/run_analysis.py"))
        with tempfile.TemporaryDirectory() as directory:
            project = Path(directory) / "uncreated_project"
            arguments = ["run_analysis.py", "--project", str(project), "--from-counts", "--rebuild-counts"]
            with patch("sys.argv", arguments), contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as caught:
                    runner["main"]()
            self.assertEqual(caught.exception.code, 2)
            self.assertFalse(project.exists())

    def test_missing_counts_stop_before_writing(self):
        runner = runpy.run_path(str(ROOT / "scripts/run_analysis.py"))
        with tempfile.TemporaryDirectory() as directory:
            project = Path(directory) / "uncreated_project"
            arguments = ["run_analysis.py", "--project", str(project), "--from-counts"]
            with patch("sys.argv", arguments), contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as caught:
                    runner["main"]()
            self.assertEqual(caught.exception.code, 2)
            self.assertFalse(project.exists())

    def test_crossref_subset_preserves_fields_and_order(self):
        subset = runpy.run_path(str(ROOT / "scripts/subset_rat_crossrefs.py"))["subset"]
        header = b"#tax_id\tGeneID\tEnsembl_gene_identifier\n"
        first = b"10116\t2\tENSRNOG00000000002\n"
        second = b"10116\t1\tENSRNOG00000000001\n"
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "source.gz"
            output = Path(directory) / "rat.gz"
            with gzip.open(source, "wb") as handle:
                handle.write(header + first + b"9606\t3\tENSG00000000003\n" + second)
            self.assertEqual(subset(source, output), 2)
            with gzip.open(output, "rb") as handle:
                self.assertEqual(handle.read(), header + first + second)
            with self.assertRaises(FileExistsError):
                subset(source, output)


if __name__ == "__main__":
    unittest.main()
