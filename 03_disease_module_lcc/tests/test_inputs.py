"""Offline regression checks for graph construction and disease-ID parsing."""

import gzip
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import canine_neuralgia_LCC_netmedpy as analysis


class InputTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)

    def graph_file(self, text):
        path = self.root / "links.txt.gz"
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            handle.write(text)
        return path

    def test_score_boundary_duplicates_and_self_loops(self):
        path = self.graph_file(
            "protein1 protein2 combined_score\n"
            "A B 899\nA C 900\nC A 950\nA C 999\nD D 901\n"
        )
        graph, qc = analysis.read_string_graph(path)
        self.assertEqual(set(graph.nodes), {"A", "C", "D"})
        self.assertEqual(graph.number_of_edges(), 2)
        self.assertEqual(graph["A"]["C"]["combined_score"], 900)
        self.assertEqual(qc["input_edge_records"], 4)
        self.assertEqual(qc["same_orientation_duplicates"], 1)
        self.assertEqual(qc["reversed_duplicates"], 1)
        self.assertEqual(qc["self_loops"], 1)
        self.assertEqual(qc["minimum_score"], 900)
        self.assertEqual(qc["maximum_score"], 999)
        self.assertEqual(qc["mean_score"], 937.5)

    def test_empty_retained_graph_is_rejected(self):
        path = self.graph_file("protein1 protein2 combined_score\nA B 899\n")
        with self.assertRaisesRegex(ValueError, "No interactome edges"):
            analysis.read_string_graph(path)

    def test_unexpected_header_is_rejected(self):
        path = self.graph_file("source target score\nA B 900\n")
        with self.assertRaisesRegex(ValueError, "Unexpected interactome header"):
            analysis.read_string_graph(path)

    def test_malformed_score_is_rejected(self):
        path = self.graph_file("protein1 protein2 combined_score\nA B unknown\n")
        with self.assertRaisesRegex(ValueError, "Invalid combined_score"):
            analysis.read_string_graph(path)

    def test_disease_ids_preserve_first_occurrence_order(self):
        path = self.root / "genes.txt"
        path.write_text("B\n\n A \nB\nC", encoding="utf-8")
        self.assertEqual(analysis.read_disease_ids(path), (["B", "A", "C"], 1))

    def test_multiple_disease_columns_are_rejected(self):
        path = self.root / "genes.txt"
        path.write_text("A B\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "expected 1"):
            analysis.read_disease_ids(path)


if __name__ == "__main__":
    unittest.main()
