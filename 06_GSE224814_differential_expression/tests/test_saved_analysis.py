"""Check archived counts, DEG selection, annotation, and sample QC offline."""

import csv
import gzip
import math
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
SAMPLES = ["A1_TG", "A2_TG", "A3_TG", "B1_TG", "B2_TG", "B3_TG"]
DE = ROOT / "results/differential_expression"


def rows(path: Path, **kwargs) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, **kwargs))


def result_rows(suffix: str) -> list[dict[str, str]]:
    return rows(DE / f"GSE224814_CCI_ION_vs_Sham_{suffix}.csv")


class SavedAnalysisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.raw = rows(ROOT / "data_processed/GSE224814_raw_count_matrix.csv")
        cls.filtered = rows(ROOT / "data_processed/GSE224814_filtered_counts.csv")
        cls.all_results = result_rows("all_results")
        cls.degs = result_rows("DEGs")

    def test_sample_design(self):
        metadata = rows(ROOT / "data_processed/sample_metadata.csv")
        self.assertEqual([row["sample_id"] for row in metadata], SAMPLES)
        self.assertEqual([row["condition"] for row in metadata], ["Sham"] * 3 + ["CCI_ION"] * 3)
        self.assertEqual(len({row["GSM"] for row in metadata}), 6)
        self.assertEqual(list(self.raw[0]), ["original_id", *SAMPLES])

    def test_raw_counts(self):
        self.assertEqual(len(self.raw), 32623)
        self.assertEqual(len({row["original_id"] for row in self.raw}), len(self.raw))
        self.assertTrue(all(int(row[sample]) >= 0 for row in self.raw for sample in SAMPLES))

    def test_featurecounts_match_matrix(self):
        for sample in SAMPLES:
            path = ROOT / "data_processed/featurecounts" / f"{sample}.featureCounts.txt"
            with path.open() as handle:
                records = list(csv.DictReader((line for line in handle if not line.startswith("#")), delimiter="\t"))
            self.assertEqual(len(records), len(self.raw))
            column = list(records[0])[-1]
            for observed, expected in zip(records, self.raw):
                self.assertEqual(observed["Geneid"], expected["original_id"])
                self.assertEqual(int(observed[column]), int(expected[sample]))

    def test_low_count_filter(self):
        expected = [row for row in self.raw if sum(int(row[sample]) >= 10 for sample in SAMPLES) >= 3]
        self.assertEqual(len(expected), 15029)
        self.assertEqual(self.filtered, expected)

    def test_complete_results(self):
        self.assertEqual(len(self.all_results), 15029)
        self.assertEqual({row["original_id"] for row in self.all_results}, {row["original_id"] for row in self.filtered})
        self.assertEqual(sum(not row["padj"] for row in self.all_results), 18)
        self.assertEqual(sum(bool(row["padj"]) and float(row["padj"]) < 0.05 for row in self.all_results), 2343)

    def test_primary_deg_rule_and_values(self):
        expected = {
            row["original_id"]: row for row in self.all_results
            if row["padj"] and float(row["padj"]) < 0.05 and abs(float(row["log2FoldChange"])) >= 1
        }
        self.assertEqual(len(expected), 431)
        self.assertEqual({row["original_id"] for row in self.degs}, set(expected))
        for row in self.degs:
            self.assertEqual({key: row[key] for key in expected[row["original_id"]]}, expected[row["original_id"]])
            self.assertEqual(row["direction"], "Up" if float(row["log2FoldChange"]) >= 1 else "Down")

    def test_directional_subsets(self):
        up, down = result_rows("upregulated"), result_rows("downregulated")
        self.assertEqual(len(up), 368)
        self.assertEqual(len(down), 63)
        self.assertEqual(up, [row for row in self.degs if row["direction"] == "Up"])
        self.assertEqual(down, [row for row in self.degs if row["direction"] == "Down"])

    def test_annotation_retains_degs_and_statistics(self):
        annotated = rows(DE / "GSE224814_DEGs_for_ortholog_mapping.csv")
        self.assertEqual(len(annotated), len(self.degs))
        self.assertEqual([row["original_id"] for row in annotated], [row["original_id"] for row in self.degs])
        for observed, original in zip(annotated, self.degs):
            for name in ("baseMean", "log2FoldChange", "pvalue", "padj", "direction"):
                self.assertEqual(observed[name], original[name])
        qc = rows(ROOT / "results/tables/gene_annotation_QC.csv")
        self.assertEqual(sum(row["mapping_status"] == "uniquely_mapped" for row in qc), 430)
        self.assertEqual(sum(row["mapping_status"] == "ambiguous" for row in qc), 1)
        self.assertEqual(sum(row["mapping_status"] == "unresolved" for row in qc), 0)

    def test_rat_annotation_snapshot(self):
        with gzip.open(ROOT / "data_raw/annotation/gene2ensembl.rat.tsv.gz", "rt") as handle:
            table = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(table), 36316)
        self.assertTrue(all(row["#tax_id"] == "10116" for row in table))

    def test_vst_sample_distances(self):
        vst = rows(ROOT / "data_processed/GSE224814_VST_expression.csv")
        self.assertEqual([row["original_id"] for row in vst], [row["original_id"] for row in self.filtered])
        vectors = {sample: [float(row[sample]) for row in vst] for sample in SAMPLES}
        saved = rows(ROOT / "results/tables/sample_distance_matrix.csv")
        self.assertEqual([row["sample_id"] for row in saved], SAMPLES)
        for row in saved:
            for sample in SAMPLES:
                observed = math.dist(vectors[row["sample_id"]], vectors[sample])
                self.assertTrue(math.isclose(observed, float(row[sample]), rel_tol=1e-12, abs_tol=1e-10))

    def test_no_samples_excluded(self):
        qc = rows(ROOT / "results/tables/PCA_outlier_diagnostics.csv")
        self.assertEqual([row["sample_id"] for row in qc], SAMPLES)
        self.assertTrue(all(row["exploratory_flag"] == "False" for row in qc))

    def test_archived_validation(self):
        checks = rows(ROOT / "results/tables/validation_checks.csv")
        self.assertTrue(checks)
        self.assertTrue(all(row["status"] == "PASS" for row in checks))


if __name__ == "__main__":
    unittest.main()
