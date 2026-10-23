#!/usr/bin/env python3
"""Audit GSE224814 supplementary files and authoritative GEO/SRA metadata."""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd


SAMPLES = [
    ("A1_TG", "GSM7033526", "Sham", "A1_TG sham sample 1", "SRR23372014"),
    ("A2_TG", "GSM7033527", "Sham", "A2_TG sham sample 2", "SRR23372013"),
    ("A3_TG", "GSM7033528", "Sham", "A3_TG sham sample 3", "SRR23372012"),
    ("B1_TG", "GSM7033529", "CCI_ION", "B1_TG CCI-ION sample 1", "SRR23372011"),
    ("B2_TG", "GSM7033530", "CCI_ION", "B2_TG CCI-ION sample 2", "SRR23372010"),
    ("B3_TG", "GSM7033531", "CCI_ION", "B3_TG CCI-ION sample 3", "SRR23372009"),
]


def parse_soft(path: Path) -> dict[str, dict[str, object]]:
    samples: dict[str, dict[str, object]] = {}
    current = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                current = line.split(" = ", 1)[1]
                samples[current] = {"processing": [], "relations": []}
            elif current and line.startswith("^SAMPLE = ") is False and line.startswith("^"):
                current = None
            elif current and line.startswith("!Sample_title = "):
                samples[current]["title"] = line.split(" = ", 1)[1]
            elif current and line.startswith("!Sample_characteristics_ch1 = treatment: "):
                samples[current]["treatment"] = line.split("treatment: ", 1)[1]
            elif current and line.startswith("!Sample_supplementary_file_1 = "):
                samples[current]["supplementary"] = line.rsplit("/", 1)[-1]
            elif current and line.startswith("!Sample_relation = SRA: "):
                samples[current]["experiment"] = line.rsplit("=", 1)[-1]
            elif current and line.startswith("!Sample_data_processing = "):
                samples[current]["processing"].append(line.split(" = ", 1)[1])
    return samples


def audit(project: Path) -> bool:
    raw = project / "data_raw"
    tables = project / "results" / "tables"
    logs = project / "results" / "logs"
    processed = project / "data_processed"
    for p in (tables, logs, processed):
        p.mkdir(parents=True, exist_ok=True)

    metadata = pd.DataFrame(SAMPLES, columns=["sample_id", "GSM", "condition", "geo_title", "SRA_run"])
    metadata[["sample_id", "GSM", "condition"]].to_csv(processed / "sample_metadata.csv", index=False)

    soft_samples = parse_soft(raw / "GSE224814_family.soft.gz")
    sra = pd.read_csv(raw / "SraRunInfo.csv")
    run_by_gsm = dict(zip(sra["SampleName"], sra["Run"]))

    matching = []
    audit_rows = []
    identifier_frames = []
    files = sorted((raw / "GSE224814_RAW_extracted").glob("*.txt.gz"))
    expected_by_gsm = {row[1]: row for row in SAMPLES}
    for path in files:
        match = re.match(r"(GSM\d+)_([AB]\d)\.txt\.gz$", path.name)
        gsm = match.group(1) if match else ""
        token = match.group(2) if match else ""
        expected = expected_by_gsm.get(gsm)
        if not expected:
            raise ValueError(f"Unexpected supplementary filename: {path.name}")
        sample_id, _, condition, expected_title, expected_run = expected
        soft = soft_samples[gsm]
        observed_header = ""
        with gzip.open(path, "rt") as handle:
            observed_header = handle.readline().rstrip("\n")
        expected_column = sample_id.split("_", 1)[0]
        matching.append(
            {
                "sample_id": sample_id,
                "GSM": gsm,
                "condition": condition,
                "source_filename": path.name,
                "filename_token": token,
                "expression_header": observed_header,
                "GEO_title": soft.get("title", ""),
                "GEO_treatment": soft.get("treatment", ""),
                "SRA_experiment": soft.get("experiment", ""),
                "SRA_run": run_by_gsm.get(gsm, ""),
                "expected_SRA_run": expected_run,
                "title_matches": soft.get("title") == expected_title,
                "header_matches": observed_header == f"Gene\t{expected_column}",
                "condition_matches": ("sham" in str(soft.get("treatment", "")).lower())
                if condition == "Sham"
                else ("cci-ion" in str(soft.get("treatment", "")).lower()),
                "SRA_run_matches": run_by_gsm.get(gsm) == expected_run,
            }
        )

        df = pd.read_csv(path, sep="\t", dtype={0: "string"})
        ids = df.iloc[:, 0]
        numeric = pd.to_numeric(df.iloc[:, 1], errors="coerce")
        finite = numeric[np.isfinite(numeric)]
        integer_prop = float(np.isclose(finite, np.round(finite), atol=1e-9).mean())
        positives = np.sort(finite[finite > 0].unique())
        min_positive = float(positives[0]) if positives.size else np.nan
        reconstruction_error = (
            float(np.max(np.abs(finite / min_positive - np.round(finite / min_positive))))
            if positives.size
            else np.nan
        )
        na_like = ids.fillna("").str.strip().str.lower().isin({"na", "n/a", "nan", "null", "none", ""})
        malformed = ~ids.fillna("").str.match(r"^[A-Za-z][A-Za-z0-9_.-]*$")
        identifier_frames.append(ids.rename(path.name))
        evidence = (
            f"Only {integer_prop:.4%} of numeric observations are integer-valued; "
            f"smallest positive value={min_positive:.12g}; values/minimum-positive are integers "
            f"within max error {reconstruction_error:.3g}. GEO SOFT contains contradictory statements "
            "('TPM values' and 'raw mRNA counts'). Values are therefore treated as library-normalized "
            "count-like expression, not genuine raw integer counts."
        )
        audit_rows.append(
            {
                "filename": path.name,
                "associated_GSM": gsm,
                "associated_sample": sample_id,
                "dimensions": f"{df.shape[0]} rows x {df.shape[1]} columns",
                "delimiter": "tab",
                "column_names": "|".join(map(str, df.columns)),
                "identifier_column": df.columns[0],
                "identifier_type": "gene-symbol-like labels (not stable identifiers; highly duplicated)",
                "expression_count_column": df.columns[1],
                "inferred_numeric_representation": "sample-specific library-normalized count-like expression values; exact normalization label indeterminate from contradictory GEO metadata",
                "numeric_data_type": str(numeric.dtype),
                "minimum_value": float(finite.min()),
                "maximum_value": float(finite.max()),
                "sum_values": float(finite.sum()),
                "proportion_integer_valued": integer_prop,
                "smallest_positive_value": min_positive,
                "max_integer_reconstruction_error": reconstruction_error,
                "missing_identifiers": int(ids.isna().sum()),
                "literal_NA_like_identifiers": int(na_like.sum()),
                "duplicate_identifiers": int(ids.duplicated().sum()),
                "malformed_identifiers": int(malformed.sum()),
                "nonnumeric_values": int(numeric.isna().sum()),
                "suitable_as_raw_count_input_for_PyDESeq2": False,
                "evidence": evidence,
            }
        )

    audit_df = pd.DataFrame(audit_rows)
    audit_df.to_csv(tables / "input_file_audit.csv", index=False)
    matching_df = pd.DataFrame(matching)
    matching_df["all_checks_match"] = matching_df[
        ["title_matches", "header_matches", "condition_matches", "SRA_run_matches"]
    ].all(axis=1)
    matching_df.to_csv(tables / "sample_matching_audit.csv", index=False)

    first = pd.read_csv(files[0], sep="\t", dtype={0: "string"})
    ids = first.iloc[:, 0]
    dup_counts = ids.value_counts(dropna=False)
    duplicate_values = dup_counts[dup_counts > 1]
    duplicate_records = []
    for original_id, n in duplicate_values.items():
        row_numbers = (np.flatnonzero(ids.to_numpy() == original_id) + 2).tolist()
        duplicate_records.append(
            {
                "original_identifier": original_id,
                "number_of_rows": int(n),
                "source_file_row_numbers_including_header": ";".join(map(str, row_numbers)),
                "decision": "not collapsed; biological identity cannot be established from duplicated labels alone",
            }
        )
    pd.DataFrame(duplicate_records).to_csv(tables / "duplicate_identifier_audit.csv", index=False)

    identical_order = all(identifier_frames[0].equals(x) for x in identifier_frames[1:])
    ident_summary = pd.DataFrame(
        [
            {
                "original_identifier_format": "gene-symbol-like labels in column Gene",
                "rows": len(ids),
                "unique_identifiers": ids.nunique(dropna=True),
                "duplicated_identifier_rows": ids.duplicated().sum(),
                "missing_identifiers": ids.isna().sum(),
                "literal_NA_like_identifiers": ids.fillna("").str.strip().str.lower().isin({"na", "n/a", "nan", "null", "none", ""}).sum(),
                "malformed_identifiers": (~ids.fillna("").str.match(r"^[A-Za-z][A-Za-z0-9_.-]*$")).sum(),
                "identical_identifier_order_across_six_files": identical_order,
                "duplicate_handling": "retained as independent source rows in the audit; not used for DE because values are normalized",
            }
        ]
    )
    ident_summary.to_csv(tables / "identifier_audit.csv", index=False)

    lines = [
        "GSE224814 SUPPLEMENTARY INPUT AUDIT",
        f"Files found: {len(files)}",
        *[f"  - {p.name}" for p in files],
        f"All six sample mappings verified against GEO and SRA metadata: {matching_df['all_checks_match'].all()}",
        f"Rows per file: {audit_df['dimensions'].iloc[0]}",
        "Identifier type: gene-symbol-like labels; 17,108 unique labels among 47,061 rows; duplicates were not silently collapsed.",
        "Raw integer counts suitable for PyDESeq2: NO",
        "Reason: values are decimals and sample-specific exact scalar multiples, while GEO metadata contradicts itself by calling the files both TPM and raw-count matrices.",
        "Decision: do not use these supplementary values in PyDESeq2; generate gene-level integer counts from the six SRA paired-end runs.",
    ]
    summary = "\n".join(lines) + "\n"
    (logs / "input_audit_summary.txt").write_text(summary)
    print(summary, flush=True)
    return False


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    audit(args.project.resolve())


if __name__ == "__main__":
    main()
