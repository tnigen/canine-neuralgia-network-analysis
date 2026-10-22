#!/usr/bin/env python3
"""Map canine gene symbols to STRING protein identifiers."""

from __future__ import annotations

import csv
import io
import sys
import time
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

try:
    import openpyxl  # noqa: F401  # Required by pandas for .xlsx files.
    import pandas as pd
    import requests
except ModuleNotFoundError as exc:
    print(f"Missing Python package: {exc.name}", file=sys.stderr)
    print(
        "Install the required packages with:\n"
        "    python -m pip install pandas requests openpyxl",
        file=sys.stderr,
    )
    raise SystemExit(1) from exc


SCRIPT_VERSION = "1.0.0"
INPUT_FILE = Path("drug_target_dog.xlsx")
OUTPUT_FILE = Path("drug_target_dog_string.xlsx")
UNMAPPED_FILE = Path("drug_target_dog_string_unmapped.xlsx")
QC_FILE = Path("drug_target_dog_string_mapping_qc.xlsx")

TARGET_COLUMN = "dog_target_protein"
STRING_ID_COLUMN = "string_protein_id"
REQUIRED_INPUT_COLUMNS = ["drug", TARGET_COLUMN]

STRING_API_URL = "https://string-db.org/api/tsv/get_string_ids"
SPECIES_ID = 9615
EXPECTED_STRING_PREFIX = "9615.ENSCAFP"
CALLER_IDENTITY = "dog_target_string_mapper"

BATCH_SIZE = 100
REQUEST_TIMEOUT = (10, 60)  # (connection timeout, response timeout), in seconds
DELAY_BETWEEN_REQUESTS = 1.0


@dataclass(frozen=True)
class MappingResult:
    string_protein_id: str
    mapping_status: str
    details: str = ""


def clean_gene_symbol(value: object) -> str | None:
    """Return a whitespace-trimmed symbol, or None for an empty Excel cell."""

    if pd.isna(value):
        return None
    symbol = str(value).strip()
    return symbol or None


def chunks(items: list[str], size: int) -> Iterator[list[str]]:
    for start in range(0, len(items), size):
        yield items[start : start + size]


def parse_string_response(
    response_text: str, batch: list[str]
) -> dict[str, list[dict[str, str]]]:
    """Group STRING response rows by submitted gene symbol."""

    grouped: dict[str, list[dict[str, str]]] = {gene: [] for gene in batch}
    if not response_text.strip():
        return grouped

    reader = csv.DictReader(io.StringIO(response_text), delimiter="\t")
    required_fields = {"queryIndex", "stringId", "ncbiTaxonId", "preferredName"}
    if reader.fieldnames is None or not required_fields.issubset(reader.fieldnames):
        raise ValueError(
            "Unexpected STRING response columns: "
            f"{reader.fieldnames if reader.fieldnames is not None else 'none'}"
        )

    for row in reader:
        gene: str | None = None
        try:
            query_index = int(row.get("queryIndex", ""))
            if 0 <= query_index < len(batch):
                gene = batch[query_index]
        except (TypeError, ValueError):
            pass

        # Fall back to the echoed identifier if queryIndex is unavailable.
        if gene is None:
            query_item = (row.get("queryItem") or "").strip()
            if query_item in grouped:
                gene = query_item

        if gene is None:
            raise ValueError("STRING returned a row that cannot be linked to an input gene")
        grouped[gene].append(row)

    return grouped


def choose_mapping(gene: str, candidates: list[dict[str, str]]) -> MappingResult:
    """Select one valid canine STRING mapping or report ambiguity."""

    if not candidates:
        return MappingResult("", "unmapped", "STRING returned no mapping")

    dog_candidates = [
        row for row in candidates if str(row.get("ncbiTaxonId", "")).strip() == str(SPECIES_ID)
    ]
    if not dog_candidates:
        returned_taxa = sorted(
            {str(row.get("ncbiTaxonId", "")).strip() for row in candidates}
        )
        return MappingResult(
            "",
            "error",
            "STRING returned no dog-taxonomy candidate; returned taxa: "
            + ", ".join(returned_taxa),
        )

    valid_by_id: dict[str, dict[str, str]] = {}
    invalid_ids: list[str] = []
    for row in dog_candidates:
        string_id = str(row.get("stringId", "")).strip()
        if string_id.startswith(EXPECTED_STRING_PREFIX):
            valid_by_id.setdefault(string_id, row)
        else:
            invalid_ids.append(string_id or "<blank>")

    valid_candidates = list(valid_by_id.values())
    if not valid_candidates:
        return MappingResult(
            "",
            "error",
            f"Dog mapping failed required prefix {EXPECTED_STRING_PREFIX}; "
            f"returned IDs: {', '.join(invalid_ids)}",
        )

    if len(valid_candidates) == 1:
        return MappingResult(valid_candidates[0]["stringId"].strip(), "mapped")

    exact_name_matches = [
        row
        for row in valid_candidates
        if str(row.get("preferredName", "")).strip() == gene
    ]
    if len(exact_name_matches) == 1:
        selected_id = exact_name_matches[0]["stringId"].strip()
        print(
            f"[multiple mappings] {gene}: selected {selected_id} because it is "
            "the only dog candidate whose preferredName exactly matches the input."
        )
        return MappingResult(selected_id, "mapped")

    candidate_ids = ", ".join(valid_by_id)
    return MappingResult(
        "",
        "error",
        "Ambiguous dog mappings; no deterministic unique selection was possible: "
        + candidate_ids,
    )


def map_unique_genes(unique_genes: list[str]) -> dict[str, MappingResult]:
    results: dict[str, MappingResult] = {}
    batches = list(chunks(unique_genes, BATCH_SIZE))

    with requests.Session() as session:
        session.headers.update(
            {"User-Agent": "dog-target-string-mapper/1.0 (STRING identifier mapping)"}
        )

        for batch_number, batch in enumerate(batches, start=1):
            params = {
                "identifiers": "\r".join(batch),
                "species": SPECIES_ID,
                "echo_query": 1,
                "caller_identity": CALLER_IDENTITY,
            }

            try:
                response = session.post(
                    STRING_API_URL,
                    data=params,
                    timeout=REQUEST_TIMEOUT,
                )

                # STRING returns 404 when no submitted identifiers resolve.
                if response.status_code == 404:
                    for gene in batch:
                        results[gene] = MappingResult(
                            "", "unmapped", "STRING returned HTTP 404 (no identifiers resolved)"
                        )
                else:
                    response.raise_for_status()
                    grouped = parse_string_response(response.text, batch)
                    for gene in batch:
                        results[gene] = choose_mapping(gene, grouped[gene])

            except (requests.RequestException, ValueError) as exc:
                message = f"STRING API request/response error: {exc}"
                print(
                    f"[error] Batch {batch_number} failed for genes: "
                    f"{', '.join(batch)}\n        {message}",
                    file=sys.stderr,
                )
                for gene in batch:
                    results[gene] = MappingResult("", "error", message)

            if batch_number < len(batches):
                time.sleep(DELAY_BETWEEN_REQUESTS)

    return results


def write_excel_outputs(
    original_df: "pd.DataFrame",
    cleaned_symbols: list[str | None],
    unique_genes: list[str],
    results: dict[str, MappingResult],
) -> None:
    output_df = original_df.copy()
    mapped_ids = [
        results[symbol].string_protein_id if symbol is not None else ""
        for symbol in cleaned_symbols
    ]

    # Preserve the original row and column order.
    target_position = output_df.columns.get_loc(TARGET_COLUMN)
    output_df.insert(target_position + 1, STRING_ID_COLUMN, mapped_ids)

    qc_rows = [
        {
            TARGET_COLUMN: gene,
            STRING_ID_COLUMN: results[gene].string_protein_id,
            "mapping_status": results[gene].mapping_status,
        }
        for gene in unique_genes
    ]
    qc_df = pd.DataFrame(
        qc_rows,
        columns=[TARGET_COLUMN, STRING_ID_COLUMN, "mapping_status"],
    )

    # Keep unresolved identifiers and request errors in one audit file.
    unmapped_rows = [
        {
            TARGET_COLUMN: gene,
            "mapping_status": results[gene].mapping_status,
            "details": results[gene].details,
        }
        for gene in unique_genes
        if results[gene].mapping_status != "mapped"
    ]
    unmapped_df = pd.DataFrame(
        unmapped_rows,
        columns=[TARGET_COLUMN, "mapping_status", "details"],
    )

    output_df.to_excel(OUTPUT_FILE, index=False, engine="openpyxl")
    unmapped_df.to_excel(UNMAPPED_FILE, index=False, engine="openpyxl")
    qc_df.to_excel(QC_FILE, index=False, engine="openpyxl")


def main() -> int:
    if not INPUT_FILE.is_file():
        print(f"Input file not found: {INPUT_FILE}", file=sys.stderr)
        return 1

    try:
        original_df = pd.read_excel(INPUT_FILE, engine="openpyxl")
    except Exception as exc:
        print(f"Could not read {INPUT_FILE}: {exc}", file=sys.stderr)
        return 1

    missing_columns = [
        column for column in REQUIRED_INPUT_COLUMNS if column not in original_df.columns
    ]
    if missing_columns:
        print(
            f"Required input column(s) not found: {missing_columns}. "
            f"Available columns: {list(original_df.columns)}",
            file=sys.stderr,
        )
        return 1
    if list(original_df.columns).count(TARGET_COLUMN) != 1:
        print(f"Column {TARGET_COLUMN!r} must occur exactly once.", file=sys.stderr)
        return 1
    if STRING_ID_COLUMN in original_df.columns:
        print(
            f"Input already contains {STRING_ID_COLUMN!r}; refusing to overwrite it.",
            file=sys.stderr,
        )
        return 1

    # Trim a separate lookup copy without changing the input values.
    cleaned_symbols = [
        clean_gene_symbol(value) for value in original_df[TARGET_COLUMN].tolist()
    ]
    unique_genes = list(
        dict.fromkeys(symbol for symbol in cleaned_symbols if symbol is not None)
    )

    print(
        f"Mapping {len(unique_genes)} unique dog gene symbols with STRING "
        f"(taxonomy ID {SPECIES_ID})..."
    )
    results = map_unique_genes(unique_genes)

    try:
        write_excel_outputs(original_df, cleaned_symbols, unique_genes, results)
    except Exception as exc:
        print(f"Could not write an output Excel file: {exc}", file=sys.stderr)
        return 1

    mapped_count = sum(
        result.mapping_status == "mapped" for result in results.values()
    )
    error_count = sum(
        result.mapping_status == "error" for result in results.values()
    )
    not_mapped_count = len(unique_genes) - mapped_count
    mapping_rate = (
        (100.0 * mapped_count / len(unique_genes)) if unique_genes else 0.0
    )

    print("\nMapping summary")
    print(f"Total input rows: {len(original_df)}")
    print(f"Unique dog gene symbols: {len(unique_genes)}")
    print(f"Successfully mapped: {mapped_count}")
    print(f"Unmapped: {not_mapped_count}")
    print(f"Errors: {error_count}")
    print(f"Mapping rate: {mapping_rate:.1f}%")
    print(f"Output file: {OUTPUT_FILE}")
    print(f"Unmapped file: {UNMAPPED_FILE}")
    print(f"QC file: {QC_FILE}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
