#!/usr/bin/env python3
"""Map human drug-target symbols to strict one-to-one dog orthologs with BioMart."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import platform
import sys
import time
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable, Sequence
from urllib.parse import urlsplit, urlunsplit

import pandas as pd
import requests


SCRIPT_VERSION = "1.0.1"
DEFAULT_BIOMART_URL = "https://mart.ensembl.org/biomart/martservice"
MART_NAME = "ENSEMBL_MART_ENSEMBL"
STRICT_ONE2ONE_VALUE = "ortholog_one2one"
SOURCE_SPECIES = "Homo sapiens"
TARGET_SPECIES = "Dog / Canis lupus familiaris"

QC_COLUMNS = [
    "drug",
    "human_target_protein",
    "human_hgnc_id",
    "human_ensembl_gene_id",
    "dog_target_protein",
    "dog_ensembl_gene_id",
    "homology_type",
    "mapping_status",
    "notes",
]
FINAL_COLUMNS = ["drug", "dog_target_protein"]
SUMMARY_COLUMNS = [
    "drug",
    "number_of_human_targets",
    "number_of_mapped_dog_targets",
    "number_of_unmapped_human_targets",
    "mapping_percentage",
]


class BioMartError(RuntimeError):
    """Raised when BioMart validation or querying cannot be completed safely."""


@dataclass(frozen=True)
class Field:
    identifier: str
    label: str
    description: str = ""
    page: str = ""


@dataclass(frozen=True)
class BioMartConfig:
    requested_service_url: str
    resolved_service_url: str
    mart: str
    mart_display_name: str
    dataset: str
    dataset_display_name: str
    symbol_filter: Field
    human_gene_id_filter: Field
    human_gene_id: Field
    human_gene_name: Field
    match_symbol: Field
    human_hgnc_id: Field
    dog_gene_id: Field
    dog_gene_name: Field
    homology_type: Field

    @property
    def resolution_attribute_ids(self) -> list[str]:
        return list(
            dict.fromkeys(
                [
                    self.human_gene_id.identifier,
                    self.human_gene_name.identifier,
                    self.match_symbol.identifier,
                    self.human_hgnc_id.identifier,
                ]
            )
        )

    @property
    def homology_attribute_ids(self) -> list[str]:
        ordered = [
            self.human_gene_id.identifier,
            self.human_gene_name.identifier,
            self.dog_gene_id.identifier,
            self.dog_gene_name.identifier,
            self.homology_type.identifier,
        ]
        return list(dict.fromkeys(ordered))


def _service_url_without_query(url: str) -> str:
    parts = urlsplit(url)
    return urlunsplit((parts.scheme, parts.netloc, parts.path, "", ""))


def _nonempty_lines(text: str) -> list[str]:
    return [line for line in text.splitlines() if line.strip()]


def _parse_fields(text: str) -> list[Field]:
    fields: list[Field] = []
    # BioMart filter rows may contain a very large enumerated-options field,
    # exceeding Python csv's conservative default field-size limit. These
    # metadata endpoints are plain tab-separated text without CSV quoting.
    for line in text.splitlines():
        row = line.split("\t")
        if not row or not row[0].strip():
            continue
        padded = row + [""] * max(0, 4 - len(row))
        fields.append(
            Field(
                identifier=padded[0].strip(),
                label=padded[1].strip(),
                description=padded[2].strip(),
                page=padded[3].strip(),
            )
        )
    return fields


def _select_field(
    fields: Sequence[Field],
    purpose: str,
    exact_labels: Sequence[str],
    preferred_ids: Sequence[str] = (),
    required_page: str | None = None,
) -> Field:
    wanted_labels = {label.casefold() for label in exact_labels}
    candidates = [f for f in fields if f.label.casefold() in wanted_labels]
    if required_page:
        page_matches = [f for f in candidates if f.page.casefold() == required_page.casefold()]
        if page_matches:
            candidates = page_matches
    for identifier in preferred_ids:
        for field in candidates:
            if field.identifier == identifier:
                return field
    unique = {f.identifier: f for f in candidates}
    if len(unique) == 1:
        return next(iter(unique.values()))
    if not candidates:
        raise BioMartError(
            f"Could not identify the current BioMart field for {purpose}. "
            f"Expected one of these live labels: {', '.join(exact_labels)}."
        )
    raise BioMartError(
        f"BioMart configuration is ambiguous for {purpose}: "
        + ", ".join(sorted(unique))
    )


class BioMartClient:
    def __init__(self, service_url: str, timeout: float, retries: int) -> None:
        self.service_url = service_url
        self.timeout = timeout
        self.retries = retries
        self.session = requests.Session()
        self.session.headers.update(
            {
                "User-Agent": (
                    f"drug-target-dog-biomart/{SCRIPT_VERSION} "
                    "(scientific orthology mapping)"
                )
            }
        )

    def _request(self, method: str, **kwargs: object) -> requests.Response:
        last_error: Exception | None = None
        for attempt in range(self.retries + 1):
            try:
                response = self.session.request(
                    method,
                    self.service_url,
                    timeout=self.timeout,
                    allow_redirects=True,
                    **kwargs,
                )
                if response.status_code in {429, 500, 502, 503, 504}:
                    raise requests.HTTPError(
                        f"temporary HTTP {response.status_code}", response=response
                    )
                response.raise_for_status()
                body_start = response.text.lstrip()[:200].casefold()
                content_type = response.headers.get("content-type", "").casefold()
                if "<html" in body_start or "text/html" in content_type:
                    raise BioMartError(
                        "BioMart returned an HTML page instead of registry/query data "
                        f"from {response.url}."
                    )
                if "query error" in body_start:
                    raise BioMartError(response.text.strip()[:1000])
                return response
            except (requests.RequestException, BioMartError) as exc:
                last_error = exc
                if attempt == self.retries:
                    break
                time.sleep(min(2**attempt, 8))
        raise BioMartError(
            f"BioMart request failed after {self.retries + 1} attempt(s): {last_error}"
        ) from last_error

    def _get_metadata(self, **params: str) -> requests.Response:
        return self._request("GET", params=params)

    def validate_configuration(self) -> BioMartConfig:
        registry_response = self._get_metadata(type="registry")
        try:
            root = ET.fromstring(registry_response.text)
        except ET.ParseError as exc:
            raise BioMartError(f"BioMart registry is not valid XML: {exc}") from exc

        mart_nodes = [
            node
            for node in root.iter("MartURLLocation")
            if node.attrib.get("name") == MART_NAME
        ]
        if not mart_nodes:
            raise BioMartError(
                f"The Ensembl Genes mart {MART_NAME!r} is absent from the live registry."
            )
        mart_node = mart_nodes[0]
        mart_display_name = mart_node.attrib.get("displayName", MART_NAME)

        datasets_response = self._get_metadata(type="datasets", mart=MART_NAME)
        dataset_rows = list(csv.reader(io.StringIO(datasets_response.text), delimiter="\t"))
        human_rows = [
            row
            for row in dataset_rows
            if len(row) >= 3
            and row[0] == "TableSet"
            and row[2].casefold().startswith("human genes")
        ]
        preferred = [row for row in human_rows if row[1] == "hsapiens_gene_ensembl"]
        selected_rows = preferred or human_rows
        if len(selected_rows) != 1:
            found = ", ".join(row[1] for row in human_rows) or "none"
            raise BioMartError(
                "Could not uniquely identify the Homo sapiens genes dataset in the "
                f"current Ensembl Genes mart (candidates: {found})."
            )
        dataset_row = selected_rows[0]
        dataset = dataset_row[1]
        dataset_display_name = dataset_row[2]

        filters_response = self._get_metadata(type="filters", dataset=dataset)
        attributes_response = self._get_metadata(type="attributes", dataset=dataset)
        filters = _parse_fields(filters_response.text)
        attributes = _parse_fields(attributes_response.text)

        symbol_filter = _select_field(
            filters,
            "human HGNC/gene-symbol input filter",
            ["HGNC symbol(s) [e.g. A1BG]", "Gene Name(s) [e.g. MT-TF]"],
            preferred_ids=["hgnc_symbol", "external_gene_name"],
        )
        human_gene_id_filter = _select_field(
            filters,
            "human Ensembl gene stable-ID filter",
            ["Gene stable ID(s) [e.g. ENSG00000000003]"],
            preferred_ids=["ensembl_gene_id"],
        )
        human_gene_id = _select_field(
            attributes,
            "human gene stable ID",
            ["Gene stable ID"],
            preferred_ids=["ensembl_gene_id"],
            required_page="homologs",
        )
        human_gene_name = _select_field(
            attributes,
            "human gene name",
            ["Gene name"],
            preferred_ids=["external_gene_name"],
            required_page="homologs",
        )
        human_hgnc_id = _select_field(
            attributes,
            "human HGNC identity",
            ["HGNC ID"],
            preferred_ids=["hgnc_id"],
            required_page="feature_page",
        )

        same_as_filter = [
            f for f in attributes if f.identifier == symbol_filter.identifier
        ]
        if same_as_filter:
            match_symbol = same_as_filter[0]
        elif symbol_filter.identifier == human_gene_name.identifier:
            match_symbol = human_gene_name
        else:
            raise BioMartError(
                f"The selected input filter {symbol_filter.identifier!r} has no matching "
                "output attribute, so returned rows cannot be audited against input symbols."
            )

        dog_gene_id = _select_field(
            attributes,
            "dog gene stable ID",
            ["Dog gene stable ID"],
            required_page="homologs",
        )
        dog_gene_name = _select_field(
            attributes,
            "dog gene name",
            ["Dog gene name"],
            required_page="homologs",
        )
        homology_type = _select_field(
            attributes,
            "dog homology type",
            ["Dog homology type"],
            required_page="homologs",
        )

        resolved_urls = [
            registry_response.url,
            datasets_response.url,
            filters_response.url,
            attributes_response.url,
        ]
        resolved_services = {_service_url_without_query(url) for url in resolved_urls}
        if len(resolved_services) != 1:
            raise BioMartError(
                "BioMart metadata requests resolved to inconsistent services: "
                + ", ".join(sorted(resolved_services))
            )
        resolved_service = next(iter(resolved_services))
        self.service_url = resolved_service

        return BioMartConfig(
            requested_service_url=_service_url_without_query(registry_response.request.url),
            resolved_service_url=resolved_service,
            mart=MART_NAME,
            mart_display_name=mart_display_name,
            dataset=dataset,
            dataset_display_name=dataset_display_name,
            symbol_filter=symbol_filter,
            human_gene_id_filter=human_gene_id_filter,
            human_gene_id=human_gene_id,
            human_gene_name=human_gene_name,
            match_symbol=match_symbol,
            human_hgnc_id=human_hgnc_id,
            dog_gene_id=dog_gene_id,
            dog_gene_name=dog_gene_name,
            homology_type=homology_type,
        )

    def query_symbols(
        self, symbols: Sequence[str], config: BioMartConfig
    ) -> pd.DataFrame:
        """Resolve input symbols to human genes using one BioMart attribute page."""
        return self._query(
            values=symbols,
            dataset_id=config.dataset,
            filter_id=config.symbol_filter.identifier,
            attribute_ids=config.resolution_attribute_ids,
            output_columns={
                "human_ensembl_gene_id": config.human_gene_id.identifier,
                "human_gene_name": config.human_gene_name.identifier,
                "match_symbol": config.match_symbol.identifier,
                "human_hgnc_id": config.human_hgnc_id.identifier,
            },
        )

    def query_gene_ids(
        self, gene_ids: Sequence[str], config: BioMartConfig
    ) -> pd.DataFrame:
        """Retrieve dog homology after human symbols have been resolved."""
        return self._query(
            values=gene_ids,
            dataset_id=config.dataset,
            filter_id=config.human_gene_id_filter.identifier,
            attribute_ids=config.homology_attribute_ids,
            output_columns={
                "human_ensembl_gene_id": config.human_gene_id.identifier,
                "human_gene_name": config.human_gene_name.identifier,
                "dog_ensembl_gene_id": config.dog_gene_id.identifier,
                "dog_target_protein": config.dog_gene_name.identifier,
                "homology_type": config.homology_type.identifier,
            },
        )

    def _query(
        self,
        values: Sequence[str],
        dataset_id: str,
        filter_id: str,
        attribute_ids: Sequence[str],
        output_columns: dict[str, str],
    ) -> pd.DataFrame:
        query = ET.Element(
            "Query",
            {
                "virtualSchemaName": "default",
                "formatter": "TSV",
                "header": "0",
                "uniqueRows": "1",
                "count": "",
                "datasetConfigVersion": "0.6",
            },
        )
        dataset = ET.SubElement(
            query, "Dataset", {"name": dataset_id, "interface": "default"}
        )
        ET.SubElement(
            dataset,
            "Filter",
            {"name": filter_id, "value": ",".join(values)},
        )
        for identifier in attribute_ids:
            ET.SubElement(dataset, "Attribute", {"name": identifier})
        xml_query = ET.tostring(query, encoding="unicode")
        response = self._request("POST", data={"query": xml_query})

        rows = list(csv.reader(io.StringIO(response.text), delimiter="\t"))
        bad_rows = [i + 1 for i, row in enumerate(rows) if len(row) != len(attribute_ids)]
        if bad_rows:
            raise BioMartError(
                "BioMart returned malformed TSV rows with unexpected column counts "
                f"(first affected response rows: {bad_rows[:5]})."
            )

        raw = pd.DataFrame(rows, columns=attribute_ids, dtype=str)
        for identifier in attribute_ids:
            raw[identifier] = raw[identifier].fillna("").astype(str).str.strip()

        def column(identifier: str) -> pd.Series:
            if identifier in raw:
                return raw[identifier]
            return pd.Series([""] * len(raw), dtype=str)

        return pd.DataFrame(
            {output_name: column(identifier) for output_name, identifier in output_columns.items()}
        ).drop_duplicates(ignore_index=True)


def read_input(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Input file does not exist: {path}")
    if not path.is_file():
        raise ValueError(f"Input path is not a regular file: {path}")
    suffix = path.suffix.casefold()
    if suffix == ".csv":
        frame = pd.read_csv(path, dtype=str, keep_default_na=False, encoding="utf-8-sig")
    elif suffix == ".xlsx":
        frame = pd.read_excel(path, dtype=str, keep_default_na=False, engine="openpyxl")
    else:
        raise ValueError(
            f"Unsupported input format {path.suffix!r}; use a .csv or .xlsx file."
        )
    required = ["drug", "human_target_protein"]
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError("Missing required input column(s): " + ", ".join(missing))
    cleaned = frame[required].copy()
    for column in required:
        cleaned[column] = cleaned[column].fillna("").astype(str).str.strip()
    return cleaned


def batched(values: Sequence[str], size: int) -> Iterable[list[str]]:
    for start in range(0, len(values), size):
        yield list(values[start : start + size])


def query_with_isolation(
    client: BioMartClient,
    config: BioMartConfig,
    values: Sequence[str],
    query_kind: str,
) -> tuple[list[pd.DataFrame], dict[str, str], int]:
    """Query a batch, recursively isolating individual failures when necessary."""
    try:
        method = client.query_symbols if query_kind == "symbols" else client.query_gene_ids
        return [method(values, config)], {}, 1
    except BioMartError as exc:
        if len(values) == 1:
            return [], {values[0]: str(exc)}, 0
        midpoint = len(values) // 2
        left_frames, left_errors, left_successes = query_with_isolation(
            client, config, values[:midpoint], query_kind
        )
        right_frames, right_errors, right_successes = query_with_isolation(
            client, config, values[midpoint:], query_kind
        )
        return (
            left_frames + right_frames,
            {**left_errors, **right_errors},
            left_successes + right_successes,
        )


def query_all_symbols(
    client: BioMartClient,
    config: BioMartConfig,
    symbols: Sequence[str],
    batch_size: int,
) -> tuple[pd.DataFrame, dict[str, str]]:
    resolution_frames: list[pd.DataFrame] = []
    symbol_errors: dict[str, str] = {}
    successful_resolution_requests = 0
    for batch in batched(symbols, batch_size):
        batch_frames, batch_errors, batch_successes = query_with_isolation(
            client, config, batch, "symbols"
        )
        resolution_frames.extend(batch_frames)
        symbol_errors.update(batch_errors)
        successful_resolution_requests += batch_successes
    if symbols and successful_resolution_requests == 0:
        first_error = next(iter(symbol_errors.values()), "unknown BioMart query failure")
        raise BioMartError(
            "Every symbol-resolution query failed; treating this as a systemic BioMart failure. "
            f"First error: {first_error}"
        )
    resolution_empty = pd.DataFrame(
        columns=[
            "human_ensembl_gene_id",
            "human_gene_name",
            "match_symbol",
            "human_hgnc_id",
        ]
    )
    resolutions = (
        pd.concat(resolution_frames, ignore_index=True).drop_duplicates()
        if resolution_frames
        else resolution_empty
    )

    gene_ids = list(
        dict.fromkeys(
            value for value in resolutions["human_ensembl_gene_id"].tolist() if value
        )
    )
    homology_frames: list[pd.DataFrame] = []
    gene_id_errors: dict[str, str] = {}
    successful_homology_requests = 0
    for batch in batched(gene_ids, batch_size):
        batch_frames, batch_errors, batch_successes = query_with_isolation(
            client, config, batch, "gene_ids"
        )
        homology_frames.extend(batch_frames)
        gene_id_errors.update(batch_errors)
        successful_homology_requests += batch_successes
    if gene_ids and successful_homology_requests == 0:
        first_error = next(iter(gene_id_errors.values()), "unknown BioMart query failure")
        raise BioMartError(
            "Every dog-homology query failed; treating this as a systemic BioMart failure. "
            f"First error: {first_error}"
        )

    # A failed gene-ID query makes every input symbol resolving to that ID
    # incomplete; record it as a query error rather than as "no ortholog".
    for row in resolutions.itertuples(index=False):
        if row.human_ensembl_gene_id in gene_id_errors:
            symbol_errors[row.match_symbol] = gene_id_errors[row.human_ensembl_gene_id]

    homology_empty = pd.DataFrame(
        columns=[
            "human_ensembl_gene_id",
            "human_gene_name",
            "dog_ensembl_gene_id",
            "dog_target_protein",
            "homology_type",
        ]
    )
    homologies = (
        pd.concat(homology_frames, ignore_index=True).drop_duplicates()
        if homology_frames
        else homology_empty
    )
    successful_resolutions = resolutions.loc[
        ~resolutions["match_symbol"].isin(symbol_errors)
    ]
    combined = successful_resolutions.merge(
        homologies,
        on="human_ensembl_gene_id",
        how="left",
        suffixes=("_resolved", "_homology"),
    )
    if "human_gene_name_resolved" in combined:
        combined["human_gene_name"] = combined["human_gene_name_resolved"].where(
            combined["human_gene_name_resolved"] != "",
            combined.get("human_gene_name_homology", ""),
        )
        combined = combined.drop(
            columns=["human_gene_name_resolved", "human_gene_name_homology"],
            errors="ignore",
        )
    for column in ["dog_ensembl_gene_id", "dog_target_protein", "homology_type"]:
        if column not in combined:
            combined[column] = ""
        combined[column] = combined[column].fillna("")
    return combined, symbol_errors


def is_usable_dog_symbol(value: str) -> bool:
    return bool(value.strip()) and value.strip().casefold() not in {
        "-",
        "na",
        "n/a",
        "none",
        "null",
    }


def build_qc(
    relationships: pd.DataFrame,
    mappings: pd.DataFrame,
    query_errors: dict[str, str],
) -> pd.DataFrame:
    groups = {symbol: group for symbol, group in mappings.groupby("match_symbol", sort=False)}
    output: list[dict[str, str]] = []

    for relation in relationships.itertuples(index=False):
        drug = relation.drug
        symbol = relation.human_target_protein
        base = {"drug": drug, "human_target_protein": symbol}
        if symbol in query_errors:
            output.append(
                {
                    **base,
                    "human_hgnc_id": "",
                    "human_ensembl_gene_id": "",
                    "dog_target_protein": "",
                    "dog_ensembl_gene_id": "",
                    "homology_type": "",
                    "mapping_status": "biomart_query_error",
                    "notes": query_errors[symbol],
                }
            )
            continue
        group = groups.get(symbol)
        if group is None or group.empty:
            note = "Input human target is blank." if not symbol else "No exact filter match returned by BioMart."
            output.append(
                {
                    **base,
                    "human_hgnc_id": "",
                    "human_ensembl_gene_id": "",
                    "dog_target_protein": "",
                    "dog_ensembl_gene_id": "",
                    "homology_type": "",
                    "mapping_status": "human_symbol_not_found",
                    "notes": note,
                }
            )
            continue

        human_ids = {value for value in group["human_ensembl_gene_id"] if value}
        hgnc_ids = {value for value in group["human_hgnc_id"] if value}
        # Ensembl may represent one HGNC gene on the primary chromosome and
        # several alternate assembly loci (for example GABBR1 in the MHC).
        # Those loci are one biological HGNC identity, not symbol ambiguity.
        ambiguous = len(hgnc_ids) > 1 or (not hgnc_ids and len(human_ids) != 1)
        has_one2one = bool((group["homology_type"] == STRICT_ONE2ONE_VALUE).any())
        strict_usable = group.loc[
            (group["homology_type"] == STRICT_ONE2ONE_VALUE)
            & group["dog_target_protein"].map(is_usable_dog_symbol)
        ]
        dog_identities = {
            row.dog_ensembl_gene_id or row.dog_target_protein
            for row in strict_usable.itertuples(index=False)
        }
        ambiguous_one2one = len(dog_identities) > 1
        for mapping in group.itertuples(index=False):
            dog_id = mapping.dog_ensembl_gene_id
            dog_symbol = mapping.dog_target_protein
            homology_type = mapping.homology_type
            notes = ""
            if ambiguous:
                status = "ambiguous_human_symbol"
                notes = (
                    "BioMart returned "
                    f"{len(hgnc_ids)} distinct HGNC IDs and {len(human_ids)} distinct "
                    "human Ensembl gene IDs for this exact symbol."
                )
            elif ambiguous_one2one:
                status = "ambiguous_one2one_mapping"
                notes = (
                    "Excluded because this HGNC identity produced multiple distinct strict "
                    "one-to-one dog orthologs across its Ensembl loci."
                )
            elif homology_type == STRICT_ONE2ONE_VALUE:
                if is_usable_dog_symbol(dog_symbol):
                    status = "mapped_one2one"
                    if not dog_id:
                        notes = "Strict one-to-one result has no dog Ensembl gene stable ID."
                else:
                    status = "dog_symbol_missing"
                    notes = "Strict one-to-one dog ortholog has no usable dog gene symbol."
            elif homology_type:
                status = "non_one2one_relationship" if has_one2one else "non_one2one_only"
                notes = "Excluded because BioMart does not classify this row as strict one-to-one."
            elif dog_id or is_usable_dog_symbol(dog_symbol):
                status = "dog_homology_type_missing"
                notes = "Dog information was returned without a homology-type classification."
            else:
                if has_one2one:
                    status = "additional_human_locus_without_dog_ortholog"
                    notes = (
                        "This Ensembl locus has no dog orthology row; another locus with the "
                        "same HGNC identity supplies the retained strict one-to-one mapping."
                    )
                else:
                    status = "no_dog_ortholog"
                    notes = (
                        "Human symbol was recognized, but no dog orthology information was returned."
                    )
            output.append(
                {
                    **base,
                    "human_hgnc_id": mapping.human_hgnc_id,
                    "human_ensembl_gene_id": mapping.human_ensembl_gene_id,
                    "dog_target_protein": dog_symbol,
                    "dog_ensembl_gene_id": dog_id,
                    "homology_type": homology_type,
                    "mapping_status": status,
                    "notes": notes,
                }
            )

    return pd.DataFrame(output, columns=QC_COLUMNS)


def build_final(qc: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    candidates = qc.loc[qc["mapping_status"] == "mapped_one2one", FINAL_COLUMNS].copy()
    before = len(candidates)
    final = candidates.drop_duplicates(FINAL_COLUMNS, keep="first", ignore_index=True)
    return final, before - len(final)


def build_per_drug_summary(
    relationships: pd.DataFrame, qc: pd.DataFrame, final: pd.DataFrame
) -> pd.DataFrame:
    mapped_pairs = set(
        qc.loc[qc["mapping_status"] == "mapped_one2one", ["drug", "human_target_protein"]]
        .drop_duplicates()
        .itertuples(index=False, name=None)
    )
    rows: list[dict[str, object]] = []
    for drug, group in relationships.groupby("drug", sort=False, dropna=False):
        human_targets = list(dict.fromkeys(group["human_target_protein"].tolist()))
        denominator = len(human_targets)
        mapped_human_count = sum((drug, symbol) in mapped_pairs for symbol in human_targets)
        mapped_dog_count = final.loc[
            final["drug"] == drug, "dog_target_protein"
        ].nunique()
        rows.append(
            {
                "drug": drug,
                "number_of_human_targets": denominator,
                "number_of_mapped_dog_targets": int(mapped_dog_count),
                "number_of_unmapped_human_targets": denominator - mapped_human_count,
                "mapping_percentage": round(
                    100.0 * mapped_human_count / denominator, 2
                )
                if denominator
                else 0.0,
            }
        )
    return pd.DataFrame(rows, columns=SUMMARY_COLUMNS)


def calculate_counts(
    total_rows: int,
    duplicate_input_rows: int,
    relationships: pd.DataFrame,
    mappings: pd.DataFrame,
    qc: pd.DataFrame,
    final: pd.DataFrame,
    duplicate_final_rows: int,
) -> dict[str, int]:
    unique_symbols = set(relationships["human_target_protein"])
    recognized = set(mappings.loc[mappings["human_ensembl_gene_id"] != "", "match_symbol"])
    any_dog = set(
        mappings.loc[
            (mappings["dog_ensembl_gene_id"] != "")
            | (mappings["dog_target_protein"] != "")
            | (mappings["homology_type"] != ""),
            "match_symbol",
        ]
    )
    strict = set(
        qc.loc[
            qc["homology_type"] == STRICT_ONE2ONE_VALUE, "human_target_protein"
        ]
    )
    usable = set(
        qc.loc[qc["mapping_status"] == "mapped_one2one", "human_target_protein"]
    )
    return {
        "total_input_rows": total_rows,
        "duplicate_input_rows_removed": duplicate_input_rows,
        "unique_drugs": int(relationships["drug"].nunique(dropna=False)),
        "unique_human_target_genes": len(unique_symbols),
        "human_symbols_successfully_recognized": len(unique_symbols & recognized),
        "human_targets_with_any_canine_orthology_information": len(unique_symbols & any_dog),
        "human_targets_with_strict_one2one_canine_orthologs": len(strict),
        "one2one_orthologs_with_usable_canine_gene_symbols": len(usable),
        "human_targets_not_retained": len(unique_symbols - usable),
        "final_unique_drug_dog_target_relationships": len(final),
        "duplicate_final_rows_removed": duplicate_final_rows,
    }


def print_summary(counts: dict[str, int]) -> None:
    total_unique = counts["unique_human_target_genes"]

    def percentage(value: int) -> str:
        return f"{(100.0 * value / total_unique):.2f}%" if total_unique else "0.00%"

    labels = [
        ("Total input rows", "total_input_rows", False),
        ("Duplicate input rows removed", "duplicate_input_rows_removed", False),
        ("Unique drugs", "unique_drugs", False),
        ("Unique human target genes", "unique_human_target_genes", False),
        ("Human symbols successfully recognized", "human_symbols_successfully_recognized", True),
        (
            "Human targets with any canine orthology information",
            "human_targets_with_any_canine_orthology_information",
            True,
        ),
        (
            "Human targets with strict one-to-one canine orthologs",
            "human_targets_with_strict_one2one_canine_orthologs",
            True,
        ),
        (
            "One-to-one orthologs with usable canine gene symbols",
            "one2one_orthologs_with_usable_canine_gene_symbols",
            True,
        ),
        ("Human targets not retained", "human_targets_not_retained", True),
        (
            "Final unique drug-dog target relationships",
            "final_unique_drug_dog_target_relationships",
            False,
        ),
        ("Duplicate final rows removed", "duplicate_final_rows_removed", False),
    ]
    print("\nMapping summary")
    print("=" * 72)
    for label, key, show_percentage in labels:
        value = counts[key]
        suffix = f" ({percentage(value)} of unique human targets)" if show_percentage else ""
        print(f"{label}: {value}{suffix}")


def write_metadata(
    path: Path,
    input_path: Path,
    config: BioMartConfig,
    counts: dict[str, int],
) -> None:
    attributes = [
        ("Human Gene stable ID", config.human_gene_id),
        ("Human Gene name", config.human_gene_name),
        ("Input-symbol audit attribute", config.match_symbol),
        ("Human HGNC identity", config.human_hgnc_id),
        ("Dog Gene stable ID", config.dog_gene_id),
        ("Dog Gene name", config.dog_gene_name),
        ("Dog homology type", config.homology_type),
    ]
    lines = [
        f"Date and time: {datetime.now().astimezone().isoformat(timespec='seconds')}",
        f"Python version: {platform.python_version()}",
        f"Script name/version: {Path(__file__).name} {SCRIPT_VERSION}",
        f"Ensembl BioMart requested service: {config.requested_service_url}",
        f"Ensembl BioMart resolved service: {config.resolved_service_url}",
        f"BioMart mart: {config.mart} ({config.mart_display_name})",
        f"Human dataset identifier: {config.dataset}",
        f"Human dataset display name: {config.dataset_display_name}",
        (
            "Human gene-symbol filter identifier: "
            f"{config.symbol_filter.identifier} ({config.symbol_filter.label})"
        ),
        (
            "Human Ensembl gene-ID filter identifier: "
            f"{config.human_gene_id_filter.identifier} ({config.human_gene_id_filter.label})"
        ),
        "Attributes used:",
    ]
    lines.extend(
        f"  - {purpose}: {field.identifier} ({field.label})"
        for purpose, field in attributes
    )
    lines.extend(
        [
            f"Source species: {SOURCE_SPECIES}",
            f"Target species: {TARGET_SPECIES}",
            "Orthology inclusion criterion: strict human-to-dog one-to-one orthologs only",
            f"Exact BioMart homology-type value used: {STRICT_ONE2ONE_VALUE}",
            f"Input filename: {input_path.name}",
            f"Input SHA-256: {sha256_file(input_path)}",
            f"Input row count: {counts['total_input_rows']}",
            f"Duplicate input rows removed: {counts['duplicate_input_rows_removed']}",
            f"Unique human target count: {counts['unique_human_target_genes']}",
            (
                "Per-drug mapping percentage definition: 100 * unique input human targets "
                "that produced at least one retained dog symbol / unique input human targets "
                "for that drug"
            ),
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_outputs(
    output_dir: Path,
    final: pd.DataFrame,
    qc: pd.DataFrame,
    per_drug: pd.DataFrame,
    input_path: Path,
    config: BioMartConfig,
    counts: dict[str, int],
) -> None:
    final.to_csv(output_dir / "drug_target_dog_one2one.csv", index=False)
    final.to_excel(output_dir / "drug_target_dog_one2one.xlsx", index=False, engine="openpyxl")
    qc.to_excel(output_dir / "drug_target_mapping_qc.xlsx", index=False, engine="openpyxl")
    per_drug.to_excel(
        output_dir / "drug_target_mapping_summary.xlsx", index=False, engine="openpyxl"
    )
    write_metadata(output_dir / "mapping_run_metadata.txt", input_path, config, counts)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prefix_chars="-–",
        description=(
            "Map human drug-target gene symbols to strict one-to-one canine orthologs "
            "using the current Ensembl Genes BioMart."
        )
    )
    parser.add_argument("--input", "–input", required=True, type=Path, help="Input .csv or .xlsx file")
    parser.add_argument(
        "--output-dir", "–output-dir", required=True, type=Path, help="Directory for result files"
    )
    parser.add_argument(
        "--biomart-url",
        default=DEFAULT_BIOMART_URL,
        help=f"BioMart martservice URL (default: {DEFAULT_BIOMART_URL})",
    )
    parser.add_argument("--batch-size", type=int, default=100, help="Unique symbols per query (default: 100)")
    parser.add_argument("--timeout", type=float, default=90.0, help="HTTP timeout in seconds (default: 90)")
    parser.add_argument("--retries", type=int, default=3, help="Retries after the first HTTP attempt (default: 3)")
    args = parser.parse_args(argv)
    if args.batch_size < 1:
        parser.error("--batch-size must be at least 1")
    if args.timeout <= 0:
        parser.error("--timeout must be greater than 0")
    if args.retries < 0:
        parser.error("--retries cannot be negative")
    return args


def run(args: argparse.Namespace) -> None:
    input_frame = read_input(args.input)
    total_rows = len(input_frame)
    relationships = input_frame.drop_duplicates(
        ["drug", "human_target_protein"], keep="first", ignore_index=True
    )
    duplicate_input_rows = total_rows - len(relationships)

    try:
        args.output_dir.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        raise OSError(f"Could not create output directory {args.output_dir}: {exc}") from exc
    if not args.output_dir.is_dir():
        raise NotADirectoryError(f"Output path is not a directory: {args.output_dir}")

    client = BioMartClient(args.biomart_url, args.timeout, args.retries)
    print(f"Validating current Ensembl BioMart configuration at {args.biomart_url} ...")
    config = client.validate_configuration()
    print(
        f"Validated {config.mart_display_name}; dataset {config.dataset} "
        f"({config.dataset_display_name})."
    )
    print(
        "Using filters/attributes: "
        + ", ".join(
            [config.symbol_filter.identifier, config.human_gene_id_filter.identifier]
            + config.resolution_attribute_ids
            + config.homology_attribute_ids
        )
    )

    unique_symbols = list(dict.fromkeys(relationships["human_target_protein"].tolist()))
    query_symbols = [symbol for symbol in unique_symbols if symbol]
    print(f"Querying {len(query_symbols)} unique nonblank human symbols in BioMart batches ...")
    mappings, query_errors = query_all_symbols(
        client, config, query_symbols, args.batch_size
    )
    qc = build_qc(relationships, mappings, query_errors)
    final, duplicate_final_rows = build_final(qc)
    per_drug = build_per_drug_summary(relationships, qc, final)
    counts = calculate_counts(
        total_rows,
        duplicate_input_rows,
        relationships,
        mappings,
        qc,
        final,
        duplicate_final_rows,
    )
    write_outputs(args.output_dir, final, qc, per_drug, args.input, config, counts)
    print_summary(counts)
    print(f"\nWrote results to: {args.output_dir}")


def main(argv: Sequence[str] | None = None) -> int:
    try:
        run(parse_args(argv))
        return 0
    except (BioMartError, FileNotFoundError, ValueError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
