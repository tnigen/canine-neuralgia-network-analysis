#!/usr/bin/env python3
"""Reproducible canine drug–disease proximity analysis using NetMedPy."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib.metadata
import math
import os
import platform
import random
import re
import shlex
import sys
from collections import OrderedDict
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable, Optional


CANINE_STRING_ID = re.compile(r"^9615\.ENSCAFP\d+$")
REQUIRED_NETWORK_COLUMNS = ["protein1", "protein2", "combined_score"]
REQUIRED_DRUG_COLUMNS = ["drug", "dog_target_protein", "string_protein_id"]
PROXIMITY_COLUMNS = [
    "rank",
    "drug",
    "analysis_status",
    "n_targets_supplied",
    "n_targets_with_string_id",
    "n_targets_in_network",
    "n_targets_used",
    "observed_dc",
    "random_mean_dc",
    "random_sd_dc",
    "z_score",
    "empirical_p_value",
    "empirical_p_value_two_sided",
    "n_randomizations",
    "null_model",
    "random_seed",
]

PATH_ARGUMENTS = {"--network", "--disease", "--drug-targets", "--output"}


def parse_score_cutoff(value: str) -> Optional[int]:
    """Parse an integer STRING score cutoff or the literal ``none``."""
    text = value.strip().lower()
    if text == "none":
        return None
    try:
        cutoff = int(text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "score cutoff must be an integer from 0 to 1000 or 'none'"
        ) from exc
    if not 0 <= cutoff <= 1000:
        raise argparse.ArgumentTypeError(
            "score cutoff must be an integer from 0 to 1000 or 'none'"
        )
    return cutoff


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("value must be at least 1")
    return parsed


def seed_int(value: str) -> int:
    parsed = int(value)
    if not 0 <= parsed <= 4_294_967_295:
        raise argparse.ArgumentTypeError("seed must be between 0 and 4294967295")
    return parsed


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Canine STRING drug–disease closest-distance analysis using "
            "NetMedPy 0.2.0.1 for distances, null sampling, and statistics."
        )
    )
    parser.add_argument("--network", required=True, type=Path)
    parser.add_argument("--disease", required=True, type=Path)
    parser.add_argument("--drug-targets", required=True, type=Path)
    parser.add_argument(
        "--score-cutoff",
        required=True,
        type=parse_score_cutoff,
        metavar="INTEGER|none",
        help="retain combined_score >= INTEGER, or disable filtering with 'none'",
    )
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--seed", type=seed_int, default=42)
    parser.add_argument(
        "--n-randomizations",
        type=positive_int,
        default=1000,
        help="NetMedPy n_iter (default: NetMedPy's 1000)",
    )
    parser.add_argument(
        "--null-model",
        choices=["degree_match", "log_binning", "uniform"],
        default="degree_match",
        help="NetMedPy null model (default: degree_match)",
    )
    parser.add_argument(
        "--bin-size",
        type=positive_int,
        default=100,
        help="NetMedPy bin_size for log_binning (default: 100)",
    )
    parser.add_argument(
        "--n-processors",
        type=positive_int,
        default=None,
        help="workers for NetMedPy all_pair_distances (default: all CPUs)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="write input/QC outputs without computing proximity",
    )
    save_group = parser.add_mutually_exclusive_group()
    save_group.add_argument(
        "--save-distance-matrix",
        dest="save_distance_matrix",
        action="store_true",
        default=True,
        help="save NetMedPy's distance matrix (default)",
    )
    save_group.add_argument(
        "--no-save-distance-matrix",
        dest="save_distance_matrix",
        action="store_false",
        help="do not save the potentially large distance matrix",
    )
    return parser


def ensure_reproducible_hash_seed(seed: int) -> None:
    """Restart once so iteration over NetMedPy sets is reproducible."""
    expected = str(seed)
    if os.environ.get("PYTHONHASHSEED") == expected:
        return
    if os.environ.get("CANINE_PROXIMITY_HASH_REEXEC") == "1":
        raise RuntimeError("could not establish the requested PYTHONHASHSEED")
    env = os.environ.copy()
    env["PYTHONHASHSEED"] = expected
    env["CANINE_PROXIMITY_HASH_REEXEC"] = "1"
    script = str(Path(sys.argv[0]).resolve())
    os.execve(sys.executable, [sys.executable, script, *sys.argv[1:]], env)


def require_input(path: Path, label: str) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        raise FileNotFoundError(f"{label} file not found: {resolved}")
    return resolved


def prepare_output_directory(path: Path) -> Path:
    resolved = path.expanduser().resolve()
    if resolved.exists() and any(resolved.iterdir()):
        raise FileExistsError(
            f"output directory is not empty; choose a new directory: {resolved}"
        )
    resolved.mkdir(parents=True, exist_ok=True)
    return resolved


def open_text_auto(path: Path):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def public_path(path: Path) -> str:
    """Return a reproducible path without user-specific directories."""
    resolved = path.expanduser().resolve()
    try:
        return str(resolved.relative_to(Path.cwd().resolve()))
    except ValueError:
        return resolved.name


def public_command(arguments: list[str]) -> str:
    """Format a command line without absolute input or output paths."""
    if not arguments:
        return ""

    sanitized = [Path(arguments[0]).name]
    path_value_expected = False
    for argument in arguments[1:]:
        if path_value_expected:
            sanitized.append(public_path(Path(argument)))
            path_value_expected = False
            continue

        option, separator, value = argument.partition("=")
        if option in PATH_ARGUMENTS:
            sanitized.append(
                f"{option}={public_path(Path(value))}" if separator else argument
            )
            path_value_expected = not separator
        else:
            sanitized.append(argument)
    return shlex.join(sanitized)


def load_network(path: Path, cutoff: Optional[int]):
    """Adapt a STRING links file to an unweighted connected NetMedPy graph."""
    import networkx as nx

    graph = nx.Graph()
    qc: OrderedDict[str, Any] = OrderedDict(
        input_file=public_path(path),
        score_filtering_mode=(
            "disabled" if cutoff is None else "combined_score_greater_than_or_equal"
        ),
        score_cutoff=("none" if cutoff is None else cutoff),
        total_input_edge_records=0,
        malformed_records=0,
        records_passing_score_filter=0,
        records_filtered_out_by_score=0,
        self_loops=0,
        duplicate_edge_records=0,
    )
    malformed_examples: list[tuple[int, str, str]] = []

    with open_text_auto(path) as handle:
        header_line = handle.readline()
        if not header_line:
            raise ValueError("network file is empty")
        header = header_line.strip().split()
        if header != REQUIRED_NETWORK_COLUMNS:
            raise ValueError(
                f"network header must be {REQUIRED_NETWORK_COLUMNS}; observed {header}"
            )
        qc["confirmed_header"] = " ".join(header)

        for line_number, raw_line in enumerate(handle, start=2):
            qc["total_input_edge_records"] += 1
            parts = raw_line.split()
            reason = ""
            if len(parts) != 3:
                reason = "wrong_number_of_fields"
            else:
                protein1, protein2, score_text = parts
                try:
                    score = int(score_text)
                except ValueError:
                    reason = "non_integer_combined_score"
                else:
                    if not 0 <= score <= 1000:
                        reason = "combined_score_out_of_range"
                    elif not protein1 or not protein2:
                        reason = "missing_protein_identifier"
            if reason:
                qc["malformed_records"] += 1
                if len(malformed_examples) < 100:
                    malformed_examples.append(
                        (line_number, raw_line.rstrip("\r\n"), reason)
                    )
                continue

            if cutoff is not None and score < cutoff:
                qc["records_filtered_out_by_score"] += 1
                continue
            qc["records_passing_score_filter"] += 1

            if protein1 == protein2:
                qc["self_loops"] += 1
                continue
            if graph.has_edge(protein1, protein2):
                qc["duplicate_edge_records"] += 1
                continue
            # Omit weights to use unweighted shortest paths.
            graph.add_edge(protein1, protein2)

    if graph.number_of_edges() == 0:
        raise ValueError("no usable network edges remain after validation/filtering")

    components = list(nx.connected_components(graph))
    components.sort(key=len, reverse=True)
    component_sizes = [len(component) for component in components]
    selected_nodes = set(graph.nodes)
    analysis_nodes = set(components[0])
    qc["retained_unique_edges"] = graph.number_of_edges()
    qc["number_of_nodes"] = graph.number_of_nodes()
    qc["number_of_edges"] = graph.number_of_edges()
    qc["number_of_connected_components"] = len(components)
    qc["largest_connected_component_size"] = len(analysis_nodes)
    qc["nodes_outside_largest_connected_component"] = (
        graph.number_of_nodes() - len(analysis_nodes)
    )
    qc["analysis_background"] = "largest_connected_component_of_selected_interactome"

    if len(components) == 1:
        analysis_graph = graph
    else:
        analysis_graph = graph.subgraph(analysis_nodes).copy()
    qc["largest_connected_component_edges"] = analysis_graph.number_of_edges()

    return (
        analysis_graph,
        selected_nodes,
        analysis_nodes,
        qc,
        component_sizes,
        malformed_examples,
    )


def read_disease_set(
    path: Path, selected_nodes: set[str], analysis_nodes: set[str]
):
    seen: set[str] = set()
    duplicated_ids: set[str] = set()
    unique_valid: list[str] = []
    excluded: list[dict[str, Any]] = []
    raw_lines = blank_lines = nonempty_rows = duplicate_occurrences = 0
    invalid_identifier_rows = 0

    with path.open("r", encoding="utf-8", newline=None) as handle:
        for line_number, line in enumerate(handle, start=1):
            raw_lines += 1
            protein = line.strip()
            if not protein:
                blank_lines += 1
                continue
            nonempty_rows += 1
            if not CANINE_STRING_ID.fullmatch(protein):
                invalid_identifier_rows += 1
                excluded.append(
                    {
                        "line_number": line_number,
                        "protein_id": protein,
                        "reason": "invalid_identifier_format",
                    }
                )
                continue
            if protein in seen:
                duplicate_occurrences += 1
                duplicated_ids.add(protein)
                excluded.append(
                    {
                        "line_number": line_number,
                        "protein_id": protein,
                        "reason": "duplicate_identifier",
                    }
                )
                continue
            seen.add(protein)
            unique_valid.append(protein)

    mapped_selected = [p for p in unique_valid if p in selected_nodes]
    absent_selected = [p for p in unique_valid if p not in selected_nodes]
    outside_analysis_component = [
        p for p in mapped_selected if p not in analysis_nodes
    ]
    used = sorted(p for p in mapped_selected if p in analysis_nodes)

    excluded.extend(
        {"line_number": "", "protein_id": p, "reason": "absent_from_selected_interactome"}
        for p in absent_selected
    )
    excluded.extend(
        {
            "line_number": "",
            "protein_id": p,
            "reason": "outside_netmedpy_background_lcc",
        }
        for p in outside_analysis_component
    )

    qc: OrderedDict[str, Any] = OrderedDict(
        input_file=public_path(path),
        raw_lines=raw_lines,
        blank_lines_ignored=blank_lines,
        disease_protein_rows_supplied=nonempty_rows,
        invalid_identifier_rows=invalid_identifier_rows,
        unique_valid_disease_proteins=len(unique_valid),
        duplicated_identifiers=len(duplicated_ids),
        duplicate_occurrences_beyond_first=duplicate_occurrences,
        mapped_to_selected_interactome=len(mapped_selected),
        absent_from_selected_interactome=len(absent_selected),
        outside_netmedpy_background_lcc=len(outside_analysis_component),
        final_number_used_in_proximity=len(used),
        disease_set_definition="complete_GeneCards_union_DEG_input_not_disease_LCC",
    )
    if not used:
        raise ValueError("no disease proteins can participate in NetMedPy proximity")
    return used, qc, excluded


def normalize_cell(value: Any) -> str:
    if value is None:
        return ""
    return str(value).strip()


def read_drug_targets(
    path: Path, selected_nodes: set[str], analysis_nodes: set[str]
):
    import openpyxl

    workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    worksheet = workbook.active
    rows = worksheet.iter_rows(values_only=True)
    try:
        header_row = next(rows)
    except StopIteration as exc:
        raise ValueError("drug-target workbook is empty") from exc
    headers = [normalize_cell(value) for value in header_row]
    missing_columns = [name for name in REQUIRED_DRUG_COLUMNS if name not in headers]
    if missing_columns:
        raise ValueError(
            f"drug-target workbook is missing required columns: {missing_columns}; "
            f"observed {headers}"
        )
    indices = {name: headers.index(name) for name in REQUIRED_DRUG_COLUMNS}

    grouped: OrderedDict[str, list[dict[str, Any]]] = OrderedDict()
    excluded: list[dict[str, Any]] = []
    blank_rows_ignored = 0
    rows_without_drug = 0
    for workbook_row, row in enumerate(rows, start=2):
        drug = normalize_cell(row[indices["drug"]])
        dog_target = normalize_cell(row[indices["dog_target_protein"]])
        string_id = normalize_cell(row[indices["string_protein_id"]])
        if not drug and not dog_target and not string_id:
            blank_rows_ignored += 1
            continue
        record = {
            "workbook_row": workbook_row,
            "drug": drug,
            "dog_target_protein": dog_target,
            "string_protein_id": string_id,
        }
        if not drug:
            rows_without_drug += 1
            excluded.append({**record, "reason": "missing_drug_name"})
            continue
        grouped.setdefault(drug, []).append(record)
    workbook.close()

    target_sets: OrderedDict[str, set[str]] = OrderedDict()
    qc_rows: list[dict[str, Any]] = []
    for drug, records in grouped.items():
        seen: set[str] = set()
        with_id_rows = missing_id_rows = invalid_id_rows = duplicate_id_rows = 0
        absent_ids: set[str] = set()
        in_selected: set[str] = set()
        outside_component: set[str] = set()
        used: set[str] = set()

        for record in records:
            string_id = record["string_protein_id"]
            if not string_id:
                missing_id_rows += 1
                excluded.append({**record, "reason": "missing_string_protein_id"})
                continue
            with_id_rows += 1
            if not CANINE_STRING_ID.fullmatch(string_id):
                invalid_id_rows += 1
                excluded.append(
                    {**record, "reason": "invalid_string_protein_id_format"}
                )
                continue
            if string_id in seen:
                duplicate_id_rows += 1
                excluded.append(
                    {**record, "reason": "duplicate_string_protein_id_within_drug"}
                )
                continue
            seen.add(string_id)
            if string_id not in selected_nodes:
                absent_ids.add(string_id)
                excluded.append(
                    {**record, "reason": "absent_from_selected_interactome"}
                )
            elif string_id not in analysis_nodes:
                in_selected.add(string_id)
                outside_component.add(string_id)
                excluded.append(
                    {**record, "reason": "outside_netmedpy_background_lcc"}
                )
            else:
                in_selected.add(string_id)
                used.add(string_id)

        target_sets[drug] = used
        coverage = len(used) / len(records) if records else 0.0
        if not used:
            coverage_flag = "zero_valid_targets"
        elif coverage < 0.5:
            coverage_flag = "low_target_coverage_below_50_percent"
        else:
            coverage_flag = "adequate"
        qc_rows.append(
            {
                "drug": drug,
                "total_target_rows_supplied": len(records),
                "targets_with_nonempty_string_id": with_id_rows,
                "unique_string_ids_after_deduplication": len(seen),
                "targets_in_selected_interactome": len(in_selected),
                "targets_absent_from_selected_interactome": len(absent_ids),
                "targets_excluded_missing_string_id": missing_id_rows,
                "targets_excluded_invalid_string_id": invalid_id_rows,
                "duplicate_target_id_rows_excluded": duplicate_id_rows,
                "targets_excluded_outside_netmedpy_background_lcc": len(
                    outside_component
                ),
                "targets_used_in_proximity": len(used),
                "target_coverage_fraction": coverage,
                "coverage_flag": coverage_flag,
            }
        )

    workbook_qc = OrderedDict(
        input_file=public_path(path),
        worksheet=worksheet.title,
        confirmed_columns="|".join(headers),
        number_of_drugs=len(grouped),
        total_target_rows=sum(len(records) for records in grouped.values()),
        blank_rows_ignored=blank_rows_ignored,
        rows_excluded_missing_drug_name=rows_without_drug,
        rows_with_missing_string_protein_id=sum(
            row["targets_excluded_missing_string_id"] for row in qc_rows
        ),
        drugs_with_zero_valid_targets=sum(
            row["targets_used_in_proximity"] == 0 for row in qc_rows
        ),
    )
    return target_sets, qc_rows, excluded, workbook_qc


def write_key_value(path: Path, values: OrderedDict[str, Any] | dict[str, Any]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for key, value in values.items():
            handle.write(f"{key}={value}\n")


def write_tsv(path: Path, rows: Iterable[dict[str, Any]], columns: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_summary_tsv(path: Path, values: OrderedDict[str, Any]) -> None:
    write_tsv(
        path,
        ({"metric": key, "value": value} for key, value in values.items()),
        ["metric", "value"],
    )


def write_qc_outputs(
    output: Path,
    network_qc: OrderedDict[str, Any],
    component_sizes: list[int],
    malformed_examples: list[tuple[int, str, str]],
    disease_qc: OrderedDict[str, Any],
    disease_used: list[str],
    excluded_disease: list[dict[str, Any]],
    drug_qc_rows: list[dict[str, Any]],
    excluded_targets: list[dict[str, Any]],
    workbook_qc: OrderedDict[str, Any],
) -> None:
    network_text = OrderedDict(network_qc)
    network_text["duplicate_definition"] = "undirected_protein_pair_after_score_filter"
    network_text["edge_weighting_for_proximity"] = "unweighted"
    network_text["malformed_examples_saved"] = len(malformed_examples)
    write_key_value(output / "network_qc.txt", network_text)
    write_summary_tsv(output / "network_qc.tsv", network_text)
    write_tsv(
        output / "network_components.tsv",
        (
            {
                "component_rank": rank,
                "size": size,
                "used_by_netmedpy": rank == 1,
            }
            for rank, size in enumerate(component_sizes, start=1)
        ),
        ["component_rank", "size", "used_by_netmedpy"],
    )
    write_tsv(
        output / "malformed_network_record_examples.tsv",
        (
            {"line_number": line, "raw_record": raw, "reason": reason}
            for line, raw, reason in malformed_examples
        ),
        ["line_number", "raw_record", "reason"],
    )
    write_summary_tsv(output / "disease_set_qc.tsv", disease_qc)
    with (output / "disease_set_used.txt").open(
        "w", encoding="utf-8", newline="\n"
    ) as handle:
        for protein in disease_used:
            handle.write(f"{protein}\n")
    write_tsv(
        output / "excluded_disease_proteins.tsv",
        excluded_disease,
        ["line_number", "protein_id", "reason"],
    )
    drug_qc_columns = list(drug_qc_rows[0].keys()) if drug_qc_rows else ["drug"]
    write_tsv(output / "drug_target_qc.tsv", drug_qc_rows, drug_qc_columns)
    write_tsv(
        output / "excluded_drug_targets.tsv",
        excluded_targets,
        [
            "workbook_row",
            "drug",
            "dog_target_protein",
            "string_protein_id",
            "reason",
        ],
    )
    write_summary_tsv(output / "drug_target_workbook_qc.tsv", workbook_qc)


def print_validation(
    network_qc: OrderedDict[str, Any],
    disease_qc: OrderedDict[str, Any],
    drug_qc_rows: list[dict[str, Any]],
    workbook_qc: OrderedDict[str, Any],
) -> None:
    print("NETWORK")
    print(f"  input file: {network_qc['input_file']}")
    print(f"  score cutoff: {network_qc['score_cutoff']}")
    print(f"  total input edges: {network_qc['total_input_edge_records']}")
    print(f"  retained unique edges: {network_qc['retained_unique_edges']}")
    print(f"  nodes: {network_qc['number_of_nodes']}")
    print(f"  edges: {network_qc['number_of_edges']}")
    print(f"  components: {network_qc['number_of_connected_components']}")
    print(f"  LCC size: {network_qc['largest_connected_component_size']}")
    matrix_gib = (
        network_qc["largest_connected_component_size"] ** 2 * 4 / 1024**3
    )
    print(f"  estimated float32 distance matrix: {matrix_gib:.3f} GiB")
    print("DISEASE")
    print(f"  supplied disease proteins: {disease_qc['disease_protein_rows_supplied']}")
    print(f"  duplicate occurrences: {disease_qc['duplicate_occurrences_beyond_first']}")
    print(f"  mapped proteins: {disease_qc['mapped_to_selected_interactome']}")
    print(f"  missing proteins: {disease_qc['absent_from_selected_interactome']}")
    print(f"  final disease-set size: {disease_qc['final_number_used_in_proximity']}")
    print("DRUG TARGETS")
    print(f"  number of drugs: {workbook_qc['number_of_drugs']}")
    print(f"  total target rows: {workbook_qc['total_target_rows']}")
    print(
        "  rows with missing string_protein_id: "
        f"{workbook_qc['rows_with_missing_string_protein_id']}"
    )
    print(
        "  targets with STRING IDs: "
        f"{sum(row['targets_with_nonempty_string_id'] for row in drug_qc_rows)}"
    )
    print(
        "  targets present in selected network: "
        f"{sum(row['targets_in_selected_interactome'] for row in drug_qc_rows)}"
    )
    print(
        "  targets absent from selected network: "
        f"{sum(row['targets_absent_from_selected_interactome'] for row in drug_qc_rows)}"
    )
    print(
        "  drugs with zero valid targets: "
        f"{workbook_qc['drugs_with_zero_valid_targets']}"
    )


def make_parameters(
    args: argparse.Namespace,
    network: Path,
    disease: Path,
    drug_targets: Path,
    network_qc: OrderedDict[str, Any],
    disease_qc: OrderedDict[str, Any],
    workbook_qc: OrderedDict[str, Any],
) -> OrderedDict[str, Any]:
    import networkx
    import numpy
    import openpyxl
    import pandas

    cutoff = "none" if args.score_cutoff is None else args.score_cutoff
    return OrderedDict(
        analysis_status="validation_only" if args.dry_run else "proximity_pending",
        analysis_datetime_started=datetime.now().astimezone().isoformat(),
        command=public_command(sys.argv),
        network_file=public_path(network),
        network_sha256=sha256_file(network),
        disease_file=public_path(disease),
        disease_sha256=sha256_file(disease),
        drug_targets_file=public_path(drug_targets),
        drug_targets_sha256=sha256_file(drug_targets),
        string_version="12.0",
        string_species_taxon="9615_Canis_lupus_familiaris",
        score_cutoff=cutoff,
        score_filtering_mode=(
            "disabled_all_valid_edges_retained"
            if args.score_cutoff is None
            else "combined_score_greater_than_or_equal"
        ),
        proximity_method="dc_NetMedPy_directed_AMSPL",
        netmedpy_function="netmedpy.proximity",
        netmedpy_symmetric=False,
        distance_method="netmedpy.all_pair_distances_shortest_path",
        edge_weighting="unweighted",
        null_model=args.null_model,
        degree_aware_randomization=(args.null_model in {"degree_match", "log_binning"}),
        degree_matching=(args.null_model == "degree_match"),
        exact_degree_matching=(args.null_model == "degree_match"),
        bin_size=args.bin_size,
        number_of_randomizations=args.n_randomizations,
        random_standard_deviation="numpy.std_ddof_0",
        empirical_p_value="NetMedPy_lower_tail_count_R_le_observed_divided_by_n_iter",
        empirical_p_value_plus_one_correction=False,
        random_seed=args.seed,
        python_hash_seed=os.environ.get("PYTHONHASHSEED"),
        netmedpy_version=importlib.metadata.version("netmedpy"),
        python_version=platform.python_version(),
        networkx_version=networkx.__version__,
        numpy_version=numpy.__version__,
        pandas_version=pandas.__version__,
        openpyxl_version=openpyxl.__version__,
        selected_interactome_nodes=network_qc["number_of_nodes"],
        selected_interactome_edges=network_qc["number_of_edges"],
        selected_interactome_components=network_qc["number_of_connected_components"],
        netmedpy_background_component="largest_connected_component_of_selected_interactome",
        netmedpy_background_nodes=network_qc["largest_connected_component_size"],
        netmedpy_background_edges=network_qc["largest_connected_component_edges"],
        disease_set_definition="complete_GeneCards_union_DEG_input_not_disease_LCC",
        disease_proteins_supplied=disease_qc["disease_protein_rows_supplied"],
        disease_proteins_used=disease_qc["final_number_used_in_proximity"],
        number_of_drugs=workbook_qc["number_of_drugs"],
        n_processors=(
            "NetMedPy_default_all_CPUs"
            if args.n_processors is None
            else args.n_processors
        ),
        distance_matrix_saved=(args.save_distance_matrix and not args.dry_run),
    )


def blank_result_rows(
    drug_qc_rows: list[dict[str, Any]],
    args: argparse.Namespace,
    status: str,
) -> list[dict[str, Any]]:
    rows = []
    for qc in drug_qc_rows:
        rows.append(
            {
                "rank": "",
                "drug": qc["drug"],
                "analysis_status": status,
                "n_targets_supplied": qc["total_target_rows_supplied"],
                "n_targets_with_string_id": qc[
                    "targets_with_nonempty_string_id"
                ],
                "n_targets_in_network": qc["targets_in_selected_interactome"],
                "n_targets_used": qc["targets_used_in_proximity"],
                "observed_dc": "",
                "random_mean_dc": "",
                "random_sd_dc": "",
                "z_score": "",
                "empirical_p_value": "",
                "empirical_p_value_two_sided": "",
                "n_randomizations": args.n_randomizations,
                "null_model": args.null_model,
                "random_seed": args.seed,
            }
        )
    return rows


def run_netmedpy(
    graph,
    disease_used: list[str],
    target_sets: OrderedDict[str, set[str]],
    drug_qc_rows: list[dict[str, Any]],
    args: argparse.Namespace,
    output: Path,
) -> list[dict[str, Any]]:
    import netmedpy
    import numpy as np

    print("Computing NetMedPy all-pairs unweighted shortest-path distances...")
    distance_kwargs: dict[str, Any] = {
        "distance": "shortest_path",
        "matrix_dtype": np.float32,
    }
    if args.n_processors is not None:
        distance_kwargs["n_processors"] = args.n_processors
        distance_kwargs["n_tasks"] = args.n_processors
    distance_matrix = netmedpy.all_pair_distances(graph, **distance_kwargs)
    if args.save_distance_matrix:
        matrix_path = output / "netmedpy_shortest_path_distances.npz"
        print(f"Saving NetMedPy distance matrix: {matrix_path}")
        netmedpy.save_distances(distance_matrix, str(matrix_path))

    random.seed(args.seed)
    np.random.seed(args.seed)
    qc_by_drug = {row["drug"]: row for row in drug_qc_rows}
    result_rows: list[dict[str, Any]] = []
    null_distributions: list[tuple[str, int, float]] = []

    for drug in sorted(target_sets):
        targets = sorted(target_sets[drug])
        qc = qc_by_drug[drug]
        base = {
            "rank": "",
            "drug": drug,
            "n_targets_supplied": qc["total_target_rows_supplied"],
            "n_targets_with_string_id": qc["targets_with_nonempty_string_id"],
            "n_targets_in_network": qc["targets_in_selected_interactome"],
            "n_targets_used": len(targets),
            "n_randomizations": args.n_randomizations,
            "null_model": args.null_model,
            "random_seed": args.seed,
        }
        if not targets:
            result_rows.append(
                {
                    **base,
                    "analysis_status": "skipped_zero_valid_targets",
                    "observed_dc": "",
                    "random_mean_dc": "",
                    "random_sd_dc": "",
                    "z_score": "",
                    "empirical_p_value": "",
                    "empirical_p_value_two_sided": "",
                }
            )
            continue

        print(
            f"Running NetMedPy proximity for {drug} "
            f"({len(targets)} targets; {args.n_randomizations} randomizations)..."
        )
        scores = netmedpy.proximity(
            graph,
            targets,
            disease_used,
            distance_matrix,
            null_model=args.null_model,
            n_iter=args.n_randomizations,
            bin_size=args.bin_size,
            symmetric=False,
        )
        for iteration, random_dc in enumerate(scores["dist"], start=1):
            null_distributions.append((drug, iteration, float(random_dc)))
        result_rows.append(
            {
                **base,
                "analysis_status": "analyzed",
                "observed_dc": float(scores["raw_amspl"]),
                "random_mean_dc": float(scores["d_mu"]),
                "random_sd_dc": float(scores["d_sigma"]),
                "z_score": float(scores["z_score"]),
                "empirical_p_value": float(scores["p_value_single_tail"]),
                "empirical_p_value_two_sided": float(
                    scores["p_value_double_tail"]
                ),
            }
        )

    def rank_key(row: dict[str, Any]):
        z = row["z_score"]
        p = row["empirical_p_value"]
        dc = row["observed_dc"]
        z_key = z if isinstance(z, float) and math.isfinite(z) else math.inf
        p_key = p if isinstance(p, float) and math.isfinite(p) else math.inf
        dc_key = dc if isinstance(dc, float) and math.isfinite(dc) else math.inf
        return (row["analysis_status"] != "analyzed", z_key, p_key, dc_key, row["drug"])

    result_rows.sort(key=rank_key)
    rank = 0
    for row in result_rows:
        if row["analysis_status"] == "analyzed":
            rank += 1
            row["rank"] = rank

    with gzip.open(
        output / "netmedpy_null_distributions.tsv.gz",
        "wt",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["drug", "iteration", "random_dc"])
        writer.writerows(null_distributions)
    return result_rows


def main(argv: Optional[list[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if argv is None:
        ensure_reproducible_hash_seed(args.seed)

    network = require_input(args.network, "network")
    disease = require_input(args.disease, "disease")
    drug_targets = require_input(args.drug_targets, "drug-target")
    output = prepare_output_directory(args.output)

    (
        graph,
        selected_nodes,
        analysis_nodes,
        network_qc,
        component_sizes,
        malformed_examples,
    ) = load_network(network, args.score_cutoff)
    disease_used, disease_qc, excluded_disease = read_disease_set(
        disease, selected_nodes, analysis_nodes
    )
    target_sets, drug_qc_rows, excluded_targets, workbook_qc = read_drug_targets(
        drug_targets, selected_nodes, analysis_nodes
    )
    write_qc_outputs(
        output,
        network_qc,
        component_sizes,
        malformed_examples,
        disease_qc,
        disease_used,
        excluded_disease,
        drug_qc_rows,
        excluded_targets,
        workbook_qc,
    )
    print_validation(network_qc, disease_qc, drug_qc_rows, workbook_qc)

    parameters = make_parameters(
        args,
        network,
        disease,
        drug_targets,
        network_qc,
        disease_qc,
        workbook_qc,
    )
    if args.dry_run:
        write_tsv(
            output / "proximity_results.tsv",
            blank_result_rows(drug_qc_rows, args, "not_run_dry_run"),
            PROXIMITY_COLUMNS,
        )
        parameters["analysis_datetime_completed"] = datetime.now().astimezone().isoformat()
        write_key_value(output / "analysis_parameters.txt", parameters)
        print(f"Dry run complete. QC outputs: {output}")
        return 0

    write_key_value(output / "analysis_parameters.txt", parameters)
    result_rows = run_netmedpy(
        graph, disease_used, target_sets, drug_qc_rows, args, output
    )
    write_tsv(output / "proximity_results.tsv", result_rows, PROXIMITY_COLUMNS)
    parameters["analysis_status"] = "complete"
    parameters["analysis_datetime_completed"] = datetime.now().astimezone().isoformat()
    write_key_value(output / "analysis_parameters.txt", parameters)
    print(f"Analysis complete. Results: {output / 'proximity_results.tsv'}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, FileExistsError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
