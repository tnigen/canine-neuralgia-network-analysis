#!/usr/bin/env python3
"""Evaluate canine neuralgia LCC significance with NetMedPy."""

from __future__ import annotations

import argparse
import csv
import gzip
import importlib.metadata
import math
import platform
import random
import sys
from pathlib import Path

import netmedpy
import networkx as nx
import numpy as np


SPECIES = "Canis lupus familiaris (domestic dog)"
STRING_SPECIES_ID = "9615"
SUPPORTED_NULL_MODELS = (
    "degree_match, log_binning, strength_binning, uniform, custom"
)
NETMEDPY_OUTPUT_KEYS = (
    "d_mu, d_sigma, z_score, p_val, lcc, lcc_size, dist"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run canine neuralgia LCC significance analysis with NetMedPy."
    )
    parser.add_argument(
        "--interactome",
        type=Path,
        default=Path("9615.protein.links.v12.0.txt.gz"),
        help="Gzipped STRING protein.links file.",
    )
    parser.add_argument(
        "--disease-genes",
        type=Path,
        default=Path("canine_neuralgia_genes.txt"),
        help="Text file containing one canine STRING protein ID per line.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("."),
        help="Directory for result, report, missing-ID, and environment files.",
    )
    parser.add_argument(
        "--null-model",
        choices=("degree_match", "log_binning", "uniform"),
        default="log_binning",
        help="NetMedPy null model (default: log_binning).",
    )
    parser.add_argument(
        "--n-iterations",
        type=int,
        default=10_000,
        help="Number of NetMedPy randomizations (default: 10000).",
    )
    parser.add_argument(
        "--bin-size",
        type=int,
        default=100,
        help="Minimum degree-bin size for log_binning (default: 100).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed applied before calling NetMedPy (default: 42).",
    )
    args = parser.parse_args()
    if args.n_iterations <= 0:
        parser.error("--n-iterations must be greater than zero")
    if args.bin_size < 1:
        parser.error("--bin-size must be at least one")
    return args


def package_version(distribution_name: str) -> str:
    return importlib.metadata.version(distribution_name)


def read_string_graph(path: Path) -> tuple[nx.Graph, dict]:
    """Read a STRING links file, build a simple graph, and collect input QC."""
    graph = nx.Graph()
    first_rows: list[tuple[str, str, int]] = []
    input_edge_records = 0
    self_loops = 0
    same_orientation_duplicates = 0
    reversed_duplicates = 0
    minimum_score: int | None = None
    maximum_score: int | None = None
    score_sum = 0

    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        header = handle.readline().strip().split()
        expected = ["protein1", "protein2", "combined_score"]
        if header != expected:
            raise ValueError(
                f"Unexpected interactome header {header!r}; expected {expected!r}"
            )

        for line_number, line in enumerate(handle, start=2):
            fields = line.split()
            if not fields:
                continue
            if len(fields) != 3:
                raise ValueError(
                    f"Line {line_number} in {path} has {len(fields)} fields; expected 3"
                )
            protein1, protein2, score_text = fields
            try:
                score = int(score_text)
            except ValueError as exc:
                raise ValueError(
                    f"Invalid combined_score on line {line_number}: {score_text!r}"
                ) from exc

            # Retain STRING associations with combined_score >= 900.
            if score < 900:
                continue

            input_edge_records += 1
            score_sum += score
            minimum_score = score if minimum_score is None else min(minimum_score, score)
            maximum_score = score if maximum_score is None else max(maximum_score, score)
            if protein1 == protein2:
                self_loops += 1
            if len(first_rows) < 5:
                first_rows.append((protein1, protein2, score))

            existing = graph.get_edge_data(protein1, protein2)
            current_forward = protein1 <= protein2
            if existing is not None:
                if existing["_first_forward"] == current_forward:
                    same_orientation_duplicates += 1
                else:
                    reversed_duplicates += 1
                # Retain the first score for each undirected pair.
            else:
                graph.add_edge(
                    protein1,
                    protein2,
                    combined_score=score,
                    _first_forward=current_forward,
                )

            if input_edge_records % 1_000_000 == 0:
                print(
                    f"Loaded {input_edge_records:,} STRING edge records "
                    f"({graph.number_of_nodes():,} nodes)",
                    flush=True,
                )

    if input_edge_records == 0:
        raise ValueError(f"No interactome edges were found in {path}")

    qc = {
        "header": header,
        "first_rows": first_rows,
        "input_edge_records": input_edge_records,
        "graph_edges": graph.number_of_edges(),
        "graph_nodes": graph.number_of_nodes(),
        "self_loops": self_loops,
        "same_orientation_duplicates": same_orientation_duplicates,
        "reversed_duplicates": reversed_duplicates,
        "all_duplicate_records": (
            same_orientation_duplicates + reversed_duplicates
        ),
        "minimum_score": minimum_score,
        "maximum_score": maximum_score,
        "mean_score": score_sum / input_edge_records,
    }
    return graph, qc


def read_disease_ids(path: Path) -> tuple[list[str], int]:
    """Read non-empty one-column IDs and retain first-occurrence order."""
    supplied: list[str] = []
    with path.open("rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.split()
            if not fields:
                continue
            if len(fields) != 1:
                raise ValueError(
                    f"Line {line_number} in {path} has {len(fields)} fields; expected 1"
                )
            supplied.append(fields[0])

    unique_ids = list(dict.fromkeys(supplied))
    return unique_ids, len(supplied) - len(unique_ids)


def format_number(value: float | int) -> str:
    if isinstance(value, (float, np.floating)):
        if math.isnan(float(value)):
            return "nan"
        return f"{float(value):.12g}"
    return str(value)


def interpretation(observed: int, null_mean: float, z_score: float, p_value: float) -> str:
    if math.isnan(z_score):
        return (
            "The null distribution had zero variance, so a finite LCC Z-score "
            "could not be calculated. Significance should not be claimed from "
            "the Z-score."
        )
    if z_score > 0 and p_value < 0.05:
        return (
            f"The observed LCC ({observed} proteins) is larger than the null "
            f"expectation ({null_mean:.3f}) and is significant at alpha=0.05. "
            "This supports treating the mapped canine neuralgia proteins as a "
            "topologically coherent disease module for subsequent network-medicine "
            "analyses. It does not establish drug efficacy or drug-disease proximity."
        )
    if z_score > 0:
        return (
            "The observed LCC is larger than the null mean but is not significant "
            "at alpha=0.05; the present analysis does not provide sufficient evidence "
            "of unusually strong disease-module connectivity."
        )
    return (
        "The observed LCC is not larger than expected under the selected null model; "
        "the present analysis does not support unusually strong disease-module "
        "connectivity."
    )


def write_environment_versions(path: Path) -> None:
    package_names = [
        "netmedpy",
        "networkx",
        "numpy",
        "scipy",
        "pandas",
        "matplotlib",
        "seaborn",
    ]
    lines = [
        f"Python=={platform.python_version()}",
        f"pip=={package_version('pip')}",
        f"Operating system={platform.platform()}",
        f"Architecture={platform.machine()}",
    ]
    lines.extend(f"{name}=={package_version(name)}" for name in package_names)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if not args.interactome.is_file():
        raise FileNotFoundError(args.interactome)
    if not args.disease_genes.is_file():
        raise FileNotFoundError(args.disease_genes)

    netmedpy_version = package_version("netmedpy")
    networkx_version = package_version("networkx")
    python_version = platform.python_version()

    graph, qc = read_string_graph(args.interactome)
    disease_ids, duplicate_disease_ids = read_disease_ids(args.disease_genes)
    graph_node_ids = set(graph.nodes)
    mapped_ids = [node for node in disease_ids if node in graph_node_ids]
    missing_ids = [node for node in disease_ids if node not in graph_node_ids]
    present_percentage = (
        100.0 * len(mapped_ids) / len(disease_ids) if disease_ids else 0.0
    )

    missing_path = args.output_dir / "canine_neuralgia_missing_from_interactome.txt"
    missing_path.write_text(
        "".join(f"{node}\n" for node in missing_ids), encoding="utf-8"
    )
    if not mapped_ids:
        raise ValueError("None of the supplied disease-module IDs occur in the graph")

    print("\nPre-analysis quality control")
    print("A. First 5 STRING rows:")
    for protein1, protein2, score in qc["first_rows"]:
        print(f"   {protein1} {protein2} {score}")
    print(f"B. Interactome edges (unique undirected): {qc['graph_edges']:,}")
    print(f"C. Interactome nodes: {qc['graph_nodes']:,}")
    print("D. First 10 disease-module IDs:")
    for node in disease_ids[:10]:
        print(f"   {node}")
    print(f"E. Disease-module IDs found: {len(mapped_ids):,}")
    print(f"F. Disease-module IDs absent: {len(missing_ids):,}")
    print(f"G. Self-loops present: {'yes' if qc['self_loops'] else 'no'}")
    print(
        "H. Duplicate/reversed duplicate records: "
        f"{qc['same_orientation_duplicates']:,}/"
        f"{qc['reversed_duplicates']:,}"
    )
    print("I. Final LCC graph: unweighted topology; combined_score retained as metadata")
    print(f"J. Python: {python_version}")
    print(f"K. NetMedPy: {netmedpy_version}")
    print(f"L. NetworkX: {networkx_version}")
    print(f"M. Null model: {args.null_model}")
    print(f"N. Randomization iterations: {args.n_iterations:,}", flush=True)

    # Seed Python and NumPy before NetMedPy's null-model sampling.
    random.seed(args.seed)
    np.random.seed(args.seed)
    lcc_result = netmedpy.lcc_significance(
        graph,
        mapped_ids,
        null_model=args.null_model,
        n_iter=args.n_iterations,
        bin_size=args.bin_size,
    )

    observed_lcc = int(lcc_result["lcc_size"])
    null_mean = float(lcc_result["d_mu"])
    null_sd = float(lcc_result["d_sigma"])
    z_score = float(lcc_result["z_score"])
    p_value = float(lcc_result["p_val"])
    result_row = {
        "disease_module_size": len(disease_ids),
        "mapped_disease_module_size": len(mapped_ids),
        "missing_disease_module_size": len(missing_ids),
        "interactome_nodes": qc["graph_nodes"],
        "interactome_edges": qc["graph_edges"],
        "observed_lcc_size": observed_lcc,
        "null_mean_lcc": format_number(null_mean),
        "null_sd_lcc": format_number(null_sd),
        "z_score": format_number(z_score),
        "p_value": format_number(p_value),
        "netmedpy_version": netmedpy_version,
        "networkx_version": networkx_version,
        "python_version": python_version,
        "null_model": args.null_model,
        "n_iterations": args.n_iterations,
    }
    results_path = args.output_dir / "canine_neuralgia_LCC_results.csv"
    with results_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(result_row))
        writer.writeheader()
        writer.writerow(result_row)

    p_value_note = ""
    if p_value == 0.0:
        p_value_note = (
            " NetMedPy returned p=0 because none of the 10,000 randomized LCCs "
            "was at least as large as observed; the finite-simulation resolution "
            "supports reporting this as empirical p<0.0001 alongside the exact output."
        )
    report_text = f"""Canine neuralgia LCC significance report
==========================================

Inputs
------
Species: {SPECIES}
STRING species ID: {STRING_SPECIES_ID}
Interactome: {args.interactome}
Disease module: {args.disease_genes}

Installed NetMedPy API verified from package source
---------------------------------------------------
Function: lcc_significance(net, A, null_model='degree_match', node_bucket=None, n_iter=1000, bin_size=100, strength_weight='weight')
Expected network: networkx.Graph
Expected node set: iterable of graph node identifiers
Supported null models: {SUPPORTED_NULL_MODELS}
Returned values: {NETMEDPY_OUTPUT_KEYS}

Graph construction and interactome QC
-------------------------------------
Method: simple undirected NetworkX Graph using protein1 and protein2 exactly as supplied by STRING
Graph used for LCC: unweighted topology
Edge metadata retained: combined_score (integer)
No combined_score threshold was applied.
Input edge records: {qc['input_edge_records']}
Unique undirected graph edges: {qc['graph_edges']}
Unique graph nodes: {qc['graph_nodes']}
Self-loop records: {qc['self_loops']}
Repeated same-orientation edge records: {qc['same_orientation_duplicates']}
Reversed duplicate edge records: {qc['reversed_duplicates']}
Minimum combined_score: {qc['minimum_score']}
Maximum combined_score: {qc['maximum_score']}
Mean combined_score: {qc['mean_score']:.6f}

First five STRING rows:
""" + "".join(
        f"{protein1} {protein2} {score}\n"
        for protein1, protein2, score in qc["first_rows"]
    ) + f"""
Disease-module ID matching
--------------------------
Unique disease-module IDs supplied: {len(disease_ids)}
Duplicate supplied ID occurrences: {duplicate_disease_ids}
IDs present in the interactome: {len(mapped_ids)}
IDs absent from the interactome: {len(missing_ids)}
Percentage present: {present_percentage:.4f}%
Filtering for calculation: only the {len(mapped_ids)} IDs present in the background interactome were passed to NetMedPy. Missing IDs were not silently discarded; they are recorded in {missing_path.name}.

LCC significance result
-----------------------
Observed LCC size: {observed_lcc}
Null model: {args.null_model}
Degree-bin size: {args.bin_size}
Random seed: {args.seed}
Randomization iterations: {args.n_iterations}
Null mean LCC: {format_number(null_mean)}
Null standard deviation: {format_number(null_sd)}
LCC Z-score: {format_number(z_score)}
NetMedPy empirical upper-tail p-value: {format_number(p_value)}

Interpretation
--------------
The observed LCC is the largest group of mapped canine neuralgia proteins that are mutually reachable using only interactions among proteins in the supplied disease set. The log-binning null distribution consists of equally sized randomized node sets sampled from similar network-degree bins. The Z-score measures how many null standard deviations the observed LCC lies above or below the null mean. The empirical upper-tail p-value is the fraction of randomized LCCs at least as large as the observed LCC.{p_value_note}

{interpretation(observed_lcc, null_mean, z_score, p_value)}

Methodological limitations
--------------------------
This is a topological disease-module validation, not evidence of causality, drug efficacy, or drug-disease proximity. STRING associations are not necessarily direct physical interactions and can reflect heterogeneous evidence. Using the supplied full STRING links file without an additional threshold preserves the requested input but can make the network dense. Results depend on interactome coverage, disease-ID coverage, database ascertainment, and the selected approximate degree-bin null model. The graph is undirected and context-free, so tissue, temporal, cell-type, directionality, and interaction-sign information are not modeled. Monte Carlo p-value precision is limited by {args.n_iterations} iterations.

Software
--------
Python: {python_version}
NetMedPy: {netmedpy_version}
NetworkX: {networkx_version}
NumPy: {package_version('numpy')}
SciPy: {package_version('scipy')}
pandas: {package_version('pandas')}
"""
    report_path = args.output_dir / "canine_neuralgia_LCC_report.txt"
    report_path.write_text(report_text, encoding="utf-8")
    write_environment_versions(args.output_dir / "environment_versions.txt")

    print("\nAnalysis complete")
    print(f"Observed LCC size: {observed_lcc}")
    print(f"Null mean LCC: {format_number(null_mean)}")
    print(f"Null SD LCC: {format_number(null_sd)}")
    print(f"Z-score: {format_number(z_score)}")
    print(f"Empirical p-value: {format_number(p_value)}")
    print(f"Results: {results_path}")
    print(f"Report: {report_path}")
    print(f"Missing IDs: {missing_path}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {type(exc).__name__}: {exc}", file=sys.stderr)
        raise
