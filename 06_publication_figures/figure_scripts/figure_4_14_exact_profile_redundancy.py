#!/usr/bin/env python3
"""Create Figure 4.14: recurrent exact drug-target profile distribution."""

from __future__ import annotations

import argparse
import collections
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from common import load_analyzed, save_figure


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-700", type=Path, required=True)
    parser.add_argument("--target-workbook", type=Path, required=True)
    parser.add_argument("--sheet", default="Sheet1")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()

    df07 = load_analyzed(args.results_700, {"drug"})
    mapping = pd.read_excel(args.target_workbook, sheet_name=args.sheet, usecols="A:C")
    targets: dict[str, set[str]] = collections.defaultdict(set)
    for drug, _gene, string_id in mapping.itertuples(index=False, name=None):
        if pd.notna(drug) and pd.notna(string_id):
            targets[str(drug).strip().lower()].add(str(string_id).strip())

    profiles: dict[tuple[str, ...], list[str]] = collections.defaultdict(list)
    for _, row in df07.iterrows():
        profile = tuple(sorted(targets[str(row["drug"]).strip().lower()]))
        profiles[profile].append(str(row["drug"]))
    counts = collections.Counter(len(values) for values in profiles.values() if len(values) >= 2)
    xs = sorted(counts)
    ys = [counts[x] for x in xs]

    fig, ax = plt.subplots(figsize=(7.8, 6.3))
    ax.set_axisbelow(True)
    ax.grid(axis="y", linewidth=0.45, alpha=0.25)
    ax.bar(xs, ys)
    ax.set_xlabel("Number of drugs sharing an exact target profile")
    ax.set_ylabel("Number of target profiles")
    ax.set_title("Distribution of Recurrent Exact Drug-Target Profiles", pad=10)
    ax.set_yscale("log")
    fig.tight_layout()
    pdf_path, _ = save_figure(fig, args.output_dir, "Figure_4_14_Exact_Profile_Redundancy")
    print(pdf_path)


if __name__ == "__main__":
    main()
