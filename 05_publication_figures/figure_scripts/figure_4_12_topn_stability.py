#!/usr/bin/env python3
"""Create Figure 4.12: Top-N stability across STRING thresholds."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt

from common import load_analyzed, save_figure


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-700", type=Path, required=True)
    parser.add_argument("--results-900", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()

    required = {"drug", "rank"}
    df07 = load_analyzed(args.results_700, required)
    df09 = load_analyzed(args.results_900, required)
    top_n_values = [50, 100, 154, 230, 307, 500]
    retention = []
    jaccard = []
    for n in top_n_values:
        top07 = set(df07.loc[df07["rank"] <= n, "drug"])
        top09 = set(df09.loc[df09["rank"] <= n, "drug"])
        shared = len(top07 & top09)
        retention.append(shared / n)
        jaccard.append(shared / len(top07 | top09))

    fig, ax = plt.subplots(figsize=(8.2, 6.2))
    ax.plot(top_n_values, retention, marker="o", linewidth=2.0, label="Retention")
    ax.plot(top_n_values, jaccard, marker="s", linewidth=2.0, label="Jaccard index")
    ax.set_xlabel("Top-N threshold")
    ax.set_ylabel("Proportion")
    ax.set_ylim(0, 1)
    ax.set_title("Top-N Stability Across STRING Confidence Thresholds", pad=10)
    ax.legend(fontsize=12)
    ax.grid(True, linewidth=0.45, alpha=0.25)
    fig.tight_layout()
    pdf_path, _ = save_figure(fig, args.output_dir, "Figure_4_12_TopN_Stability")
    print(pdf_path)


if __name__ == "__main__":
    main()
