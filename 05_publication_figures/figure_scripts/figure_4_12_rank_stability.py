#!/usr/bin/env python3
"""Create Figure 4.12: rank stability across STRING thresholds."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

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
    common = df07[["drug", "rank"]].merge(
        df09[["drug", "rank"]], on="drug", suffixes=("_07", "_09")
    )
    rank07 = common["rank_07"].to_numpy()
    rank09 = common["rank_09"].to_numpy()
    rho = pd.Series(rank07).rank().corr(pd.Series(rank09).rank())

    fig, ax = plt.subplots(figsize=(7.8, 7.0))
    ax.scatter(rank07, rank09, s=12, alpha=0.45)
    limit = max(rank07.max(), rank09.max())
    ax.plot([1, limit], [1, limit], linestyle="--", linewidth=1.4)
    ax.set_xlabel("Rank at STRING ≥0.7")
    ax.set_ylabel("Rank at STRING ≥0.9")
    ax.set_title("Rank Stability Across STRING Confidence Thresholds", pad=10)
    ax.text(
        0.04,
        0.94,
        f"n = {len(common):,}\nSpearman ρ = {rho:.3f}",
        transform=ax.transAxes,
        va="top",
        fontsize=13,
    )
    ax.grid(True, linewidth=0.45, alpha=0.25)
    fig.tight_layout()
    pdf_path, _ = save_figure(fig, args.output_dir, "Figure_4_12_Rank_Stability")
    print(pdf_path)


if __name__ == "__main__":
    main()
