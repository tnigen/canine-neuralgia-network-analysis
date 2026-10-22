#!/usr/bin/env python3
"""Create Figure 4.14: drug-target count versus proximity Z-score."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from common import load_analyzed, save_figure


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-700", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()

    df07 = load_analyzed(args.results_700, {"n_targets_used", "z_score"})
    n_targets = df07["n_targets_used"].astype(int).to_numpy()
    z_scores = df07["z_score"].astype(float).to_numpy()
    rho = pd.Series(n_targets).rank().corr(pd.Series(z_scores).rank())

    fig, ax = plt.subplots(figsize=(7.8, 6.3))
    ax.set_axisbelow(True)
    ax.grid(True, linewidth=0.45, alpha=0.25)
    ax.scatter(n_targets, z_scores, s=15, alpha=0.35)
    ax.set_xlabel("Number of drug targets used in proximity analysis")
    ax.set_ylabel("Network-proximity Z-score")
    ax.set_title("Association Between Drug Target Count and Network Proximity", pad=10)
    ax.text(
        0.97,
        0.05,
        f"Spearman ρ = {rho:.3f}\nn = {len(df07):,}",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=13,
    )
    fig.tight_layout()
    pdf_path, _ = save_figure(fig, args.output_dir, "Figure_4_14_Target_Count_vs_Z")
    print(pdf_path)


if __name__ == "__main__":
    main()
