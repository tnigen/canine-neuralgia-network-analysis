#!/usr/bin/env python3
"""Create the network-proximity ranking figure with the top 20 highlighted."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from common import load_analyzed, save_figure


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-700", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()

    df07 = load_analyzed(args.results_700, {"drug", "z_score"})
    ranked = df07.sort_values(["z_score", "drug"], ascending=[True, True]).reset_index(drop=True)
    ranked["plot_rank"] = np.arange(1, len(ranked) + 1)
    top20 = ranked.head(20)

    fig, ax = plt.subplots(figsize=(14.2, 8.6))
    ax.plot(ranked["plot_rank"], ranked["z_score"], linewidth=1.8, label="_nolegend_")
    ax.axhline(0, linestyle="--", linewidth=1.2, alpha=0.6)
    for _, row in top20.iterrows():
        ax.scatter(
            row["plot_rank"],
            row["z_score"],
            s=58,
            zorder=3,
            label=f"{str(row['drug']).capitalize()} (rank {int(row['plot_rank'])})",
        )
    ax.set_title(f"Network Proximity Ranking of {len(ranked):,} Analyzable Drugs", pad=12)
    ax.set_xlabel("Network proximity rank")
    ax.set_ylabel("Proximity Z-score")
    ax.legend(
        title="Top 20 drugs",
        bbox_to_anchor=(1.02, 0.5),
        loc="center left",
        frameon=True,
        fontsize=11,
        title_fontsize=13,
        labelspacing=0.55,
    )
    fig.tight_layout()
    pdf_path, _ = save_figure(fig, args.output_dir, "network_proximity_ranking_top20")
    print(pdf_path)


if __name__ == "__main__":
    main()
