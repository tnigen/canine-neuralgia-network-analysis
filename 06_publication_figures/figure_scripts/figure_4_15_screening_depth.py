#!/usr/bin/env python3
"""Create Figure 4.15: clinical-reference recovery by screening depth."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from common import CT_POSITIVE, load_analyzed, save_figure


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-700", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()

    df07 = load_analyzed(args.results_700, {"drug", "rank"})
    ranked = df07.sort_values("rank").reset_index(drop=True)
    labels = ranked["drug"].str.strip().str.lower().isin(CT_POSITIVE).to_numpy()
    total_positive = int(labels.sum())
    if total_positive == 0:
        raise ValueError("No clinical-reference drugs were found in the analyzed rows.")

    depths = np.array([10, 15, 20])
    ns = np.floor(len(ranked) * depths / 100 + 0.5).astype(int)
    recovered = np.array([int(labels[:n].sum()) for n in ns])
    percentages = recovered / total_positive * 100

    fig, ax = plt.subplots(figsize=(7.6, 6.1))
    ax.set_axisbelow(True)
    ax.grid(True, linewidth=0.45, alpha=0.25)
    ax.plot(depths, percentages, marker="o", linewidth=2.0)
    for x, y, n in zip(depths, percentages, ns):
        ax.annotate(
            f"n={n}",
            (x, y),
            xytext=(0, 10),
            textcoords="offset points",
            ha="center",
            fontsize=12,
        )
    ax.set_xlabel("Primary ranking retained (%)")
    ax.set_ylabel("Clinical reference drugs recovered (%)")
    ax.set_title("Clinical-Reference Recovery Across Operational Screening Depths", pad=10)
    ax.set_xticks([10, 15, 20])
    ax.margins(y=0.12)
    fig.tight_layout()
    pdf_path, _ = save_figure(fig, args.output_dir, "Figure_4_15_Screening_Depth")
    print(pdf_path)


if __name__ == "__main__":
    main()
