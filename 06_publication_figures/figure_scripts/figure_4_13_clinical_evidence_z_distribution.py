#!/usr/bin/env python3
"""Create Figure 4.13: Z-score distributions by clinical evidence."""

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

    df07 = load_analyzed(args.results_700, {"drug", "z_score"})
    clinical_mask = df07["drug"].str.strip().str.lower().isin(CT_POSITIVE)
    z_clinical = df07.loc[clinical_mask, "z_score"].astype(float).to_numpy()
    z_background = df07.loc[~clinical_mask, "z_score"].astype(float).to_numpy()

    fig, ax = plt.subplots(figsize=(8.0, 6.4))
    ax.boxplot(
        [z_background, z_clinical],
        tick_labels=["Background", "ClinicalTrials-positive"],
        showfliers=False,
        patch_artist=False,
    )
    rng = np.random.default_rng(42)
    ax.scatter(rng.normal(1, 0.035, len(z_background)), z_background, s=10, alpha=0.14)
    ax.scatter(rng.normal(2, 0.035, len(z_clinical)), z_clinical, s=26, alpha=0.72)
    ax.axhline(0, linestyle="--", linewidth=1.2)
    ax.set_ylabel("Network-proximity Z-score")
    ax.set_title("Network Proximity by Independent Clinical Evidence", pad=10)
    ax.grid(axis="both", linewidth=0.45, alpha=0.25)
    fig.tight_layout()
    pdf_path, _ = save_figure(
        fig,
        args.output_dir,
        "Figure_4_13_Clinical_Evidence_Z_Distribution",
    )
    print(pdf_path)


if __name__ == "__main__":
    main()
