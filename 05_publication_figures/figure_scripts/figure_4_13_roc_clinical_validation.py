#!/usr/bin/env python3
"""Create Figure 4.13: clinical-validation ROC curve without a grid."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from common import CT_POSITIVE, load_analyzed, save_figure


def roc_points(labels: np.ndarray, scores: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    order = np.argsort(-scores, kind="mergesort")
    labels = labels[order]
    scores = scores[order]
    true_positives = np.cumsum(labels)
    false_positives = np.cumsum(1 - labels)
    threshold_ends = np.r_[np.where(np.diff(scores))[0], labels.size - 1]
    tps = np.r_[0, true_positives[threshold_ends]]
    fps = np.r_[0, false_positives[threshold_ends]]
    tpr = tps / tps[-1]
    fpr = fps / fps[-1]
    return fpr, tpr, float(np.trapezoid(tpr, fpr))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-700", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()

    df07 = load_analyzed(args.results_700, {"drug", "z_score"})
    labels = df07["drug"].str.strip().str.lower().isin(CT_POSITIVE).astype(int).to_numpy()
    scores = -df07["z_score"].astype(float).to_numpy()
    fpr, tpr, auc = roc_points(labels, scores)

    fig, ax = plt.subplots(figsize=(7.0, 6.4))
    ax.plot(fpr, tpr, linewidth=2.2, label=f"AUC = {auc:.3f}")
    ax.plot([0, 1], [0, 1], linestyle="--", linewidth=1.4)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("False-positive rate")
    ax.set_ylabel("True-positive rate")
    ax.set_title("Clinical Validation of Network-Proximity Ranking", pad=10)
    ax.legend(loc="lower right", frameon=True, fontsize=12)
    fig.tight_layout()
    pdf_path, _ = save_figure(fig, args.output_dir, "Figure_4_13_ROC_Clinical_Validation")
    print(pdf_path)


if __name__ == "__main__":
    main()
