"""Shared data-loading and Times New Roman plotting utilities."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib import font_manager
import pandas as pd


CT_POSITIVE = {
    "adenosine",
    "capsaicin",
    "clonidine",
    "dextromethorphan",
    "diclofenac",
    "dronabinol",
    "duloxetine",
    "eslicarbazepine acetate",
    "gabapentin",
    "gabapentin enacarbil",
    "galcanezumab",
    "guanfacine",
    "hydromorphone",
    "lacosamide",
    "lidocaine",
    "lx9211",
    "mirogabalin",
    "nabiximols",
    "botulinum toxin type a",
    "perampanel",
    "pregabalin",
    "tapentadol",
    "tariquidar",
    "tramadol",
}


def configure_times_new_roman() -> None:
    """Configure Matplotlib for publication figures using Times New Roman."""
    candidate_fonts = (
        Path("/System/Library/Fonts/Supplemental/Times New Roman.ttf"),
        Path("/System/Library/Fonts/Supplemental/Times New Roman Bold.ttf"),
        Path("/System/Library/Fonts/Supplemental/Times New Roman Italic.ttf"),
        Path("C:/Windows/Fonts/times.ttf"),
        Path("C:/Windows/Fonts/timesbd.ttf"),
        Path("C:/Windows/Fonts/timesi.ttf"),
    )
    for font_path in candidate_fonts:
        if font_path.exists():
            font_manager.fontManager.addfont(str(font_path))

    try:
        font_manager.findfont(
            font_manager.FontProperties(family="Times New Roman"),
            fallback_to_default=False,
        )
    except ValueError as exc:
        raise RuntimeError(
            "Times New Roman is not installed. Install it before generating figures."
        ) from exc

    plt.rcParams.update(
        {
            "font.family": "Times New Roman",
            "font.size": 13,
            "axes.titlesize": 17,
            "axes.labelsize": 15,
            "legend.fontsize": 12,
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "mathtext.fontset": "custom",
            "mathtext.rm": "Times New Roman",
            "mathtext.it": "Times New Roman:italic",
            "mathtext.bf": "Times New Roman:bold",
        }
    )


def load_analyzed(path: Path, required_columns: Iterable[str]) -> pd.DataFrame:
    """Load a proximity TSV and keep rows marked as analyzed."""
    frame = pd.read_csv(path, sep="\t")
    required = {"analysis_status", *required_columns}
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"Missing required columns in {path}: {', '.join(missing)}")
    return frame.loc[frame["analysis_status"].eq("analyzed")].copy()


def save_figure(fig, output_dir: Path, stem: str) -> tuple[Path, Path]:
    """Save vector PDF and 600-dpi PNG versions of a figure."""
    output_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = output_dir / f"{stem}.pdf"
    png_path = output_dir / f"{stem}.png"
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=600, bbox_inches="tight")
    plt.close(fig)
    return pdf_path, png_path


configure_times_new_roman()
