# Plotting scripts

This directory contains nine plotting scripts and shared input/output and typography helpers in `common.py`. The supplied Python source files and dependency ranges are preserved unchanged.

Use the [repository README](../README.md) for installation, input descriptions, and commands that run from the repository root. Inputs are in `../data/`, and verified PDF and PNG outputs are in `../outputs/`.

Every script accepts `--results-700`. The two stability scripts also require `--results-900`, and exact-profile redundancy additionally requires `--target-workbook`. Use `--output-dir` to select the output directory; its default is `outputs/` relative to the working directory.

The default output format is vector PDF plus PNG at 600 dpi. Times New Roman is required. See [METHODS.md](../METHODS.md) for calculations and [REPRODUCIBILITY.md](../REPRODUCIBILITY.md) for input provenance and verification limits.
