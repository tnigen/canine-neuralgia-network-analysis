# Canine neuralgia network proximity figures

Plotting scripts, analysis inputs, and rendered figures for the canine neuralgia network proximity analysis. The primary ranking uses STRING combined_score >= 700; the stability panels compare it with the >= 900 ranking.

The Python scripts retain their supplied calculations, plotting parameters, labels, clinical-reference list, and output filenames. Both ROC variants are included: one with a grid and one without.

## Setup

Use Python 3.10 or later and install the dependencies from this directory:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

On Windows, use `.venv\Scripts\activate`. Times New Roman must be installed and discoverable by Matplotlib. The scripts register standard macOS and Windows font locations and stop if the requested font is unavailable. Font files are not distributed here.

The original dependency ranges are preserved in `figure_scripts/requirements.txt`. The rendering environment used during publication preparation is recorded in `REPRODUCIBILITY.md`.

## Inputs

| File | Purpose |
| --- | --- |
| `data/proximity_results_700.tsv` | Primary proximity ranking |
| `data/proximity_results_900.tsv` | Higher-confidence ranking for stability comparisons |
| `data/drug_target_dog_string.xlsx` | Drug-target STRING mapping used for exact-profile redundancy |
| `data/analysis_parameters_700.txt` | Analysis configuration and input checksums for the 700 analysis |
| `data/analysis_parameters_900.txt` | Analysis configuration and input checksums for the 900 analysis |

Only rows with `analysis_status` equal to `analyzed` are used. The input copies come from `../04_network_proximity/`. The mapping workbook is identical to the final workbook in `../02_canine_string_mapping/`. The input and both ranking tables contain the same 1,546 drugs; analyzed rows are selected separately for each panel.

## Generate the figures

Run these commands from this directory:

```bash
python figure_scripts/figure_4_12_rank_stability.py \
  --results-700 data/proximity_results_700.tsv \
  --results-900 data/proximity_results_900.tsv

python figure_scripts/figure_4_12_topn_stability.py \
  --results-700 data/proximity_results_700.tsv \
  --results-900 data/proximity_results_900.tsv

python figure_scripts/network_proximity_ranking_top20.py \
  --results-700 data/proximity_results_700.tsv

python figure_scripts/figure_4_13_roc_clinical_validation.py \
  --results-700 data/proximity_results_700.tsv

python figure_scripts/figure_4_13_roc_clinical_validation_grid.py \
  --results-700 data/proximity_results_700.tsv

python figure_scripts/figure_4_13_clinical_evidence_z_distribution.py \
  --results-700 data/proximity_results_700.tsv

python figure_scripts/figure_4_14_target_count_vs_z.py \
  --results-700 data/proximity_results_700.tsv

python figure_scripts/figure_4_14_exact_profile_redundancy.py \
  --results-700 data/proximity_results_700.tsv \
  --target-workbook data/drug_target_dog_string.xlsx

python figure_scripts/figure_4_15_screening_depth.py \
  --results-700 data/proximity_results_700.tsv
```

Each script writes a vector PDF and a 600-dpi PNG to `outputs/`. Existing files with the same names are replaced. Use `--output-dir` to render to another directory. The nine rendered pairs included here were generated during publication preparation; original article figure files were not supplied for comparison.

## Methods and verification

See [METHODS.md](METHODS.md) for the exact plot calculations and clinical-label limitations. [REPRODUCIBILITY.md](REPRODUCIBILITY.md) records input provenance and verification. Check the archived files with:

```bash
shasum -a 256 -c INPUT_SHA256SUMS
shasum -a 256 -c FIGURE_SHA256SUMS
shasum -a 256 -c SCRIPT_SHA256SUMS
```

The hard-coded clinical-reference names are reproduced as supplied. Trial identifiers, source citations, eligibility criteria, and evidence-curation records were not included. The clinical plot labels therefore should not be taken as a separately verified assessment of clinical evidence or its independence.

## Citation and license

Software citation metadata are in [CITATION.cff](CITATION.cff). The plotting code is released under the [MIT License](LICENSE). Third-party datasets, software, and font files retain their own terms and attribution requirements.
