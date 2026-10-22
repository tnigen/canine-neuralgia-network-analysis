# Canine neuralgia LCC significance analysis

This repository evaluates the largest connected component (LCC) induced by canine neuralgia proteins in a STRING association network. The script uses NetMedPy to compare the observed LCC size with randomized node sets under a selected null model.

## Installation

The archived analyses record Python 3.13.15 and NetMedPy 0.2.0.1. Create an environment with Python 3.13 and install the recorded scientific package versions:

```bash
python3.13 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

On Windows, activate the environment with `.venv\Scripts\activate`. The original environment records are retained alongside the results. `requirements-lock.txt` captures all package versions found in the supplied local environment; it does not replace the historical environment records.

## Inputs

- `9615.protein.links.v12.0.txt.gz`: gzipped STRING v12.0 dog association table with columns `protein1`, `protein2`, and `combined_score`.
- `canine_neuralgia_genes.txt`: one canine STRING protein ID per non-empty line. The input contains 3,363 unique IDs.

The disease-ID file is supplied as an analysis input. Its upstream selection and biological curation procedure are not included in this bundle.

## Run the supplied analysis

From the repository directory:

```bash
python canine_neuralgia_LCC_netmedpy.py \
  --interactome 9615.protein.links.v12.0.txt.gz \
  --disease-genes canine_neuralgia_genes.txt \
  --output-dir runs/score900 \
  --null-model log_binning \
  --n-iterations 10000 \
  --bin-size 100 \
  --seed 42
```

The score cutoff is fixed at **combined_score >= 900** in the supplied script. There is no cutoff option. The graph is simple, undirected, and unweighted for LCC calculations; scores remain edge metadata. Run `python canine_neuralgia_LCC_netmedpy.py --help` for the supported options.

Use a separate output directory as shown above. Files with the generated output names are replaced if they already exist in that directory.

## Archived analyses

| Archived setting | Location | Corresponding script supplied? |
| --- | --- | --- |
| Score >= 900 | `results/score_900/` | Yes |
| Score >= 700 | `results/score_700/` | No |
| No additional score cutoff | `results/unfiltered/` | No |

The 700 and unfiltered archives are preserved for comparison. The supplied script directly implements only the 900 setting. Archives use the script's standard output filenames within their score-specific directories. File contents are unchanged by this reorganization.

**Read [REPORT_NOTES.md](REPORT_NOTES.md) before quoting the reports.** The 900 and 700 reports incorrectly state that no score threshold was applied. The supplied code and recorded score minima show the threshold used. Archived files are preserved unchanged.

## Outputs

- `canine_neuralgia_LCC_results.csv`: module coverage, graph size, observed LCC, null statistics, empirical p-value, and software versions.
- `canine_neuralgia_LCC_report.txt`: run parameters, input checks, results, and interpretation.
- `canine_neuralgia_missing_from_interactome.txt`: disease IDs absent from the thresholded graph.
- `environment_versions.txt`: Python, platform, and scientific package versions.

The program does not export the null distribution or LCC membership, although NetMedPy returns them in memory. See [METHODS.md](METHODS.md) for the exact procedure and [REPOSITORY_CONTENTS.md](REPOSITORY_CONTENTS.md) for the file inventory.

## Verification and reproducibility

Verify the supplied inputs and archived outputs on macOS or Linux:

```bash
shasum -a 256 -c INPUT_SHA256SUMS
shasum -a 256 -c RESULT_SHA256SUMS
python -m unittest discover -s tests -v
```

The script seeds Python and NumPy with 42 by default. The supplied NetMedPy implementation samples while iterating over Python sets, so the seed alone does not guarantee identical Monte Carlo draws across separate processes. The historical `PYTHONHASHSEED` was not recorded. The archived numerical results and checksums remain the reference for the original analyses.

For a new analysis, record an explicit `PYTHONHASHSEED` before starting Python, together with the package versions, input hashes, and command. This controls an additional source of variation but cannot reconstruct the unrecorded historical hash seed. Full 10,000-iteration analyses can take substantial time, especially on the unfiltered graph.

## Citation and license

Software citation metadata are in [CITATION.cff](CITATION.cff). Release version `1.0.0` identifies this publication package. The Python code is released under the [MIT License](LICENSE). STRING data and NetMedPy retain their own source terms and attribution requirements.
