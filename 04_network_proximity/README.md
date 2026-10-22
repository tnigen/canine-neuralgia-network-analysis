# Canine drug–disease network proximity

This repository implements a reproducible canine STRING v12 network-proximity
analysis for the neuralgia `GeneCards ∪ DEG` disease set. It uses NetMedPy
0.2.0.1 for shortest-path distances, degree-matched randomization, Z-scores,
and empirical P-values.

The analysis prioritizes drug targets by network proximity. It does not
establish efficacy, dose suitability, or safety in dogs.

## Inputs

- `9615.protein.links.v12.0.txt.gz`: STRING v12 canine protein links with
  columns `protein1`, `protein2`, and `combined_score`.
- `canine_neuralgia_union.txt`: one canine STRING protein identifier per line.
- `drug_target_dog_string.xlsx`: worksheet `Sheet1` with columns `drug`,
  `dog_target_protein`, and `string_protein_id`.

The code validates all three inputs. `dog_target_protein` is retained as
provenance; `string_protein_id` is used for network mapping. Input hashes are
listed in `INPUT_SHA256SUMS`.

## Installation

Python 3.9 or newer is recommended.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

## Analysis

The primary analysis retains STRING interactions with
`combined_score >= 700`:

```bash
python run_proximity.py \
    --network 9615.protein.links.v12.0.txt.gz \
    --disease canine_neuralgia_union.txt \
    --drug-targets drug_target_dog_string.xlsx \
    --score-cutoff 700 \
    --output results/score_700
```

The higher-confidence sensitivity analysis uses a cutoff of 900:

```bash
python run_proximity.py \
    --network 9615.protein.links.v12.0.txt.gz \
    --disease canine_neuralgia_union.txt \
    --drug-targets drug_target_dog_string.xlsx \
    --score-cutoff 900 \
    --output results/score_900
```

To retain every valid interaction, use `--score-cutoff none`. This disables
score filtering; it is not equivalent to converting the cutoff to zero.

```bash
python run_proximity.py \
    --network 9615.protein.links.v12.0.txt.gz \
    --disease canine_neuralgia_union.txt \
    --drug-targets drug_target_dog_string.xlsx \
    --score-cutoff none \
    --output results/score_none
```

Output directories must be empty. Use `--dry-run` to generate validation and
coverage files without calculating the distance matrix or proximity scores.
See `AWS_RUNBOOK.md` for the resource-intensive unfiltered analysis.

## Reproducibility settings

The defaults reproduce the reported workflow:

- 1,000 randomizations (`--n-randomizations 1000`)
- exact degree matching (`--null-model degree_match`)
- random seed 42 (`--seed 42`)
- unweighted shortest paths
- directed closest distance from each drug target to the disease set

The wrapper fixes both Python and NumPy random seeds and restarts once with a
matching `PYTHONHASHSEED`. `--n-processors` changes only distance-matrix
parallelism. The optional `--no-save-distance-matrix` flag omits the large
derived matrix.

NetMedPy requires a connected graph. The analysis therefore uses the largest
connected component of the score-filtered background interactome. It does not
replace the disease set with a disease-induced component. All excluded disease
proteins and drug targets are recorded in the QC outputs.

## Results

Publication results for cutoffs 700 and 900 are included under `results/`.
The publication input and both ranking tables cover the same 1,546 drugs.
The primary table contains 1,531 analyzed drugs; the sensitivity table contains
1,420. Other drugs retain their explicit skip status.
The derived distance matrices are omitted because they are large and can be
regenerated from the recorded inputs and parameters. The included result
tables have not been recalculated or rounded during repository preparation.

Each completed run writes:

- `proximity_results.tsv`: ranked drug-level results and explicit skip states.
- `netmedpy_null_distributions.tsv.gz`: randomized closest-distance values.
- `netmedpy_shortest_path_distances.npz`: distance matrix, unless disabled.
- `analysis_parameters.txt`: hashes, package versions, parameters, and status.
- network, disease-set, workbook, and drug-target QC tables.
- exact included and excluded disease proteins and drug targets.

Results are ordered by ascending Z-score, followed by lower-tail empirical
P-value and observed closest distance. More negative Z-scores indicate
stronger-than-null proximity.

## Validation

```bash
sha256sum -c INPUT_SHA256SUMS
sha256sum -c RESULT_SHA256SUMS
python -m unittest discover -s tests -v
```

`METHODS.md` documents the equation, null model, graph preprocessing, and the
division of responsibility between NetMedPy and this wrapper.

## Citation

Citation metadata are provided in `CITATION.cff`. Cite the immutable `v1.0.0`
release rather than the moving default branch, and record the release tag or
commit identifier in the manuscript. Add the repository URL and release date
to `CITATION.cff` when the GitHub repository and release are created.

## License

The original software and documentation are available under the MIT License;
see `LICENSE`. Third-party data retain their original terms. The STRING
interaction file is distributed under CC BY 4.0 and requires attribution.

## References

- Ma M, Huang M, He Y, et al. Network Medicine: A Potential Approach for
  Virtual Drug Screening. *Pharmaceuticals*. 2024;17:899.
  https://doi.org/10.3390/ph17070899
- NetMedPy: https://github.com/menicgiulia/NetMedPy
- NetMedPy 0.2.0.1: https://pypi.org/project/netmedpy/0.2.0.1/
