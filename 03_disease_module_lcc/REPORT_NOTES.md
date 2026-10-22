# Notes on the supplied reports

The following issues were identified during publication preparation. Scientific calculations, archived numerical results, missing-ID lists, and original reports have been preserved.

## Score thresholds

The supplied script contains `if score < 900: continue`. Both the 900 and 700 archived reports nevertheless state `No combined_score threshold was applied.` Their limitations paragraphs repeat the same unfiltered-network description. These statements do not describe those thresholded analyses.

| Archive | Recorded minimum score | Recorded graph nodes | Recorded graph edges |
| --- | --- | --- | --- |
| 900 | 900 | 8,674 | 49,361 |
| 700 | 700 | 15,037 | 157,242 |
| No additional cutoff | 150 | 19,921 | 5,384,258 |

The 900 cutoff is established directly by the supplied code. The 700 setting is indicated by its archive name and recorded score minimum; its original executable was not supplied. No additional cutoff means no filtering beyond the records available in the supplied file, whose reported minimum is 150.

The report-generating strings have also been preserved, so a future run of this script will repeat the inherited threshold wording. Use the code, README, and these notes when describing the method.

## Seed and Monte Carlo variation

The script seeds Python and NumPy, but the supplied NetMedPy source converts the disease IDs to a set and samples while iterating over that set. Set order can change with Python's process-level hash randomization. The archived hash seed was not recorded. Identical input files and `--seed 42` alone therefore do not establish bitwise reproducibility of the historical Monte Carlo statistics.

The null distribution itself was not saved. The archived CSVs and their checksums identify the available historical results. No full randomization analysis was rerun during editorial preparation.

## Interpretation text for non-default options

The generated explanatory paragraph always describes a log-binning null model even when another supported CLI null model is selected. The p=0 explanatory note always refers to 10,000 draws and p<0.0001, even if `--n-iterations` is changed. These fixed strings fit the archived default configuration but should not be quoted for other configurations. Numerical calculations use the selected options.

An empirical p-value of zero records no exceedances in the sampled null distribution; it does not establish zero underlying probability. The original numerical output is retained without changing the estimator or adding a pseudocount.

## Scope of publication cleanup

Only comments, the module docstring, and whitespace in the supplied Python file were edited. Executable statements, parameters, report strings, CSV formatting, and data files were preserved. The environment and operating-system records contain software provenance rather than personal identifiers and remain intact.
