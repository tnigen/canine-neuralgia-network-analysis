# Computational methods

## Graph construction

The supplied script reads the gzipped STRING v12.0 links file and requires the header `protein1 protein2 combined_score`. Empty lines are skipped. Every other record must have exactly three fields, and scores must parse as integers.

Records with `combined_score < 900` are excluded. Retained protein pairs form a simple undirected NetworkX graph. Repeated and reversed pairs do not create additional edges. The first retained score for a pair is stored as `combined_score`; the first orientation is stored for duplicate counting. Self-loop records are counted and retained if present.

LCC extraction and the default degree-based null model use topology. STRING scores are not used as path lengths or LCC weights.

## Disease-set matching

The disease file must contain one identifier per non-empty line. IDs are deduplicated in first-occurrence order without symbol conversion or identifier substitution. IDs present in the retained graph are passed to NetMedPy. Absent IDs are written to a separate file. The analysis stops if none of the disease IDs are present.

## Null model and statistics

The default call is:

```python
netmedpy.lcc_significance(
    graph,
    mapped_ids,
    null_model="log_binning",
    n_iter=10000,
    bin_size=100,
)
```

Python's `random` and NumPy's random generator are seeded with 42 immediately before this call. CLI options also allow `degree_match` and `uniform`; the archived results all record `log_binning`.

In the supplied NetMedPy 0.2.0.1 source, the observed LCC is extracted from the disease-induced graph. Randomized node sets yield a distribution of LCC sizes. The null mean and standard deviation use `numpy.mean` and `numpy.std` with its default population standard deviation. The Z-score is `(observed_lcc_size - null_mean) / null_sd`; zero null variance returns `NaN`. The empirical upper-tail p-value is the number of randomized LCC sizes at least as large as observed divided by the number of iterations. No pseudocount is added.

For default log-binning, the supplied implementation constructs degree bins and samples unique replacement nodes from the bin associated with each observed node's degree. It does not globally exclude observed disease nodes from those bins. The uniform option excludes the observed disease set from its candidate pool.

## Results and provenance

The CSV records disease coverage, graph counts, observed LCC size, null mean, null standard deviation, Z-score, p-value, package versions, selected null model, and iteration count. Floating-point statistics are formatted to 12 significant digits. Seed and bin size are recorded in the text report.

The archived 700 and unfiltered analyses have their own reports, CSVs, missing-ID lists, and environment records. Their corresponding script versions were not supplied. The executable code in this repository implements the 900 cutoff only. [REPORT_NOTES.md](REPORT_NOTES.md) documents inherited report inconsistencies and Monte Carlo reproducibility limits.
