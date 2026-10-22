# Methods

## Study design and inputs

The background interactome is the supplied STRING v12 canine (`Canis lupus familiaris`, taxon 9615) protein-links file. Each valid undirected pair is represented once. `combined_score` is used only for the command-line-selected inclusion threshold; it is not transformed and is not attached as a shortest-path weight. With an integer cutoff, records satisfy `combined_score >= cutoff`. With `--score-cutoff none`, all valid interaction records pass the score filter.

The disease set is the complete set of unique, valid canine STRING identifiers in `canine_neuralgia_union.txt`, defined as the `GeneCards ∪ DEG` union. This union, after documented network-coverage exclusions, is passed to proximity. A disease-induced LCC is not extracted or substituted.

Drug targets are grouped by `drug` from `drug_target_dog_string.xlsx`. `string_protein_id` supplies graph identifiers; `dog_target_protein` is retained for provenance. The workbook contains the canine target-to-STRING mapping. No ortholog inference, alternative-identifier search, or artificial identifier construction is performed. Empty `string_protein_id` cells are excluded as `missing_string_protein_id`, not interpreted as biological absence. Nonempty identifiers absent from the selected network are separately excluded as `absent_from_selected_interactome`.

## What comes from Ma et al. (2024)

The single-drug closest distance is

\[
d_c(S,T)=\frac{1}{|T|}\sum_{t\in T}\min_{s\in S} d(s,t),
\]

where `S` is the disease-protein set, `T` is one drug's complete valid target set, and `d(s,t)` is shortest-path distance in the background interactome. The paper describes standardizing observed proximity against repeated paired random node-set distances and ranking drugs in ascending Z-score order; negative Z-scores represent closer-than-null proximity. It discusses LCC analysis as disease-module validation, which is methodologically distinct from replacing the disease union in the proximity calculation.

Reference: Ma M, Huang M, He Y, et al. *Network Medicine: A Potential Approach for Virtual Drug Screening.* Pharmaceuticals. 2024;17:899. https://doi.org/10.3390/ph17070899

## What comes from NetMedPy 0.2.0.1

The graph is a `networkx.Graph` and identifiers are arbitrary hashable node names; NetMedPy contains no human-identifier requirement in this workflow. The wrapper calls:

1. `netmedpy.all_pair_distances(graph, distance="shortest_path")` on an unweighted graph.
2. `netmedpy.proximity(graph, T, S, distance_matrix, null_model="degree_match", n_iter=1000, symmetric=False)` for each drug.

Passing drug targets as NetMedPy's first set `T` and disease proteins as its second set `S`, with `symmetric=False`, implements Equation (1): the function averages each target's minimum distance to the disease set. It does not collapse a multi-target drug to a representative target.

NetMedPy's default `degree_match` null model constructs exact-degree buckets. In every iteration it independently replaces every drug target and every disease protein with a node having the same graph degree, sampling without replacement within each randomized set. Randomized target and disease sets may overlap each other. This preserves both set sizes and their exact degree distributions. The current function default of 1000 iterations is used unless `--n-randomizations` explicitly changes it.

For the randomized distances `R`, NetMedPy calculates the population standard deviation (`numpy.std`, `ddof=0`) and

\[
Z=\frac{d_c-\operatorname{mean}(R)}{\operatorname{sd}(R)}.
\]

The reported `empirical_p_value` is NetMedPy's lower-tail fraction `count(R <= observed_dc) / n_iter`; no plus-one correction is added. Its two-sided empirical value, `count(|R-mean(R)| >= |observed_dc-mean(R)|) / n_iter`, is retained separately. The full returned null values are saved in `netmedpy_null_distributions.tsv.gz`.

NetMedPy does not expose a random-seed argument for `proximity`; internally it uses Python's `random` module. The wrapper sets `random.seed(--seed)` and `numpy.random.seed(--seed)`. It also re-executes once with `PYTHONHASHSEED` equal to the requested seed so NetMedPy's internal set iteration is reproducible. Calls are sequential and drugs are processed in sorted name order.

NetMedPy source and documentation: https://github.com/menicgiulia/NetMedPy and https://pypi.org/project/netmedpy/0.2.0.1/

## Necessary canine-specific preprocessing

The STRING gzip file is read directly and its exact header is validated as `protein1 protein2 combined_score`. Malformed rows, score-filtered records, self-loops, and orientation-independent duplicates are counted separately. Because shortest paths are topological here, the resulting NetworkX graph has no `weight` edge attribute.

NetMedPy 0.2.0.1 explicitly raises an error for a disconnected graph and its distance matrix contains no policy for converting unreachable distances to finite values. Components are therefore calculated on the complete selected background interactome, and NetMedPy receives the largest connected component of that background. This is not a disease-induced LCC. Disease proteins and targets that occur in the selected interactome but outside this background component are excluded with reason `outside_netmedpy_background_lcc`; identifiers absent from the selected interactome remain a different exclusion class. No infinite distance is imputed.

Only drugs with at least one valid target in the NetMedPy background component are analyzed. Drugs with zero valid targets remain in QC and the main results with an explicit skip status. Coverage below 50% of supplied target rows is flagged as low coverage but is not otherwise reweighted or imputed.

## Interpretation

More negative Z-scores indicate closer-than-expected placement of a drug's targets relative to the canine neuralgia union under NetMedPy's degree-matched null. Empirical P-values quantify this network result at the requested randomization resolution. These outputs prioritize candidates for follow-up; they do not demonstrate canine efficacy, dose suitability, or safety.
