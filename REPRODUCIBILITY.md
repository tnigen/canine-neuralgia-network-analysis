# Reproducibility

## Final dataset and workflow

The final mapping workbook has 4,919 drug-target rows for 1,546 distinct drugs. The copies in `02_canine_string_mapping`, `04_network_proximity`, and `05_publication_figures/data` are byte-identical. The ortholog output is identical to the STRING-mapping input.

Both proximity rankings and their supporting drug-level QC tables contain the same 1,546 drugs. The 700 ranking includes 1,531 analyzed drugs and 15 skipped drugs. The 900 ranking includes 1,420 analyzed drugs and 126 skipped drugs. Null-distribution rows correspond to the analyzed drugs, with 1,000 draws per drug.

Scientific statistics and plotting parameters are retained. The configuration files in the archived result directories describe the analysis settings and final publication inputs; their `results_source` identifies archived outputs. No new execution dates are assigned to those statistics. The root manifests identify the files in this release snapshot.

## Environments

Use each stage's own requirements file and environment. The historical LCC environment differs from the proximity environment, and a single combined installation is not provided. The plotting dependency ranges also differ from the pinned proximity installation.

The original stage-level environment records are retained where supplied. A new analysis should record its command, versions, input checksums, random seed, hash seed, and generated statistics. Exact Monte Carlo results depend on software and sampling order; archived values should not be described as a newly executed run.

## Available and missing source material

The human-target input for the first stage is not distributed; its checksum remains in the ortholog metadata. BioMart and STRING mapping query live services, so exact historical responses may not be recoverable by repeating those calls.

The LCC folder includes executable code for the 900 setting and archives for the 700 and unfiltered settings. The latter two executable versions and the historical LCC hash seed were not supplied. The inherited report wording is explained in `03_disease_module_lcc/REPORT_NOTES.md`.

The disease-ID file is supplied as the GeneCards/DEG union. Its upstream curation and selection records are not part of this collection. The plotting code retains its supplied clinical-reference names, but clinical trial identifiers, source citations, and eligibility records were not supplied.

The plotting folder includes regenerated PDFs and PNGs from the supplied rankings. Original thesis image files were not supplied for visual identity comparison. The expensive full historical analyses were not rerun during repository preparation.

## Verification environment

The final review checked all saved random-distance means and standard deviations against the retained null distributions. All six LCC regression tests passed. Four proximity checks passed; the full proximity integration test could not complete because the review environment denied operating-system shared-memory allocation. Run the complete proximity suite in a normal local or server environment before a new production analysis. The test and scientific implementation were left unchanged.
