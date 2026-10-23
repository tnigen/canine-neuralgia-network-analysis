# Reproducibility record

## Source inputs

The two ranking TSVs and their configuration files match `../04_network_proximity/results/score_700/` and `../04_network_proximity/results/score_900/`. The mapping workbook matches the input in that analysis's root directory and the final output of `../02_canine_string_mapping/`. Its SHA-256 is `1b1f1f5f8ca85b5c3bddbdbd29b6373f26c3b425e852aee307eecb838fc9031c`.

The configuration records describe the upstream network analysis and publication inputs; their relative paths refer to that analysis. Figure-generation commands and copied input paths are documented separately in this repository's README.

All ten supplied Python source files are preserved byte-for-byte, including shared utilities. No algorithms, annotations, clinical labels, statistical parameters, plot dimensions, colors, font settings, or output names were changed. The duplicate supplied ZIP is backed up outside the publication directory to avoid distributing a stale second copy.

## Rendering environment

The nine plotting scripts were executed on 22 October 2026 with Python 3.13.15, Matplotlib 3.11.1, NumPy 2.5.2, pandas 3.0.5, and openpyxl 3.1.5. Times New Roman was available from the operating system. Each script produced a PDF and a 600-dpi PNG.

The original package constraints remain unchanged. Different library versions or font installations can change rendering. PDF metadata can include a new creation timestamp on each run, so the figure checksums identify this rendered snapshot rather than promising byte-identical future PDFs.

## Verification scope

Verification checked input columns, one analyzed row per drug, finite Z-scores, non-empty clinical label classes, workbook coverage of the primary analyzed drugs, and agreement between the workbook hash and the configuration records. The ROC AUC was independently checked against pairwise positive-versus-background score comparisons, including ties.

All nine scripts completed successfully. The PNGs and PDFs were checked for expected counts and readable output. The original article image files were not supplied, so equivalence to previously published or selected figure versions cannot be established from this bundle alone.

The clinical-reference list is preserved, but its upstream evidence provenance was not supplied. These scripts reproduce that labeling choice without independently verifying clinical claims.

Both archived rankings and the input workbook contain the same 1,546 drugs. The primary figure cohort consists of the 1,531 analyzed rows in the 700 ranking; 1,420 analyzed rows enter the 900 comparisons.
