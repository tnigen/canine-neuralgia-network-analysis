# Methods

## Source data and count generation

GSE224814 contains six *Rattus norvegicus* trigeminal-ganglion samples: three
sham and three chronic constriction injury of the infraorbital nerve (CCI-ION),
collected on day 14. The source metadata and sample-to-run assignments are bundled.

The GEO supplementary files contain 47,061 rows each, with 17,108 unique
gene-symbol-like labels and 29,953 duplicated rows. Their expression values are
predominantly non-integer, and the SOFT record describes both TPM and raw counts.
These files were audited but not rounded, aggregated, reverse-normalized, or
used as the PyDESeq2 count input. The audit tables retain the numerical evidence.

Integer counts were generated from the six paired FASTQ runs. ENA URLs, byte
sizes, and MD5 checksums are in `data_raw/ENA_fastq_manifest.tsv`. Each downloaded
pair is checksum-verified before alignment.

HISAT2 2.2.1 uses the official UCSC rn6 index with `--dta --seed 20260904`.
samtools 1.22.1 converts the paired SAM stream to a temporary BAM. featureCounts
2.0.8 counts gene-level fragments using:

```text
-p --countReadPairs -B -C -s 0 -t exon -g gene_id
```

Counting is unstranded. No additional read-trimming step is applied. Ensembl
release 104 Rnor_6.0 supplies the GTF. For counting, chromosome labels
`1..20,X,Y,MT` are converted to `chr1..chr20,chrX,chrY,chrM`; coordinates and gene
identifiers are unchanged. Nonchromosomal contigs are excluded from the counting
GTF. The original GTF is retained for annotation. FASTQs and BAMs are temporary;
the count matrix, featureCounts files, and QC summaries are retained.

## Differential expression

The raw matrix contains 32,623 features and six samples. A gene is retained when
its count is at least 10 in at least three samples: 15,029 features pass and
17,594 are removed. No normalized expression values are supplied to the model.

PyDESeq2 0.5.4 uses design `~condition`, ordered conditions `Sham`, `CCI_ION`,
and the contrast `["condition", "CCI_ION", "Sham"]`. `refit_cooks=True`, Cook's
filtering, independent filtering, and `alpha=0.05` are retained. Positive
log2 fold changes indicate higher expression in CCI-ION. The primary criterion
is nonmissing BH-adjusted `padj < 0.05` and `abs(log2FoldChange) >= 1`.
No fold-change shrinkage step is added. Missing adjusted p values remain missing.

The complete table contains 15,029 rows, including 18 with missing adjusted
p values. There are 2,343 features with `padj < 0.05` irrespective of effect
size, and 431 primary DEGs: 368 upregulated and 63 downregulated. Descriptive
threshold-sensitivity results are stored separately and do not replace this rule.

## Quality control and plots

QC includes assigned library size, detected genes, alignment and assignment
rates, PyDESeq2 VST expression (`use_design=False`), two-component PCA, and
Euclidean sample distances with average-linkage clustering. The NumPy and PCA
seed is 20260904. The archived PCA explains 57.77% and 12.62% of variance on
PC1 and PC2, respectively.

Distances from each sample to its condition centroid in PC1/PC2 are compared
with the pooled upper Tukey fence `Q3 + 1.5 * IQR`. This is exploratory QC;
no samples were excluded. Heatmaps use gene-wise z-scores of VST expression.
Volcano, MA, PCA, heatmap, and library-size plotting settings are unchanged.

## Gene annotation

DEGs are annotated with the original Ensembl-104 GTF and the supplied NCBI
`gene2ensembl` and rat `gene_info` snapshots. Multiple Entrez IDs produce an
ambiguous mapping; a missing GTF gene record produces an unresolved mapping.
Missing secondary identifiers do not remove DEGs. Of 431 features, 430 are
uniquely mapped, one is ambiguous, none is unresolved, and 354 have Entrez IDs.
The ortholog-ready output contains rat identifiers only.

The count-based rerun uses a compact `gene2ensembl` snapshot containing all
original taxon-10116 rows. Only reference-file selection and provenance logging
differ between the two entry-point modes; the statistical analysis is shared.

## Scope

This is a reanalysis of three biological replicates per condition, not a
reconstruction of the submitter's undocumented processing. Differential
expression does not establish causality or cross-species transferability.
Archive integrity and internal consistency checks do not substitute for a fresh
FASTQ-to-results run. Numerical details of refits and rendered figure bytes may
vary across platforms and library builds; the original results remain bundled.
