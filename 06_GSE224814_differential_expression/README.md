# GSE224814 differential-expression analysis

Analysis of six rat trigeminal-ganglion RNA-seq samples from
[GSE224814](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE224814)
(BioProject PRJNA932751): three CCI-ION and three sham samples collected 14 days
after surgery. This folder includes raw integer counts, source metadata, code,
archived results, and figures. It does not perform rat-to-dog ortholog mapping.

The archived analysis identified **431 DEGs: 368 upregulated and 63 downregulated**
in CCI-ION versus sham, using `padj < 0.05` and `abs(log2FoldChange) >= 1`.

## Samples

| Sample | GEO sample | Condition | SRA run |
| --- | --- | --- | --- |
| A1_TG | GSM7033526 | Sham | SRR23372014 |
| A2_TG | GSM7033527 | Sham | SRR23372013 |
| A3_TG | GSM7033528 | Sham | SRR23372012 |
| B1_TG | GSM7033529 | CCI_ION | SRR23372011 |
| B2_TG | GSM7033530 | CCI_ION | SRR23372010 |
| B3_TG | GSM7033531 | CCI_ION | SRR23372009 |

## Environment

Use this stage's environment separately from the network-analysis environments:

```bash
conda env create -f environment.yml
conda activate gse224814-analysis
```

The Python packages in `requirements.txt` and command-line tools in
`environment.yml` retain the supplied version pins. The archived run used Python
3.11.16, PyDESeq2 0.5.4, HISAT2 2.2.1, samtools 1.22.1, featureCounts 2.0.8,
and four threads. See `results/tables/software_versions.csv`.

## Rerun from the supplied counts

Run in a separate copy of this folder: the pipeline writes to `results/` and
`data_processed/`, so running in the publication copy would replace archived files.

```bash
python scripts/run_analysis.py --from-counts --threads 4
```

This audits the bundled GEO files, uses the supplied integer count matrix, expands
the bundled Ensembl GTF, and runs the original filtering, PyDESeq2 model, annotation,
QC, and plotting steps. It does not download FASTQs or an alignment index. Allow
about 280 MB for the expanded GTF, in addition to the environment and results.

The bundled `gene2ensembl.rat.tsv.gz` contains every rat row from the original
NCBI snapshot, in its original order. The annotation parser already selects rat
taxonomy ID 10116; the compact file therefore provides the same cross-references.
The full rat `gene_info` snapshot is also included. See [DATA.md](DATA.md).

## Rebuild counts from FASTQs

In a separate working copy, move the existing raw-count matrix and the six
`data_processed/featurecounts/*.featureCounts.txt` files and their summaries out
of that copy. The pipeline reuses any existing count files; removing only the
matrix would merge the saved counts rather than realign reads. Then run:

```bash
python scripts/run_analysis.py --threads 4
```

The runner downloads missing inputs, aligns each paired run to the official
prebuilt rn6 index, counts fragments, and performs the analysis. FASTQs and BAMs
are removed after successful counting; `--keep-fastq` retains the FASTQs. Allow
substantial disk space: the compressed reads total about 20 GB, the index is
about 4 GB, and temporary BAMs require additional space. The original genome
FASTA is restored for provenance, although alignment uses the prebuilt index.

The full workflow retrieves live NCBI annotation tables, which may differ from
the archived snapshot. Use the bundled-count route for the supplied rat annotation
snapshot. Verify downloaded references against the original checksums before
claiming an exact historical reconstruction.

## Files and verification

| Location | Contents |
| --- | --- |
| `scripts/` | Input audit, FASTQ counting, PyDESeq2 analysis, snapshot utility |
| `data_raw/` | GEO archive and extracted files, metadata, reference and annotation snapshots |
| `data_processed/` | Raw/filtered counts, VST expression, metadata, featureCounts outputs |
| `results/differential_expression/` | Complete statistics, DEGs, directional subsets, rat identifiers |
| `results/tables/` | QC, audits, annotation evidence, versions, original source checksums |
| `results/figures/` | PCA, volcano, MA, heatmaps, library sizes; original PNGs and PDFs |
| `results/qc/` | Read-alignment and fragment-assignment summaries |

```bash
python -m unittest discover -s tests -v
```

These offline tests check the saved data and results; they do not refit PyDESeq2
or realign reads. Run `python3 tools/verify_repository.py` from the repository
root to check all recorded file hashes, including this stage's manifests.

See [METHODS.md](METHODS.md) for the model and QC settings, and [DATA.md](DATA.md)
for data provenance, large-file exclusions, and reproduction limits. The software
is covered by [MIT](LICENSE); upstream data and third-party software retain their
own terms. Citation metadata are in [CITATION.cff](CITATION.cff).

[VALIDATION.md](VALIDATION.md) records publication checks and the scope of the
count-based and original-code reruns.
