# Data and provenance

## Bundled files

The repository includes the original GEO supplementary archive and its six
compressed expression files, GEO SOFT metadata, SRA run information, the ENA
paired-read manifest, and the downloaded publication XML. These source files
are unchanged. Public source records retain their original attribution.

The integer raw-count matrix, filtered matrix, sample metadata, VST expression,
featureCounts gene coordinates, all six featureCounts outputs, DE tables,
annotation tables, QC summaries, and rendered figures are included. Scientific
values are unchanged. Personal machine prefixes in featureCounts headers and
the HISAT2 version record are replaced with relative paths; temporary-directory
prefixes are represented by `<temporary>/`.

Two original annotation snapshots are bundled:

- `data_raw/reference/Rattus_norvegicus.Rnor_6.0.104.gtf.gz`:
  the complete Ensembl release-104 GTF, expanded on demand.
- `data_raw/annotation/Rattus_norvegicus.gene_info.gz`:
  the original rat NCBI Gene information table.

`data_raw/annotation/gene2ensembl.rat.tsv.gz` is a compact derivative of the
original global `gene2ensembl.gz`: it preserves the header and all 36,316 rows
whose first field is taxonomy ID 10116, in their original order. No gene IDs,
relationships, or fields are altered. Regenerate it from the original global
snapshot with:

```bash
python scripts/subset_rat_crossrefs.py --source /path/to/gene2ensembl.gz \
  --output /path/to/gene2ensembl.rat.tsv.gz
```

The two `/path/to/` arguments are placeholders for files on the reader's system.
Re-fetching today's global file may produce a different rat snapshot.

## Files not distributed in Git

Installed environments, caches, session history, and machine-specific execution
logs are not part of the publication bundle. Large reference FASTA files,
expanded GTFs, the prebuilt HISAT2 index, FASTQs, and BAMs are also excluded.
They are not required to inspect the saved results or run `--from-counts`.

Reference retrieval for the full raw-read workflow is implemented in the scripts:

| Reference | Source |
| --- | --- |
| rn6 HISAT2 index | [Official prebuilt index](https://genome-idx.s3.amazonaws.com/hisat/rn6_genome.tar.gz) |
| Rnor_6.0 FASTA | [Ensembl release 104](https://ftp.ensembl.org/pub/release-104/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz) |
| Rnor_6.0 GTF | [Ensembl release 104](https://ftp.ensembl.org/pub/release-104/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.104.gtf.gz) |
| Global gene cross-references | [NCBI gene2ensembl](https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz) |
| Rat gene information | [NCBI Gene](https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Rattus_norvegicus.gene_info.gz) |

The global cross-reference table and genome/index files can exceed GitHub's
[100 MiB regular-file limit](https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github).
The original local project is retained outside the repository; it is not a
public data archive. Original FASTQs are retrieved from the public ENA locations
and verified against the MD5 values in the bundled manifest.

## Checksums and reproduction

`INPUT_SHA256SUMS` covers the bundled source/reference files, raw counts, sample
metadata, and featureCounts records. `RESULT_SHA256SUMS` covers saved processed
matrices and results. `SCRIPT_SHA256SUMS` covers scripts and tests. The repository's
standard-library verifier checks these manifests without scientific dependencies.

`results/tables/source_file_checksums.csv` is the unchanged provenance record
from the original run. It includes larger files not bundled in Git; it is not
a list of files expected to be present in a fresh clone. A rerun replaces that
table with checksums for the sources actually used. Run in a working copy to
preserve the archived table and figures.

The bundled-count route preserves the annotation snapshot and all model settings.
The full route may retrieve newer NCBI annotations. Compare their original
checksums before asserting exact historical equivalence; no immutable public
archive of the excluded global NCBI snapshot is claimed here.

The downloaded SOFT record links a publication title different from its Series
title. The original XML and SOFT are retained; sample assignments are based on
GEO/SRA sample metadata, not inferred from the linked paper.

## Attribution

The MIT software license does not relicense GEO/SRA/ENA, Ensembl, NCBI, or
HISAT2 content. Cite GSE224814 and the relevant original dataset publication
alongside this analysis software. Consult upstream terms for redistributed data
and third-party tools.
