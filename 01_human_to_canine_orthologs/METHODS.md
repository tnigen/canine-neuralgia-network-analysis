# Computational methods

## Input processing

The program accepts a CSV or Excel workbook containing `drug` and `human_target_protein` columns. Leading and trailing whitespace is removed, empty values are retained for quality-control reporting, and exact duplicate drug-target relationships are removed before querying.

## BioMart configuration

Each run queries the live Ensembl BioMart registry and selects the current Ensembl Genes mart and Homo sapiens genes dataset. The required human symbol filter, human gene identifiers, HGNC identifier, canine gene identifiers, canine gene symbol, and canine homology-type attribute are identified from the live metadata. The program stops if a required field cannot be identified unambiguously.

## Mapping procedure

Mapping is performed in two stages:

1. Unique non-empty human symbols are resolved to human Ensembl gene IDs, human gene names, and HGNC IDs.
2. Resolved human Ensembl gene IDs are queried for canine orthology records.

Queries are batched. A failed batch is divided recursively to isolate individual values that cannot be queried. A failed gene-ID query is reported as a query error rather than as absence of an ortholog.

A human symbol is considered unambiguous when its returned records correspond to one HGNC identity. Multiple Ensembl loci sharing that HGNC identifier are treated as representations of the same biological gene.

## Inclusion criteria

A canine mapping is retained only when the returned BioMart homology type is exactly `ortholog_one2one` and the canine gene symbol is non-empty and usable. Percentage identity is not queried or used. Symbols are not corrected, inferred, or substituted.

The final network table contains unique `drug` and `dog_target_protein` pairs.

## Quality control and summary

The quality-control workbook records the human and canine identifiers, homology type, mapping status, and notes for each deduplicated input relationship. Excluded homology classes and unresolved or ambiguous human symbols remain visible in this audit output.

For each drug, the mapping percentage is:

```text
100 × (unique input human targets producing at least one retained canine symbol)
      / (unique input human targets for that drug)
```

Run metadata record the software version, Python version, resolved BioMart configuration, inclusion criterion, input filename and SHA-256 checksum, and aggregate input counts.

## Reproducibility

The included results were generated with script version 1.0.1 and Ensembl Genes 116. BioMart is a live resource, so rerunning the program against another Ensembl release can produce different mappings. The saved metadata and checksums identify the exact input and outputs used for the archived run.
