# Human drug targets to canine one-to-one orthologs

This repository maps human drug-target gene symbols to canine gene symbols using the current Ensembl Genes BioMart. It retains only strict one-to-one orthologs and produces network-ready mappings, quality-control records, per-drug summaries, and run metadata.

## Requirements

- Python 3.9 or later
- Internet access to the Ensembl BioMart service

Create an isolated environment and install the pinned dependency ranges:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

On Windows, activate the environment with `.venv\Scripts\activate`.

## Input

Provide a CSV or Excel workbook with these columns:

| Column | Description |
| --- | --- |
| `drug` | Drug identifier or name |
| `human_target_protein` | Human HGNC gene symbol |

Values are trimmed, and exact duplicate drug-target relationships are removed. Symbols are not corrected, inferred, or substituted.

The source dataset used for the included results is not distributed in this repository. Its SHA-256 checksum is recorded in `results/mapping_run_metadata.txt`.

## Usage

```bash
python map_drug_targets_to_dog_biomart.py \
  --input path/to/drug_target_human.xlsx \
  --output-dir results
```

CSV input is also supported. Run `python map_drug_targets_to_dog_biomart.py --help` for query and retry options.

## Method

The script inspects the live BioMart registry, human dataset, filters, and attributes on each run. It first resolves human symbols to Ensembl gene and HGNC identifiers, then queries canine homology by human Ensembl gene ID. This two-stage query avoids mixing BioMart feature and homology attribute pages.

A relationship is retained only when:

1. the human symbol resolves to one HGNC identity;
2. the BioMart homology type is exactly `ortholog_one2one`; and
3. the canine record has a usable gene symbol.

Alternative Ensembl loci that share one HGNC identifier are treated as one biological gene rather than as symbol ambiguity. See [METHODS.md](METHODS.md) for the complete computational procedure.

Because BioMart is a live service, mappings can change across Ensembl releases. Each run records the resolved service, dataset, attributes, software version, input checksum, and summary counts.

## Outputs

- `drug_target_dog_one2one.csv` and `.xlsx`: deduplicated `drug` and `dog_target_protein` relationships.
- `drug_target_mapping_qc.xlsx`: audit records for retained and excluded input relationships.
- `drug_target_mapping_summary.xlsx`: input, mapped, unmapped, and percentage counts by drug.
- `mapping_run_metadata.txt`: run provenance and mapping criteria.

The `results/` directory contains a frozen run generated with script version 1.0.1 and Ensembl Genes 116. File checksums are provided in `RESULT_SHA256SUMS`.

## Citation

Citation metadata are provided in [CITATION.cff](CITATION.cff). GitHub displays these details through its **Cite this repository** feature.

## License

This software is released under the [MIT License](LICENSE).
