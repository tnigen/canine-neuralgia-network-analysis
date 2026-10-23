# Canine neuralgia network medicine analysis

Code, inputs, archived results, and plotting scripts used for canine neuralgia network analysis. The folders are numbered by workflow stage. Each stage retains its own dependency specification because the archived environments differ.

The publication dataset contains 1,546 drugs. The primary ranking uses STRING score >= 700; the sensitivity ranking uses >= 900. The archived Z-scores, p-values, and ranks are retained.

## Workflow

| Folder | Purpose | Main entry point |
| --- | --- | --- |
| [01_human_to_canine_orthologs](01_human_to_canine_orthologs/README.md) | Human drug-target symbols to canine one-to-one orthologs | `map_drug_targets_to_dog_biomart.py` |
| [02_canine_string_mapping](02_canine_string_mapping/README.md) | Canine gene symbols to STRING protein IDs | `map_dog_targets_to_string.py` |
| [03_disease_module_lcc](03_disease_module_lcc/README.md) | Disease-module connectivity and LCC significance | `canine_neuralgia_LCC_netmedpy.py` |
| [04_network_proximity](04_network_proximity/README.md) | Drug-target proximity to the disease set | `run_proximity.py` |
| [05_GSE224814_differential_expression](05_GSE224814_differential_expression/README.md) | Rat trigeminal-ganglion RNA-seq reanalysis, counts, DEGs, and QC | `scripts/run_analysis.py` |
| [06_publication_figures](06_publication_figures/README.md) | Ranking, stability, clinical-reference, and target-profile plots | `figure_scripts/` |

The ortholog output matches the STRING-mapping input exactly. The final STRING-mapping workbook is used by both the proximity analysis and the figure scripts. It contains 4,919 target rows for 1,546 drugs.

The LCC analysis and proximity analysis answer different questions. LCC evaluates connectivity among disease proteins. Proximity uses the complete supplied disease union, subject to background-network coverage exclusions; it does not substitute the disease-induced LCC. Figure generation uses the archived primary 700 ranking, with the 900 ranking used for stability comparisons.

The GSE224814 folder contains a separate rat differential-expression analysis. It
includes six-sample raw counts, source metadata, archived results, and figures.
Its 431 DEGs use `padj < 0.05` and `abs(log2FoldChange) >= 1`. It does not perform
rat-to-dog mapping or establish the provenance of the complete disease-ID union.

## Setup and execution

Enter the relevant numbered directory before following its README. Use a separate virtual environment for each stage; do not combine the different requirements into one installation. The root directory has no global requirements file for this reason.

The original human input to the ortholog stage was intentionally excluded from the supplied folder. A matching private copy is required to rerun that first stage; its checksum is recorded in the archived metadata. All supplied archived data and figures are retained.

Use fresh output directories when rerunning analyses. Live BioMart and STRING queries can change mappings, and the historical LCC hash seed is unknown. The stage-specific documentation explains these limits and identifies missing historical scripts.

## Verification

From the repository root, check all recorded file hashes with the standard-library verifier:

```bash
python3 tools/verify_repository.py
```

Run the existing regression tests in the appropriate scientific environments:

```bash
cd 03_disease_module_lcc
python -m unittest discover -s tests -v
```

```bash
cd 04_network_proximity
python -m unittest discover -s tests -v
```

These two examples start from the repository root. File-integrity verification confirms the files match the recorded snapshots. [DATA_SHA256SUMS](DATA_SHA256SUMS) covers the supplied scientific data and rendered outputs. [CODE_SHA256SUMS](CODE_SHA256SUMS) identifies the Python and shell files. See [REPRODUCIBILITY.md](REPRODUCIBILITY.md) for the scope of reproduction supported by the supplied files.

The GSE224814 saved-data checks use the Python standard library:

```bash
cd 05_GSE224814_differential_expression
python -m unittest discover -s tests -v
```

## Naming, citation, and release

[RENAMING.md](RENAMING.md) maps the original folder and LCC archive names to their new locations. Main scientific script filenames were retained so their imports, recorded commands, and output conventions remain recognizable.

Citation metadata for the collection are provided in [CITATION.cff](CITATION.cff); individual stages retain their own metadata. The public repository is [tnigen/canine-neuralgia-network-analysis](https://github.com/tnigen/canine-neuralgia-network-analysis). Version `1.0.0` is the configured publication-package version; a release tag has not yet been created. Cite an immutable release or commit in the thesis.

The original software is covered by the [MIT License](LICENSE). Third-party software, database content, and fonts retain their own terms and attribution requirements. Clinical-reference source records and the upstream disease-set curation are not supplied in this collection.
