# Naming changes

| Original folder | Current folder |
| --- | --- |
| `target ortholog/` | `01_human_to_canine_orthologs/` |
| `drug target/` | `02_canine_string_mapping/` |
| `Neuralgia_LCC/` | `03_disease_module_lcc/` |
| `proximity/` | `04_network_proximity/` |
| `figures/` | `05_publication_figures/` |
| `GSE224814/GSE224814_analysis/` | `06_GSE224814_differential_expression/` |

LCC archives were also reorganized:

| Original location within the LCC folder | Current location |
| --- | --- |
| Root report, results, and missing-ID files ending in ` 900`; root `environment_versions.txt` | `results/score_900/` |
| `LCC 0.7 cut off/` | `results/score_700/` |
| `LCC no cut off/` | `results/unfiltered/` |

Each LCC archive now uses `canine_neuralgia_LCC_report.txt`, `canine_neuralgia_LCC_results.csv`, `canine_neuralgia_missing_from_interactome.txt`, and `environment_versions.txt`. Numeric suffixes and spaces were removed from those filenames because the containing directory identifies the analysis setting. File bytes are unchanged, and the checksum manifest was updated to the new paths.

Other scientific input, output, and script basenames are retained. Recorded commands use stage-relative paths. Configuration files identify archived statistics without assigning new execution dates.
