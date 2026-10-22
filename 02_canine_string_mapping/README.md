# Canine drug targets to STRING protein identifiers

This repository converts canine drug-target gene symbols to STRING protein identifiers using the STRING identifier-mapping API. Original input rows and columns are retained, with a `string_protein_id` column inserted after `dog_target_protein`.

## Installation

The supplied environment used Python 3.9.6. Create a virtual environment and install the direct dependencies:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

On Windows, activate the environment with `.venv\Scripts\activate`. The complete package versions captured from the supplied environment are in `requirements-lock.txt`; use that file instead of `requirements.txt` to reproduce those package versions. The captured environment uses urllib3 2.x and requires a compatible OpenSSL-backed Python build.

## Input and execution

Run the program from the directory containing `drug_target_dog.xlsx`:

```bash
python map_dog_targets_to_string.py
```

The first worksheet must contain `drug` and `dog_target_protein`. Additional columns are retained. The input must not already contain `string_protein_id`. File names are configured by constants at the top of the script; there are no command-line options.

The input column contains canine gene symbols. This step performs identifier mapping after ortholog conversion; it does not perform sequence alignment or ortholog inference.

Running the script writes the three output files listed below and replaces existing files with those names. To retain the archived results, run a copy of the script and input workbook in a separate directory.

## Mapping rules

Whitespace is trimmed in a separate lookup representation, while the original input values are retained. Each unique non-empty symbol is submitted once, in batches of 100, with dog taxonomy ID `9615`. Requests use connection and response timeouts of 10 and 60 seconds, respectively, with a one-second pause between batches.

Candidates must have taxonomy ID `9615` and an identifier beginning with `9615.ENSCAFP`. One unique valid identifier is accepted. If multiple identifiers remain, the script accepts a single candidate whose preferred name exactly matches the trimmed input symbol. Remaining ambiguities are reported as errors with a blank identifier. See [METHODS.md](METHODS.md) for response handling and summary definitions.

## Included files

| File | Purpose |
| --- | --- |
| `map_dog_targets_to_string.py` | Identifier-mapping pipeline |
| `drug_target_dog.xlsx` | Input drug-target relationships |
| `drug_target_dog_string.xlsx` | Input relationships with STRING identifiers |
| `drug_target_dog_string_mapping_qc.xlsx` | One mapping status per unique non-empty symbol |
| `drug_target_dog_string_unmapped.xlsx` | Unresolved symbols and errors, with details |
| `requirements.txt` | Direct dependency versions |
| `requirements-lock.txt` | Captured complete environment package versions |
| `INPUT_SHA256SUMS`, `RESULT_SHA256SUMS` | Checksums of the archived workbooks |
| `RUN_PROVENANCE.md` | Archived data counts and provenance limitations |

On macOS or Linux, verify the archived files with:

```bash
shasum -a 256 -c INPUT_SHA256SUMS
shasum -a 256 -c RESULT_SHA256SUMS
```

## Reproducibility

The script retains the original unversioned STRING endpoint. Live database updates can change mappings. The STRING release and raw responses used for the archived run were not recorded, so the archived files and their checksums are the reference for the reported results. The captured environment describes the supplied installation, not a verified historical environment for the archived run.

STRING recommends version-specific API addresses for reproducible future analyses; changing the endpoint should be treated as a new analysis. Refer to the [STRING API documentation](https://string-db.org/help/api/).

## Citation and license

Software citation metadata are provided in [CITATION.cff](CITATION.cff), using release version `1.0.0`. The Python code is available under the [MIT License](LICENSE). STRING-derived data remain subject to the source database's terms and attribution requirements; the software license does not replace those terms.
