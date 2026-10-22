# Repository contents

## Analysis

- `run_proximity.py`: NetMedPy wrapper, input validation, and QC pipeline.
- `run_aws_none.sh`: launcher for the unfiltered analysis.
- `requirements.txt`: pinned Python dependencies.
- `tests/test_run_proximity.py`: parser, QC, integration, and reproducibility
  tests.

## Inputs

- `9615.protein.links.v12.0.txt.gz`: canine STRING v12 interaction network.
- `canine_neuralgia_union.txt`: canine neuralgia disease-protein union.
- `drug_target_dog_string.xlsx`: canine drug-target mapping.
- `INPUT_SHA256SUMS`: SHA-256 checksums for the three inputs.

## Documentation and results

- `README.md`: installation, analysis commands, outputs, and citation guidance.
- `METHODS.md`: scientific and implementation methods.
- `AWS_RUNBOOK.md`: unfiltered analysis on AWS EC2.
- `CITATION.cff`: software citation metadata.
- `LICENSE`: MIT license for the original software and documentation.
- `results/score_700` and `results/score_900`: published tables, null
  distributions, parameters, and QC outputs.

Large derived distance matrices, local virtual environments, validation runs,
logs, and deployment archives are not part of the repository.
