#!/usr/bin/env bash
set -euo pipefail

if [[ ! -x .venv/bin/python ]]; then
  echo "ERROR: .venv is missing. Follow AWS_RUNBOOK.md section 3 first." >&2
  exit 2
fi

for required_file in \
  9615.protein.links.v12.0.txt.gz \
  canine_neuralgia_union.txt \
  drug_target_dog_string.xlsx; do
  if [[ ! -f "$required_file" ]]; then
    echo "ERROR: required input is missing: $required_file" >&2
    exit 2
  fi
done

processors="${N_PROCESSORS:-28}"

exec .venv/bin/python run_proximity.py \
  --network 9615.protein.links.v12.0.txt.gz \
  --disease canine_neuralgia_union.txt \
  --drug-targets drug_target_dog_string.xlsx \
  --score-cutoff none \
  --output results/score_none \
  --seed 42 \
  --n-processors "$processors"
