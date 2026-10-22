# Running the unfiltered analysis on AWS EC2

This procedure runs the unfiltered network with NetMedPy 0.2.0.1, exact-degree
matching, 1,000 randomizations, and seed 42.

## 1. Provision an instance

Recommended configuration:

- Ubuntu Server 22.04 LTS, x86-64
- `m7i.8xlarge` (32 vCPUs, 128 GiB RAM)
- 50 GiB or more of gp3 EBS storage
- SSH access restricted to the operator's IP address

An `r7i.8xlarge` provides additional memory. On a 64-GiB `c7i.8xlarge`, limit
NetMedPy to 12–14 workers because multiprocessing creates a graph copy for
each worker.

## 2. Obtain the repository

Clone the repository and enter the network-proximity directory:

```bash
git clone <repository-url>
cd <repository-directory>
cd 04_network_proximity
```

Verify the inputs before running the analysis:

```bash
sha256sum -c INPUT_SHA256SUMS
ls -lh 9615.protein.links.v12.0.txt.gz \
    canine_neuralgia_union.txt \
    drug_target_dog_string.xlsx
```

Do not continue if a checksum fails.

## 3. Create the environment

```bash
sudo apt-get update
sudo apt-get install -y python3-venv tmux
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python -m pip check
python -c "import importlib.metadata; print(importlib.metadata.version('netmedpy'))"
python -m unittest discover -s tests -v
```

The version command must print `0.2.0.1`, and all tests must pass.

## 4. Validate the unfiltered inputs

```bash
python run_proximity.py \
    --network 9615.protein.links.v12.0.txt.gz \
    --disease canine_neuralgia_union.txt \
    --drug-targets drug_target_dog_string.xlsx \
    --score-cutoff none \
    --output results/validation_none \
    --seed 42 \
    --dry-run
```

Inspect the network summary:

```bash
grep -E 'score_cutoff|records_passing_score_filter|retained_unique_edges|number_of_nodes|number_of_connected_components' \
    results/validation_none/network_qc.txt
```

Expected network values:

```text
score_cutoff=none
records_passing_score_filter=10768516
retained_unique_edges=5384258
number_of_nodes=19921
number_of_connected_components=1
```

## 5. Run the analysis

Use `tmux` so the process survives an SSH disconnect:

```bash
tmux new -s canine-proximity
source .venv/bin/activate
N_PROCESSORS=28 bash run_aws_none.sh 2>&1 | tee run_score_none.log
```

Detach with `Ctrl-b`, then `d`. Reattach with:

```bash
tmux attach -t canine-proximity
```

For a 64-GiB `c7i.8xlarge`, use:

```bash
N_PROCESSORS=14 bash run_aws_none.sh 2>&1 | tee run_score_none.log
```

The all-pairs shortest-path calculation is the dominant cost. Monitor memory,
storage, and the run log from a separate SSH session:

```bash
tail -f run_score_none.log
free -h
df -h
ps -eo pid,ppid,%cpu,%mem,rss,cmd --sort=-rss | head -25
```

If the process is killed for insufficient memory, retain the partial output
for diagnosis and retry with fewer workers and a new output directory.

## 6. Verify completion

```bash
grep '^analysis_status=complete$' results/score_none/analysis_parameters.txt
test -s results/score_none/proximity_results.tsv
test -s results/score_none/netmedpy_null_distributions.tsv.gz
test -s results/score_none/netmedpy_shortest_path_distances.npz
```

The result table is ranked by ascending Z-score, not raw distance alone.

## 7. Retrieve results

Create an archive on the instance:

```bash
tar -czf score_none_results.tar.gz results/score_none run_score_none.log
sha256sum score_none_results.tar.gz
```

Download the archive, verify its checksum, and then stop or terminate the EC2
instance to avoid further charges.

## Scientific constraints

- Use the supplied disease union; do not substitute a disease-induced LCC.
- Do not infer missing STRING identifiers.
- Treat `none` as disabled filtering, not a numeric cutoff.
- Keep the exact-degree null model and Z-score sign convention.
- Do not interpret network prioritization as evidence of efficacy or safety.
