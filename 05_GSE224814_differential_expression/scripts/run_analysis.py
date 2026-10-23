#!/usr/bin/env python3
"""Run the GSE224814 input audit, counting, and differential-expression analysis."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tarfile
from pathlib import Path


URLS = {
    "GSE224814_RAW.tar": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE224nnn/GSE224814/suppl/GSE224814_RAW.tar",
    "GSE224814_family.soft.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE224nnn/GSE224814/soft/GSE224814_family.soft.gz",
    "SraRunInfo.csv": "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=PRJNA932751",
    "ENA_fastq_manifest.tsv": "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA932751&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes&format=tsv",
}

DEFERRED_URLS = {
    "data_raw/reference/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz": "https://ftp.ensembl.org/pub/release-104/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz",
    "data_raw/annotation/gene2ensembl.gz": "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz",
    "data_raw/annotation/Rattus_norvegicus.gene_info.gz": "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Rattus_norvegicus.gene_info.gz",
}


def run(cmd: list[str], project: Path, log_name: str) -> None:
    log_path = project / "results" / "logs" / log_name
    log_path.parent.mkdir(parents=True, exist_ok=True)
    print("RUN:", " ".join(cmd), flush=True)
    with log_path.open("w") as log:
        process = subprocess.Popen(cmd, cwd=project, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        assert process.stdout is not None
        for line in process.stdout:
            print(line, end="", flush=True)
            log.write(line)
        rc = process.wait()
    if rc:
        raise subprocess.CalledProcessError(rc, cmd)


def download_sources(project: Path) -> None:
    raw = project / "data_raw"
    raw.mkdir(parents=True, exist_ok=True)
    for filename, url in URLS.items():
        target = raw / filename
        if target.exists() and target.stat().st_size > 0:
            continue
        partial = target.with_suffix(target.suffix + ".partial")
        subprocess.run(
            [
                "curl", "-L", "--fail", "--retry", "20", "--retry-all-errors",
                "--retry-delay", "5", "--connect-timeout", "30", "--continue-at", "-",
                "-o", str(partial), url,
            ],
            check=True,
        )
        partial.replace(target)
    extracted = raw / "GSE224814_RAW_extracted"
    expected = {f"GSM70335{i}_{label}.txt.gz" for i, label in zip(range(26, 32), ["A1", "A2", "A3", "B1", "B2", "B3"], strict=True)}
    if not extracted.exists() or {p.name for p in extracted.glob("*.gz")} != expected:
        extracted.mkdir(parents=True, exist_ok=True)
        with tarfile.open(raw / "GSE224814_RAW.tar") as archive:
            members = archive.getmembers()
            names = {Path(m.name).name for m in members}
            if names != expected or any(not m.isfile() for m in members):
                raise ValueError(f"Unexpected archive contents: {sorted(names)}")
            archive.extractall(extracted, filter="data")


def download_deferred_sources(project: Path) -> None:
    """Restore archival reference/annotation inputs after disk-heavy counting."""
    for relative_path, url in DEFERRED_URLS.items():
        target = project / relative_path
        if target.exists() and target.stat().st_size > 0:
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        partial = target.with_suffix(target.suffix + ".partial")
        subprocess.run(
            [
                "curl", "-L", "--fail", "--retry", "20", "--retry-all-errors",
                "--retry-delay", "5", "--connect-timeout", "30", "--continue-at", "-",
                "-o", str(partial), url,
            ],
            check=True,
        )
        partial.replace(target)


def prepare_bundled_references(project: Path) -> None:
    """Expand the bundled GTF and check the rat cross-reference snapshot."""
    from build_raw_counts import decompress_once

    decompress_once(
        project / "data_raw/reference/Rattus_norvegicus.Rnor_6.0.104.gtf.gz",
        project / "data_processed/reference/Rattus_norvegicus.Rnor_6.0.104.gtf",
    )
    for name in ("gene2ensembl.rat.tsv.gz", "Rattus_norvegicus.gene_info.gz"):
        if not (project / "data_raw/annotation" / name).is_file():
            raise FileNotFoundError(f"Missing bundled annotation: {name}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--threads", type=int, default=min(8, os.cpu_count() or 1))
    parser.add_argument("--keep-fastq", action="store_true", help="Keep ~20 GB of verified gzipped FASTQs after counting")
    parser.add_argument("--rebuild-counts", action="store_true")
    parser.add_argument(
        "--from-counts", action="store_true",
        help="Use bundled raw counts and rat annotations without downloading reads or an alignment index",
    )
    args = parser.parse_args()
    if args.from_counts and args.rebuild_counts:
        parser.error("--from-counts and --rebuild-counts cannot be combined")
    project = args.project.resolve()
    count_matrix = project / "data_processed" / "GSE224814_raw_count_matrix.csv"
    if args.from_counts and not count_matrix.is_file():
        parser.error("--from-counts requires the supplied raw-count matrix")
    for subdir in [
        "data_raw", "data_processed", "scripts", "results/differential_expression", "results/qc",
        "results/figures", "results/tables", "results/logs",
    ]:
        (project / subdir).mkdir(parents=True, exist_ok=True)
    download_sources(project)
    python = sys.executable
    run([python, str(project / "scripts" / "audit_inputs.py"), "--project", str(project)], project, "01_audit_inputs.log")
    if args.rebuild_counts and count_matrix.exists():
        raise RuntimeError("Refusing to overwrite existing raw-count matrix; move it aside explicitly before rebuilding")
    if not count_matrix.exists():
        cmd = [python, str(project / "scripts" / "build_raw_counts.py"), "--project", str(project), "--threads", str(args.threads)]
        if args.keep_fastq:
            cmd.append("--keep-fastq")
        run(cmd, project, "02_build_raw_counts.log")
    if args.from_counts:
        prepare_bundled_references(project)
    else:
        from build_raw_counts import prepare_reference

        prepare_reference(project, args.threads)
        download_deferred_sources(project)
    analysis_cmd = [
        python, str(project / "scripts" / "analyze_counts.py"),
        "--project", str(project), "--threads", str(args.threads),
    ]
    if args.from_counts:
        analysis_cmd.append("--from-counts")
    run(
        analysis_cmd,
        project,
        "03_analyze_counts.log",
    )


if __name__ == "__main__":
    main()
