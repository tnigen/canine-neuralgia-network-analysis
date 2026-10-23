#!/usr/bin/env python3
"""Align paired GSE224814 FASTQs to rn6 and count gene-level fragments."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import gzip
import hashlib
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

import pandas as pd


RUN_TO_SAMPLE = {
    "SRR23372014": "A1_TG",
    "SRR23372013": "A2_TG",
    "SRR23372012": "A3_TG",
    "SRR23372011": "B1_TG",
    "SRR23372010": "B2_TG",
    "SRR23372009": "B3_TG",
}
SAMPLE_TO_GSM = {
    "A1_TG": "GSM7033526", "A2_TG": "GSM7033527", "A3_TG": "GSM7033528",
    "B1_TG": "GSM7033529", "B2_TG": "GSM7033530", "B3_TG": "GSM7033531",
}
SAMPLE_ORDER = ["A1_TG", "A2_TG", "A3_TG", "B1_TG", "B2_TG", "B3_TG"]
FASTA_URL = (
    "https://ftp.ensembl.org/pub/release-104/fasta/rattus_norvegicus/dna/"
    "Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz"
)
GTF_URL = (
    "https://ftp.ensembl.org/pub/release-104/gtf/rattus_norvegicus/"
    "Rattus_norvegicus.Rnor_6.0.104.gtf.gz"
)
GENE2ENSEMBL_URL = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz"
RAT_GENE_INFO_URL = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Rattus_norvegicus.gene_info.gz"
HISAT2_INDEX_URL = "https://genome-idx.s3.amazonaws.com/hisat/rn6_genome.tar.gz"


def run(cmd: list[str], log: Path | None = None) -> None:
    print("RUN:", " ".join(cmd), flush=True)
    if log:
        log.parent.mkdir(parents=True, exist_ok=True)
        with log.open("w") as handle:
            subprocess.run(cmd, check=True, stdout=handle, stderr=subprocess.STDOUT)
    else:
        subprocess.run(cmd, check=True)


def download(url: str, target: Path) -> None:
    if target.exists() and target.stat().st_size > 0:
        print(f"Using existing download: {target}", flush=True)
        return
    target.parent.mkdir(parents=True, exist_ok=True)
    partial = target.with_suffix(target.suffix + ".partial")
    run([
        "curl", "-L", "--fail", "--retry", "20", "--retry-all-errors",
        "--retry-delay", "5", "--connect-timeout", "30", "--continue-at", "-",
        "-o", str(partial), url,
    ])
    partial.replace(target)


def download_ranged(url: str, target: Path, expected_size: int, connections: int = 16) -> None:
    """Resume downloads using fixed byte ranges encoded in each piece's name."""
    if target.exists():
        if target.stat().st_size != expected_size:
            raise RuntimeError(f"Unexpected size for {target}: {target.stat().st_size} != {expected_size}")
        print(f"Using existing download: {target}", flush=True)
        return
    target.parent.mkdir(parents=True, exist_ok=True)
    piece_size = 32 * 1024 * 1024
    n_pieces = (expected_size + piece_size - 1) // piece_size
    parts_dir = target.with_name(target.name + ".parts")
    parts_dir.mkdir(parents=True, exist_ok=True)
    jobs: list[tuple[int, int, Path]] = []
    for i in range(n_pieces):
        start = i * piece_size
        end = min(expected_size - 1, start + piece_size - 1)
        piece = parts_dir / f"{i:05d}.{start}-{end}.part"
        required = end - start + 1
        if piece.exists() and piece.stat().st_size != required:
            piece.unlink()
        jobs.append((start, end, piece))

    def fetch_range(job: tuple[int, int, Path]) -> Path:
        start, end, piece = job
        required = end - start + 1
        if piece.exists() and piece.stat().st_size == required:
            return piece
        cmd = [
            "curl", "-L", "--fail", "--retry", "30", "--retry-all-errors",
            "--retry-delay", "3", "--connect-timeout", "30", "--speed-limit", "10240",
            "--speed-time", "45", "--range", f"{start}-{end}", "-o", str(piece), url,
        ]
        print("RUN:", " ".join(cmd), flush=True)
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if piece.stat().st_size != required:
            piece.unlink(missing_ok=True)
            raise RuntimeError(f"Byte-range size mismatch for {piece}")
        return piece

    with ThreadPoolExecutor(max_workers=min(connections, n_pieces)) as pool:
        pieces = list(pool.map(fetch_range, jobs))
    assembling = target.with_suffix(target.suffix + ".assembling")
    with assembling.open("wb") as out:
        for piece in pieces:
            with piece.open("rb") as src:
                shutil.copyfileobj(src, out, length=8 * 1024 * 1024)
    if assembling.stat().st_size != expected_size:
        raise RuntimeError(
            f"Assembled size mismatch for {assembling}: {assembling.stat().st_size} != {expected_size}"
        )
    assembling.replace(target)
    for piece in pieces:
        piece.unlink()
    parts_dir.rmdir()


def clear_ranged_cache(target: Path) -> None:
    """Remove only transfer artifacts associated with one explicit target."""
    target.unlink(missing_ok=True)
    target.with_suffix(target.suffix + ".assembling").unlink(missing_ok=True)
    parts_dir = target.with_name(target.name + ".parts")
    if parts_dir.exists():
        shutil.rmtree(parts_dir)
    # Clean artifacts written by the superseded offset-relative implementation.
    target.with_suffix(target.suffix + ".partial").unlink(missing_ok=True)
    for artifact in target.parent.glob(target.name + ".range*.partial*"):
        artifact.unlink()


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def stream_fastq_to_fifo(
    url: str,
    expected_size: int,
    expected_md5: str,
    fifo: Path,
    temp_dir: Path,
    label: str,
) -> str:
    """Stream one remote gzipped FASTQ through a FIFO in retryable 64 MiB ranges."""
    digest = hashlib.md5()
    bytes_streamed = 0
    chunk_bytes = 64 * 1024 * 1024
    with fifo.open("wb", buffering=0) as out:
        for part_number, start in enumerate(range(0, expected_size, chunk_bytes)):
            end = min(expected_size - 1, start + chunk_bytes - 1)
            required = end - start + 1
            chunk = temp_dir / f"{label}.{part_number:04d}.range.partial"
            if not chunk.exists() or chunk.stat().st_size != required:
                if chunk.exists():
                    chunk.unlink()
                cmd = [
                    "curl", "-L", "--fail", "--retry", "30", "--retry-all-errors",
                    "--retry-delay", "3", "--connect-timeout", "30", "--speed-limit", "10240",
                    "--speed-time", "45", "--range", f"{start}-{end}", "-o", str(chunk), url,
                ]
                subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            if chunk.stat().st_size != required:
                raise RuntimeError(f"FASTQ range size mismatch for {chunk}")
            with chunk.open("rb") as src:
                for block in iter(lambda: src.read(8 * 1024 * 1024), b""):
                    digest.update(block)
                    out.write(block)
                    bytes_streamed += len(block)
            chunk.unlink()
    observed = digest.hexdigest()
    if bytes_streamed != expected_size or observed.lower() != expected_md5.lower():
        raise RuntimeError(
            f"Stream validation failed for {label}: bytes={bytes_streamed}/{expected_size}, "
            f"md5={observed}/{expected_md5}"
        )
    return observed


def decompress_once(source: Path, target: Path) -> None:
    if target.exists() and target.stat().st_size > 0:
        return
    target.parent.mkdir(parents=True, exist_ok=True)
    partial = target.with_suffix(target.suffix + ".partial")
    with gzip.open(source, "rb") as src, partial.open("wb") as dst:
        shutil.copyfileobj(src, dst, length=8 * 1024 * 1024)
    partial.replace(target)


def prepare_reference(project: Path, threads: int) -> tuple[Path, Path, Path]:
    raw_ref = project / "data_raw" / "reference"
    proc_ref = project / "data_processed" / "reference"
    fasta_gz = raw_ref / "Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz"
    gtf_gz = raw_ref / "Rattus_norvegicus.Rnor_6.0.104.gtf.gz"
    download(GTF_URL, gtf_gz)
    gtf = proc_ref / gtf_gz.name.removesuffix(".gz")
    decompress_once(gtf_gz, gtf)
    # Use the official rn6 index rather than rebuilding it locally.
    index_tar = raw_ref / "rn6_genome.hisat2_index.tar.gz"
    index_dir = proc_ref / "hisat2_prebuilt_rn6"
    index_prefix = index_dir / "genome"
    candidates = [p for p in index_dir.rglob("genome.1.ht2")] if index_dir.exists() else []
    valid_candidates = [
        p.parent / "genome" for p in candidates
        if all((p.parent / f"genome.{i}.ht2").exists() for i in range(1, 9))
    ]
    if not valid_candidates:
        download(HISAT2_INDEX_URL, index_tar)
        index_dir.mkdir(parents=True, exist_ok=True)
        with tarfile.open(index_tar) as archive:
            archive.extractall(index_dir, filter="data")
        candidates = list(index_dir.rglob("genome.1.ht2"))
        if len(candidates) != 1:
            raise RuntimeError(f"Could not identify one HISAT2 index prefix in {index_dir}: {candidates}")
        actual_parent = candidates[0].parent
        index_prefix = actual_parent / "genome"
    elif len(valid_candidates) == 1:
        index_prefix = valid_candidates[0]
    else:
        raise RuntimeError(f"Multiple complete HISAT2 index prefixes found: {valid_candidates}")

    # UCSC uses chr-prefixed names while Ensembl release 104 uses 1..20,X,Y,MT.
    # Make a derived chromosome-only GTF; coordinates and biological IDs are unchanged.
    count_gtf = proc_ref / "Rattus_norvegicus.Rnor_6.0.104.ucsc_chrnames.gtf"
    if not count_gtf.exists() or count_gtf.stat().st_size == 0:
        standard = {str(i) for i in range(1, 21)} | {"X", "Y", "MT"}
        with gtf.open("rt", encoding="utf-8") as src, count_gtf.open("w", encoding="utf-8") as dst:
            for line in src:
                if line.startswith("#"):
                    dst.write(line)
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 9 or fields[0] not in standard:
                    continue
                fields[0] = "chrM" if fields[0] == "MT" else "chr" + fields[0]
                dst.write("\t".join(fields) + "\n")
    return fasta_gz, count_gtf, index_prefix


def count_one_sample(
    project: Path,
    run_accession: str,
    sample_id: str,
    urls: list[str],
    checksums: list[str],
    sizes: list[int],
    index_prefix: Path,
    gtf: Path,
    threads: int,
    keep_fastq: bool,
) -> Path:
    count_dir = project / "data_processed" / "featurecounts"
    count_dir.mkdir(parents=True, exist_ok=True)
    count_file = count_dir / f"{sample_id}.featureCounts.txt"
    if count_file.exists() and count_file.stat().st_size > 0:
        print(f"Using existing count result: {count_file}", flush=True)
        return count_file

    logs = project / "results" / "logs"
    qc = project / "results" / "qc"
    logs.mkdir(parents=True, exist_ok=True)
    qc.mkdir(parents=True, exist_ok=True)
    align_log = logs / f"{sample_id}_{run_accession}_hisat2.log"
    fc_log = logs / f"{sample_id}_{run_accession}_featureCounts.log"
    fc_tmp = Path(tempfile.gettempdir()) / "GSE224814_featureCounts_tmp"
    fc_tmp.mkdir(parents=True, exist_ok=True)
    fastq_dir = project / "data_raw" / "sra_fastq_cache" / run_accession
    fastq_dir.mkdir(parents=True, exist_ok=True)
    remote_urls = ["https://" + u.removeprefix("ftp://").removeprefix("https://") for u in urls]
    def fetch_and_verify(item: tuple[str, str, int]) -> Path:
        url, expected_md5, expected_size = item
        target = fastq_dir / Path(url).name
        download_ranged(url, target, expected_size)
        observed = md5(target)
        if observed.lower() != expected_md5.lower():
            print(
                f"MD5 mismatch for {target}; discarding transfer artifacts and retrying once",
                flush=True,
            )
            clear_ranged_cache(target)
            download_ranged(url, target, expected_size)
            observed = md5(target)
        if observed.lower() != expected_md5.lower():
            clear_ranged_cache(target)
            raise RuntimeError(f"MD5 mismatch for {target}: {observed} != {expected_md5}")
        return target

    with ThreadPoolExecutor(max_workers=2) as download_pool:
        fastqs = list(
            download_pool.map(fetch_and_verify, zip(remote_urls, checksums, sizes, strict=True))
        )
    bam_dir = project / "data_processed" / "alignment_cache"
    bam_dir.mkdir(parents=True, exist_ok=True)
    bam_file = bam_dir / f"{sample_id}.pair_adjacent.bam"
    if bam_file.exists():
        bam_file.unlink()
    summary_file = qc / f"{sample_id}_{run_accession}_hisat2_summary.txt"
    hisat_cmd = [
        "hisat2",
        "-p",
        str(threads),
        "--dta",
        "--seed",
        "20260904",
        "--summary-file",
        str(summary_file),
        "-x",
        str(index_prefix),
        "-1",
        str(fastqs[0]),
        "-2",
        str(fastqs[1]),
    ]
    view_cmd = [
        "samtools", "view", "-@", str(max(1, threads // 2)), "-b",
        "-o", str(bam_file), "-",
    ]
    fc_cmd = [
        "featureCounts",
        "-T",
        str(threads),
        "-p",
        "--countReadPairs",
        "-B",
        "-C",
        "-s",
        "0",
        "-t",
        "exon",
        "-g",
        "gene_id",
        "--tmpDir",
        str(fc_tmp),
        "-a",
        str(gtf),
        "-o",
        str(count_file),
        str(bam_file),
    ]
    print("RUN PIPE:", " ".join(hisat_cmd), "|", " ".join(view_cmd), flush=True)
    with align_log.open("w") as align_err:
        align = subprocess.Popen(hisat_cmd, stdout=subprocess.PIPE, stderr=align_err)
        assert align.stdout is not None
        view = subprocess.Popen(view_cmd, stdin=align.stdout, stdout=subprocess.DEVNULL, stderr=align_err)
        align.stdout.close()
        view_rc = view.wait()
        align_rc = align.wait()
    if align_rc or view_rc or not bam_file.exists():
        if bam_file.exists():
            bam_file.unlink()
        raise RuntimeError(
            f"Alignment failed for {sample_id}: hisat2={align_rc}, "
            f"samtools_view={view_rc}"
        )
    print("RUN:", " ".join(fc_cmd), flush=True)
    with fc_log.open("w") as fc_err:
        fc_rc = subprocess.run(fc_cmd, stdout=fc_err, stderr=subprocess.STDOUT).returncode
    if fc_rc:
        if count_file.exists():
            count_file.unlink()
        if bam_file.exists():
            bam_file.unlink()
        raise RuntimeError(f"featureCounts failed for {sample_id}: return code {fc_rc}")
    if not count_file.exists() or count_file.stat().st_size == 0:
        raise RuntimeError(f"No count file produced for {sample_id}")
    fc_summary = Path(str(count_file) + ".summary")
    if fc_summary.exists():
        shutil.copy2(fc_summary, qc / f"{sample_id}_{run_accession}_featureCounts_summary.txt")
    bam_file.unlink()
    if not keep_fastq:
        for fastq in fastqs:
            fastq.unlink()
        try:
            fastq_dir.rmdir()
        except OSError:
            pass
    return count_file


def merge_counts(project: Path, count_files: dict[str, Path]) -> Path:
    merged = None
    annotation = None
    for sample in SAMPLE_ORDER:
        df = pd.read_csv(count_files[sample], sep="\t", comment="#")
        if "Geneid" not in df.columns:
            raise ValueError(f"Unexpected featureCounts header in {count_files[sample]}")
        count_column = df.columns[-1]
        current = df[["Geneid", count_column]].rename(columns={count_column: sample}).set_index("Geneid")
        if merged is None:
            merged = current
            annotation = df[["Geneid", "Chr", "Start", "End", "Strand", "Length"]].copy()
        else:
            if not merged.index.equals(current.index):
                raise ValueError(f"Gene-row structure differs for {sample}")
            merged = merged.join(current, how="inner")
    assert merged is not None and annotation is not None
    if merged.index.duplicated().any():
        raise ValueError("featureCounts returned duplicated Ensembl gene IDs")
    if (merged.to_numpy() < 0).any() or not (merged.to_numpy() == merged.to_numpy().astype(int)).all():
        raise ValueError("Generated featureCounts matrix is not non-negative integer data")
    merged = merged.astype(int)
    output = project / "data_processed" / "GSE224814_raw_count_matrix.csv"
    merged.to_csv(output, index_label="original_id")
    annotation.to_csv(project / "data_processed" / "featurecounts_gene_coordinates.csv", index=False)
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--threads", type=int, default=min(8, os.cpu_count() or 1))
    parser.add_argument("--keep-fastq", action="store_true")
    args = parser.parse_args()
    project = args.project.resolve()
    for tool in ("curl", "hisat2", "hisat2-build", "samtools", "featureCounts"):
        if shutil.which(tool) is None:
            sys.exit(f"Required executable not found: {tool}")
    _, gtf, index_prefix = prepare_reference(project, args.threads)
    manifest = pd.read_csv(project / "data_raw" / "ENA_fastq_manifest.tsv", sep="\t")
    source_records = []
    for _, row in manifest.iterrows():
        run_accession = row["run_accession"]
        sample_id = RUN_TO_SAMPLE.get(run_accession, "")
        source_records.append(
            {
                "sample_id": sample_id,
                "GSM": SAMPLE_TO_GSM.get(sample_id, ""),
                "SRA_run": run_accession,
                "paired_FASTQ_URLs": row["fastq_ftp"],
                "paired_FASTQ_MD5": row["fastq_md5"],
                "paired_FASTQ_bytes": row["fastq_bytes"],
            }
        )
    pd.DataFrame(source_records).sort_values("sample_id").to_csv(
        project / "results" / "tables" / "raw_read_source_manifest.csv", index=False
    )
    manifest = manifest.set_index("run_accession")
    count_files = {}
    for run_accession, sample_id in RUN_TO_SAMPLE.items():
        row = manifest.loc[run_accession]
        urls = str(row["fastq_ftp"]).split(";")
        checksums = str(row["fastq_md5"]).split(";")
        sizes = [int(x) for x in str(row["fastq_bytes"]).split(";")]
        if len(urls) != 2 or len(checksums) != 2 or len(sizes) != 2:
            raise ValueError(f"Expected exactly two paired FASTQs for {run_accession}")
        count_files[sample_id] = count_one_sample(
            project,
            run_accession,
            sample_id,
            urls,
            checksums,
            sizes,
            index_prefix,
            gtf,
            args.threads,
            args.keep_fastq,
        )
    output = merge_counts(project, count_files)
    print(f"Raw integer count matrix written to {output}", flush=True)


if __name__ == "__main__":
    main()
