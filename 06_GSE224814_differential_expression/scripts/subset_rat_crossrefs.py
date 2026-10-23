#!/usr/bin/env python3
"""Extract rat rows from an NCBI gene2ensembl snapshot without altering fields."""

import argparse
import gzip
from pathlib import Path


def subset(source: Path, output: Path) -> int:
    count = 0
    with gzip.open(source, "rb") as reader, output.open("xb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as writer:
            header = reader.readline()
            if not header.startswith(b"#tax_id\t"):
                raise ValueError("Unexpected gene2ensembl header")
            writer.write(header)
            for line in reader:
                if line.split(b"\t", 1)[0] == b"10116":
                    writer.write(line)
                    count += 1
    return count


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    count = subset(args.source, args.output)
    print(f"Retained {count} rat cross-reference rows.")


if __name__ == "__main__":
    main()
