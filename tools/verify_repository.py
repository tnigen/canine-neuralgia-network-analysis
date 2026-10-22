#!/usr/bin/env python3
"""Verify recorded SHA-256 snapshots without third-party dependencies."""

import hashlib
from pathlib import Path


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    manifests = sorted(root.rglob('*SHA256SUMS'))
    failures = []
    checked = 0
    for manifest in manifests:
        for number, line in enumerate(manifest.read_text(encoding='utf-8').splitlines(), 1):
            if not line.strip() or line.startswith('#'):
                continue
            try:
                expected, relative = line.split('  ', 1)
                target = (manifest.parent / relative).resolve()
                target.relative_to(root)
                if len(expected) != 64 or not target.is_file():
                    raise ValueError('invalid checksum or missing file')
                digest = hashlib.sha256()
                with target.open('rb') as handle:
                    for block in iter(lambda: handle.read(1024 * 1024), b''):
                        digest.update(block)
                if digest.hexdigest() != expected:
                    raise ValueError('checksum mismatch')
                checked += 1
            except (ValueError, OSError) as exc:
                failures.append(f'{manifest.relative_to(root)}:{number}: {exc}')
    for failure in failures:
        print(f'FAIL: {failure}')
    print(f'Checked {checked} file references across {len(manifests)} manifests; {len(failures)} failure(s).')
    print('These are file-integrity checks. See REPRODUCIBILITY.md for reproduction requirements.')
    return 1 if failures else 0


if __name__ == '__main__':
    raise SystemExit(main())
