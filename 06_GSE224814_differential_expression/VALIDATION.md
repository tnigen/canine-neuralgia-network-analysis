# Publication checks

Checked on 2026-10-23; the archived scientific outputs were not replaced.

- All 15 offline tests passed, including agreement of all six featureCounts
  files with the raw matrix, low-count filtering, DEG selection, directional
  subsets, annotation membership, and VST sample distances.
- A count-based rerun with four threads completed and passed all 18 built-in
  analysis checks. It recovered the same 431 DEG identifiers and directions
  (368 up, 63 down), with the same annotation values.
- An independent rerun of the unmodified original analysis script in the same
  environment produced byte-identical differential-expression CSVs to the
  prepared code. The compact rat reference did not change annotation results.
- Refitted statistics are not all numerically identical to the archived run,
  including with the unmodified original script. Numerical identity to that
  historical run is not claimed. The archived numerical tables are unchanged.
- Dependency specifications and every archived PNG/PDF are byte-identical to
  the supplied files. Changes to scientific data files are limited to removing
  local machine prefixes from headers and one software-version record.
- All 36,316 rat cross-reference rows match the original global snapshot byte
  for byte, including their order. The complete rat gene-information snapshot
  and compressed Ensembl GTF are unchanged.

The raw FASTQs were not downloaded or realigned during publication preparation.
The stage checksum manifests identify the distributed snapshots; the repository
verifier checks their integrity separately from scientific reruns.
