# Keep the decontX result in memory for HDF5 input (#76, PR 4 of 5)

Agreed in chat on 2026-11-06. ADR-0004
(`dev/adr/0004-container-agnostic-core.md`, "HDF5 input"), audit finding
A5-06 (`dev/audits/2026-10-deep-review.md`), decided on 2026-10-19.
Version 1.11.17. Part of #76.

## Today

`.asInputClass()` wraps the result as a DelayedMatrix for DelayedArray
input. When the input is backed by an HDF5 file, it first writes the
whole result to a new HDF5 file with `as(, "HDF5Matrix")`. On pbmc4k
that write took 63% of the run time (10.8 of 17.2 s, audit), and the file
lands in R's temporary folder, so a saved SingleCellExperiment points at
a file that is gone in the next session.

## Changes

- `.asInputClass()` drops its HDF5 branch: any DelayedArray input gets
  `DelayedArray(decontXcounts)`, an in-memory dgCMatrix wrapped as a
  DelayedMatrix. The `seedApply()` check for `HDF5ArraySeed` and the
  `as(, "HDF5Matrix")` coercion go.
- The ".. Converting decontaminated matrix to a DelayedMatrix" log line
  stays.
- The class stays DelayedMatrix and the values are identical.
- HDF5Array is not added to DESCRIPTION (ADR-0004).

## Tests

- Every DelayedArray input now takes one path. The existing test "a
  DelayedMatrix over a dense matrix gives a DelayedMatrix" checks the
  class and the dgCMatrix seed, and the batched DelayedArray log
  snapshot covers the log line.
- The suite can't test HDF5 input without HDF5Array in Suggests, which
  ADR-0004 rules out. `dev/compare-commits.R` gets an `hdf5` path
  instead: it writes the simulated counts to a temporary HDF5 file, runs
  `decontX()` on it, and returns the values as a dgCMatrix, the
  contamination, and whether the result's seed is in memory. The values
  and contamination should be identical between commits and the
  in-memory flag should differ. The path is skipped without HDF5Array,
  as `legacyInit` is without scater.

## Verification

- `make test`, `make coverage` (should rise: the uncovered HDF5 branch
  goes), `make lint`.
- `Rscript dev/compare-commits.R <devel> HEAD`: the 13 existing paths
  identical; `hdf5` differs only in the in-memory flag.
- `make check-full` and `make bioccheck`.
- `dev/bench/bench-decontx.R real` (pbmc4k, `z_estimated` and
  `sce_hdf5`) before and after, numbers in the PR.

## Docs and records

- `decontX()` `@return`: the result is a DelayedMatrix held in memory,
  HDF5 input included; save it to disk explicitly, for example with
  `HDF5Array::saveHDF5SummarizedExperiment()` (as `\code{}`, not a
  link, since HDF5Array isn't declared). `make docs`.
- NEWS for 1.11.17: HDF5-backed input returns an in-memory DelayedMatrix
  instead of one backed by a temporary HDF5 file; faster, needs memory
  for the sparse result; save to disk explicitly (#76).
- DESCRIPTION: 1.11.17.
- No ADR change.

## Related packages

celda and singleCellTK still get a DelayedMatrix for DelayedArray input.
Code that relied on the result being an `HDF5Matrix` (for example
calling `path()` on it) sees the change; that file was temporary anyway.
