# Unchecked indexing in the EM loop (#80, part 1)

Approved in chat 2026-11-08. Audit finding A5-03.

## Problem

`decontXEM()` in `src/DecontX.cpp` loops over every nonzero count on every
EM iteration. Two of its accessors are expensive in Rcpp 1.1.2:

- `new_phi(i, k) += px`: `Matrix::operator()(i, j)` calls `offset()`, which
  runs `Rf_isMatrix()`, reads the `dim` attribute and checks both indices
  for every nonzero.
- `native_total(j) += px`: `Vector::operator()(j)` calls `Rf_xlength()` and
  checks the index for every nonzero.
- The `[]` accesses (`phi[...]`, `eta[...]`, `theta[j]`, `z[j]`) do a size
  comparison that can warn.

The checks protect nothing there. The guards at the top of the function
already validate the sizes of `theta`, `z`, `phi` and `eta` and the range of
`z`, and Eigen's sparse iterator yields only row indices below
`counts.rows()`. `decontXLogLik()` has the same loop with only `[]` checks;
it runs every `iterLogLik` iterations.

## Change

`src/DecontX.cpp` only:

- In the loops of `decontXEM()` and `decontXLogLik()`, read and write
  through raw pointers taken once before the loop.
- Per column, load `k = z[j] - 1` and `theta[j]` once, and accumulate
  `native_total[j]` in a local `double` stored after the column (same
  terms, same order, starting from 0).
- Every floating-point expression stays as it is, so the results are
  bit-identical (precomputed `theta + pseudocount` terms or ratio matrices
  are not, per the audit). Guards and R code are unchanged.
- Out of scope: `decontXInitialize()` (runs once) and
  `calculateNativeMatrix()` (#80 part 2, A5-04).

## Tests

No new tests: no behaviour changes, so no assertion could fail first. The
regression net covers it: the committed references in the suite, and
`Rscript dev/compare-commits.R ddaed80 HEAD` must be identical on every
path. The existing C++ guard tests stay.

## Benchmark

Installed builds, devel ddaed80 against the branch, run type
`z_supplied_noembed` (so the embedding doesn't hide the EM): sim mode at
10k, 50k and 100k cells, and real mode on PBMC4k. Results below.

## Other

Version 1.11.20; one NEWS line. No vignette or doc changes. Results don't
change, so celda and singleCellTK aren't affected.

## Results

(filled in after the benchmark)
