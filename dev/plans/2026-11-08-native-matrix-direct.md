# Direct native-count matrix (#80, part 2)

Approved in chat 2026-11-08. Audit findings A5-04 and D2-09; the 1e-14
tolerance was agreed in the 2026-10 audit.

## Problem

`calculateNativeMatrix()` in `src/DecontX.cpp` builds the decontaminated
counts once per batch, after the EM. It:

1. computes the native probability through an exp/log round trip, where
   the EM uses the direct product (the audit measured a max relative
   difference of 7.7e-15 between the two);
2. writes every nonzero with `coeffRef(i, j)`, a binary search in the
   column;
3. copies the counts into an Eigen matrix, and Rcpp copies the result
   again (`x`, `i` and `p`) when it returns it, although only `x` changes.

## Change

- C++: `calculateNativeMatrix()` becomes `calculateNativeCounts()`, which
  returns only the new `x` values (length nnz, in the order of the counts'
  `x` slot). It uses the EM's direct formula and raw pointers, walks the
  columns in storage order (no `coeffRef()`), and keeps the same guards.
- R (`.decontXbatch()`): copy the counts, replace `@x` with the result and
  set `@Dimnames` to `list(NULL, NULL)`. `@<-` copies the object shallowly,
  so `i` and `p` are shared, not copied. Clearing the dimnames keeps the
  result as before (unnamed; `.inCellOrder()` names it).
- `R/RcppExports.R` and `src/RcppExports.cpp` regenerated; the profile
  phase label in `dev/bench/bench-decontx.R` follows the new name.
- Out of scope: `decontXInitialize()` and a shared pointer helper.

## Results

Only `decontXcounts` changes, by at most about 1e-14 relative.
Contamination, estimates, `z` and the sparsity pattern (`i`, `p`, explicit
zeros) stay identical. Checked with `dev/compare-commits.R` (DIFFERS
expected only in `decontXcounts` and the `hdf5` path) and a scratch script
that reports the max relative difference of `decontXcounts` between
installed devel and the branch on the reference paths, sim 100k and
PBMC4k. `tests/testthat/fixtures/decontx-refs.rds` is regenerated with
`fixtures/make-refs.R` (the tests compare at 1e-10, so they pass either
way).

## Tests

- `test-decon.R` tests that call `calculateNativeMatrix()` directly switch
  to `calculateNativeCounts()`, including the C++ guard test.
- New: `decontXcounts` equals a dense R oracle of the direct formula at
  1e-12. It passes on devel too: this is a speedup within tolerance, with
  no behaviour change for a test to catch first.

## Benchmark

`dev/bench/bench-decontx.R`, run `z_supplied_noembed`, sim 10k/50k/100k
and PBMC4k; time, R-allocated memory and peak RSS. Baseline: the part 1
install (6d28421, same code as devel 3df858e), rerun in the same session.

## Other

Version 1.11.21; one NEWS line. No doc or vignette changes. Values move in
the last digits only, so celda and singleCellTK aren't affected in
practice.

## Measurements

(filled in after the benchmark)
