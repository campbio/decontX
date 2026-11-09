# Remove dead compiled code and migrate off MappedSparseMatrix (#78)

Approved in chat 2026-11-09. Audit findings D2-01, D2-04 to D2-07, D2-10,
D2-11, D2-13, D2-14, A2-02 and A7-16. D2-15 (the eta loop bound) is already
fixed on devel: the loop runs over `new_eta.ncol()`.

## Problem

- `src/matrixSumsSparse.cpp` (`colSumByGroupSparse`, `rowSumByGroupSparse`
  and the two `*ChangeSparse` functions, which change their input in place
  and don't check dimensions), `fastNormPropLog`, `fastNormPropSqrt` and
  `nonzero` in `src/matrixNorm.cpp`, and the numeric
  `_colSumByGroup_numeric` in `src/matrixSums.c` (segfaults on NA) have no
  caller in decontX, celda or singleCellTK.
- `.colSumByGroup()` has one caller, `simulateContamination()`, on an
  integer matrix. Its C code uses `int` offsets that overflow at
  G x C >= 2^31.
- `Eigen::MappedSparseMatrix` is deprecated and removed at Eigen HEAD.
- The C++ `delta` message says "positive integers"; the check is
  `delta < 0`.
- Unexported helpers have man pages: `fastNormProp`, `fastNormPropLog`,
  `fastNormPropSqrt`, `nonzero`, `.call_stan_vb`, `.process_stan_vb_out`.
- `inst/stan/include/license.stan` (the rstantools skeleton default) puts a
  GPL-3 notice for decontX into `src/stanExports_shrinkage.h`; DESCRIPTION
  says MIT.
- The `RcppParallel` import has no comment saying why it must stay.

## Change

- Delete `src/matrixSumsSparse.cpp`, `src/matrixSums.c`, `R/matrixSums.R`,
  and the three dead functions in `src/matrixNorm.cpp`.
- `simulateContamination()`: `t(rowsum(t(cellRmat), z))` replaces
  `.colSumByGroup()`. `z` is always 1..K with every group present, so the
  columns line up; `fastNormProp()` drops the dimnames `rowsum()` adds, so
  `eta` is identical.
- `src/DecontX.cpp`: `Eigen::Map<Eigen::SparseMatrix<double> >` replaces
  `Eigen::MappedSparseMatrix<double>`.
- `delta` message: "non-negative values".
- `@noRd` on the remaining C++ helper and the two Stan helpers.
- `license.stan`: an MIT notice with the holders from `LICENSE`; the
  regenerated Stan header then matches DESCRIPTION.
- A "do not remove" comment on `importFrom(RcppParallel, ...)`.
- Regenerate `R/RcppExports.R`, `src/RcppExports.cpp`, `man/`, `NAMESPACE`
  and `src/stanExports_shrinkage.h`.

## Tests

- Test first: the existing C++ `delta` check in `test-decon.R` matches
  "non-negative"; it fails on the old message.
- `test-simulateContamination.R` pins `eta` to the formula written out by
  hand.
- `test-matrixSums.R` is deleted: it only tests the deleted code.

## Docs and checks

- NEWS (1.11.25): unused compiled code and internal help pages removed;
  results unchanged.
- AGENTS.md Layout and Tests sections; ROADMAP marks #78 done.
- No ADR: no DESCRIPTION dependency changes.
- `dev/compare-commits.R` against devel: all paths identical. Then
  `make test`, `make coverage`, `make check-full`, `make bioccheck`.
