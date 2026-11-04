# Test coverage to target (#75)

Approved in chat, 2026-11-03. Audit findings A3-02 to A3-20
(`dev/audits/2026-10-deep-review.md`). Version 1.11.13.

## Starting point

On devel 1.11.12, `R/` is at 90.9% (1551/1707 lines) and the overall
figure on codecov's denominator is 85.1% (1740/2045), so the numeric
targets (85% and 80%) are already met. This change covers the rest of
A3-02 (every exported method and `stop()` path tested) and the test
quality findings.

## Changes

1. **`.covrignore`** (A3-03): list the generated files that `codecov.yml`
   ignores (`R/RcppExports.R`, `R/stanmodels.R`, `src/RcppExports.cpp`,
   `src/stanExports_shrinkage.cc`, `src/stanExports_shrinkage.h`), so
   `make coverage` and codecov agree. Add `^\.covrignore$` to
   `.Rbuildignore`.
2. **`.hasScater()`** (A3-14): a one-line internal helper wrapping
   `requireNamespace("scater", quietly = TRUE)`, called by
   `.decontxInitializeZLegacy()`, so the scater-missing error can be
   tested with a mock. No behavior change.
3. **New tests for paths reachable from exported functions:**
   - `retrieveFeatureIndex()`: NULL rownames; `by` as a `rowData` column
     (found and missing) and as a data.frame column (found and missing);
     a partial match hitting several rows (warning, first match kept); no
     match with either `exactMatch` setting (errors); `removeNA = TRUE`.
   - `plotDecontXContamination()`: input without estimates (error);
     `batch` naming a batch, and an unknown batch (error).
   - Marker plots: `markers` or `groupClusters` not a named list;
     duplicated `groupClusters`; SCE without `decontX_clusters`; a `z`
     column not in `colData`; `z` of the wrong length; `log1p = TRUE`
     (values and axis label); `plotDots = TRUE` (a jitter layer);
     `plotDecontXMarkerExpression()` without groups (faceted by cluster).
   - `plotBoxByCluster()` with `file = NULL` returns the plot; assert its
     data.
   - `decontX()`: `seed = NULL` gives the result of `seed = s` run inside
     `withr::with_seed(s, ...)`; `logfile` writes the log to the file and
     nothing to the console; background without column names warns; an
     SCE whose assays have no names gives "It has no named assays";
     `varGenes = NULL` uses 5000.
   - `simulateContamination()`: `seed = NULL` matches `seed = s` inside
     `withr::with_seed(s, ...)`; a single `delta` value.
4. **Internal guards no exported call reaches**, tested by direct calls:
   `.logMessages()` with a non-string `logfile`, `.checkCountsDecon()`
   without dimensions, `.processCellLabels()` with a `z` of the wrong
   length, and the 22 C++ guards in `DecontX.cpp` through
   `decontXEM()`, `decontXLogLik()`, `decontXInitialize()` and
   `calculateNativeMatrix()` (A3-18).
5. **Hollow tests replaced** (A3-04):
   - "Testing DecontX on counts matrix" and "Testing DecontX on SCE":
     build each plot and assert on its data (for example, the
     contamination plot's colour values are the contamination).
   - Background runs: assert each `eta` column equals the normalized row
     sums of the background (plus the pseudocount), not only dimensions.
   - decontPro: a Stan-free test feeds `.process_stan_vb_out()` a
     hand-built fit whose rates split each count in a hand-checked ratio,
     and checks the parameter names and the `delta` transpose. This, not
     the sum-to-input identity, shows the decomposition is right.
   - `legacyInit`: with `z` given, `legacyInit = TRUE` gives the same
     contamination and decontaminated counts as the default, since it
     only changes the UMAP; the no-`z` run still checks the labels and
     ranges. Test names say what they check.
   - The unused `res` in the Seurat test goes.
6. **Stan fits merged** (A3-09): the sum-to-input, retained-draws and
   Pareto-k tests share one fit (7 real fits become 5).
7. **Warnings and noise** (A3-08, A3-15): Seurat fixtures are built from a
   `dgCMatrix`, removing the suite's one warning. The redundant
   `skip_if_not_installed("SeuratObject")` calls go (SeuratObject is in
   Imports) and AGENTS.md's Tests section says so. `decontX()` calls pass
   `verbose = FALSE` unless the log is under test. A helper in
   `helper-decontx.R` captures Stan's console output around real fits.
   The package's own output doesn't change.
8. **Fixtures** (A3-16, withr part): the top-level `deconSim` in
   `test-decon.R` becomes a helper function, so no test depends on state
   from another. No new snapshots: `regexp =` already pins the wording,
   and snapshots are skipped on Bioconductor builders.
9. **`maxIter`** (A3-20): tests that assert results drop `maxIter = 2`
   and run the EM to convergence; tests of errors and argument passing
   keep it.

## Left uncovered on purpose

Code that I-d deletes: the `GET_FUN`/`SET_FUN` closures (A3-13), the
HDF5 branch (A5-06; HDF5Array isn't in Suggests), the old scrapper
`runUmap` branch (A3-11), and the compiled code going in #78
(`matrixNorm.cpp`, `matrixSumsSparse.cpp`, `matrixSums.c`,
`R/matrixSums.R`). `test-matrixSums.R` is left as it is for the same
reason.

## Amendments during implementation

- The no-groups branch of `plotDecontXMarkerExpression()`
  (`plot_decontx.R:449-457`) can't be reached: the input processing
  always fills in `groupClusters`, one group per cluster. That is A3-12,
  removed in I-d. The new test checks the plot without groups (one group
  per cluster) and the branch stays uncovered.
- The sparse and DelayedArray equivalence tests keep `maxIter = 2`, as
  their comment explains: the EM is deterministic, so two iterations
  are enough to compare input types.
- The C++ guards were checked by coverage (each guard line runs only when
  its condition holds, and each test asserts its message), not by
  removing them, since a missing guard can segfault.

## Checks

Each new test of an untested path must fail on its assertion when the
path is broken (checked by temporarily breaking it, not committed).
`make test` passes with no warnings; `make coverage` reported per file,
before and after; `make check-full` and `make bioccheck` before the PR.

## Results

Results, plots and messages don't change; confirm with
`dev/compare-commits.R`. Version 1.11.13, with no NEWS section, as for
the test-only 1.11.4.

## Related packages

None affected: tests and one internal helper only.
