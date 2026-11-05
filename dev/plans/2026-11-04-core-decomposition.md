# Decompose the decontX core (#76, PR 2 of 5)

Agreed in chat on 2026-11-04. ADR-0004
(`dev/adr/0004-container-agnostic-core.md`), audit findings A4-03 and
A4-06 to A4-08 (`dev/audits/2026-10-deep-review.md`), and the
`bgBatch`-without-`background` error the maintainer agreed to during PR 1
(#100). Version 1.11.15. Part of #76.

## Already done: A7-22

"`z` is duplicated inside `estimates`" was fixed in c908af7, which
removed the second `z =` key from each batch's estimates. What remains,
`estimates[[b]]$z`, is that batch's own labels, which name the `phi` and
`eta` columns, so it is not a duplicate and stays.

## Changes in this PR

Results stay `identical()`. The only behavior change is the new error
below.

### `.decontX()` split into helpers

`.decontX(counts, z, batch, background, bgBatch, params)` keeps its
signature and its order of steps: data checks, background checks, the
"Starting DecontX" banner, background rows, batches, per-batch runs,
assembly, the DelayedMatrix conversion, the closing banner. Its body
becomes calls to these helpers, all in `R/decon.R`, none with defaults:

- `.decontXBatches(batch, bgBatch, background, numCells)`: today's
  batch block (the `"all_cells"` default, the `bgBatch` pairing checks,
  `as.character()`, the `bgBatch` values checks) plus the new error.
  Returns `list(batch, bgBatch, given)`.
- `.checkBatchClusters(z, batch, given, background)`: the per-batch
  one-cluster loop.
- `.decontXLabels(z, batch)`: the user labels without unused levels and
  whether batched labels get the batch prefix. Returns
  `list(labels, prefix)`.
- `.runBatch(counts, z, background, bgBatch, batch, bat, params)`: the
  "Analyzing" log line, the batch subset, the sparse conversion, the
  seeded `.decontXoneBatch()` call and `calculateNativeMatrix()`.
  Returns `list(res, decontXcounts)`.
- `.batchEstimates(res, labels)`: the `estimates[[b]]` list, labelled
  by `.labelEstimates()`, with today's elements in today's order.
- `.setBatchColumns(estRmat, batchCounts, cells)`: today's per-batch
  column insertion (rename, drop, `cbind`, reorder), unchanged. PR 3
  replaces it.
- `.asInputClass(decontXcounts, counts, params)`: today's DelayedMatrix
  and HDF5 conversion, unchanged. PR 4 replaces it.

The batch loop stays a `for` loop in `.decontX()`, so batches run, and
draw random numbers, in today's order. Every function in the core ends
up at 50 lines or fewer.

### `.decontXoneBatch()` split in two

`.decontXoneBatch(counts, z, background, params)` keeps its signature
and its `@importFrom MCMCprecision fit_dirichlet` note, and calls:

- `.initializeBatch(counts, z, params)`: the UMAP log line,
  `.processvarGenes()`, `.processdbscanEps()`, `.decontxInitializeZ()`,
  the UMAP names, and `.processCellLabels()`. Returns `list(z, umap)`.
- `.runEM(counts, z, background, params)`: the `theta` draw,
  `decontXInitialize()`, the background `eta`, and the EM loop with its
  log lines. Returns the log-likelihoods, contamination, `theta`,
  `delta`, `phi`, `eta` and the iteration count.

Initialization runs before the `theta` draw, as today.

### One seed block (A4-06)

`.withSeed(seed, code)` evaluates `code` as is when `seed` is `NULL`,
and inside `withr::with_seed(seed, code)` otherwise, returning its value.
It replaces the three copies: in `.decontX()` (now in `.runBatch()`), in
`.decontxInitializeZLegacy()`, and in `simulateContamination()`.

### Dead code (A4-07)

- `deconMethod`, always `"clustering"`, and its `if`.
- `stopIter` and `numIterWithoutImprovement`, which is never
  incremented, so its loop condition is always true.
- `umap <- NULL`, which is always overwritten.
- The two `decontXEM()` calls, which differ only in `estimate_eta`:
  one call with `estimate_eta = is.null(background)`.
- `totalCells`, the same as `nC`.
- The five `1e-20` literals: one internal constant,
  `.decontXPseudocount <- 1e-20`.
- The checks `.decontXoneBatch()` repeats for each batch:
  `.checkCountsDecon(counts)`, `.checkDelta(delta)` and
  `.checkOneCluster(z)`. `.decontX()` and `.decontXParams()` already run
  them on the whole input.

### Logging (A4-08)

- `.logDecontX(..., params, sep = " ")` calls `.logMessages()` with
  `logfile = params$logfile`, `append = TRUE` and
  `verbose = params$verbose`, and replaces those repeated arguments in
  the core's log calls.
- `.logBanner(..., params)` writes a line of 50 dashes, the text, and
  another line of dashes: the "Starting DecontX" and "Completed DecontX"
  banners, which are six of the core's 14 log calls today.
- `.checkBackground()` and `.matchBackgroundRows()` keep their
  `logfile` and `verbose` arguments, because decontPro calls
  `.matchBackgroundRows()` too. The core passes them from `params`.
- `.logMessages()` writes to the log file with `write()` instead of
  `cat()`, with the same text: `paste(..., "", sep = sep)` and the
  newline `write()` adds give the trailing space before the newline
  that `cat(paste(..., "\n", sep = sep))` gives today. Its logfile check
  and console path don't change. This clears BiocCheck's `cat()` NOTE.
  `.logMessages()` is a copy from celda (`R/celda_functions.R`); its
  output stays the same. At PR time, the maintainer decides whether to
  file a celda issue suggesting the same change.

### `bgBatch` without `background`

`.decontXBatches()` stops with "'bgBatch' was given without
'background'; give both or neither." when `batch` and `bgBatch` are given
and `background` isn't. Today `bgBatch` is silently ignored then. When
`batch` is also missing, the existing "'bgBatch' was given without
'batch'" error still comes first, so no existing error changes. The SCE
and matrix methods both get it from the core.

## Tests

Written first.

- **New, failing on devel:** `bgBatch` with `batch` and without
  `background` stops with the new message, for a matrix and for a
  SingleCellExperiment.
- **New, passing on devel, guarding the refactor:**
  - The log file text, line by line, with dates, the total time and
    the per-iteration numbers masked: one unbatched run, and one batched
    run on a DelayedMatrix with background that has a barcode also in
    the counts and a gene not in them, so every core log line and both
    background log lines are written.
  - The console messages of the unbatched run, masked the same way.
  - `.withSeed()`: `NULL` evaluates the code and draws from the current
    stream; a seed gives the same value as `withr::with_seed()` and
    leaves the caller's stream unchanged.
- **Ported, approved by the maintainer on 2026-11-04:** the two
  `.decontXoneBatch()` calls in "Testing .decontXoneBatch" expect the
  one-cluster and `NA` counts errors that leave `.decontXoneBatch()`.
  They call `decontX()` instead, with the same counts and `z`,
  `verbose = FALSE`, and the same messages, which `.decontX()` raises.
- Unchanged: `test-references.R` and every other test.

## Verification

- `make test-one FILTER=decon` while developing; `make test`,
  `make coverage` (must not drop) and `make lint` before hand-off.
- `Rscript dev/compare-commits.R devel HEAD`: all 13 paths identical.
- `make check-full` and `make bioccheck`. The `cat()` NOTE goes, and
  the count of functions over 50 lines falls from 15.

## Docs and records

- `decontX()` help: `@param bgBatch` says it needs `background`.
  `make docs`.
- NEWS for 1.11.15: the new `bgBatch` error (#76).
- DESCRIPTION: 1.11.15.
- No ADR change: ADR-0004 already records these decisions.

## Related packages

celda forwards to `decontX::decontX()` and reads nothing back.
singleCellTK's `runParams` and results don't change. Code that passes
`bgBatch` with `batch` but without `background` now gets an error
instead of having `bgBatch` ignored.

## Amendments

Made while implementing, on 2026-11-04:

- `.decontXoneBatch()` calls three helpers, not two: `.initializeBatch()`,
  `.initializeEM()` (the `theta` draw, `decontXInitialize()` and the
  background `eta`) and `.runEM()` (the loop). `.runEM()` hands each
  log-likelihood step, with its log line, to `.emLogLik()`. The batch
  loop moved from `.decontX()` into `.runBatches()`, which gathers the
  per-cell results. This keeps every core function at 50 lines or fewer.
- `.logMessages()` builds the text exactly as `cat()` joined it (pieces
  ending in `" \n"`, joined by a space), drops the last newline, and
  hands it to `write()`, which adds it back. `paste(..., "", sep = sep)`
  alone would match `cat()` only for single strings. `write()` treats
  `logfile` as `cat()` does, so `logfile = ""` still logs to standard
  output; a new test pins that.
- The `.withSeed()` test can't pass on devel, where the helper doesn't
  exist; it was written with the helper.
- The snapshot of the "'bgBatch' was given without 'batch'" error (from
  #74) now shows "Error in `.decontXBatches()`" instead of
  "Error in `.decontX()`": R names the function that calls `stop()`. The
  message doesn't change. The other batch errors, now raised in
  `.decontXBatches()`, show that name in their header too.
- The new `bgBatch` error also comes before the two `bgBatch` value
  errors ("values that aren't in 'batch'", "Every batch needs droplets")
  when there is no background; NEWS says `bgBatch` used to be ignored
  unless its values didn't match `batch`.
- A third log snapshot, without `z`, pins the "estimating cell types"
  line. It was recorded on the branch; that line's text is unchanged.
- `.logBanner()` reads the total time after writing the first line of
  dashes, microseconds later than before; only the logged time moves.
- `@importFrom withr with_seed` is removed with its last unqualified use;
  `.withSeed()` calls `withr::with_seed()`.
