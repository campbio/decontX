# Skip the embedding when `z` is given, and `numThreads` (#79)

Agreed in chat on 2026-11-07. Audit findings A5-01, A5-02 and A1-16
(`dev/audits/2026-10-deep-review.md`), decided on 2026-10-19; ADR-0004
("Decided 2026-11-03") puts the new arguments in this spec rather than
an ADR. Phase I-e, version 1.11.19, Fixes #79.

## Today

`.initializeBatch()` always runs the initialization (normalization,
variable genes, PCA, UMAP via scrapper, or scater with
`legacyInit = TRUE`). With a user `z`, its only product is the UMAP
stored in the results, yet it is about 90% of the run time at 50k to
100k cells and about 55% on pbmc4k. Every scrapper call and the legacy
`calculateUMAP()` run with one thread.

## Changes

### New arguments

Both go at the end of the signatures of the `SingleCellExperiment` and
`ANY` methods, after `verbose`, so positional calls don't shift. Both go
through `.decontXParams()` (also at the end), so they are recorded in
`runParams`, and the formals-agreement test covers them.

- `embedding = TRUE`: whether to compute the low-dimensional embedding of
  the cells (today the UMAP). The name is generic on purpose: the method
  may change, the argument doesn't.
- `numThreads = 1`: threads for the initialization.

### `embedding`

- Checked as TRUE or FALSE with the other flags in `.decontXParams()`.
- `.decontX()` stops before the run when `embedding` is FALSE and `z` is
  NULL: "'embedding = FALSE' needs 'z': the cell clusters are estimated
  from the embedding."
- With `embedding = FALSE`, `.initializeBatch()` doesn't call
  `.decontxInitializeZ()` and skips the ".... Generating UMAP" log line.
  It still runs the "All cells must have at least one count" check,
  moved into a small helper that the scrapper initialization calls too,
  so the opt-out doesn't let empty cells reach the EM. It returns
  `umap = NULL`.
- `.checkBatchSizes()` lowers its minimum from 3 cells per batch (15
  with `legacyInit = TRUE`), which only the PCA and UMAP need, to 2. With
  one cell, the default `estimateDelta = TRUE` can't fit the
  contamination prior (MCMCprecision warns "Maximum number of iterations
  reached" each iteration), so the error says the 2 cells are for that.
  The EM's other need (two clusters per batch, or one with background)
  is still checked by `.checkBatchClusters()`. Added after review,
  approved 2026-11-07.
- Each batch's `estimates$UMAP` is `NULL`. `.decontXUMAPs()` returns an
  empty list, so the SCE method writes no `decontX_UMAP` or
  `decontX_<batch>_UMAP`. A decontX UMAP from an earlier run stays in
  `reducedDims` (the 2026-11-01 decision for re-runs) and the docs say
  so.
- The decontamination is bit-identical: the scrapper initialization
  doesn't use R's random numbers, and the legacy one runs inside its own
  `withr::with_seed()`, so the EM sees the same random stream either
  way. Only `estimates$UMAP`, the reducedDims and `runParams` differ.

### `numThreads`

- A whole number of at least 1 (`.isWholeNumber(numThreads, 1)`, finite).
- Passed as `num.threads` to scrapper's `modelGeneVariances()`,
  `runPca()` and `runUmap()`, and as `n_threads` to the legacy
  `scater::calculateUMAP()`. Nothing else in decontX is threaded; it has
  no effect with `embedding = FALSE`.
- The audit measured the scrapper UMAP as bit-identical across thread
  counts and the PCA as differing by about 1e-11. The docs say that
  `numThreads > 1` may change the UMAP coordinates slightly and, rarely,
  the estimated clusters. The legacy path keeps uwot's
  `n_sgd_threads = 0`; whether its UMAP is identical across thread counts
  is checked while implementing and stated in the docs and the PR.

### `plotDecontXContamination()`

- New argument `reducedDimName = NULL`, used only when `x` is a
  SingleCellExperiment (an error otherwise).
- When given: coordinates are the first two columns of
  `reducedDim(x, reducedDimName)`, contamination is
  `colData(x)$decontX_contamination`. With `batch = NULL` every cell is
  plotted; with a batch name, only the cells of that batch
  (`metadata(x)$decontX$runParams$batch`). A name not in
  `reducedDimNames(x)` is an error listing the available names. Axis
  labels come from the reducedDim's column names, falling back to
  `Dim_1`/`Dim_2`.
- When NULL: the decontX UMAP of the batch, as today. If that batch has
  no UMAP (a run with `embedding = FALSE`), the error says so and points
  to `reducedDimName`.
- The `batch` documentation states the two meanings of NULL.

### Docs and NEWS

- Roxygen for `embedding` and `numThreads`. The "If z is provided, it is
  used only for the UMAP stored in the results" sentences in `varGenes`
  and `legacyInit` mention `embedding = FALSE`. `@return` says
  `estimates$UMAP` is NULL without the embedding.
- NEWS for 1.11.19: both arguments, the plot argument, and the
  `numThreads` caveat.
- No vignette change (the rewrite is #83).

## Tests

Written first, in `test-decon.R` and the plot tests:

- `embedding = FALSE` gives identical counts, contamination, z,
  estimates (minus UMAP) to `embedding = TRUE`, single batch and batched,
  matrix and SCE.
- `estimates$UMAP` is NULL, the SCE gets no decontX reducedDim, and an
  earlier `decontX_UMAP` is kept.
- The `z = NULL` error; non-flag `embedding` and invalid `numThreads`
  (0, 1.5, NA, "2", c(1, 2)) errors.
- Zero-count cells still error with `embedding = FALSE`.
- `numThreads = 2` runs on a small fixture (scrapper path; legacy path
  when scater is installed) and gives identical decontamination with
  user `z`.
- `runParams` records both arguments; the formals test passes.
- `plotDecontXContamination()`: `reducedDimName` with all cells, one
  batch, a missing name, a non-SCE `x`; the error when the UMAP is
  missing.
- Log snapshots: the `embedding = FALSE` run's log, if it differs from
  existing snapshots.

`dev/compare-commits.R`: identical on all paths, apart from `runParams`
gaining `embedding` and `numThreads`.

## Benchmark

`dev/bench/bench-decontx.R` gains `DECONTX_BENCH_EMBEDDING` and
`DECONTX_BENCH_THREADS` (or equivalent options) so the PR reports, in
sim mode (10k, 50k, 100k) and real mode (pbmc4k): `z` given with and
without the embedding, and estimated `z` with `numThreads` 1 and 4.

## Checks

`make test`, `make check-full`, `make bioccheck`; no new lints.
