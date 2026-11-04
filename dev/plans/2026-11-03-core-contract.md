# Container-agnostic core contract for decontX (#76, PR 1 of 5)

Agreed in chat on 2026-11-03; approved by the maintainer on 2026-11-04.
ADR-0004 (`dev/adr/0004-container-agnostic-core.md`), audit findings
A4-01, A4-02, A4-04 and A4-05 (`dev/audits/2026-10-deep-review.md`), and
#96. Version 1.11.14.

## The five PRs of #76

Each branches from devel, targets devel and says "Part of #76"; the last
one closes it. Every one is held to the references from #66.

1. **This PR:** the core contract, `.decontXParams()`, the formals test,
   the logic moved out of the SCE method, and `runParams` rebuilt
   (closing #96). ADR-0004 is finalized here.
2. Decomposition into named helpers, one seed block, dead code removed,
   `.logMessages` tidied and `cat()` removed, `z` out of `estimates`
   (A4-03, A4-06 to A4-08, A7-22).
3. Batched assembly with one `cbind` and one reorder (A5-05).
4. HDF5 input returns an in-memory DelayedMatrix (A5-06, NEWS).
5. `decontXcounts` and `simulateContamination` move to their own files
   (A4-17).

No `refactor/` base branch: each piece gives the same results and leaves
devel releasable.

## Changes in this PR

### Three layers

1. **Methods read and write the container, nothing else.**
   - The SCE method checks `assayName` and takes the assay. If
     `background` is given, it checks that `background` is a
     SummarizedExperiment, checks `bgAssayName` (default `assayName`)
     and takes that assay. It checks each matrix's class with
     `.checkCountsMatrix()`.
   - It calls `params <- .decontXParams(...)`, then
     `result <- .decontX(counts, z, batch, background, bgBatch, params)`.
   - It writes `colData(x)$decontX_contamination` and
     `$decontX_clusters`, each UMAP from `.decontXUMAPs()` with
     `reducedDim<-`, the assay with `decontXcounts<-`, and
     `metadata(x)$decontX` as `result` without `decontXcounts` and
     `batch`.
   - The `ANY` method does the same without the container steps (it
     keeps today's class check for `x` and, when `background` has
     dimensions, for `background`) and returns `result`.
   - Each method calls `.decontX()` once; `.checkDecontXArgs()` is
     removed.
2. **`.decontXParams(maxIter, delta, estimateDelta, convergence,
   iterLogLik, varGenes, dbscanEps, seed, legacyInit, logfile,
   verbose)`**, with no defaults. It runs today's parameter checks (the
   part of `.checkDecontXArgs()` that isn't about the data, including
   `.checkEMArgs()`, `.checkDelta()`, `.processvarGenes()` and
   `.processdbscanEps()` as validators) and returns a plain named list
   of the values as given, in that order. `varGenes = NULL` and
   `dbscanEps = NULL` stay `NULL`; `.decontXoneBatch()` still resolves
   them, so `runParams` keeps today's values.
3. **`.decontX(counts, z, batch, background, bgBatch, params)`**, with
   no defaults.
   - It checks the data first: `.checkCountsDecon(counts)`, then
     `.checkCellValues()` for `z` and `batch`.
   - It then handles `background`, which moves here from the methods:
     `.checkBackground()` against `colnames(counts)` (shape, `bgBatch`
     length and `NA`, removing barcodes that are also in `counts`),
     `.checkCountsDecon(background, "background")`, and
     `.matchBackgroundRows()`.
   - Batches, labels and the per-batch runs follow as today, reading
     every parameter from `params`.
   - `.decontXoneBatch(counts, z, background, params)` loses its unused
     `batch` argument and its defaults.
   - `.decontxInitializeZ()` loses its defaults (`varGenes = 2000`,
     `dbscanEps = 1`, `estimateCellTypes = TRUE`, `seed = 12345`,
     `legacyInit = FALSE`); its one caller passes all of them.
   - The core's background arguments are named `background` and
     `bgBatch`, like the user arguments and the error messages, not
     `countsBackground` and `batchBackground`.

The loop bodies of `.decontX()` and `.decontXoneBatch()` keep their
logic; splitting them is PR 2.

### Results

- **`.decontX()` returns** `runParams`, `estimates`, `decontXcounts`,
  `contamination`, `z` (today's names, in today's order) and then
  `batch`: the per-cell character vector used for the run,
  `"all_cells"` for every cell when no batch was given. The `ANY`
  method returns this list, so matrix users get the new `batch`
  element.
- **`runParams`** is built once, with today's elements in today's order
  and today's values: `z` as given, `batch` and `batchBackground` as the
  processed per-cell vectors, and the parameters as given.
  `batchBackground` is always present, `NULL` when there's no
  background (#96). Methods no longer read `runParams`; they read
  `result$batch`.
- **UMAPs** stay in `estimates[[b]]$UMAP`, with each batch's cells only.
  `.decontXUMAPs(result, cellNames)` returns the matrices to store, one
  row per cell: `list(decontX_UMAP = <the single batch's UMAP, as it
  is>)`, or one `decontX_<b>_UMAP` per batch, padded with `NA` for the
  other batches' cells, with columns `UMAP_1`, `UMAP_2` and rows named by
  cell. The SCE method uses it now; the Phase II AnnData method can use
  it too. The single-batch columns stay `DecontX_UMAP_1/2`, so the
  stored values don't change; the naming mismatch is left for a
  separate issue.
- **SCE objects** are identical to today's, except that
  `metadata(x)$decontX$runParams` always has `batchBackground`.
- **Related packages:** celda forwards to `decontX::decontX()` and
  reads nothing back. singleCellTK reads `runParams$batch` and
  `$batchBackground` as per-cell vectors, which don't change.

### Error order

The parameter checks now run before the data checks (counts values,
`z` and `batch` lengths), and the background checks run in the core,
after the method has taken the assays. So when an input has more than
one problem a different error may come first (for example, a bad
`bgAssayName` is now reported before the background-barcode checks).
No message changes. If an
existing test depends on the old order, the maintainer decides; tests
are not changed to pass.

## Tests

Written first. Each fails on devel on its assertion, except the method
half of the formals test, which already passes (the two methods agree
today) and guards against future drift.

- **No defaults in the core:** every formal of `.decontX()` (now) and
  `.decontXParams()` (once it exists) is an empty symbol. On devel this
  fails on `.decontX()`'s `maxIter = 200`.
- **Formals agree:** the SCE and `ANY` methods have the same arguments,
  defaults and order, apart from `x`, `assayName` and `bgAssayName`;
  their parameter arguments (all but `x`, `assayName`, `z`, `batch`,
  `background`, `bgAssayName` and `bgBatch`) are exactly
  `names(formals(.decontXParams))`, in order. The method comparison is
  written against devel first; the `.decontXParams()` part is added
  with the constructor.
- **`runParams` always has `batchBackground`**, with and without
  background (fails on devel without background, #96).
- **The `ANY` result has `batch`:** `"all_cells"` per cell without
  `batch`, and `as.character(batch)` with it.
- **`.decontXUMAPs()`:** one batch gives one unpadded `decontX_UMAP`;
  two batches give two padded matrices, `NA` outside their batch, rows
  named by cell, columns `UMAP_1`/`UMAP_2`.
- **SCE metadata** has no `decontXcounts` or `batch` element.
- Unchanged: `test-references.R` and every existing test.

## Verification

- `make test-one FILTER=decon` while developing; `make test`,
  `make coverage` (must not drop) and `make lint` before hand-off.
- `Rscript dev/compare-commits.R devel HEAD`. Expected: every path
  identical except `noZ`, `noZSce` and `legacyInit`, which differ only
  by the new `batch` element (`noZ`, `legacyInit`) and
  `runParams$batchBackground`. A scratch comparison that drops those
  two fields confirms it; the PR reports the result.
- `make check-full` and `make bioccheck`.

## Docs and records

- `decontX()` help: `@return` names `batch` and says `runParams` always
  has `batchBackground`. `make docs`.
- NEWS for 1.11.14: the `batch` element in the matrix result (#76) and
  `runParams$batchBackground` always present (#96).
- DESCRIPTION: 1.11.14.
- ADR-0004: the decisions of 2026-11-03 (below), status Accepted on the
  maintainer's approval.
- `dev/ROADMAP.md`: the five-PR order under I-d.
- `dev/agent-log.md`: an entry for this PR.

## Decisions of 2026-11-03

- PR order and base as above.
- `runParams`: option A, keep its contents and stop methods reading it;
  fix #96 here. Option B (store the raw `batch`) was rejected: it would
  change saved metadata and the per-cell vectors singleCellTK splits by
  sample.
- UMAPs: option A, no separate UMAP element; they stay in `estimates`,
  and `.decontXUMAPs()` pads them for containers that need one row per
  cell. Option B (a padded `UMAP` element in the result) was rejected:
  it duplicates every UMAP in the matrix result.
- New arguments from I-e (A5-01, A5-02) are recorded in #79's spec, not
  an ADR. ADR-0004 states the rule: new decontX arguments go into
  `.decontXParams()` and every method signature, and the formals test
  covers them.
- decontPro: no change in #76. `.decontPro()` is already its core and
  keeps today's list; ADR-0005 decides what containers write and
  ADR-0006 the argument names. A decontPro params constructor waits for
  ADR-0006.

## Amendments

- 2026-11-04, approved by the maintainer: two existing tests call
  internals whose signatures change, so they are ported with the same
  inputs and assertions. "decontX still accepts integer, factor and Inf
  arguments" calls `.decontXParams()` with the same parameter values, and
  `.checkCellValues()` for the same `z` and `batch`. "Testing
  .decontXoneBatch" passes `background = NULL` and
  `params = defaultParams()`, a new helper in `helper-decontx.R` that
  builds `.decontXParams()` from the `ANY` method's defaults.
