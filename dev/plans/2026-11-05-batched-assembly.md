# Assemble batched results with one cbind (#76, PR 3 of 5)

Agreed in chat on 2026-11-05. ADR-0004
(`dev/adr/0004-container-agnostic-core.md`, "Batched assembly"), audit
finding A5-05 (D1-27, `dev/audits/2026-10-deep-review.md`). Version
1.11.16. Part of #76.

## Today

`.runBatches()` starts from an all-zero sparse matrix the size of the
counts. For each batch, `.setBatchColumns()` renames every column, drops
the batch's columns, `cbind`s the batch's result, reorders every column
and restores the names: about four full copies per batch, one batch
included.

## Changes

Results stay `identical()`.

- `.runBatches()` keeps each batch's decontaminated counts in a list,
  named by batch, and no longer builds the zero matrix.
- `.bindBatches(batchCounts, batch, dimNames)`, new, with no defaults,
  replaces `.setBatchColumns()`:
  - one batch: that batch's matrix as it is;
  - more: one `cbind` in batch order (`unique(batch)`), then one
    column reorder back to cell order, skipped when the cells are
    already in that order;
  - then the counts' dimnames. Subsetting uses `drop = FALSE`.
- The per-cell contamination and `z`, `.asInputClass()` (PR 4) and the
  rest of the core don't change.

Prototype at 2,000 genes x 40,000 cells (about 92 MB per copy): one
batch 0.030 s and 200 MB extra peak before, 0.002 s and 5 MB after; four
batches 0.132 s and 268 MB before, 0.024 s and 239 MB after. The audit's
"7.5 copies to 2" counts allocations; the peak falls less.

## Tests

- **New, passing on devel, guarding the refactor:** a batched run with
  interleaved batches (`rep(c("A", "B"), 50)`; the references use
  contiguous blocks) and user `z`. Each batch's columns of
  `decontXcounts`, and its contamination, are `identical()` to a run on
  that batch alone, and the dimnames are the counts'.
- **New, for the helper:** `.bindBatches()` with one batch, interleaved
  batches, and no column names. It can't run on devel, where the helper
  doesn't exist.
- Unchanged: `test-references.R` and every other test.

## Verification

- `make test`, `make coverage` (must not drop), `make lint`.
- `Rscript dev/compare-commits.R devel HEAD`: all 13 paths identical.
- `make check-full` and `make bioccheck`.
- `dev/bench/bench-decontx.R` before and after, `z_supplied` and
  `batched4` at 50,000 cells, numbers in the PR.

## Docs and records

- NEWS for 1.11.16: batched runs put the decontaminated counts together
  faster and with less memory; results don't change (#76).
- DESCRIPTION: 1.11.16.
- No ADR change.

## Related packages

celda and singleCellTK get identical results.

## Amendments

Made while implementing and reviewing, on 2026-11-05:

- The spec review found that `dimnames(counts)` keeps the names of the
  dimnames list (for example from `table()`), which the old zero matrix
  dropped. `.runBatches()` passes `list(rownames(counts),
  colnames(counts))`, and a test pins the unnamed dimnames, with and
  without batches.
- `/code-review` found that looking batches up by name loses a batch
  labelled `""`: batched runs with that label stopped with "subscript
  out of bounds", where devel assembled the counts correctly. The
  per-batch counts are a positional list in `unique(batch)` order, and
  a test checks that batch's counts and contamination against a run on
  it alone. On devel, the cells of a `""` batch already get the cluster
  label `"-"` and its estimates can't be read by name; that older bug is
  unchanged here and is filed as #102.
- The helper is split in two, so `.runBatches()` can drop the per-batch
  list before the reorder: `.bindBatches(batchCounts)` binds (one batch
  as it is), and `.inCellOrder(bound, batch, dimNames)` reorders, only
  when batches interleave, and sets the dimnames. The cell order comes
  from one `order(match(batch, unique(batch)))`.
- Measured at 2,000 genes x 40,000 cells (about 94 MB per copy), each
  batch's result freshly allocated as `calculateNativeMatrix()` does:
  one batch 0.10 s and 289 MB extra peak before, 0.02 s and 94 MB (the
  result itself) after; four interleaved batches 0.14 s and 448 MB
  before, 0.05 s and 394 MB after. Matrix's `cbind` alone peaks at about
  1.75 copies, so the batched peak stays near 2.75 copies, not the
  ADR's "about 2". Going lower would mean building the sparse slots by
  hand, which this PR doesn't do.
