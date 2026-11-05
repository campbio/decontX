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
