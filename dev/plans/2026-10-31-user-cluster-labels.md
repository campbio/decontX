# Return user cluster labels, and fix the marker plots (#70)

Status: approved by the maintainer, 2026-10-31. Roadmap Phase I-b; audit
findings A1-09, A1-22, D2-08, A1-07, A1-20.

## Why

Reproduced on devel at 61a109c.

- **User labels come back as codes (A1-09).** `.decontX()` fills
  `returnZ <- rep(NA, nC)` (logical), so assigning a batch's factor
  keeps its integer codes; batched runs use `paste0(bat, "-", res$z)`,
  which gives labels. Character `z` `T/B/Mono` comes back as `1 2 1 3`
  unbatched and `"A-1" "A-2"` batched; a factor comes back as codes
  unbatched and `"A-T"` batched. `.processCellLabels()` maps labels to
  first-appearance numbers and then `as.factor()`, which sorts character
  numbers as text: with 12 clusters, `"c2"` gets code 5. `phi` and `eta`
  have no column names. (NA in `z` is already rejected since #67.)
- **Unused factor levels give NaN columns (A1-22, D2-08).** A level with
  no cells in a batch, or an unused level before a used one, gives an
  all-NaN `phi` column: `decontXEM` divides 0/0 for an empty cluster.
- **Marker percentages are wrong for labels that aren't 1..K (A1-07).**
  `.calculateDecontXBarplotPercent()` groups by label value but sizes
  clusters with `tabulate(z)` (by position), and the plot indexes the
  label names by value. 0-based labels error ("differing number of rows:
  9, 6"), 10/11/12 give NA labels, 2/5/9 give two NA labels, and a factor
  with an unused level gives percentages of `Inf` and 168.
- **Marker plot crashes (A1-20).** One marker gene drops the matrix to a
  vector (`x[geneMarkerIndex, ]`), and `DelayedArray()` fails.
  `plotDecontXMarkerExpression()` looks clusters up with `z[df$Cell]`,
  but `reshape2::melt` turns numeric-like cell names into integers, so
  it indexes by position: with permuted names 576 of 900 rows got the
  wrong cluster, and names like `"1001"` gave NA for every row.

## What changes

1. **`result$z` is always a factor.**
   - User `z`: the user's labels, unprefixed even when batched. A factor
     keeps its level order minus unused levels; any other `z` gets
     sorted levels, as `factor()` gives.
   - Estimated `z`: the same values as today, now as a factor. Unbatched
     levels are `1..K`, so `as.integer()` gives today's codes; batched
     labels stay `"<batch>-<k>"`, ordered by batch, then numerically.
   - `colData(x)$decontX_clusters` gets the same factor.
2. **Labelled estimates, same fit.** The internal codes and the EM are
   unchanged. On output, each batch's `estimates[[batch]]$z` is a factor
   of the labels in that batch (estimated: `1..K`), and the `phi` and
   `eta` columns are named by label and ordered as its levels, so
   `phi[, as.integer(z)]` picks each cell's column. Values stay
   bit-identical.
3. **Unused levels are dropped per batch** before the fit, so
   `decontXEM` never gets an empty cluster. No C++ guard: `decontXEM`
   isn't exported, and nothing can reach it with an empty cluster after
   this change.
4. **Marker plots.**
   - The shared input helper drops unused levels, and the percentage
     calculation numbers the levels 1..K, with the level names as the
     cluster names, so the groups, sizes and labels line up.
   - `drop = FALSE` when subsetting to the marker genes.
   - `plotDecontXMarkerExpression()` attaches each cell's cluster by
     position, not by cell name.
5. **Docs:** `z` and the returned `z` in `?decontX`; `groupClusters` in
   both marker plot pages says it takes the labels in `z`.

## Results

- A factor `z` with unused levels loses its NaN `phi` column;
  contamination moves by about 1e-16. Nothing else changes numerically.
- `z` changes type and values (labels instead of codes); `phi`/`eta` gain
  column names, and their columns move when the first-appearance order
  differs from the sorted order.
- `groupClusters` with a user `z` now takes the labels, not codes.
- Marker plot labels change only where they are wrong today.
- `decontx-refs.rds` is regenerated (`z` type, `phi`/`eta` names). Before
  that, a scratch comparison confirms every number matches the stored
  references with names and the type of `z` set aside.

## Tests

Written first; each fails on its assertion on devel.

- Character and factor `z`, unbatched and batched: `result$z` is a factor
  with the user's labels and levels; the SCE `decontX_clusters` matches.
- 12 character clusters: each cell gets its own label back.
- `phi`/`eta` column names equal `levels(estimates$z)`; numeric `z`
  given out of sorted order gives the same values as sorted `z`, with
  columns reordered.
- Estimated `z`: unbatched factor with `as.integer()` equal to today's
  codes; batched `"<batch>-<k>"`.
- An unused factor level, and a level missing from one batch: no NaN in
  `phi`, and `estimates$z` has no unused levels.
- Marker percentages for 0-based, 10/11/12 and 2/5/9 labels match the
  1..K percentages with the right labels; a factor with an unused level
  gives no `Inf`; one marker gene works; the expression plot with
  permuted numeric-like or `"1001"`-style cell names gives each cell its
  own cluster.

## Not changing

- The release branch: the result type changes, so this is devel-only.
- `.processCellLabels()`'s internal code order (keeps the fit
  bit-identical).
- reshape2 (A4-10 replaces it separately).

## Checklist

Version 1.11.8; NEWS entries (`as.integer(result$z)` may change;
`groupClusters` takes labels for a user `z`); `make docs`; `make test`,
`make coverage`; `make check-full`, `make bioccheck` before the PR;
`dev/compare-commits.R` differs only in `z` type and `phi`/`eta` names.
