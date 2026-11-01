# Allow one cluster with background, and run the k-means fallback on PCA (#69)

Status: approved by the maintainer, 2026-10-31. Roadmap Phase I-b; audit
findings A8-01 and A8-02. Refs #43.

## Why

Reproduced on devel at 0040d2c.

- **One cluster is refused even with background (A8-01).**
  `.processCellLabels()` stops with "No need to decontaminate when only
  one cluster is in the dataset" for any one-cluster `z`, per batch. The
  model is valid when the cells really are one cell type and the
  contamination profile comes from `background`. On 323 simulated cells
  of one type, with the check patched out:
  - one cluster + background: correlation 0.9978, MAE 0.0042
  - `z = rep(1:2, ...)` + background (the #43 workaround): 0.9978, 0.0042
  - estimated z + background (4 spurious clusters): 0.983, 0.015
  - estimated z, no background: 0.449, 0.22

  Without background, one cluster crashes inside the EM ("missing value
  where TRUE/FALSE needed"), so the error stays there. The error doesn't
  name the batch. Forcing three true cell types into one cluster gives a
  correlation of 0.10, so the docs must say one cluster fits only cells
  of one type.
- **The k-means fallback runs on a dense all-genes matrix (A8-02).** When
  dbscan finds one cluster, `stats::kmeans(as.matrix(t(normed)), 2)`
  builds a dense cells x genes matrix: 1.36 GiB at 5,000 cells x 36,601
  genes, which matches #43's 1.4 GiB warning. On simulated data the same
  k-means on the PCA components is as accurate (correlation 0.80 vs 0.78
  over three seeds), faster, and uses about a third less peak memory.

## What changes

1. **One cluster with background.** `.decontXoneBatch()` allows a `z`
   with one cluster when the batch has background (since #68 every batch
   has background droplets when `background` is given). Without
   background it stops; the message names the batch when `batch` is
   given and says that `background` allows one cluster.
2. **k-means fallback on PCA.** `.decontxInitializeZScrapper()` also
   returns its PCA components, and the fallback runs k-means on them
   (cells x at most 50 PCs). `legacyInit = TRUE` keeps the dense
   all-genes matrix, transposed with `Matrix::t()`: base `t()` failed
   on the sparse log counts, so this path had always stopped with
   "argument is not a matrix" (added after review, with the
   maintainer's approval).
3. **Docs:** `z` and `background` say that one cluster is allowed with
   `background`, and fits only when the cells are one cell type.

## Tests

Written first; each fails on its assertion on devel.

- One-cluster `z` with background runs, unbatched and with one batch down
  to one cluster; one cluster comes back.
- One cluster without background errors; the batched message names the
  batch.
- With `stats::kmeans` mocked and `dbscanEps = 0` (which forces the
  fallback), k-means gets one row per cell and at most 50 columns.
  `legacyInit = TRUE` runs, gets a dense cells x all genes matrix, and
  returns two clusters (skipped without scater; on devel it fails on the
  `t()` error).
- No change elsewhere: the stored references still match, and
  `dev/compare-commits.R devel HEAD` reports every path identical, unless
  a path hits the fallback, which is then named in the PR.

## NEWS, version, release

- NEWS: one cluster allowed with background; the fallback's results
  change (default initialization only, when dbscan finds one cluster);
  the legacy fallback no longer fails.
- Version 1.11.7.
- Devel only, with the other Phase I-b fixes.

## Effect on related packages

celda forwards to `decontX()`. A one-cluster `z` with background now runs
instead of stopping; runs that hit the k-means fallback get different
clusters and contamination estimates.

## Out of scope

- using one cluster instead of k-means when z is estimated, background is
  given and dbscan finds one cluster
- a batch with exactly one cell (A1-08)
- the class of the returned `z` (#70)
