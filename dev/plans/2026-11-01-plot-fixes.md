# Plot fixes: shared bandwidth, `file`, `groupClusters` labels, unnamed UMAPs (#73)

Status: approved by the maintainer, 2026-11-01 (pooled bandwidth, after
comparing it with per-curve and smaller-of-two bandwidths on the
decontPro vignette data). Roadmap Phase I-b; audit findings
A1-24, A1-25, A7-17 and A6-04.

## Why

Reproduced on devel at d15a881 (1.11.10).

- **Bandwidth (A1-24).** `plotDensity()` (`plot_decontPro.R:65`) lets
  `geom_density()` choose a bandwidth per curve, so the original and
  decontaminated curves of one ADT are smoothed differently; the audit
  measured 0.204 vs 0.390 and peaks inflated 1.18-1.59x.
- **`file` (A7-17).** `plotDensity()` and `plotBoxByCluster()` always
  append ".pdf", so `file = "x.pdf"` writes `x.pdf.pdf`. The docs say
  they return "a pdf file"; they return the path, invisibly (from
  `ggsave()`).
- **`groupClusters` (A1-25).** After a batched run the labels are
  prefixed (`"b1-1"`), so `groupClusters = list(T = c(1, 2))` stops with
  `'groupClusters' not found in 'z': 1,2`, which doesn't show the labels
  to use. The docs don't mention the prefix.
- **Stale UMAPs (A1-25).** A re-run with a different batch layout leaves
  the earlier run's `decontX_UMAP` / `decontX_<b>_UMAP` in
  `reducedDims`. Maintainer decision (2026-11-01): keep them, since
  someone may run decontX several ways on purpose, and document it.
- **Unnamed UMAP (A6-04).** `plotDecontXContamination()` reads the UMAP
  columns by name (`plot_decontx.R:67-68`) and fails with
  `.data[[NULL]]` when they have none. decontX always names them
  (`decon.R:275, 858`); names are lost outside decontX, as in an h5ad
  round trip through `uns`. Maintainer decision: keep a small positional
  fix now; the Phase II AnnData getter restores the names.

## What changes

1. **`plotDensity()`.** For each feature, one bandwidth from both curves
   pooled: `stats::bw.nrd0()` (geom_density's default rule) on the
   pseudo-log scale the density is computed on, `asinh(x / 2)` (what
   `transform = "pseudo_log"` applies with its defaults), passed as
   `geom_density(bw = )`. No new import.
2. **`file`.** A shared internal helper for both functions: append
   ".pdf" only when `file` doesn't already end in ".pdf" (any case), save
   with `ggsave()` as now (same size), return the path invisibly.
   `@param file` and `@return` say so.
3. **`groupClusters`.** The not-found error in
   `.processPlotDecontXMarkerInupt()` also lists the labels in `z` (the
   first 10, then "..."). `@param groupClusters` of
   `plotDecontXMarkerPercentage()` and `plotDecontXMarkerExpression()`
   says batched decontX runs can prefix labels with the batch
   (`"<batch>-<label>"`) and to use the labels as they appear in `z`.
4. **Stale UMAPs.** No code change. The `decontX()` `@return` says a
   re-run replaces the assay, the colData columns and
   `metadata(x)$decontX`, while a UMAP from an earlier run under a
   different batch name stays in `reducedDims` and doesn't match the
   current results.
5. **`plotDecontXContamination()`.** Take UMAP columns 1 and 2 by
   position. Axis titles are the column names when present (so named
   UMAPs plot as before), else `UMAP_1` / `UMAP_2`.
6. **Phase II note.** ADR-0005 (AnnData) and the ROADMAP Phase II
   AnnData bullet record that the AnnData getter must restore the names
   `uns` drops when it rebuilds the estimates (A6-04); the #73 fix is a
   stopgap. The ROADMAP #73 line drops "clear stale per-batch UMAPs".

7. **Added at review.** `plotDensity()` leaves NA and Inf values out of
   the bandwidth, as `geom_density()` drops them, and keeps the default
   bandwidth with fewer than 2 values left (the pooled bandwidth stopped
   with "need at least 2 data points"). The transform comes from the
   plot's own x scale, not a separate `asinh(x / 2)`. The UMAP axis
   titles also fall back for empty or NA names. The `groupClusters`
   message uses the existing `.firstFew()`. `@return` names the
   patchwork object. Not taken: removing old UMAPs (maintainer decision
   above), a shape check on the UMAP (decontX always writes two columns;
   Phase II), and asserting on built density grids instead of the
   layer's `bw` (more fragile).

8. **Added at hand-off** (maintainer request, 2026-11-02). The decontPro
   vignette passes `Idents(adt_seurat)` to `decontPro()` instead of
   `as.integer(Idents(adt_seurat))`, which turned Seurat's 0-based
   clusters into codes 1..K, so the `plotBoxByCluster()` numbers were one
   higher than the Seurat plots. decontPro has accepted any labels since
   #72; on the vignette data the decontaminated counts are identical, and
   the `r_est` rows are named "0" to "8".

## Tests

Written first; each fails on its assertion on devel.

- `plotDensity()`: the density layer's `bw` equals `bw.nrd0()` of the
  pooled pseudo-log values, for each panel.
- `plotDensity()` and `plotBoxByCluster()` with `file = "<tmp>/a.pdf"`
  write `a.pdf` (no `a.pdf.pdf`); with `"<tmp>/b"` write `b.pdf`; the
  return is invisible and equals the path.
- `plotDecontXMarkerPercentage()` with prefixed labels in `z` and
  unprefixed `groupClusters` stops with a message listing `"b1-1"`.
- `plotDecontXContamination()` on a result whose UMAP has no column
  names builds without error, with axis titles `UMAP_1` / `UMAP_2`;
  named UMAPs keep their names as axis titles.

## NEWS, version, release

- Version 1.11.11. NEWS entries: shared bandwidth, the `file` argument,
  the `groupClusters` message and docs, unnamed UMAP columns.
- Devel only: plotting fixes, no change to results.

## Effect on results and related packages

No change to decontX or decontPro results, so `dev/compare-commits.R`
isn't needed. `plotDensity()` curves look different (both now smoothed
alike). celda and singleCellTK call `decontX()`, not these plots.
