# decontX 1.11.1

* `decontX` no longer depends on the `celda` package (ADR-0003 Stage 1).
The one internal use of `celda::normalizeCounts()` (proportion
normalization inside `simulateContamination()`) is replaced by the
package's own `fastNormProp` routine, which produces identical results,
and the example and vignette now plot cluster labels and marker-gene
expression on the UMAP with `scater::plotReducedDim()` instead of celda's
plotting functions. This drops `celda` (and its transitive dependencies)
entirely and removes the `R CMD check` deprecation warnings that came from
celda's masked, stale copy of `decontX` being invoked after
`library(celda)` in the examples.
* The cell cluster initialization used by `decontX()` when `z` is not
supplied (log-normalization, variable gene selection, PCA, and UMAP) has
been migrated from scater/scuttle to the scrapper package, because
scuttle deprecated `logNormCounts()`/`normalizeCounts()` in Bioconductor
3.24. This removes the deprecation warnings from `R CMD check`.
* Results of `decontX()` without `z` may differ slightly from previous
versions (same workflow and parameters; different variable-gene
criterion and UMAP implementation). A new `legacyInit = TRUE` argument
reproduces the previous behavior; it requires the scater package (now in
Suggests) and is kept for backwards compatibility for as long as the
upstream functions remain available. Users supplying their own `z` are
unaffected.
* With the default initialization, the UMAP step itself is deterministic
even when `seed = NULL` (scrapper seeds it explicitly), but overall
`decontX()` results with `seed = NULL` remain non-reproducible because
the EM initialization is unseeded — set a seed for full reproducibility.
* The plotting functions (`plotDecontXContamination()`,
`plotDecontXMarkerPercentage()`, `plotDecontXMarkerExpression()`,
`plotDensity()`, `plotBoxByCluster()`) no longer use the deprecated
`ggplot2::aes_string()` or `trans =` scale argument, which removes the
ggplot2 deprecation warnings they emitted. decontX now requires
ggplot2 (>= 3.5.0).
* `plotDecontXMarkerPercentage()` and `plotDecontXMarkerExpression()` now
work when `x` is a matrix (they previously failed with a ggplot2
aesthetics error when the plot was drawn). The facets of
`plotDecontXMarkerExpression()` now follow the order of `groupClusters`,
matching the x-axis.

# decontX 1.10.0

* No significant changes (Bioconductor 3.23 release).

# decontX 1.8.0

* No significant changes (Bioconductor 3.22 release).

# decontX 1.6.0

* No significant changes (Bioconductor 3.21 release). Includes the
DecontPro changes released in 1.4.1 on the Bioconductor 3.20 branch.

# decontX 1.4.1 (2025-02-16)

* `decontPro()` gained an `ambient_counts` argument, so an empty-droplet
count matrix can be supplied to estimate the ambient contamination
profile empirically instead of inferring it from the filtered droplets.
* **Breaking:** the first argument of `decontPro()` was renamed from
`object` to `filtered_counts` for consistency with the new
`ambient_counts` argument. Code that passed the count matrix by name
(`decontPro(object = x, ...)`) must be updated; positional calls are
unaffected.
* The estimated ambient contamination profile is now returned in the
`parameters` element of the `decontPro()` output (`p_est`).
* Updated deprecated array syntax in the Stan model, fixed a C compiler
warning in `matrixSums`, and corrected a `deconPro` typo in the
DecontPro vignette.

# decontX 1.4.0

* No significant changes (Bioconductor 3.20 release).

# decontX 1.2.0

* No significant changes (Bioconductor 3.19 release).

# decontX 1.0.0

* First Bioconductor release (Bioconductor 3.18), providing `decontX()`
for ambient RNA decontamination of single-cell RNA-seq counts and
`decontPro()` for single-cell protein expression (CITE-seq/TotalSeq)
data.

# decontX 0.99.5 (2023-10-19)

* First submission to Bioconductor after review. 
