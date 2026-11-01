# decontX 1.11.7

* `decontX()` allows a `z` with one cluster when `background` is
supplied, since the contamination profile then comes from the empty
droplets. This fits only when the cells are a single cell type. Without
`background`, one cluster still stops, and the error now names the batch
(#69).
* When dbscan finds only one cluster, the k-means fallback now runs on the
PCA components instead of a dense matrix of all genes, which could need
gigabytes of memory (#43). This changes the clusters and contamination
estimates on that path only; `legacyInit = TRUE` still clusters all
genes.
* With `legacyInit = TRUE`, the k-means fallback no longer stops with
"argument is not a matrix". It had failed on this path since
`legacyInit` was added, so no earlier results change.

# decontX 1.11.6

* `decontX()` now matches the rows of `background` to the genes of `x` by
name, and `decontPro()` matches the rows of `ambient_counts` to the ADTs
of `filtered_counts`. Before, rows were matched by position: rows in a
different order ran without an error and gave wrong estimates, extra or
missing genes gave a cryptic error (#43), and names that didn't overlap
ran. Now extra rows are dropped with a log message, missing genes stop
with an error that names them, and rows are matched by position only
when either matrix has no row names. Results don't change when the rows
already match (#68).
* `bgBatch` must have one value per background droplet, with no `NA`,
and its values must be the same as those of `batch`. Before, a short
`bgBatch` was recycled, and a batch without matching background droplets
silently ran without background. `decontX()` also stops when no
background droplets are left after removing barcodes that are in `x`,
and a batch with a single background droplet no longer fails.

# decontX 1.11.5

* `decontX()` now checks its inputs before running and stops with an error
that names the argument. Before, a `z` or `batch` shorter than the number
of cells was padded with `NA` and gave wrong results; `z` or `batch` with
`NA`, invalid `maxIter`, `iterLogLik`, `convergence`, `seed`, `varGenes`,
`dbscanEps` or `logfile` values, non-logical `estimateDelta`,
`legacyInit` or `verbose`, a `delta` that isn't two finite positive
numbers, and negative counts in `x` or `background` ran silently or
failed with cryptic errors. `maxIter`, `iterLogLik` and `varGenes` must
be whole numbers (or `Inf`, which keeps working), and `maxIter = Inf`
needs `convergence` above 0. An `assayName` or `bgAssayName` that isn't
an assay of the input now stops before decontX runs; an assay index
still works. Results don't change for valid input (#67).
* The `varGenes` and `legacyInit` docs now say they still set the UMAP
stored in the results when `z` is given, and the error for cells with no
counts no longer suggests supplying `z`.

# decontX 1.11.3

* The "Working with Seurat" section of the decontX vignette now runs as
pasted with SeuratObject 5: it reads counts with `LayerData()` instead of
`GetAssayData(slot = )` (defunct since SeuratObject 5.3), attaches
SingleCellExperiment, and uses one name for the Seurat object (#65).

# decontX 1.11.2

* The `decontX-package` help page now has a Value section, which fixes a
BiocCheck warning.
* Documentation is regenerated with roxygen2 8.1.0, and the Stan model
exports with rstantools 2.7.1. Results don't change.

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
work when `x` is a dense or sparse matrix (they previously failed
with a ggplot2 aesthetics error when the plot was drawn). The facets of
`plotDecontXMarkerExpression()` now follow the order of `groupClusters`,
matching the x-axis.
* `decontPro()` no longer draws 1000 unused samples from the variational
posterior. Only the mean of the approximation is used, but the draws
scaled with the number of ADTs x droplets and could take hundreds of GB
of disk and memory on large datasets (#46). On a 50 x 400 test dataset
the temporary file shrank from 444 MB to 2 MB and peak memory from 3.0 GB
to 1.2 GB. Results are unchanged. rstan's "Pareto k diagnostic" warning
is no longer shown: it is estimated from the draws and would always
report `Inf` with the minimal number of draws now requested.

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
