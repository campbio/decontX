# Plots expression of marker genes before and after decontamination

Generates a violin plot that shows the counts of marker genes in cells
across specific clusters or cell types. Can be used to view the
expression of marker genes in different cell types before and after
decontamination with
[`decontX`](https://campbio.github.io/decontX/dev/reference/decontX.md).

## Usage

``` r
plotDecontXMarkerExpression(
  x,
  markers,
  groupClusters = NULL,
  assayName = c("counts", "decontXcounts"),
  z = NULL,
  exactMatch = TRUE,
  by = "rownames",
  log1p = FALSE,
  ncol = NULL,
  plotDots = FALSE,
  dotSize = 0.1
)
```

## Arguments

- x:

  Either a
  [SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  or a matrix-like object of counts.

- markers:

  Character Vector or List. A character vector or list of character
  vectors with the names of the marker genes of interest.

- groupClusters:

  List. A named list that allows cell clusters labels coded in `z` to be
  regrouped and renamed on the fly. For example,
  `list(Tcells=c(1, 2), Bcells=7)` would recode clusters 1 and 2 to
  "Tcells" and cluster 7 to "Bcells". The values are the labels in `z`,
  such as the labels of a `z` supplied to
  [`decontX`](https://campbio.github.io/decontX/dev/reference/decontX.md).
  After a batched
  [`decontX`](https://campbio.github.io/decontX/dev/reference/decontX.md)
  run, the labels can have the batch as a prefix (`"<batch>-<label>"`,
  such as `"b1-1"`); use the labels as they appear in `z`. Note that if
  this is used, clusters in `z` not found in `groupClusters` will be
  excluded. Default `NULL`.

- assayName:

  Character vector. Name(s) of the assay(s) to plot if `x` is a
  [SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html).
  If more than one assay is listed, then side-by-side violin plots will
  be generated. Default `c("counts", "decontXcounts")`.

- z:

  Character, Integer, or Vector. Indicates the cluster labels for each
  cell. If `x` is a
  [SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  and `z = NULL`, then the cluster labels from
  [`decontX`](https://campbio.github.io/decontX/dev/reference/decontX.md)
  will be retrieved from the `colData` of `x` (i.e.
  `colData(x)$decontX_clusters`). If `z` is a single character or
  integer, then that column will be retrieved from `colData` of `x`.
  (i.e. `colData(x)[,z]`). If `x` is a counts matrix, then `z` will need
  to be a vector the same length as the number of columns in `x` that
  indicate the cluster to which each cell belongs. Default `NULL`.

- exactMatch:

  Boolean. Whether to only identify exact matches for the markers or to
  identify partial matches using
  [`grep`](https://rdrr.io/r/base/grep.html). See
  [`retrieveFeatureIndex`](https://campbio.github.io/decontX/dev/reference/retrieveFeatureIndex.md)
  for more details. Default `TRUE`.

- by:

  Character. Where to search for the markers if `x` is a
  [SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html).
  See
  [`retrieveFeatureIndex`](https://campbio.github.io/decontX/dev/reference/retrieveFeatureIndex.md)
  for more details. If `x` is a matrix, then this must be set to
  `"rownames"`. Default `"rownames"`.

- log1p:

  Boolean. Whether to apply the function `log1p` to the data before
  plotting. This function will add a pseudocount of 1 and then log
  transform the expression values. Default `FALSE`.

- ncol:

  Integer. Number of columns to make in the plot. Default `NULL`.

- plotDots:

  Boolean. If `TRUE`, the expression of features will be plotted as
  points in addition to the violin curve. Default `FALSE`.

- dotSize:

  Numeric. Size of points if `plotDots = TRUE`. Default `0.1`.

## Value

Returns a `ggplot` object.

## See also

See
[`decontX`](https://campbio.github.io/decontX/dev/reference/decontX.md)
for a full example of how to estimate and plot contamination.

## Author

Shiyi Yang, Joshua Campbell

## Examples

``` r
# Generate matrix with contamination
s <- simulateContamination(seed = 12345)

library(SingleCellExperiment)
library(scater)
sce <- SingleCellExperiment(list(counts = s$observedCounts))
sce <- decontX(sce)
#> --------------------------------------------------
#> Starting DecontX
#> --------------------------------------------------
#> Mon Nov 23 16:40:21 2026 .. Analyzing all cells
#> Mon Nov 23 16:40:21 2026 .... Converting to sparse matrix
#> Mon Nov 23 16:40:21 2026 .... Generating UMAP and estimating cell types
#> Mon Nov 23 16:40:22 2026 .... Estimating contamination
#> Mon Nov 23 16:40:22 2026 ...... Completed iteration: 9 | converge: 0.0009154
#> Mon Nov 23 16:40:22 2026 .. Calculating final decontaminated matrix
#> --------------------------------------------------
#> Completed DecontX. Total time: 0.2385576 secs
#> --------------------------------------------------

# Plot contamination on UMAP
plotDecontXContamination(sce)


# Plot decontX cluster labels on the UMAP
plotReducedDim(sce, dimred = "decontX_UMAP",
               colour_by = "decontX_clusters",
               text_by = "decontX_clusters")


# Plot percentage of marker genes detected
# in each cell cluster before decontamination
s$markers
#> $CellType_1_Markers
#> [1] "Gene_47" "Gene_32" "Gene_86"
#> 
#> $CellType_2_Markers
#> [1] "Gene_70" "Gene_33" "Gene_48"
#> 
#> $CellType_3_Markers
#> [1] "Gene_74" "Gene_26" "Gene_20"
#> 
plotDecontXMarkerPercentage(sce, markers = s$markers, assayName = "counts")


# Plot percentage of marker genes detected
# in each cell cluster after contamination
plotDecontXMarkerPercentage(sce, markers = s$markers,
                            assayName = "decontXcounts")


# Plot percentage of marker genes detected in each cell
# comparing original and decontaminated counts side-by-side
plotDecontXMarkerPercentage(sce, markers = s$markers,
                            assayName = c("counts", "decontXcounts"))


# Plot raw counts of indiviual markers genes before
# and after decontamination
plotDecontXMarkerExpression(sce, unlist(s$markers))
```
