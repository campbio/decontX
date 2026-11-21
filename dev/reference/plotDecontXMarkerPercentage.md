# Plots percentage of cells cell types expressing markers

Generates a barplot that shows the percentage of cells within clusters
or cell types that have detectable levels of given marker genes. Can be
used to view the expression of marker genes in different cell types
before and after decontamination with
[`decontX`](https://campbio.github.io/decontX/dev/reference/decontX.md).

## Usage

``` r
plotDecontXMarkerPercentage(
  x,
  markers,
  groupClusters = NULL,
  assayName = c("counts", "decontXcounts"),
  z = NULL,
  threshold = 1,
  exactMatch = TRUE,
  by = "rownames",
  ncol = round(sqrt(length(markers))),
  labelBars = TRUE,
  labelSize = 3
)
```

## Arguments

- x:

  Either a
  [SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  or a matrix-like object of counts.

- markers:

  List. A named list indicating the marker genes for each cell type of
  interest. Multiple markers can be supplied for each cell type. For
  example,
  `list(Tcell_Markers=c("CD3E", "CD3D"), Bcell_Markers=c("CD79A", "CD79B", "MS4A1")`
  would specify markers for human T-cells and B-cells. A cell will be
  considered "positive" for a cell type if it has a count greater than
  `threshold` for at least one of the marker genes in the list.

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
  excluded from the barplot. Default `NULL`.

- assayName:

  Character vector. Name(s) of the assay(s) to plot if `x` is a
  [SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html).
  If more than one assay is listed, then side-by-side barplots will be
  generated. Default `c("counts", "decontXcounts")`.

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

- threshold:

  Numeric. Markers greater than or equal to this value will be
  considered detected in a cell. Default 1.

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
  `"rownames"`.Default `"rownames"`.

- ncol:

  Integer. Number of columns to make in the plot. Default
  `round(sqrt(length(markers))`.

- labelBars:

  Boolean. Whether to display percentages above each bar Default `TRUE`.

- labelSize:

  Numeric. Size of the percentage labels in the barplot. Default 3.

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
#> Sat Nov 21 17:43:14 2026 .. Analyzing all cells
#> Sat Nov 21 17:43:14 2026 .... Converting to sparse matrix
#> Sat Nov 21 17:43:14 2026 .... Generating UMAP and estimating cell types
#> Sat Nov 21 17:43:14 2026 .... Estimating contamination
#> Sat Nov 21 17:43:14 2026 ...... Completed iteration: 9 | converge: 0.0009154
#> Sat Nov 21 17:43:14 2026 .. Calculating final decontaminated matrix
#> --------------------------------------------------
#> Completed DecontX. Total time: 0.229171 secs
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
