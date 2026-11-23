# Get or set decontaminated counts matrix

Gets or sets the decontaminated counts matrix from a a
[SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
object.

## Usage

``` r
decontXcounts(object, ...)

decontXcounts(object, ...) <- value

# S4 method for class 'SingleCellExperiment'
decontXcounts(object, ...)

# S4 method for class 'SingleCellExperiment'
decontXcounts(object, ...) <- value
```

## Arguments

- object:

  A
  [SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  object.

- ...:

  For the generic, further arguments to pass to each method.

- value:

  A matrix to save as an assay called `decontXcounts`

## Value

If getting, the assay from `object` with the name `decontXcounts` will
be returned. If setting, a
[SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
object will be returned with `decontXcounts` listed in the `assay` slot.

## See also

[`assay`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
and `assay<-`

## Examples

``` r
# Simulate contaminated counts and store them in a SingleCellExperiment
library(SingleCellExperiment)
sim <- simulateContamination(seed = 12345)
sce <- SingleCellExperiment(list(counts = sim$observedCounts))

# Run decontX with the known cell populations, skipping the UMAP
sce <- decontX(sce, z = sim$z, embedding = FALSE)
#> --------------------------------------------------
#> Starting DecontX
#> --------------------------------------------------
#> Mon Nov 23 16:40:15 2026 .. Analyzing all cells
#> Mon Nov 23 16:40:15 2026 .... Converting to sparse matrix
#> Mon Nov 23 16:40:15 2026 .... Estimating contamination
#> Mon Nov 23 16:40:15 2026 ...... Completed iteration: 9 | converge: 0.0009154
#> Mon Nov 23 16:40:15 2026 .. Calculating final decontaminated matrix
#> --------------------------------------------------
#> Completed DecontX. Total time: 0.004794598 secs
#> --------------------------------------------------

# Get the decontaminated counts
decontaminated <- decontXcounts(sce)
decontaminated[1:5, 1:5]
#> 5 x 5 sparse Matrix of class "dgCMatrix"
#>        Cell_1 Cell_2     Cell_3 Cell_4     Cell_5
#> Gene_1      .      . 0.32539945      .  0.9963283
#> Gene_2      .      . 0.06889239      . 23.9966630
#> Gene_3      .      . .               .  .        
#> Gene_4      .      . .               .  .        
#> Gene_5      .      . 0.77355527      .  .        

# Replace them, here with the counts rounded to whole numbers
decontXcounts(sce) <- round(decontaminated)
```
