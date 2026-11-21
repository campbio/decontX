# Retrieve row index for a set of features

This will return indices of features among the rownames or rowData of a
data.frame, matrix, or a
[SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
object including a
[SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html).
Partial matching (i.e. grepping) can be used by setting
`exactMatch = FALSE`.

## Usage

``` r
retrieveFeatureIndex(
  features,
  x,
  by = "rownames",
  exactMatch = TRUE,
  removeNA = FALSE
)
```

## Arguments

- features:

  Character vector of feature names to find in the rows of `x`. Features
  that are `NA` never match.

- x:

  A data.frame, matrix, or
  [SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  object to search.

- by:

  A single character string. Where to search for features in `x`. If set
  to `"rownames"` then the features will be searched for among
  `rownames(x)`. If `x` inherits from class
  [SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html),
  then `by` can be one of the fields in the row annotation data.frame
  (i.e. one of `colnames(rowData(x))`). Default `"rownames"`.

- exactMatch:

  Boolean. Whether to only identify exact matches or to identify partial
  matches using [`grep`](https://rdrr.io/r/base/grep.html). Default
  `TRUE`.

- removeNA:

  Boolean. If set to `FALSE`, features not found in `x` will be given
  `NA` and the returned vector will be the same length as `features`. If
  set to `TRUE`, then the `NA` values will be removed from the returned
  vector. Default `FALSE`.

## Value

A vector of row indices for the matching features in `x`.

## See also

'[retrieveFeatureInfo](https://rdrr.io/pkg/scater/man/retrieveFeatureInfo.html)'
from package `'scater'` and [`regex`](https://rdrr.io/r/base/regex.html)
for how to use regular expressions when `exactMatch = FALSE`.

## Author

Yusuke Koga, Joshua Campbell

## Examples

``` r
# A small count matrix with 20 genes and 10 cells
counts <- matrix(seq_len(200), nrow = 20, ncol = 10,
                 dimnames = list(paste0("Gene_", 1:20),
                                 paste0("Cell_", 1:10)))

# Find genes by their exact names
retrieveFeatureIndex(c("Gene_1", "Gene_5"), counts)
#> [1] 1 5

# Find genes by part of their names
retrieveFeatureIndex(c("ne_5", "ne_7"), counts, exactMatch = FALSE)
#> [1] 5 7
```
