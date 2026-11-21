# Density of each ADT, raw counts overlapped with decontaminated counts

Both curves of an ADT are smoothed with the same bandwidth, computed
from the raw and decontaminated counts together, so their heights can be
compared.

## Usage

``` r
plotDensity(counts, decontaminated_counts, features, file = NULL)
```

## Arguments

- counts:

  original count matrix of nADT x nDroplet.

- decontaminated_counts:

  decontaminated count matrix, with the droplets of `counts` in the same
  order. If both have column names, they must be the same.

- features:

  names of ADT to plot. Each must be a rowname of both `counts` and
  `decontaminated_counts`.

- file:

  file name to save the plot into as a pdf. ".pdf" is added unless
  `file` already ends in ".pdf". If omitted, return the plot.

## Value

A `patchwork` object with one `ggplot` per feature, or, if `file` is
given, the path of the saved pdf, invisibly.

## Examples

``` r
# Simulate a dataset with 3 cells and 2 ADTs
counts <- matrix(c(60, 72, 52, 49, 89, 112),
  nrow = 2,
  dimnames = list(
    c("CD3", "CD4"),
    c(
      "CTGTTTACACCGCTAG",
      "CTCTACGGTGTGGCTC",
      "AGCAGCCAGGCTCATT"
    )
  )
)

decontaminated_counts <- matrix(c(58, 36, 26, 45, 88, 110),
  nrow = 2,
  dimnames = list(
    c("CD3", "CD4"),
    c(
      "CTGTTTACACCGCTAG",
      "CTCTACGGTGTGGCTC",
      "AGCAGCCAGGCTCATT"
    )
  )
)

plotDensity(
  counts,
  decontaminated_counts,
  c("CD3", "CD4")
)
```
