# Simulate contaminated count matrix

This function generates a list containing two count matrices – one for
real expression, the other one for contamination, as well as other
parameters used in the simulation which can be useful for running
decontamination.

## Usage

``` r
simulateContamination(
  C = 300,
  G = 100,
  K = 3,
  NRange = c(500, 1000),
  beta = 0.1,
  delta = c(1, 10),
  numMarkers = 3,
  seed = 12345
)
```

## Arguments

- C:

  Integer. Number of cells to be simulated. Default `300`.

- G:

  Integer. Number of genes to be simulated. Default `100`.

- K:

  Integer. Number of cell populations to be simulated, at least 2, since
  the contamination comes from the other populations. With few cells,
  fewer populations may be drawn; at least 2 are still needed. Default
  `3`.

- NRange:

  Integer vector. A vector of length 2 that specifies the lower and
  upper bounds of the number of counts generated for each cell. Default
  `c(500, 1000)`.

- beta:

  Numeric. Concentration parameter for Phi. Default `0.1`.

- delta:

  Numeric or Numeric vector. Concentration parameter for Theta. If input
  as a single numeric value, symmetric values for beta distribution are
  specified; if input as a vector of length 2, the two values will be
  the shape1 and shape2 parameters of the beta distribution
  respectively. Default `c(1, 10)`.

- numMarkers:

  Integer. Number of markers for each cell population. Default `3`.

- seed:

  Integer. Passed to
  [`with_seed`](https://withr.r-lib.org/reference/with_seed.html). For
  reproducibility, a default value of 12345 is used. If NULL, no calls
  to [`with_seed`](https://withr.r-lib.org/reference/with_seed.html) are
  made.

## Value

A list with the following items:

- `nativeCounts`::

  Matrix of the native (real) counts, genes by cells.

- `observedCounts`::

  Matrix of the observed counts: the native counts plus the
  contamination. Both count matrices have gene and cell names.

- `NByC`::

  Integer vector. The total observed counts of each cell.

- `z`::

  Integer vector. The cell population of each cell.

- `eta`::

  Matrix of the contamination distribution of each population, genes by
  populations.

- `phi`::

  Matrix of the native expression distribution of each population, genes
  by populations.

- `markers`::

  List of the marker gene names of each population.

- `numMarkers`::

  The number of markers for each population.

- `contamination`::

  Numeric vector. The true proportion of contamination in each cell.

## Author

Shiyi Yang, Yuan Yin, Joshua Campbell

## Examples

``` r
contaminationSim <- simulateContamination(K = 3, delta = c(1, 10))
```
