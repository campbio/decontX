
<!-- README.md is generated from README.Rmd. Please edit that file -->

# decontX

<!-- badges: start -->

[![R-CMD-check](https://github.com/campbio/decontX/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/campbio/decontX/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/campbio/decontX/graph/badge.svg?token=z5YnsXuWqh)](https://codecov.io/gh/campbio/decontX)
<!-- badges: end -->

Methods for decontamination of single cell data. This package implements
both DecontX (Yang et al.,
[2020](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1950-6))
for single-cell RNA-seq data and DecontPro (Yin et al.,
[2023](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkad1032/7420100)) for
single-cell protein expression data.

## Installation Instructions

You can install the package through [Bioconductor](https://bioconductor.org/packages/decontX) with:

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("decontX")
```

Alternatively, install the development version through
[GitHub](https://github.com/campbio/decontX) using `devtools` package:

``` r
# install.packages("devtools")
devtools::install_github("campbio/decontX")
```

## Vignettes

Vignettes are available on [Bioconductor](https://bioconductor.org/packages/decontX).

To build vignette when installing from GitHub, add the flag `build_vignettes = TRUE`:

``` r
library(devtools)
install_github("campbio/decontX", build_vignettes = TRUE)
```

Then vignettes can be accessed through:

``` r
vignette('decontX', package = 'decontX')
vignette('decontPro', package = 'decontX')
```
