#' The 'decontX' package.
#'
#' @description DecontX is a decontamination algorithm for single-cell
#' RNA-seq data. DecontPro is a decontamination and background removal
#' algorithm for single cell protein expression data such as CITE-seq or
#' Total-seq.
#'
#' @return None; this page documents the package.
#' @name decontX-package
#' @aliases decontX-package
#' @useDynLib decontX, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import rstantools
#' @importFrom rstan sampling
#' @importFrom ggplot2 .data
#'
#' @references
#' Stan Development Team (2022). RStan: the R interface to Stan. R package
#' version 2.21.7. https://mc-stan.org
#'
NULL

# Do not remove this import. The compiled Stan model links against the TBB
# library that RcppParallel ships, and importing from RcppParallel loads that
# library before decontX's own. rstantools sets packages up this way.
#' @importFrom RcppParallel RcppParallelLibs
NULL
