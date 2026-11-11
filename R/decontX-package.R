#' The 'decontX' package.
#'
#' @description DecontX is a decontamination algorithm for single-cell
#' RNA-seq data. DecontPro is a decontamination and background removal
#' algorithm for single cell protein expression data such as CITE-seq or
#' Total-seq.
#'
#' See \code{vignette("decontX", package = "decontX")} and
#' \code{vignette("decontPro", package = "decontX")} to get started, and
#' \code{citation("decontX")} to cite the methods.
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
#' Yang S, Corbett SE, Koga Y, Wang Z, Johnson WE, Yajima M, Campbell JD
#' (2020). Decontamination of ambient RNA in single-cell RNA-seq with
#' DecontX. Genome Biology, 21(1), 57.
#' \doi{10.1186/s13059-020-1950-6}
#'
#' Yin Y, Yajima M, Campbell JD (2024). Characterization and
#' decontamination of background noise in droplet-based single-cell
#' protein expression data with DecontPro. Nucleic Acids Research, 52(1),
#' e4. \doi{10.1093/nar/gkad1032}
#'
#' Stan Development Team (2022). RStan: the R interface to Stan. R package
#' version 2.21.7. https://mc-stan.org
#'
NULL

# Do not remove this import. The compiled Stan model links against the TBB
# library that RcppParallel ships, and importing from RcppParallel loads that
# library before decontX's own. rstantools sets packages up this way.
#' @importFrom RcppParallel RcppParallelLibs
NULL
