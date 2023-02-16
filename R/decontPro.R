#' Decontaminate using decontPro
#'
#' @name decontPro
#'
#' @param object Data matrix NxM (feature x droplet).
#' @param cell_type 1xM vector of cell_type. 1-based.
#' @param delta_sd Prior for delta. Set to 2e-5 for now.
#' @param background_sd Prior for background. Set to 2e-6 for now.
#' @param ... Additional arguments for generics.
#'
#' @return A list containing decontaminated counts, and estimated parameters.
#'
#' @examples
#'
NULL
#' NULL


setGeneric("decontPro", function(object,
                                 cell_type,
                                 delta_sd,
                                 background_sd,
                                 ...)
  standardGeneric("decontPro"))


#' @export
#' @rdname decontPro
setMethod("decontPro", "SingleCellExperiment", function(object,
                                                        cell_type,
                                                        delta_sd,
                                                        background_sd,
                                                        ...) {
  counts <- SummarizedExperiment::assay(object, 'counts')
  output <- .decontPro(counts,
                       cell_type,
                       delta_sd,
                       background_sd)

})



#' @export
#' @rdname decontPro
setMethod("decontPro", "Seurat", function(object,
                                          cell_type,
                                          delta_sd,
                                          background_sd,
                                          ...) {
  counts <- Seurat::GetAssayData(object, slot = 'counts')
  output <- .decontPro(counts,
                       cell_type,
                       delta_sd,
                       background_sd)

})



#' @export
#' @rdname decontPro
setMethod("decontPro", "ANY", function(object,
                                                  cell_type,
                                                  delta_sd,
                                                  background_sd,
                                                  ...) {
  output <- .decontPro(object,
                                  cell_type,
                                  delta_sd,
                                  background_sd)

})





.decontPro <- function(counts,
                                  cell_type,
                                  delta_sd,
                                  background_sd) {
  ## Prep data
  N <- nrow(counts)
  M <- ncol(counts)

  p <- rowSums(counts)
  p <- p / sum(p)

  OC <- colSums(counts)

  counts <- as.matrix(counts)

  dat <- list(
    N = N,
    M = M,
    K = length(unique(cell_type)),
    cell_type = cell_type,
    counts = counts,
    OC = OC,
    p = p,
    run_estimation = 1,

    delta_sd = delta_sd,
    background_sd = background_sd
  )


  init <- list(
    delta = matrix(rep(1e-4, N * M),
                   nrow = M,
                   ncol = N),
    background = matrix(rep(1e-2, N * M),
                        nrow = N,
                        ncol = M)
  )




  ## Call inference
  out <- .call_stan_vb(dat, init)



  ## Process Stan output
  re <- .process_stan_vb_out(out, dat)


  return(re)

}
