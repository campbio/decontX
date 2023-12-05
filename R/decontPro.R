#' Decontaminate using decontPro
#'
#' @name decontPro
#'
#' @param filtered_counts Count matrix NxM (feature x droplet) with only 
#' filtered droplets after cell calling. If this is a
#' \linkS4class{SingleCellExperiment} or a \linkS4class{Seurat} object, 
#' expect \code{counts} in the assay slot.
#' @param cell_type 1xM 1-based integer vector indicating cell type of each 
#' droplet.
#' @param delta_sd Prior variance for ambient contamination level.
#' Default to 2e-5.
#' @param background_sd Prior variance for background contamination level.
#' Default to 2e-6.
#' @param ambient_counts Count matrix NxM (feature x droplet) with only ambient
#' droplets. Similar to \code{filtered_counts} param, if it is a wrapper object, 
#' expect \code{counts} in the assay slot. Default to NULL.
#' @param ... Additional arguments for generics.
#'
#' @return A list containing decontaminated counts, and estimated parameters.
#'
#' @examples
#' # Simulated count matrix with 100 features x 10 droplets
#' counts <- matrix(sample(1:10,
#'                         1000,
#'                         replace = TRUE),
#'                  ncol = 10)
#'
#' # Cell type indicator
#' k <- c(1, 1, 2, 2, 2, 3, 3, 4, 4, 4)
#' 
#' # Simulated ambient count matrix (optional input)
#' ambient_counts <- matrix(sample(1:2,
#'                         1000,
#'                         replace = TRUE),
#'                  ncol = 10)
#'
#' # Decontamination
#' out <- decontPro(counts, k, 1e-2, 1e-2, ambient_counts)
#'
#' # Decontaminated counts
#' decontaminated_counts <- out$decontaminated_counts
NULL
#' NULL


setGeneric("decontPro", function(filtered_counts,
                                 cell_type,
                                 ...)
  standardGeneric("decontPro"))


#' @export
#' @rdname decontPro
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
setMethod("decontPro", "SingleCellExperiment", function(filtered_counts,
                                                        cell_type,
                                                        delta_sd = 2e-5,
                                                        background_sd = 2e-6,
                                                        ambient_counts = NULL,
                                                        ...) {
  counts <- SummarizedExperiment::assay(filtered_counts, 'counts')
  if (!is.null(ambient_counts)) {
    ambient_counts <- SummarizedExperiment::assay(ambient_counts, 'counts')
  }
  output <- .decontPro(counts = counts,
                       cell_type = cell_type,
                       delta_sd = delta_sd,
                       background_sd = background_sd,
                       ambient_counts = ambient_counts)

})



#' @export
#' @rdname decontPro
#' @importClassesFrom Seurat Seurat
setMethod("decontPro", "Seurat", function(filtered_counts,
                                          cell_type,
                                          delta_sd = 2e-5,
                                          background_sd = 2e-6,
                                          ambient_counts = NULL,
                                          ...) {
  counts <- Seurat::GetAssayData(filtered_counts, slot = 'counts')
  if (!is.null(ambient_counts)){
    ambient_counts <- Seurat::GetAssayData(ambient_counts, slot = 'counts')
  }
  output <- .decontPro(counts = counts,
                       cell_type = cell_type,
                       delta_sd = delta_sd,
                       background_sd = background_sd,
                       ambient_counts = ambient_counts)

})



#' @export
#' @rdname decontPro
setMethod("decontPro", "ANY", function(filtered_counts,
                                       cell_type,
                                       delta_sd = 2e-5,
                                       background_sd = 2e-6,
                                       ambient_counts = NULL,
                                       ...) {
  output <- .decontPro(counts = filtered_counts,
                       cell_type = cell_type,
                       delta_sd = delta_sd,
                       background_sd = background_sd,
                       ambient_counts = ambient_counts)

})





.decontPro <- function(counts,
                       cell_type,
                       delta_sd,
                       background_sd,
                       ambient_counts) {
  ## Prep data
  N <- nrow(counts)
  M <- ncol(counts)
  
  if (is.null(ambient_counts)){
    temp <- counts
    
  } else {
    
    # Quick check ambient counts have same num of ADTs as counts
    if (nrow(counts) != nrow(ambient_counts)){
      stop(
        "Input ambient_counts and counts have different number of ADTs/rows."
      )
    }
    temp <- ambient_counts
  }
  
  p <- Matrix::rowSums(temp)
  p <- p / sum(p)
  

  OC <- Matrix::colSums(counts)

  if (any(OC == 0)) {
    stop(
      "Droplets with 0 total counts detected.",
      " Remove them and try again."
    )
  }

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
