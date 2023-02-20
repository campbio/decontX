#' Decontaminate using decontPro
#'
#' @name decontPro
#'
#' @param object Data matrix NxM (feature x droplet).
#' @param cell_type 1xM vector of cell type. 1-based.
#' @param delta_sd Prior variance for ambient contamination level. Default to 2e-5.
#' @param background_sd Prior variance for background contamination level. Default to 2e-6.
#' @param ... Additional arguments for generics.
#'
#' @return A list containing decontaminated counts, and estimated parameters.
#'
#' @examples
#'# # Load CITE-seq dataset
#'# library(DropletUtils)
#'# sce_filtered <- read10xCounts('filtered_feature_bc_matrix')
#'# adt_filtered <- sce_filtered[rowData(sce_filtered)$Type == 'Antibody Capture',]
#'#
#'# # Generate cell clusters
#'# library(Seurat)
#'# library(dplyr)
#'# adt_seurat <- CreateSeuratObject(counts(adt_filtered), assay = 'ADT')
#'# adt_seurat <- NormalizeData(adt_seurat, normalization.method = "CLR",
#'#   margin = 2) %>%
#'#   ScaleData(assay = "ADT") %>%
#'#   RunPCA(assay = "ADT", features = rownames(adt_seurat),
#'#   reduction.name = "pca_adt", npcs = 10) %>%
#'#   FindNeighbors(dims = 1:10, assay = "ADT", reduction = "pca_adt") %>%
#'#   FindClusters(resolution = 0.2)
#'#
#'# clusters = as.integer(Idents(adt_seurat))
#'#
#'# # Run DecontPro
#'# counts <- as.matrix(counts(adt_filtered))
#'# out <- decontPro(counts, clusters)
#'#
#'# decontaminated_counts <- out$decontaminated_counts
#'
#'
NULL
#' NULL


setGeneric("decontPro", function(object,
                                 cell_type,
                                 ...)
  standardGeneric("decontPro"))


#' @export
#' @rdname decontPro
setMethod("decontPro", "SingleCellExperiment", function(object,
                                                        cell_type,
                                                        delta_sd = 2e-5,
                                                        background_sd = 2e-6,
                                                        ...) {
  counts <- SummarizedExperiment::assay(object, 'counts')
  output <- .decontPro(counts = counts,
                       cell_type = cell_type,
                       delta_sd = delta_sd,
                       background_sd = background_sd)

})



#' @export
#' @rdname decontPro
setMethod("decontPro", "Seurat", function(object,
                                          cell_type,
                                          delta_sd = 2e-5,
                                          background_sd = 2e-6,
                                          ...) {
  counts <- Seurat::GetAssayData(object, slot = 'counts')
  output <- .decontPro(counts = counts,
                       cell_type = cell_type,
                       delta_sd = delta_sd,
                       background_sd = background_sd)

})



#' @export
#' @rdname decontPro
setMethod("decontPro", "ANY", function(object,
                                       cell_type,
                                       delta_sd = 2e-5,
                                       background_sd = 2e-6,
                                       ...) {
  output <- .decontPro(counts = object,
                       cell_type = cell_type,
                       delta_sd = delta_sd,
                       background_sd = background_sd)

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
