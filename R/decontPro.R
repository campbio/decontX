#' Decontaminate using decontPro
#'
#' @name decontPro
#'
#' @param filtered_counts Count matrix NxM (feature x droplet) with only
#' filtered droplets after cell calling. The counts must be non-negative
#' whole numbers with no NA. If this is a
#' \linkS4class[SingleCellExperiment]{SingleCellExperiment} or a
#' \linkS4class[SeuratObject]{Seurat} object,
#' expect \code{counts} in the assay slot. Split \code{counts} layers of a
#' Seurat v5 assay are joined first, so the droplets follow
#' \code{Cells()}.
#' @param cell_type Vector of length M giving the cell type of each droplet,
#' with no NA. Any labels can be used, such as integers, character strings
#' or a factor. The rows of \code{r_est} in the returned parameters follow
#' the sorted labels, or a factor's levels with unused levels dropped.
#' @param delta_sd Prior variance for ambient contamination level.
#' Default to 2e-5.
#' @param background_sd Prior variance for background contamination level.
#' Default to 2e-6.
#' @param ambient_counts Count matrix NxM (feature x droplet) with only ambient
#' droplets. Similar to \code{filtered_counts} param, if it is a wrapper object,
#' expect \code{counts} in the assay slot. When both have row names, its
#' rows are matched to the ADTs of \code{filtered_counts} by name: ADTs
#' that aren't in \code{filtered_counts} are dropped, and every ADT of
#' \code{filtered_counts} must be present. Otherwise the rows are matched
#' by position. The ADTs used must have no NA or negative counts.
#' Default to NULL.
#' @param ... Not used. Any other argument stops with an error.
#'
#' @return A list containing decontaminated counts, and estimated parameters.
#' The row names of \code{parameters$r_est}, the native ADT profile of each
#' cell type, are the \code{cell_type} labels.
#'
#' @examples
#' # Simulated count matrix with 100 features x 10 droplets
#' counts <- matrix(
#'   sample(1:10,
#'     1000,
#'     replace = TRUE
#'   ),
#'   ncol = 10
#' )
#'
#' # Cell type indicator
#' k <- c(1, 1, 2, 2, 2, 3, 3, 4, 4, 4)
#'
#' # Simulated ambient count matrix (optional input)
#' ambient_counts <- matrix(
#'   sample(1:2,
#'     1000,
#'     replace = TRUE
#'   ),
#'   ncol = 10
#' )
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
                                 ...) {
  standardGeneric("decontPro")
})


#' @export
#' @rdname decontPro
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
setMethod("decontPro", "SingleCellExperiment", function(filtered_counts,
                                                        cell_type,
                                                        delta_sd = 2e-5,
                                                        background_sd = 2e-6,
                                                        ambient_counts = NULL,
                                                        ...) {
  .checkDecontProDots(...)
  counts <- SummarizedExperiment::assay(filtered_counts, "counts")
  if (!is.null(ambient_counts)) {
    ambient_counts <- SummarizedExperiment::assay(ambient_counts, "counts")
  }
  .decontPro(
    counts = counts,
    cell_type = cell_type,
    delta_sd = delta_sd,
    background_sd = background_sd,
    ambient_counts = ambient_counts
  )
})


#' @export
#' @rdname decontPro
#' @importClassesFrom SeuratObject Seurat
setMethod("decontPro", "Seurat", function(filtered_counts,
                                          cell_type,
                                          delta_sd = 2e-5,
                                          background_sd = 2e-6,
                                          ambient_counts = NULL,
                                          ...) {
  .checkDecontProDots(...)
  counts <- .seuratCounts(filtered_counts)
  if (!is.null(ambient_counts)) {
    ambient_counts <- .seuratCounts(ambient_counts)
  }
  .decontPro(
    counts = counts,
    cell_type = cell_type,
    delta_sd = delta_sd,
    background_sd = background_sd,
    ambient_counts = ambient_counts
  )
})


#' @export
#' @rdname decontPro
setMethod("decontPro", "ANY", function(filtered_counts,
                                       cell_type,
                                       delta_sd = 2e-5,
                                       background_sd = 2e-6,
                                       ambient_counts = NULL,
                                       ...) {
  .checkDecontProDots(...)
  .decontPro(
    counts = filtered_counts,
    cell_type = cell_type,
    delta_sd = delta_sd,
    background_sd = background_sd,
    ambient_counts = ambient_counts
  )
})


.decontPro <- function(counts,
                       cell_type,
                       delta_sd,
                       background_sd,
                       ambient_counts) {
  ## Prep data
  counts <- as.matrix(counts)
  .checkDecontProCounts(counts)
  N <- nrow(counts)
  M <- ncol(counts)
  cellType <- .recodeCellType(cell_type, M)

  if (is.null(ambient_counts)) {
    temp <- counts
  } else {
    temp <- .matchBackgroundRows(
      counts = counts,
      background = ambient_counts,
      countsArg = "filtered_counts",
      bgArg = "ambient_counts",
      feature = "ADTs"
    )
    .checkCountsDecon(temp, "ambient_counts")
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

  dat <- list(
    N = N,
    M = M,
    K = length(cellType$labels),
    cell_type = cellType$codes,
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
      ncol = N
    ),
    background = matrix(rep(1e-2, N * M),
      nrow = N,
      ncol = M
    )
  )


  ## Call inference
  out <- .call_stan_vb(dat, init)


  ## Process Stan output
  dat$cell_type_labels <- cellType$labels
  re <- .process_stan_vb_out(out, dat)


  re
}


# Reads the counts of the default assay. Split Seurat v5 layers are joined
# on the extracted assay, so the columns follow Cells() and the user's
# object is left as it is.
.seuratCounts <- function(x) {
  assay <- x[[SeuratObject::DefaultAssay(x)]]
  countLayers <- SeuratObject::Layers(assay, search = "counts")
  if (methods::is(assay, "Assay5") && length(countLayers) > 1) {
    assay <- SeuratObject::JoinLayers(assay)
  }
  SeuratObject::LayerData(assay, layer = "counts")
}


# The Stan model declares the counts as non-negative integers
.checkDecontProCounts <- function(counts) {
  if (!is.numeric(counts)) {
    stop("'filtered_counts' must be numeric, not ", typeof(counts), ".")
  }
  if (anyNA(counts)) {
    stop("'filtered_counts' has NA values.")
  }
  if (any(counts < 0)) {
    stop("'filtered_counts' has negative values.")
  }
  # Integer matrices hold whole numbers already
  wholeNumbers <- is.integer(counts) ||
    all(is.finite(counts) & counts == round(counts))
  if (!wholeNumbers) {
    stop("'filtered_counts' must hold whole-number counts.")
  }
}


# Stan needs cell types coded 1..K with every code used. Labels are sorted,
# or kept in a factor's level order with unused levels dropped.
.recodeCellType <- function(cell_type, M) {
  if (length(cell_type) != M) {
    stop(
      "'cell_type' has ", length(cell_type), " entries but ",
      "'filtered_counts' has ", M, " droplets (columns)."
    )
  }
  if (anyNA(cell_type)) {
    stop("'cell_type' has NA values. Give every droplet a cell type.")
  }
  if (is.factor(cell_type)) {
    cell_type <- droplevels(cell_type)
  } else {
    cell_type <- factor(cell_type)
  }
  list(
    codes = as.integer(cell_type),
    labels = levels(cell_type)
  )
}


.checkDecontProDots <- function(...) {
  if (...length() == 0) {
    return(invisible(NULL))
  }
  argNames <- names(substitute(list(...)))[-1]
  if (is.null(argNames)) {
    argNames <- rep("", ...length())
  }
  argNames <- ifelse(argNames == "", "an unnamed value",
    paste0("'", argNames, "'")
  )
  stop(
    "Unknown arguments to decontPro(): ",
    paste(argNames, collapse = ", "), ". See ?decontPro."
  )
}
