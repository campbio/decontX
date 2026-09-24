#' @title Contamination estimation with decontX
#'
#' @description Identifies contamination from factors such as ambient RNA
#' in single cell genomic datasets.
#'
#' @name decontX
#'
#' @param x A numeric matrix of counts or a \linkS4class{SingleCellExperiment}
#' with the matrix located in the assay slot under \code{assayName}.
#' Cells in each batch will be subsetted and converted to a sparse matrix
#' of class \code{dgCMatrix} from package \link{Matrix} before analysis. This
#' object should only contain filtered cells after cell calling. Empty
#' cell barcodes (low expression droplets before cell calling) are not needed
#' to run DecontX.
#' @param assayName Character. Name of the assay to use if \code{x} is a
#' \linkS4class{SingleCellExperiment}.
#' @param z Numeric or character vector. Cell cluster labels. If NULL,
#' PCA will be used to reduce the dimensionality of the dataset initially,
#' '\link[scrapper]{runUmap}' from the 'scrapper' package
#' will be used to further reduce the dataset to 2 dimensions and
#' the '\link[dbscan]{dbscan}' function from the 'dbscan' package
#' will be used to identify clusters of broad cell types. Default NULL.
#' @param batch Numeric or character vector. Batch labels for cells.
#' If batch labels are supplied, DecontX is run on cells from each
#' batch separately. Cells run in different channels or assays
#' should be considered different batches. Default NULL.
#' @param background A numeric matrix of counts or a
#' \linkS4class{SingleCellExperiment} with the matrix located in the assay
#' slot under \code{assayName}. It should have the same data format as \code{x}
#' except it contains the empty droplets instead of cells. When supplied,
#' empirical distribution of transcripts from these empty droplets
#' will be used as the contamination distribution. Default NULL.
#' @param bgAssayName Character. Name of the assay to use if \code{background}
#' is a \linkS4class{SingleCellExperiment}. Default to same as
#' \code{assayName}.
#' @param bgBatch Numeric or character vector. Batch labels for
#' \code{background}. Its unique values should be the same as those in
#' \code{batch}, such that each batch of cells have their corresponding batch
#' of empty droplets as background, pointed by this parameter. Default to NULL.
#' @param maxIter Integer. Maximum iterations of the EM algorithm. Default 500.
#' @param convergence Numeric. The EM algorithm will be stopped if the maximum
#' difference in the contamination estimates between the previous and
#' current iterations is less than this. Default 0.001.
#' @param iterLogLik Integer. Calculate log likelihood every \code{iterLogLik}
#' iteration. Default 10.
#' @param delta Numeric Vector of length 2. Concentration parameters for
#' the Dirichlet prior for the contamination in each cell. The first element
#' is the prior for the native counts while the second element is the prior for
#' the contamination counts. These essentially act as pseudocounts for the
#' native and contamination in each cell. If \code{estimateDelta = TRUE},
#' this is only used to produce a random sample of proportions for an initial
#' value of contamination in each cell. Then
#' \code{\link[MCMCprecision]{fit_dirichlet}} is used to update
#' \code{delta} in each iteration.
#' If \code{estimateDelta = FALSE}, then \code{delta} is fixed with these
#' values for the entire inference procedure. Fixing \code{delta} and
#' setting a high number in the second element will force \code{decontX}
#' to be more aggressive and estimate higher levels of contamination at
#' the expense of potentially removing native expression.
#' Default \code{c(10, 10)}.
#' @param estimateDelta Boolean. Whether to update \code{delta} at each
#' iteration.
#' @param varGenes Integer. The number of variable genes to use in
#' dimensionality reduction before clustering. Variability is calculated
#' using \code{\link[scrapper]{modelGeneVariances}} from the 'scrapper'
#' package (or by the variance of the log-normalized counts when
#' \code{legacyInit = TRUE}). Used only when z is not provided.
#' Default 5000.
#' @param dbscanEps Numeric. The clustering resolution parameter
#' used in '\link[dbscan]{dbscan}' to estimate broad cell clusters.
#' Used only when z is not provided. Default 1.
#' @param seed Integer. Seed for reproducibility: the estimation for each
#'  batch (including the EM algorithm) is wrapped in
#'  \link[withr]{with_seed}, and the 'scrapper' UMAP used for cell
#'  cluster initialization is seeded explicitly. A default value of 12345
#'  is used. If NULL, no seeding is performed and results are not
#'  reproducible between runs (only the scrapper UMAP itself, which has
#'  fixed internal default seeds, stays deterministic).
#' @param legacyInit Logical. Use the original scater/scuttle-based
#'  initialization (log-normalization and UMAP) to reproduce results from
#'  decontX 1.11.0 and earlier. This option requires the 'scater' package
#'  and may emit deprecation warnings from scater/scuttle, whose
#'  normalization functions were deprecated in favor of the 'scrapper'
#'  package. It is kept for backwards compatibility and will remain
#'  available as long as the upstream functions are still provided by
#'  scater/scuttle. Used only when z is not provided. Default FALSE.
#' @param logfile Character. Messages will be redirected to a file named
#'  `logfile`. If NULL, messages will be printed to stdout.  Default NULL.
#' @param verbose Logical. Whether to print log messages. Default TRUE.
#' @param ... For the generic, further arguments to pass to each method.
#'
#' @return If \code{x} is a matrix-like object, a list will be returned
#' with the following items:
#' \describe{
#' \item{\code{decontXcounts}:}{The decontaminated matrix. Values obtained
#' from the variational inference procedure may be non-integer. However,
#' integer counts can be obtained by rounding,
#' e.g. \code{round(decontXcounts)}.}
#' \item{\code{contamination}:}{Percentage of contamination in each cell.}
#' \item{\code{estimates}:}{List of estimated parameters for each batch. If z
#' was not supplied, then the UMAP coordinates used to generated cell
#' cluster labels will also be stored here.}
#' \item{\code{z}:}{Cell population/cluster labels used for analysis.}
#' \item{\code{runParams}:}{List of arguments used in the function call.}
#' }
#'
#' If \code{x} is a \linkS4class{SingleCellExperiment}, then the decontaminated
#' counts will be stored as an assay and can be accessed with
#' \code{decontXcounts(x)}. The contamination values and cluster labels
#' will be stored in \code{colData(x)}. \code{estimates} and \code{runParams}
#' will be stored in \code{metadata(x)$decontX}. The UMAPs used to generated
#' cell cluster labels will be stored in
#' \code{reducedDims} slot in \code{x}.
#'
#' @author Shiyi Yang, Yuan Yin, Joshua Campbell
#'
#' @example man/examples/decontX.R
#'
NULL

#' @export
#' @rdname decontX
setGeneric("decontX", function(x, ...) standardGeneric("decontX"))


#########################
# Setting up S4 methods #
#########################


#' @export
#' @rdname decontX
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importClassesFrom Matrix dgCMatrix
setMethod("decontX", "SingleCellExperiment", function(x,
                                                      assayName = "counts",
                                                      z = NULL,
                                                      batch = NULL,
                                                      background = NULL,
                                                      bgAssayName = NULL,
                                                      bgBatch = NULL,
                                                      maxIter = 500,
                                                      delta = c(10, 10),
                                                      estimateDelta = TRUE,
                                                      convergence = 0.001,
                                                      iterLogLik = 10,
                                                      varGenes = 5000,
                                                      dbscanEps = 1,
                                                      seed = 12345,
                                                      legacyInit = FALSE,
                                                      logfile = NULL,
                                                      verbose = TRUE) {
  countsBackground <- NULL
  if (!is.null(background)) {
    # Remove cells with the same ID between x and the background matrix
    # Also update bgBatch when background is updated and bgBatch is not null
    temp <- .checkBackground(
      x = x,
      background = background,
      bgBatch = bgBatch,
      logfile = logfile,
      verbose = verbose
    )

    background <- temp$background
    bgBatch <- temp$bgBatch

    if (is.null(bgAssayName)) {
      bgAssayName <- assayName
    }
    countsBackground <- SummarizedExperiment::assay(background, i = bgAssayName)
  }

  mat <- SummarizedExperiment::assay(x, i = assayName)

  result <- .decontX(
    counts = mat,
    z = z,
    batch = batch,
    countsBackground = countsBackground,
    batchBackground = bgBatch,
    maxIter = maxIter,
    convergence = convergence,
    iterLogLik = iterLogLik,
    delta = delta,
    estimateDelta = estimateDelta,
    varGenes = varGenes,
    dbscanEps = dbscanEps,
    seed = seed,
    legacyInit = legacyInit,
    logfile = logfile,
    verbose = verbose
  )

  ## Add results into column annotation
  SummarizedExperiment::colData(x)$decontX_contamination <- result$contamination
  SummarizedExperiment::colData(x)$decontX_clusters <- as.factor(result$z)

  ## Put estimated UMAPs into SCE
  batchIndex <- unique(result$runParams$batch)
  if (length(batchIndex) > 1) {
    for (i in batchIndex) {
      ## Each individual UMAP will only be for one batch so need
      ## to put NAs in for cells in other batches
      tempUMAP <- matrix(NA, ncol = 2, nrow = ncol(mat))
      tempUMAP[result$runParams$batch == i, ] <- result$estimates[[i]]$UMAP
      colnames(tempUMAP) <- c("UMAP_1", "UMAP_2")
      rownames(tempUMAP) <- colnames(mat)

      SingleCellExperiment::reducedDim(
        x,
        paste("decontX", i, "UMAP", sep = "_")
      ) <- tempUMAP
    }
  } else {
    SingleCellExperiment::reducedDim(x, "decontX_UMAP") <-
      result$estimates[[batchIndex]]$UMAP
  }

  ## Save the rest of the result object into metadata
  decontXcounts(x) <- result$decontXcounts
  result$decontXcounts <- NULL
  S4Vectors::metadata(x)$decontX <- result

  x
})

#' @export
#' @rdname decontX
setMethod("decontX", "ANY", function(x,
                                     z = NULL,
                                     batch = NULL,
                                     background = NULL,
                                     bgBatch = NULL,
                                     maxIter = 500,
                                     delta = c(10, 10),
                                     estimateDelta = TRUE,
                                     convergence = 0.001,
                                     iterLogLik = 10,
                                     varGenes = 5000,
                                     dbscanEps = 1,
                                     seed = 12345,
                                     legacyInit = FALSE,
                                     logfile = NULL,
                                     verbose = TRUE) {
  countsBackground <- NULL
  if (!is.null(background)) {
    # Remove cells with the same ID between x and the background matrix
    # Also update bgBatch when background is updated and bgBatch is not null
    temp <- .checkBackground(
      x = x,
      background = background,
      bgBatch = bgBatch,
      logfile = logfile,
      verbose = verbose
    )

    background <- temp$background
    countsBackground <- background

    bgBatch <- temp$bgBatch
  }

  .decontX(
    counts = x,
    z = z,
    batch = batch,
    countsBackground = countsBackground,
    batchBackground = bgBatch,
    maxIter = maxIter,
    convergence = convergence,
    iterLogLik = iterLogLik,
    delta = delta,
    estimateDelta = estimateDelta,
    varGenes = varGenes,
    dbscanEps = dbscanEps,
    seed = seed,
    legacyInit = legacyInit,
    logfile = logfile,
    verbose = verbose
  )
})


## Copied from SingleCellExperiment Package

GET_FUN <- function(exprs_values, ...) {
  (exprs_values) # To ensure evaluation
  function(object, ...) {
    SummarizedExperiment::assay(object, i = exprs_values, ...)
  }
}

SET_FUN <- function(exprs_values, ...) {
  (exprs_values) # To ensure evaluation
  function(object, ..., value) {
    SummarizedExperiment::assay(object, i = exprs_values, ...) <- value
    object
  }
}


#' @title Get or set decontaminated counts matrix
#'
#' @description Gets or sets the decontaminated counts matrix from a
#' a \linkS4class{SingleCellExperiment} object.
#' @name decontXcounts
#' @param object A \linkS4class{SingleCellExperiment} object.
#' @param value A matrix to save as an assay called \code{decontXcounts}
#' @param ... For the generic, further arguments to pass to each method.
#' @return If getting, the assay from \code{object} with the name
#' \code{decontXcounts} will be returned. If setting, a
#' \linkS4class{SingleCellExperiment} object will be returned with
#' \code{decontXcounts} listed in the \code{assay} slot.
#' @seealso \code{\link{assay}} and \code{\link{assay<-}}
NULL

#' @export
#' @rdname decontXcounts
setGeneric("decontXcounts", function(object, ...) {
  standardGeneric("decontXcounts")
})

#' @export
#' @rdname decontXcounts
setGeneric("decontXcounts<-", function(object, ..., value) {
  standardGeneric("decontXcounts<-")
})

#' @export
#' @rdname decontXcounts
setMethod("decontXcounts", "SingleCellExperiment", GET_FUN("decontXcounts"))

#' @export
#' @rdname decontXcounts
setMethod(
  "decontXcounts<-", c("SingleCellExperiment", "ANY"),
  SET_FUN("decontXcounts")
)


##########################
# Core Decontx Functions #
##########################

.decontX <- function(counts,
                     z = NULL,
                     batch = NULL,
                     countsBackground = NULL,
                     batchBackground = NULL,
                     maxIter = 200,
                     convergence = 0.001,
                     iterLogLik = 10,
                     delta = c(10, 10),
                     estimateDelta = TRUE,
                     varGenes = NULL,
                     dbscanEps = NULL,
                     seed = 12345,
                     legacyInit = FALSE,
                     logfile = NULL,
                     verbose = TRUE) {
  startTime <- Sys.time()
  .logMessages(paste(rep("-", 50), collapse = ""),
    logfile = logfile,
    append = TRUE,
    verbose = verbose
  )
  .logMessages("Starting DecontX",
    logfile = logfile,
    append = TRUE,
    verbose = verbose
  )
  .logMessages(paste(rep("-", 50), collapse = ""),
    logfile = logfile,
    append = TRUE,
    verbose = verbose
  )

  runParams <- list(
    z = z,
    batch = batch,
    batchBackground = batchBackground,
    maxIter = maxIter,
    delta = delta,
    estimateDelta = estimateDelta,
    convergence = convergence,
    varGenes = varGenes,
    dbscanEps = dbscanEps,
    legacyInit = legacyInit,
    logfile = logfile,
    verbose = verbose
  )


  totalGenes <- nrow(counts)
  totalCells <- ncol(counts)
  geneNames <- rownames(counts)
  nC <- ncol(counts)
  allCellNames <- colnames(counts)

  ## Set up final decontaminated matrix
  estRmat <- Matrix::Matrix(
    data = 0,
    ncol = totalCells,
    nrow = totalGenes,
    sparse = TRUE,
    dimnames = list(geneNames, allCellNames)
  )

  ## Generate batch labels if none were supplied
  if (is.null(batch)) {
    batch <- rep("all_cells", nC)

    # If batch null, bgBatch has to be null
    if (!is.null(batchBackground)) {
      stop(
        "When experiment default to no bacth, background should ",
        "also default to no batch."
      )
    }

    if (!is.null(countsBackground)) {
      batchBackground <- rep("all_cells", ncol(countsBackground))
    }
  } else {
    # If batch not null and countsBackground supplied,
    # user has to supply batchBackground as well
    if (!is.null(countsBackground) && is.null(batchBackground)) {
      stop(
        "Cell batch, and background are supplied. Please also ",
        "supply background batch."
      )
    }
  }
  batch <- as.character(batch)
  if (!is.null(batchBackground)) {
    batchBackground <- as.character(batchBackground)
  }
  runParams$batch <- batch
  runParams$batchBackground <- batchBackground
  batchIndex <- unique(batch)

  ## Set result lists upfront for all cells from different batches
  estConp <- rep(NA, nC)
  returnZ <- rep(NA, nC)
  resBatch <- list()

  ## Cycle through each sample/batch and run DecontX
  for (bat in batchIndex) {
    if (length(batchIndex) == 1) {
      .logMessages(
        date(),
        ".. Analyzing all cells",
        logfile = logfile,
        append = TRUE,
        verbose = verbose
      )
    } else {
      .logMessages(
        date(),
        " .. Analyzing cells in batch '",
        bat, "'",
        sep = "",
        logfile = logfile,
        append = TRUE,
        verbose = verbose
      )
    }

    zBat <- NULL
    countsBat <- counts[, batch == bat]
    bgBat <- countsBackground[, batchBackground == bat]

    ## Convert to sparse matrix
    if (!inherits(countsBat, "dgCMatrix")) {
      .logMessages(
        date(),
        ".... Converting to sparse matrix",
        logfile = logfile,
        append = TRUE,
        verbose = verbose
      )
      countsBat <- methods::as(countsBat, "dgCMatrix")
    }
    if (!is.null(bgBat)) {
      if (!inherits(bgBat, "dgCMatrix")) {
        bgBat <- methods::as(bgBat, "dgCMatrix")
      }
    }

    if (!is.null(z)) {
      zBat <- z[batch == bat]
    }
    if (is.null(seed)) {
      res <- .decontXoneBatch(
        counts = countsBat,
        z = zBat,
        batch = bat,
        countsBackground = bgBat,
        maxIter = maxIter,
        delta = delta,
        estimateDelta = estimateDelta,
        convergence = convergence,
        iterLogLik = iterLogLik,
        logfile = logfile,
        verbose = verbose,
        varGenes = varGenes,
        dbscanEps = dbscanEps,
        seed = seed,
        legacyInit = legacyInit
      )
    } else {
      withr::with_seed(
        seed,
        res <- .decontXoneBatch(
          counts = countsBat,
          z = zBat,
          batch = bat,
          countsBackground = bgBat,
          maxIter = maxIter,
          delta = delta,
          estimateDelta = estimateDelta,
          convergence = convergence,
          iterLogLik = iterLogLik,
          logfile = logfile,
          verbose = verbose,
          varGenes = varGenes,
          dbscanEps = dbscanEps,
          seed = seed,
          legacyInit = legacyInit
        )
      )
    }

    ## Try to convert class of new matrix to class of original matrix

    .logMessages(
      date(),
      ".. Calculating final decontaminated matrix",
      logfile = logfile,
      append = TRUE,
      verbose = verbose
    )

    estRmat.temp <- calculateNativeMatrix(
      counts = countsBat,
      theta = res$theta,
      eta = res$eta,
      phi = res$phi,
      z = as.integer(res$z),
      pseudocount = 1e-20
    )

    # Speed up sparse matrix value assignment by cbind -> order recovery
    allCol <- paste0("col_", seq_len(ncol(estRmat)))
    colnames(estRmat) <- allCol

    subCol <- paste0("col_", which(batch == bat))
    colnames(estRmat.temp) <- subCol

    estRmat <- estRmat[, !(allCol %in% subCol)]
    estRmat <- cbind(estRmat, estRmat.temp)

    # Recover order and set names
    estRmat <- estRmat[, allCol]
    dimnames(estRmat) <- list(geneNames, allCellNames)

    resBatch[[bat]] <- list(
      z = res$z,
      phi = res$phi,
      eta = res$eta,
      delta = res$delta,
      theta = res$theta,
      contamination = res$contamination,
      logLikelihood = res$logLikelihood,
      UMAP = res$UMAP,
      z = res$z,
      iteration = res$iteration
    )

    estConp[batch == bat] <- res$contamination
    if (length(batchIndex) > 1) {
      returnZ[batch == bat] <- paste0(bat, "-", res$z)
    } else {
      returnZ[batch == bat] <- res$z
    }
  }
  names(resBatch) <- batchIndex

  returnResult <- list(
    "runParams" = runParams,
    "estimates" = resBatch,
    "decontXcounts" = estRmat,
    "contamination" = estConp,
    "z" = returnZ
  )


  if (inherits(counts, c("DelayedMatrix", "DelayedArray"))) {
    .logMessages(
      date(),
      ".. Converting decontaminated matrix to", class(counts),
      logfile = logfile,
      append = TRUE,
      verbose = verbose
    )

    ## Determine class of seed in DelayedArray
    seed.class <- unique(DelayedArray::seedApply(counts, class))[[1]]
    if (seed.class == "HDF5ArraySeed") {
      returnResult$decontXcounts <-
        methods::as(returnResult$decontXcounts, "HDF5Matrix")
    } else {
      if (isTRUE(methods::canCoerce(returnResult$decontXcounts, seed.class))) {
        returnResult$decontXcounts <-
          methods::as(returnResult$decontXcounts, seed.class)
      }
    }
    returnResult$decontXcounts <-
      DelayedArray::DelayedArray(returnResult$decontXcounts)
  } else {
    try(
      {
        if (methods::canCoerce(returnResult$decontXcounts, class(counts))) {
          returnResult$decontXcounts <-
            methods::as(returnResult$decontXcounts, class(counts))
        }
      },
      silent = TRUE
    )
  }

  endTime <- Sys.time()
  .logMessages(paste(rep("-", 50), collapse = ""),
    logfile = logfile,
    append = TRUE,
    verbose = verbose
  )
  .logMessages("Completed DecontX. Total time:",
    format(difftime(endTime, startTime)),
    logfile = logfile,
    append = TRUE,
    verbose = verbose
  )
  .logMessages(paste(rep("-", 50), collapse = ""),
    logfile = logfile,
    append = TRUE,
    verbose = verbose
  )

  returnResult
}


# This function updates decontamination for one batch
# seed passed to this function is to be furhter passed to
# function .decontxInitializeZ()
#
# Do not remove the @importFrom below: fit_dirichlet is fetched from the
# MCMCprecision namespace by the C++ EM step (src/DecontX.cpp, via
# Environment::namespace_env), which R's static analysis cannot see.
# Dropping the import triggers an R CMD check NOTE and can leave
# MCMCprecision unloaded when decontXEM() runs.
#' @importFrom MCMCprecision fit_dirichlet
#' @noRd
.decontXoneBatch <- function(counts,
                             z = NULL,
                             batch = NULL,
                             countsBackground = NULL,
                             maxIter = 200,
                             delta = c(10, 10),
                             estimateDelta = TRUE,
                             convergence = 0.01,
                             iterLogLik = 10,
                             logfile = NULL,
                             verbose = TRUE,
                             varGenes = NULL,
                             dbscanEps = NULL,
                             seed = 12345,
                             legacyInit = FALSE) {
  .checkCountsDecon(counts)
  .checkDelta(delta)

  nC <- ncol(counts)
  deconMethod <- "clustering"

  ## Generating UMAP and cell cluster labels if none are provided
  umap <- NULL
  if (is.null(z)) {
    m <- ".... Generating UMAP and estimating cell types"
    estimateCellTypes <- TRUE
  } else {
    m <- ".... Generating UMAP"
    estimateCellTypes <- FALSE
  }
  .logMessages(
    date(),
    m,
    logfile = logfile,
    append = TRUE,
    verbose = verbose
  )

  varGenes <- .processvarGenes(varGenes)
  dbscanEps <- .processdbscanEps(dbscanEps)

  celda.init <- .decontxInitializeZ(
    counts = counts,
    varGenes = varGenes,
    dbscanEps = dbscanEps,
    estimateCellTypes = estimateCellTypes,
    seed = seed,
    legacyInit = legacyInit
  )
  if (is.null(z)) {
    z <- celda.init$z
  }
  umap <- celda.init$umap
  colnames(umap) <- c(
    "DecontX_UMAP_1",
    "DecontX_UMAP_2"
  )
  rownames(umap) <- colnames(counts)

  z <- .processCellLabels(z, numCells = nC)

  iter <- 1L
  numIterWithoutImprovement <- 0L
  stopIter <- 3L

  .logMessages(
    date(),
    ".... Estimating contamination",
    logfile = logfile,
    append = TRUE,
    verbose = verbose
  )

  if (deconMethod == "clustering") {
    ## Initialization
    theta <- stats::rbeta(
      n = nC,
      shape1 = delta[1],
      shape2 = delta[2]
    )

    nextDecon <- decontXInitialize(
      counts = counts,
      theta = theta,
      z = z,
      pseudocount = 1e-20
    )
    phi <- nextDecon$phi
    eta <- nextDecon$eta

    # if countsBackground is not null, use empirical dist. to replace eta
    if (!is.null(countsBackground)) {
      # Add pseudocount to each gene in eta
      eta_tilda <- Matrix::rowSums(countsBackground) + 1e-20
      eta <- eta_tilda / sum(eta_tilda)
      # Make eta a matrix same dimension as phi
      eta <- matrix(eta, length(eta), dim(phi)[2])
    }

    ll <- c()

    ## EM updates
    theta.previous <- theta
    converged <- FALSE
    counts.colsums <- Matrix::colSums(counts)
    while (iter <= maxIter && !isTRUE(converged) &&
             numIterWithoutImprovement <= stopIter) {
      if (is.null(countsBackground)) {
        nextDecon <- decontXEM(
          counts = counts,
          counts_colsums = counts.colsums,
          phi = phi,
          estimate_eta = TRUE,
          eta = eta,
          theta = theta,
          z = z,
          estimate_delta = isTRUE(estimateDelta),
          delta = delta,
          pseudocount = 1e-20
        )
      } else {
        nextDecon <- decontXEM(
          counts = counts,
          counts_colsums = counts.colsums,
          phi = phi,
          estimate_eta = FALSE,
          eta = eta,
          theta = theta,
          z = z,
          estimate_delta = isTRUE(estimateDelta),
          delta = delta,
          pseudocount = 1e-20
        )
      }


      theta <- nextDecon$theta
      phi <- nextDecon$phi
      eta <- nextDecon$eta
      delta <- nextDecon$delta

      max.divergence <- max(abs(theta.previous - theta))
      if (max.divergence < convergence) {
        converged <- TRUE
      }
      theta.previous <- theta

      ## Calculate likelihood and check for convergence
      if (iter %% iterLogLik == 0 || converged) {
        llTemp <- decontXLogLik(
          counts = counts,
          z = z,
          phi = phi,
          eta = eta,
          theta = theta,
          pseudocount = 1e-20
        )

        ll <- c(ll, llTemp)

        .logMessages(date(),
          "...... Completed iteration:",
          iter,
          "| converge:",
          signif(max.divergence, 4),
          logfile = logfile,
          append = TRUE,
          verbose = verbose
        )
      }

      iter <- iter + 1L
    }
  }

  resConp <- nextDecon$contamination
  names(resConp) <- colnames(counts)

  list(
    "logLikelihood" = ll,
    "contamination" = resConp,
    "theta" = theta,
    "delta" = delta,
    "phi" = phi,
    "eta" = eta,
    "UMAP" = umap,
    "iteration" = iter - 1L,
    "z" = z
  )
}


## Make sure provided count matrix is the right type
.checkCountsDecon <- function(counts) {
  if (sum(is.na(counts)) > 0) {
    stop("Missing value in 'counts' matrix.")
  }
  if (is.null(dim(counts))) {
    stop("At least 2 genes need to have non-zero expressions.")
  }
}


## Make sure provided cell labels are the right type
#' @importFrom plyr mapvalues
.processCellLabels <- function(z, numCells) {
  if (length(z) != numCells) {
    stop(
      "'z' must be of the same length as the number of cells in the",
      " 'counts' matrix."
    )
  }
  if (length(unique(z)) < 2) {
    stop(
      "No need to decontaminate when only one cluster",
      " is in the dataset."
    ) # Even though
    # everything runs smoothly when length(unique(z)) == 1, result is not
    # trustful
  }
  if (!is.factor(z)) {
    z <- plyr::mapvalues(z, unique(z), seq_along(unique(z)))
    z <- as.factor(z)
  }
  z
}


.decontxInitializeZ <- function(counts,
                                varGenes = 2000,
                                dbscanEps = 1,
                                estimateCellTypes = TRUE,
                                seed = 12345,
                                legacyInit = FALSE) {
  if (isTRUE(legacyInit)) {
    init <- .decontxInitializeZLegacy(counts,
      varGenes = varGenes,
      seed = seed
    )
  } else {
    init <- .decontxInitializeZScrapper(counts,
      varGenes = varGenes,
      seed = seed
    )
  }
  normed <- init$normed
  resUmap <- init$umap

  z <- NULL
  if (isTRUE(estimateCellTypes)) {
    # Find clusters with dbSCAN
    totalClusters <- 1
    iter <- 1
    while (totalClusters <= 1 && dbscanEps > 0 && iter < 10) {
      resDbscan <- dbscan::dbscan(resUmap, dbscanEps)
      dbscanEps <- dbscanEps - (0.25 * dbscanEps)
      totalClusters <- length(unique(resDbscan$cluster))
      iter <- iter + 1
    }

    # If dbscan was not able to get more than 2 clusters,
    # use kmeans to force 2 clusters as a last resort
    if (totalClusters == 1) {
      ## transpose before realizing so only one dense copy is allocated
      cl <- stats::kmeans(as.matrix(t(normed)), 2)
      z <- cl$cluster
    } else {
      z <- resDbscan$cluster
    }
  }

  list(
    "z" = z,
    "umap" = resUmap
  )
}

## Default initialization: log-normalization, feature selection, PCA, and
## UMAP via the 'scrapper' package. scrapper does not use R's RNG; the
## UMAP is seeded explicitly, so results are deterministic for a given
## seed (and with scrapper's own default seeds when seed = NULL).
.decontxInitializeZScrapper <- function(counts, varGenes, seed) {
  libSizes <- Matrix::colSums(counts)
  if (any(libSizes == 0)) {
    stop(
      "All cells must have at least one count to estimate cell ",
      "clusters. Remove empty cells or droplets before running ",
      "decontX, or supply cluster labels with the 'z' parameter."
    )
  }
  sf <- scrapper::centerSizeFactors(libSizes)
  normed <- scrapper::normalizeCounts(counts, size.factors = sf)

  ## scrapper::modelGeneVariances errors on duplicated gene names (common
  ## in 10x data labeled with gene symbols); names are not needed here
  ## because variable genes are selected by index
  rownames(normed) <- NULL

  ## mean.filter = FALSE ranks all genes like the previous scater-based
  ## selection did; the default abundance filter (min.mean = 0.1) can
  ## remove every gene in very sparse datasets and abort
  geneVar <- scrapper::modelGeneVariances(normed,
    mean.filter = FALSE,
    num.threads = 1
  )
  hvg <- scrapper::chooseHighlyVariableGenes(geneVar$statistics$residuals,
    top = varGenes
  )

  ## 50 PCs matches the internal default of scater::calculateUMAP,
  ## which this pipeline replaces
  pca <- scrapper::runPca(normed[hvg, , drop = FALSE],
    number = 50,
    num.threads = 1
  )

  ## scrapper >= 1.5 split runUmap's 'seed' into 'initialize.seed' and
  ## 'optimize.seed'; support both APIs
  umapArgs <- list(pca$components, num.threads = 1)
  if (!is.null(seed)) {
    if ("seed" %in% names(formals(scrapper::runUmap))) {
      umapArgs$seed <- seed
    } else {
      umapArgs$initialize.seed <- seed
      umapArgs$optimize.seed <- seed
    }
  }
  resUmap <- do.call(scrapper::runUmap, umapArgs)

  list(normed = normed, umap = resUmap)
}

## Original scater/scuttle-based initialization, kept for backwards
## compatibility so results from decontX 1.11.0 and earlier can be
## reproduced with legacyInit = TRUE. logNormCounts/normalizeCounts are
## deprecated upstream in favor of 'scrapper'; this path keeps working
## for as long as scater/scuttle continue to export them.
.decontxInitializeZLegacy <- function(counts, varGenes, seed) {
  if (!requireNamespace("scater", quietly = TRUE)) {
    stop(
      "'legacyInit = TRUE' requires the 'scater' package. Install it ",
      "with BiocManager::install(\"scater\") or use the default ",
      "initialization (legacyInit = FALSE)."
    )
  }
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts)
  )
  sce <- scater::logNormCounts(sce, log = TRUE)

  if (!is.null(seed)) {
    with_seed(
      seed,
      resUmap <- scater::calculateUMAP(sce,
        ntop = varGenes,
        n_threads = 1,
        exprs_values = "logcounts"
      )
    )
  } else {
    resUmap <- scater::calculateUMAP(sce,
      ntop = varGenes,
      n_threads = 1,
      exprs_values = "logcounts"
    )
  }

  list(normed = SingleCellExperiment::logcounts(sce), umap = resUmap)
}


## process varGenes
.processvarGenes <- function(varGenes) {
  if (is.null(varGenes)) {
    varGenes <- 5000
  } else {
    if (varGenes < 2 || length(varGenes) > 1) {
      stop("Parameter 'varGenes' must be an integer larger than 1.")
    }
  }
  varGenes
}

## process dbscanEps for resolusion threshold using DBSCAN
.processdbscanEps <- function(dbscanEps) {
  if (is.null(dbscanEps)) {
    dbscanEps <- 1
  } else {
    if (dbscanEps < 0) {
      stop("Parameter 'dbscanEps' needs to be non-negative.")
    }
  }
  dbscanEps
}

.checkDelta <- function(delta) {
  if (!is.numeric(delta) || length(delta) != 2 || any(delta < 0)) {
    stop(
      "'delta' needs to be a numeric vector of length 2",
      " containing positive values."
    )
  }
  delta
}


#########################
# Simulating Data       #
#########################

#' @title Simulate contaminated count matrix
#' @description This function generates a list containing two count matrices --
#'  one for real expression, the other one for contamination, as well as other
#'  parameters used in the simulation which can be useful for running
#'  decontamination.
#' @param C Integer. Number of cells to be simulated. Default \code{300}.
#' @param G Integer. Number of genes to be simulated. Default \code{100}.
#' @param K Integer. Number of cell populations to be simulated.
#' Default \code{3}.
#' @param NRange Integer vector. A vector of length 2 that specifies the lower
#'  and upper bounds of the number of counts generated for each cell. Default
#'  \code{c(500, 1000)}.
#' @param beta Numeric. Concentration parameter for Phi. Default \code{0.1}.
#' @param delta Numeric or Numeric vector. Concentration parameter for Theta.
#'  If input as a single numeric value, symmetric values for beta
#'  distribution are specified; if input as a vector of lenght 2, the two
#'  values will be the shape1 and shape2 paramters of the beta distribution
#'  respectively. Default \code{c(1, 5)}.
#' @param numMarkers Integer. Number of markers for each cell population.
#' Default \code{3}.
#' @param seed Integer. Passed to \code{\link[withr]{with_seed}}.
#' For reproducibility, a default value of 12345 is used. If NULL, no calls to
#'  \code{\link[withr]{with_seed}} are made.
#' @return A list containing the \code{nativeMatirx} (real expression),
#' \code{observedMatrix} (real expression + contamination), as well as other
#' parameters used in the simulation.
#' @author Shiyi Yang, Yuan Yin, Joshua Campbell
#' @examples
#' contaminationSim <- simulateContamination(K = 3, delta = c(1, 10))
#' @importFrom withr with_seed
#' @export
simulateContamination <- function(C = 300,
                                  G = 100,
                                  K = 3,
                                  NRange = c(500, 1000),
                                  beta = 0.1,
                                  delta = c(1, 10),
                                  numMarkers = 3,
                                  seed = 12345) {
  if (is.null(seed)) {
    res <- .simulateContaminatedMatrix(
      C = C,
      G = G,
      K = K,
      NRange = NRange,
      beta = beta,
      delta = delta,
      numMarkers = numMarkers
    )
  } else {
    with_seed(
      seed,
      res <- .simulateContaminatedMatrix(
        C = C,
        G = G,
        K = K,
        NRange = NRange,
        beta = beta,
        delta = delta,
        numMarkers = numMarkers
      )
    )
  }

  res
}


.simulateContaminatedMatrix <- function(C = 300,
                                        G = 100,
                                        K = 3,
                                        NRange = c(500, 1000),
                                        beta = 0.5,
                                        delta = c(1, 2),
                                        numMarkers = 3) {
  if (length(delta) == 1) {
    cpByC <- stats::rbeta(
      n = C,
      shape1 = delta,
      shape2 = delta
    )
  } else {
    cpByC <- stats::rbeta(
      n = C,
      shape1 = delta[1],
      shape2 = delta[2]
    )
  }

  z <- sample(seq(K), size = C, replace = TRUE)
  if (length(unique(z)) < K) {
    warning(
      "Only ",
      length(unique(z)),
      " clusters are simulated. Try to increase numebr of cells 'C' if",
      " more clusters are needed"
    )
    K <- length(unique(z))
    z <- plyr::mapvalues(z, unique(z), seq_along(unique(z)))
  }

  NbyC <- sample(seq(min(NRange), max(NRange)),
    size = C,
    replace = TRUE
  )
  cNbyC <- vapply(seq(C), function(i) {
    stats::rbinom(
      n = 1,
      size = NbyC[i],
      p = cpByC[i]
    )
  }, integer(1))
  rNbyC <- NbyC - cNbyC

  phi <- .rdirichlet(K, rep(beta, G))

  ## Select random genes to be markers in each cell population
  ## by setting their values to zero.
  if (K * numMarkers > G) {
    stop(
      "The number of markers ('numMarkers') times the number of cell",
      " populations ('K') cannot be greater than the number of",
      " genes ('G')."
    )
  }
  markerKIndex <- rep(seq(K), each = numMarkers)
  markerRowIndex <- sample(seq(G), numMarkers * K)
  for (i in seq(K)) {
    ix <- markerRowIndex[markerKIndex == i]
    phi[i, ix] <- max(phi[i, ])
    for (j in setdiff(seq(K), i)) {
      phi[j, ix] <- 0
    }
  }
  phi <- prop.table(phi, margin = 1)

  ## sample real expressed count matrix
  cellRmat <- vapply(seq(C), function(i) {
    stats::rmultinom(1, size = rNbyC[i], prob = phi[z[i], ])
  }, integer(G))

  rownames(cellRmat) <- paste0("Gene_", seq(G))
  colnames(cellRmat) <- paste0("Cell_", seq(C))

  ## Get list of marker names
  markerNames <- list()
  for (i in seq(K)) {
    markerNames[[i]] <- rownames(cellRmat)[markerRowIndex[markerKIndex == i]]
  }
  names(markerNames) <- paste0("CellType_", seq(K), "_Markers")

  ## sample contamination count matrix
  nGByK <-
    rowSums(cellRmat) - .colSumByGroup(cellRmat, group = z, K = K)
  ## Column-wise proportion normalization. This reproduces
  ## celda::normalizeCounts(normalize = "proportion") exactly (verified to
  ## the bit on the simulation path) via the package's own fastNormProp
  ## C++ routine with a zero pseudocount, removing the celda dependency
  ## (ADR-0003 Stage 1).
  eta <- fastNormProp(nGByK, 0)

  cellCmat <- vapply(seq(C), function(i) {
    stats::rmultinom(1, size = cNbyC[i], prob = eta[, z[i]])
  }, integer(G))
  cellOmat <- cellRmat + cellCmat
  contamination <- colSums(cellCmat) / colSums(cellOmat)

  rownames(cellOmat) <- paste0("Gene_", seq(G))
  colnames(cellOmat) <- paste0("Cell_", seq(C))

  list(
    "nativeCounts" = cellRmat,
    "observedCounts" = cellOmat,
    "NByC" = NbyC,
    "z" = z,
    "eta" = eta,
    "phi" = t(phi),
    "markers" = markerNames,
    "numMarkers" = numMarkers,
    "contamination" = contamination
  )
}


.checkBackground <- function(x, background, bgBatch,
                             logfile = NULL, verbose = FALSE) {
  # Remove background barcodes that have already appeared in x
  # If bgBatch param is supplied, also remove duplicate bgBatch
  if (!is.null(colnames(background))) {
    dupBarcode <- colnames(background) %in% colnames(x)
  } else {
    dupBarcode <- FALSE
    warning(
      "No column names were found for the 'background' matrix. ",
      "No checking was performed between the ids in the 'backgroud' ",
      "matrix and 'x'.",
      " Please ensure that no true cells are included in the background ",
      "matrix. Otherwise, results will be incorrect."
    )
  }

  if (any(dupBarcode)) {
    .logMessages(
      date(),
      ".. ",
      sum(dupBarcode),
      " cells in the background matrix were removed as they were found in",
      " the filtered matrix.",
      logfile = logfile,
      append = TRUE,
      verbose = verbose
    )
    background <- background[, !(dupBarcode), drop = FALSE]

    if (!is.null(bgBatch)) {
      if (length(bgBatch) != length(dupBarcode)) {
        stop(
          "Length of bgBatch must be equal to the number of columns",
          "of background matrix."
        )
      }
      bgBatch <- bgBatch[!(dupBarcode)]
    }
  }

  re <- list(
    background = background,
    bgBatch = bgBatch
  )

  re
}
