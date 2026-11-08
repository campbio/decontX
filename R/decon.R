#' @title Contamination estimation with decontX
#'
#' @description Identifies contamination from factors such as ambient RNA
#' in single cell genomic datasets.
#'
#' @name decontX
#'
#' @param x A matrix of counts or a
#' \linkS4class[SingleCellExperiment]{SingleCellExperiment} with the matrix
#' located in the assay slot under \code{assayName}. The matrix can be a
#' numeric matrix, a sparse matrix from package \link[Matrix]{Matrix} or
#' 'SparseArray', or a \linkS4class[DelayedArray]{DelayedMatrix}.
#' Cells in each batch will be subsetted and converted to a sparse matrix
#' of class \code{dgCMatrix} from package \link[Matrix]{Matrix} before
#' analysis. Each batch needs at least 3 cells (15 with
#' \code{legacyInit = TRUE}) for the embedding, or 2 with
#' \code{embedding = FALSE}. This
#' object should only contain filtered cells after cell calling. Empty
#' cell barcodes (low expression droplets before cell calling) are not needed
#' to run DecontX.
#' @param assayName Character. Name of the assay to use if \code{x} is a
#' \linkS4class[SingleCellExperiment]{SingleCellExperiment}.
#' @param z Numeric, character or factor vector. Cell cluster labels,
#' which the result's \code{z} returns as a factor of the same labels.
#' If NULL,
#' PCA will be used to reduce the dimensionality of the dataset initially,
#' '\link[scrapper]{runUmap}' from the 'scrapper' package
#' will be used to further reduce the dataset to 2 dimensions and
#' the '\link[dbscan]{dbscan}' function from the 'dbscan' package
#' will be used to identify clusters of broad cell types. If dbscan finds
#' only one cluster, k-means on the PCA components (on all genes with
#' \code{legacyInit = TRUE}) splits the cells into two. \code{z} needs at
#' least two clusters in each batch, unless \code{background} is
#' supplied: with \code{background}, one cluster is allowed, which fits
#' only when the cells are a single cell type.
#' Default NULL.
#' @param batch Numeric or character vector. Batch labels for cells.
#' If batch labels are supplied, DecontX is run on cells from each
#' batch separately. Cells run in different channels or assays
#' should be considered different batches. Default NULL.
#' @param background A matrix of counts or a
#' \linkS4class[SingleCellExperiment]{SingleCellExperiment} with the matrix
#' located in the assay slot under \code{assayName}. It must be a
#' \linkS4class[SingleCellExperiment]{SingleCellExperiment} when \code{x} is
#' one, and otherwise one of the matrix classes \code{x} accepts. It
#' contains the empty droplets instead of cells. When supplied,
#' empirical distribution of transcripts from these empty droplets
#' will be used as the contamination distribution. When both have row
#' names, its rows are matched to the genes of \code{x} by name: genes
#' that aren't in \code{x} are dropped, and every gene of \code{x} must be
#' present. Otherwise the rows are matched by position. Barcodes that are
#' also in \code{x} are removed. With \code{background}, \code{z} may have
#' a single cluster. Default NULL.
#' @param bgAssayName Character. Name of the assay to use if \code{background}
#' is a \linkS4class[SingleCellExperiment]{SingleCellExperiment}. Default to
#' same as
#' \code{assayName}.
#' @param bgBatch Numeric or character vector. Batch labels for
#' \code{background}, one per column. Its unique values must be the same as
#' those in \code{batch}, such that each batch of cells have their
#' corresponding batch of empty droplets as background, pointed by this
#' parameter: every batch needs at least one droplet, after barcodes that are
#' also in \code{x} are removed. It needs \code{batch} and \code{background}.
#' Default to NULL.
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
#' \code{legacyInit = TRUE}). If \code{z} is provided, it is used only for
#' the UMAP stored in the results, and not at all with
#' \code{embedding = FALSE}. Default 5000.
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
#'  scater/scuttle. If \code{z} is provided, it is used only for the UMAP
#'  stored in the results, and not at all with \code{embedding = FALSE}.
#'  Default FALSE.
#' @param logfile Character. Messages will be redirected to a file named
#'  `logfile`. If NULL, messages will be printed to stdout.  Default NULL.
#' @param verbose Logical. Whether to print log messages. Default TRUE.
#' @param embedding Logical. Whether to compute the low-dimensional
#'  embedding of the cells (currently a UMAP) that decontX uses to estimate
#'  the cell clusters and stores in the results. When \code{z} is provided,
#'  the embedding is used only for the stored results, so
#'  \code{embedding = FALSE} skips it, which saves most of the run time on
#'  large data sets; the decontamination is the same either way (except
#'  with \code{legacyInit = TRUE} and \code{seed = NULL}, where the legacy
#'  UMAP draws from the session's random numbers before the EM). Then no
#'  UMAP is stored, and \code{\link{plotDecontXContamination}} can plot the
#'  contamination on another embedding with \code{reducedDimName}.
#'  \code{embedding = FALSE} needs \code{z}. Default TRUE.
#' @param numThreads Integer. Number of threads for computing the
#'  embedding (variable genes, PCA and UMAP). With more than one thread the
#'  PCA can differ in the last digits on large data sets, which can move
#'  the UMAP coordinates and change estimated clusters, so use the same
#'  \code{numThreads} to reproduce a run.
#'  With \code{legacyInit = TRUE}, only the UMAP uses the threads; it was
#'  the same with one and four threads on 8,000 simulated cells. Has no
#'  effect with \code{embedding = FALSE}. Default 1.
#' @param ... For the generic, further arguments to pass to each method.
#'
#' @return If \code{x} is a matrix-like object, a list will be returned
#' with the following items:
#' \describe{
#' \item{\code{decontXcounts}:}{The decontaminated matrix, a sparse
#' \code{dgCMatrix} whatever the class of \code{x}, or a
#' \linkS4class[DelayedArray]{DelayedMatrix} when \code{x} is one. The
#' DelayedMatrix holds the sparse result in memory, also when \code{x} is
#' on disk (for example an \code{HDF5Matrix}); to keep it on disk, save it
#' explicitly, for example with \code{HDF5Array::writeHDF5Array()}. Values
#' obtained from the variational inference procedure may be non-integer.
#' However, integer counts can be obtained by rounding,
#' e.g. \code{round(decontXcounts)}.}
#' \item{\code{contamination}:}{Percentage of contamination in each cell.}
#' \item{\code{estimates}:}{List of estimated parameters for each batch,
#' including the UMAP coordinates (\code{UMAP}), which are used to generate
#' the cell cluster labels when z is not supplied. \code{UMAP} is
#' \code{NULL} with \code{embedding = FALSE}. The columns of \code{phi} and
#' \code{eta} are named by cluster label, in the order of the levels of the
#' batch's \code{z}.}
#' \item{\code{z}:}{Factor of the cell cluster labels used for analysis. A
#' supplied \code{z} comes back with its own labels: a factor keeps its
#' level order, without unused levels, and other vectors get sorted levels.
#' Estimated clusters are numbered 1 to K in each batch. When \code{batch}
#' is supplied, the labels are prefixed with the batch
#' (\code{"<batch>-<label>"}) for estimated clusters, and for a supplied
#' \code{z} that has a label in more than one batch, since the same label
#' may then name different clusters. Prefixed labels are ordered by the
#' batch's first cell, then by label. The \code{z} in each batch's
#' \code{estimates} has no prefix.}
#' \item{\code{runParams}:}{List of arguments used in the function call.
#' Its \code{batch} and \code{batchBackground} are the batch of each cell
#' and of each background droplet used for the run. It always has
#' \code{batchBackground}, which is \code{NULL} when \code{background} and
#' \code{bgBatch} are both NULL.}
#' \item{\code{batch}:}{Character vector of the batch of each cell used for
#' the run, \code{"all_cells"} for every cell when \code{batch} is NULL.}
#' }
#'
#' If \code{x} is a \linkS4class[SingleCellExperiment]{SingleCellExperiment},
#' then the decontaminated
#' counts will be stored as an assay and can be accessed with
#' \code{decontXcounts(x)}. The contamination values and cluster labels
#' will be stored in \code{colData(x)}. \code{estimates} and \code{runParams}
#' will be stored in \code{metadata(x)$decontX}. The UMAPs used to generated
#' cell cluster labels will be stored in
#' \code{reducedDims} slot in \code{x}, except with
#' \code{embedding = FALSE}. Running \code{decontX} again on
#' \code{x} replaces the assay, the \code{colData} columns and
#' \code{metadata(x)$decontX}. A UMAP from an earlier run
#' (\code{decontX_UMAP} or \code{decontX_<batch>_UMAP}) stays in
#' \code{reducedDims} when the new run doesn't replace it, because the
#' batches differ or the new run has \code{embedding = FALSE}, and it
#' doesn't match the new results. When the counts
#' are on disk, the \code{decontXcounts} assay is held in memory; save the
#' object with, for example,
#' \code{HDF5Array::saveHDF5SummarizedExperiment()} to keep it on disk.
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
                                                      verbose = TRUE,
                                                      embedding = TRUE,
                                                      numThreads = 1) {
  .checkAssayName(x, assayName, "assayName", "x")
  mat <- SummarizedExperiment::assay(x, i = assayName)
  .checkCountsMatrix(mat, "The 'x' assay")
  params <- .decontXParams(
    maxIter = maxIter,
    delta = delta,
    estimateDelta = estimateDelta,
    convergence = convergence,
    iterLogLik = iterLogLik,
    varGenes = varGenes,
    dbscanEps = dbscanEps,
    seed = seed,
    legacyInit = legacyInit,
    logfile = logfile,
    verbose = verbose,
    embedding = embedding,
    numThreads = numThreads
  )

  bgCounts <- NULL
  if (!is.null(background)) {
    if (!inherits(background, "SummarizedExperiment")) {
      stop(
        "'background' must be a SingleCellExperiment when 'x' is one; ",
        "it is a ", class(background)[1], "."
      )
    }
    bgArg <- "bgAssayName"
    if (is.null(bgAssayName)) {
      bgAssayName <- assayName
      bgArg <- "assayName"
    }
    .checkAssayName(background, bgAssayName, bgArg, "background")
    bgCounts <- SummarizedExperiment::assay(background, i = bgAssayName)
    .checkCountsMatrix(bgCounts, "The 'background' assay")
  }

  result <- .decontX(
    counts = mat,
    z = z,
    batch = batch,
    background = bgCounts,
    bgBatch = bgBatch,
    params = params
  )

  ## Add results into column annotation
  SummarizedExperiment::colData(x)$decontX_contamination <- result$contamination
  SummarizedExperiment::colData(x)$decontX_clusters <- result$z

  ## Put estimated UMAPs into SCE
  umaps <- .decontXUMAPs(result, colnames(mat))
  for (name in names(umaps)) {
    SingleCellExperiment::reducedDim(x, name) <- umaps[[name]]
  }

  ## Save the rest of the result object into metadata
  decontXcounts(x) <- result$decontXcounts
  result$decontXcounts <- NULL
  result$batch <- NULL
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
                                     verbose = TRUE,
                                     embedding = TRUE,
                                     numThreads = 1) {
  hint <- ""
  if (inherits(x, "SummarizedExperiment")) {
    hint <- " Convert it with as(x, \"SingleCellExperiment\")."
  }
  .checkCountsMatrix(x, "'x'",
    choices = "a SingleCellExperiment, ",
    hint = hint
  )
  params <- .decontXParams(
    maxIter = maxIter,
    delta = delta,
    estimateDelta = estimateDelta,
    convergence = convergence,
    iterLogLik = iterLogLik,
    varGenes = varGenes,
    dbscanEps = dbscanEps,
    seed = seed,
    legacyInit = legacyInit,
    logfile = logfile,
    verbose = verbose,
    embedding = embedding,
    numThreads = numThreads
  )

  # Objects without dimensions get .checkBackground()'s message
  if (!is.null(dim(background))) {
    .checkCountsMatrix(background, "'background'",
      when = " when 'x' is a matrix"
    )
  }

  .decontX(
    counts = x,
    z = z,
    batch = batch,
    background = background,
    bgBatch = bgBatch,
    params = params
  )
})


##########################
# Core Decontx Functions #
##########################

## The pseudocount the EM and the background profile add to avoid zeros
.decontXPseudocount <- 1e-20

## The container-agnostic core: counts, labels, background and the
## parameters from .decontXParams() in, a plain list of results out. Every
## method calls it, passing every argument
.decontX <- function(counts, z, batch, background, bgBatch, params) {
  .checkCountsDecon(counts)
  .checkCellValues(z, "z", ncol(counts))
  .checkCellValues(batch, "batch", ncol(counts))
  .checkEmbeddingZ(z, params$embedding)

  if (!is.null(background)) {
    # Drop background barcodes that are also cells, with their bgBatch
    temp <- .checkBackground(counts, background, bgBatch,
      logfile = params$logfile, verbose = params$verbose
    )
    background <- temp$background
    bgBatch <- temp$bgBatch
    .checkCountsDecon(background, "background")
  }

  startTime <- Sys.time()
  .logBanner("Starting DecontX", params = params)
  if (!is.null(background)) {
    background <- .matchBackgroundRows(counts, background,
      logfile = params$logfile, verbose = params$verbose
    )
  }

  batches <- .decontXBatches(batch, bgBatch, background, ncol(counts))
  # Checked for every batch before any of them runs
  .checkBatchSizes(batches$batch, batches$given, params)
  .checkBatchClusters(z, batches$batch, batches$given, background)

  res <- .runBatches(counts, z, background, batches, params)
  returnResult <- list(
    ## The batches used for the run, then the parameters as given
    "runParams" = c(
      list(z = z, batch = batches$batch, batchBackground = batches$bgBatch),
      params
    ),
    "estimates" = res$estimates,
    "decontXcounts" = .asInputClass(res$decontXcounts, counts, params),
    "contamination" = res$contamination,
    "z" = res$z,
    "batch" = batches$batch
  )

  .logBanner("Completed DecontX. Total time:",
    format(difftime(Sys.time(), startTime)),
    params = params
  )

  returnResult
}

## The batch of each cell and of each background droplet used for the run,
## and whether 'batch' was given. Without 'batch', every cell and droplet
## is in batch "all_cells"
.decontXBatches <- function(batch, bgBatch, background, numCells) {
  given <- !is.null(batch)
  if (is.null(batch)) {
    batch <- rep("all_cells", numCells)

    # If batch null, bgBatch has to be null
    if (!is.null(bgBatch)) {
      stop("'bgBatch' was given without 'batch'; give both or neither.")
    }

    if (!is.null(background)) {
      bgBatch <- rep("all_cells", ncol(background))
    }
  } else if (is.null(background) && !is.null(bgBatch)) {
    stop("'bgBatch' was given without 'background'; give both or neither.")
  } else if (!is.null(background) && is.null(bgBatch)) {
    # If batch not null and background supplied,
    # user has to supply bgBatch as well
    stop(
      "Cell batch, and background are supplied. Please also ",
      "supply background batch."
    )
  }
  batch <- as.character(batch)
  if (!is.null(bgBatch)) {
    bgBatch <- as.character(bgBatch)
    notInBatch <- setdiff(bgBatch, batch)
    if (length(notInBatch) > 0) {
      stop(
        "'bgBatch' has values that aren't in 'batch': ",
        .firstFew(notInBatch), "."
      )
    }
    # Barcodes that are also in x were already removed from the background
    noBackground <- setdiff(batch, bgBatch)
    if (length(noBackground) > 0) {
      stop(
        "Every batch needs droplets in 'background'; none for: ",
        .firstFew(noBackground), "."
      )
    }
  }
  list(batch = batch, bgBatch = bgBatch, given = given)
}

## Without background, the labels of each batch need two clusters or more
.checkBatchClusters <- function(z, batch, given, background) {
  if (is.null(z) || !is.null(background)) {
    return(invisible(NULL))
  }
  for (bat in unique(batch)) {
    .checkOneCluster(z[batch == bat], if (given) bat)
  }
  invisible(NULL)
}

## User labels come back as given, without unused levels. Batched results
## get the batch prefix when a label could name different clusters:
## estimated clusters, which are numbered within each batch, and user
## labels found in more than one batch
.decontXLabels <- function(z, batch) {
  labels <- NULL
  if (!is.null(z)) {
    labels <- if (is.factor(z)) droplevels(z) else factor(z)
  }
  prefix <- length(unique(batch)) > 1
  if (prefix && !is.null(z)) {
    batchesPerLabel <- tapply(batch, labels, function(b) length(unique(b)))
    prefix <- any(batchesPerLabel > 1)
  }
  list(labels = labels, prefix = prefix)
}

## Run decontX on each batch in turn and gather the results of every cell
.runBatches <- function(counts, z, background, batches, params) {
  batch <- batches$batch
  labels <- .decontXLabels(z, batch)

  ## Set result lists upfront for all cells from different batches
  batchCounts <- list()
  estConp <- rep(NA, ncol(counts))
  returnZ <- character(ncol(counts))
  zLevels <- character()
  resBatch <- list()

  for (bat in unique(batch)) {
    cells <- batch == bat
    run <- .runBatch(counts, z, background, batches$bgBatch, batch, bat,
                     params)
    batchCounts <- c(batchCounts, list(run$decontXcounts))
    batchLabels <- if (is.null(z)) run$res$z else labels$labels[cells]
    resBatch[[bat]] <- .batchEstimates(run$res, batchLabels)

    estConp[cells] <- run$res$contamination
    prefix <- if (labels$prefix) paste0(bat, "-") else ""
    returnZ[cells] <- paste0(prefix, resBatch[[bat]]$z)
    zLevels <- c(zLevels, paste0(prefix, levels(resBatch[[bat]]$z)))
  }
  names(resBatch) <- unique(batch)

  # Unprefixed user labels keep the level order of the whole 'z'
  if (is.null(z) || labels$prefix) {
    returnZ <- factor(returnZ, levels = zLevels)
  } else {
    returnZ <- labels$labels
  }
  decontXcounts <- .bindBatches(batchCounts)
  # Free the per-batch counts before the reorder copies the result
  rm(batchCounts)
  # The dimnames without the names of the list, as decontX always gave them
  dimNames <- list(rownames(counts), colnames(counts))
  list(
    estimates = resBatch,
    decontXcounts = .inCellOrder(decontXcounts, batch, dimNames),
    contamination = estConp,
    z = returnZ
  )
}

## Run decontX on the cells and background droplets of one batch, and
## compute their decontaminated counts
.runBatch <- function(counts, z, background, bgBatch, batch, bat, params) {
  if (length(unique(batch)) == 1) {
    .logDecontX(date(), ".. Analyzing all cells", params = params)
  } else {
    .logDecontX(date(), " .. Analyzing cells in batch '", bat, "'",
      sep = "", params = params
    )
  }

  countsBat <- counts[, batch == bat, drop = FALSE]
  bgBat <- background[, bgBatch == bat, drop = FALSE]
  zBat <- if (is.null(z)) NULL else z[batch == bat]

  ## Convert to sparse matrix
  if (!inherits(countsBat, "dgCMatrix")) {
    .logDecontX(date(), ".... Converting to sparse matrix", params = params)
    countsBat <- .asDgCMatrix(countsBat)
  }
  if (!is.null(bgBat)) {
    bgBat <- .asDgCMatrix(bgBat)
  }

  res <- .withSeed(
    params$seed,
    .decontXoneBatch(
      counts = countsBat,
      z = zBat,
      background = bgBat,
      params = params
    )
  )

  .logDecontX(date(), ".. Calculating final decontaminated matrix",
    params = params
  )
  decontXcounts <- .nativeCounts(countsBat, res)
  list(res = res, decontXcounts = decontXcounts)
}

## The decontaminated counts of one batch: each count times its native
## probability. Only the values change, so the result is a plain dgCMatrix
## copy of the counts with new values (a subclass is dropped). Slot
## assignment copies shallowly, so the row indices and column pointers are
## shared, not copied. Unnamed, so the batches bind without names, as they
## always have; .inCellOrder() names the bound result
.nativeCounts <- function(counts, res) {
  native <- methods::as(counts, "dgCMatrix")
  methods::slot(native, "x") <- calculateNativeCounts(
    counts = native,
    theta = res$theta,
    eta = res$eta,
    phi = res$phi,
    z = as.integer(res$z),
    pseudocount = .decontXPseudocount
  )
  dimnames(native) <- NULL
  native
}

## The decontaminated counts of each batch, in a list in batch order, bound
## into one matrix with a single cbind. One batch's matrix is not copied
.bindBatches <- function(batchCounts) {
  if (length(batchCounts) == 1) {
    return(batchCounts[[1]])
  }
  do.call(cbind, batchCounts)
}

## Put the columns of .bindBatches() back in the cells' order and name them.
## Its columns hold the cells batch by batch, in order of first appearance
## of each batch, so one reorder is needed only when batches interleave
.inCellOrder <- function(bound, batch, dimNames) {
  byBatch <- order(match(batch, unique(batch)))
  if (is.unsorted(byBatch)) {
    bound <- bound[, order(byBatch), drop = FALSE]
  }
  dimnames(bound) <- dimNames
  bound
}

## The estimates of one batch, with z, phi and eta labelled
.batchEstimates <- function(res, labels) {
  labelled <- .labelEstimates(res, labels)
  list(
    z = labelled$z,
    phi = labelled$phi,
    eta = labelled$eta,
    delta = res$delta,
    theta = res$theta,
    contamination = res$contamination,
    logLikelihood = res$logLikelihood,
    UMAP = res$UMAP,
    iteration = res$iteration
  )
}

## DelayedMatrix counts give a DelayedMatrix result over the sparse result,
## held in memory even when the counts are on disk (for example in HDF5)
.asInputClass <- function(decontXcounts, counts, params) {
  if (!inherits(counts, c("DelayedMatrix", "DelayedArray"))) {
    return(decontXcounts)
  }
  .logDecontX(date(),
    ".. Converting decontaminated matrix to a DelayedMatrix",
    params = params
  )
  DelayedArray::DelayedArray(decontXcounts)
}

## Evaluate code with the random number generator seeded by
## withr::with_seed(), or in the current stream when seed is NULL
.withSeed <- function(seed, code) {
  if (is.null(seed)) {
    return(code)
  }
  withr::with_seed(seed, code)
}

## Log a message to params$logfile, or to the console without one
.logDecontX <- function(..., params, sep = " ") {
  .logMessages(...,
    sep = sep,
    logfile = params$logfile,
    append = TRUE,
    verbose = params$verbose
  )
}

## Log a message between two lines of dashes
.logBanner <- function(..., params) {
  rule <- paste(rep("-", 50), collapse = "")
  .logDecontX(rule, params = params)
  .logDecontX(..., params = params)
  .logDecontX(rule, params = params)
}


# This function updates decontamination for one batch
# params$seed is passed on to .decontxInitializeZ()
# .decontX() checks the counts, delta and the clusters of each batch
# before calling it
#
# Do not remove the @importFrom below: fit_dirichlet is fetched from the
# MCMCprecision namespace by the C++ EM step (src/DecontX.cpp, via
# Environment::namespace_env), which R's static analysis cannot see.
# Dropping the import triggers an R CMD check NOTE and can leave
# MCMCprecision unloaded when decontXEM() runs.
#' @importFrom MCMCprecision fit_dirichlet
#' @noRd
.decontXoneBatch <- function(counts,
                             z,
                             background,
                             params) {
  init <- .initializeBatch(counts, z, params)

  .logDecontX(date(), ".... Estimating contamination", params = params)
  start <- .initializeEM(counts, init$z, background, params$delta)
  em <- .runEM(counts, init$z, background, start, params)

  list(
    "logLikelihood" = em$logLikelihood,
    "contamination" = em$contamination,
    "theta" = em$theta,
    "delta" = em$delta,
    "phi" = em$phi,
    "eta" = em$eta,
    "UMAP" = init$umap,
    "iteration" = em$iteration,
    "z" = init$z
  )
}

## The UMAP of one batch, and its cluster labels as integer codes:
## estimated when z is NULL, otherwise the given ones. Without the
## embedding (z is then given) the UMAP is NULL
.initializeBatch <- function(counts, z, params) {
  if (!params$embedding) {
    .checkLibSizes(counts)
    return(list(z = .processCellLabels(z, numCells = ncol(counts)),
                umap = NULL))
  }
  estimateCellTypes <- is.null(z)
  if (estimateCellTypes) {
    m <- ".... Generating UMAP and estimating cell types"
  } else {
    m <- ".... Generating UMAP"
  }
  .logDecontX(date(), m, params = params)

  varGenes <- .processvarGenes(params$varGenes)
  dbscanEps <- .processdbscanEps(params$dbscanEps)

  celda.init <- .decontxInitializeZ(
    counts = counts,
    varGenes = varGenes,
    dbscanEps = dbscanEps,
    estimateCellTypes = estimateCellTypes,
    seed = params$seed,
    legacyInit = params$legacyInit,
    numThreads = params$numThreads
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

  list(z = .processCellLabels(z, numCells = ncol(counts)), umap = umap)
}

## The starting point of the EM: contamination drawn from the delta prior,
## and phi and eta from it. With background, eta is the background's gene
## profile
.initializeEM <- function(counts, z, background, delta) {
  theta <- stats::rbeta(
    n = ncol(counts),
    shape1 = delta[1],
    shape2 = delta[2]
  )

  nextDecon <- decontXInitialize(
    counts = counts,
    theta = theta,
    z = z,
    pseudocount = .decontXPseudocount
  )
  phi <- nextDecon$phi
  eta <- nextDecon$eta

  # if background is not null, use empirical dist. to replace eta
  if (!is.null(background)) {
    # Add pseudocount to each gene in eta
    eta_tilda <- Matrix::rowSums(background) + .decontXPseudocount
    eta <- eta_tilda / sum(eta_tilda)
    # Make eta a matrix same dimension as phi
    eta <- matrix(eta, length(eta), dim(phi)[2])
  }
  list(theta = theta, phi = phi, eta = eta)
}

## The EM updates, until theta converges or maxIter. With background, eta
## stays the background's profile
.runEM <- function(counts, z, background, start, params) {
  theta <- start$theta
  phi <- start$phi
  eta <- start$eta
  delta <- params$delta
  ll <- c()
  iter <- 1L
  converged <- FALSE
  counts.colsums <- Matrix::colSums(counts)
  while (iter <= params$maxIter && !isTRUE(converged)) {
    theta.previous <- theta
    nextDecon <- decontXEM(
      counts = counts,
      counts_colsums = counts.colsums,
      phi = phi,
      estimate_eta = is.null(background),
      eta = eta,
      theta = theta,
      z = z,
      estimate_delta = isTRUE(params$estimateDelta),
      delta = delta,
      pseudocount = .decontXPseudocount
    )
    theta <- nextDecon$theta
    phi <- nextDecon$phi
    eta <- nextDecon$eta
    delta <- nextDecon$delta

    max.divergence <- max(abs(theta.previous - theta))
    if (max.divergence < params$convergence) {
      converged <- TRUE
    }

    ## Calculate likelihood and check for convergence
    if (iter %% params$iterLogLik == 0 || converged) {
      ll <- c(ll, .emLogLik(counts, z, phi, eta, theta, iter, max.divergence,
                            params))
    }
    iter <- iter + 1L
  }

  contamination <- nextDecon$contamination
  names(contamination) <- colnames(counts)
  list(
    logLikelihood = ll, contamination = contamination, theta = theta,
    delta = delta, phi = phi, eta = eta, iteration = iter - 1L
  )
}

## The log-likelihood after an EM update, logged with the iteration and the
## largest change in theta
.emLogLik <- function(counts, z, phi, eta, theta, iter, divergence, params) {
  ll <- decontXLogLik(
    counts = counts,
    z = z,
    phi = phi,
    eta = eta,
    theta = theta,
    pseudocount = .decontXPseudocount
  )
  .logDecontX(date(), "...... Completed iteration:", iter,
    "| converge:", signif(divergence, 4),
    params = params
  )
  ll
}


## Check the decontX parameters and gather them in one list, in the order
## runParams records them. No defaults: the methods' signatures hold them,
## and a test checks that they agree
.decontXParams <- function(maxIter,
                           delta,
                           estimateDelta,
                           convergence,
                           iterLogLik,
                           varGenes,
                           dbscanEps,
                           seed,
                           legacyInit,
                           logfile,
                           verbose,
                           embedding,
                           numThreads) {
  .checkEMArgs(maxIter, convergence, iterLogLik)
  .checkDelta(delta)
  .processvarGenes(varGenes)
  .processdbscanEps(dbscanEps)
  .checkFlags(list(
    estimateDelta = estimateDelta,
    legacyInit = legacyInit,
    verbose = verbose,
    embedding = embedding
  ))
  .checkRunArgs(seed, logfile, numThreads)
  list(
    maxIter = maxIter,
    delta = delta,
    estimateDelta = estimateDelta,
    convergence = convergence,
    iterLogLik = iterLogLik,
    varGenes = varGenes,
    dbscanEps = dbscanEps,
    seed = seed,
    legacyInit = legacyInit,
    logfile = logfile,
    verbose = verbose,
    embedding = embedding,
    numThreads = numThreads
  )
}

## The UMAPs to store with one row per cell, named as the SCE reducedDims:
## a single batch's UMAP as it is, or each batch's padded with NA for the
## cells of the other batches. The single UMAP keeps its DecontX_UMAP_1/2
## columns and the padded ones get UMAP_1/2, as decontX always stored them.
## None when the run had no embedding
.decontXUMAPs <- function(result, cellNames) {
  batches <- names(result$estimates)
  # Without the embedding, no batch has a UMAP
  if (is.null(result$estimates[[1]]$UMAP)) {
    return(list())
  }
  if (length(batches) == 1) {
    return(list(decontX_UMAP = result$estimates[[1]]$UMAP))
  }
  umaps <- lapply(batches, function(b) {
    umap <- matrix(NA, ncol = 2, nrow = length(cellNames))
    umap[result$batch == b, ] <- result$estimates[[b]]$UMAP
    colnames(umap) <- c("UMAP_1", "UMAP_2")
    rownames(umap) <- cellNames
    umap
  })
  names(umaps) <- paste("decontX", batches, "UMAP", sep = "_")
  umaps
}

## Arguments that must be TRUE or FALSE, in a named list
.checkFlags <- function(flags) {
  for (flag in names(flags)) {
    if (!(isTRUE(flags[[flag]]) || isFALSE(flags[[flag]]))) {
      stop("'", flag, "' must be TRUE or FALSE.")
    }
  }
}

## The seed, the log file and the number of threads
.checkRunArgs <- function(seed, logfile, numThreads) {
  if (!is.null(seed) &&
        !(.isNumber(seed) && abs(seed) <= .Machine$integer.max)) {
    stop("'seed' must be NULL or a single number in R's integer range.")
  }
  if (!is.null(logfile) &&
        !(is.character(logfile) && length(logfile) == 1 && !is.na(logfile))) {
    stop("'logfile' must be NULL or a single character string.")
  }
  if (!.isWholeNumber(numThreads, min = 1)) {
    stop("'numThreads' must be a single whole number of at least 1.")
  }
}

## Without the embedding, the clusters can't be estimated
.checkEmbeddingZ <- function(z, embedding) {
  if (!embedding && is.null(z)) {
    stop(
      "'embedding = FALSE' needs 'z': the cell clusters are estimated ",
      "from the embedding."
    )
  }
}

## The EM loop arguments
.checkEMArgs <- function(maxIter, convergence, iterLogLik) {
  # Inf runs until convergence
  if (!(.isWholeNumber(maxIter, min = 1) || .isInf(maxIter))) {
    stop("'maxIter' must be a single whole number of at least 1, or Inf.")
  }
  # Inf calculates the log likelihood only at convergence
  if (!(.isWholeNumber(iterLogLik, min = 1) || .isInf(iterLogLik))) {
    stop("'iterLogLik' must be a single whole number of at least 1, or Inf.")
  }
  if (!.isNumber(convergence, min = 0)) {
    stop("'convergence' must be a single non-negative number.")
  }
  if (.isInf(maxIter) && convergence == 0) {
    stop(
      "'convergence' must be above 0 when 'maxIter' is Inf, or the EM ",
      "never stops."
    )
  }
  invisible(NULL)
}

## z and batch need one value per cell and no NA
.checkCellValues <- function(values, arg, numCells) {
  if (is.null(values)) {
    return(invisible(NULL))
  }
  if (length(values) != numCells) {
    stop(
      "'", arg, "' must be of the same length as the number of cells in ",
      "the 'counts' matrix."
    )
  }
  if (anyNA(values)) {
    stop("'", arg, "' must not contain NA.")
  }
  invisible(NULL)
}

## A single finite number of at least 'min'
.isNumber <- function(x, min = -Inf) {
  is.numeric(x) && length(x) == 1 && is.finite(x) && x >= min
}

## A single finite whole number of at least 'min'
.isWholeNumber <- function(x, min = -Inf) {
  .isNumber(x, min = min) && x == round(x)
}

## Positive infinity, which some arguments accept to mean "no limit"
.isInf <- function(x) {
  is.numeric(x) && length(x) == 1 && isTRUE(x == Inf)
}

## The assay must exist, so a wrong name stops before decontX runs
.checkAssayName <- function(object, assayName, arg, objectArg) {
  available <- SummarizedExperiment::assayNames(object)
  numAssays <- length(SummarizedExperiment::assays(object))
  isName <- is.character(assayName) && length(assayName) == 1 &&
    assayName %in% available
  isIndex <- .isWholeNumber(assayName, min = 1) && assayName <= numAssays
  if (!(isName || isIndex)) {
    found <- if (length(available) > 0) {
      paste0("Its assays are: ", paste(available, collapse = ", "), ".")
    } else {
      "It has no named assays."
    }
    stop(
      "'", arg, "' (", deparse(assayName), ") is not an assay of '",
      objectArg, "'. ", found
    )
  }
  invisible(NULL)
}

## Make sure provided count matrix is the right type
.checkCountsDecon <- function(counts, name = "counts") {
  if (is.null(dim(counts))) {
    stop("At least 2 genes need to have non-zero expressions.")
  }
  if (anyNA(counts)) {
    stop("Missing value in '", name, "' matrix.")
  }
  if (min(counts) < 0) {
    stop("Negative value in '", name, "' matrix.")
  }
}

## The numeric matrix classes that .asDgCMatrix() converts
.isCountsMatrix <- function(x) {
  if (is.matrix(x)) {
    return(is.numeric(x))
  }
  if (methods::is(x, "Matrix")) {
    # Logical and pattern Matrix classes aren't counts
    return(methods::is(x, "dMatrix"))
  }
  isDelayed <- methods::is(x, "DelayedMatrix") ||
    (methods::is(x, "SparseArray") && length(dim(x)) == 2)
  isDelayed && DelayedArray::type(x) %in% c("double", "integer")
}

## Stops unless x is one of those classes; 'what' names it in the message
.checkCountsMatrix <- function(x, what, choices = "", when = "", hint = "") {
  if (.isCountsMatrix(x)) {
    return(invisible(NULL))
  }
  # A matrix of the wrong type is named by its type, e.g. "character matrix"
  found <- class(x)[1]
  if (is.matrix(x)) {
    found <- paste(typeof(x), "matrix")
  } else if (methods::is(x, "DelayedArray") ||
               methods::is(x, "SparseArray")) {
    found <- paste(DelayedArray::type(x), found)
  }
  stop(
    what, " must be ", choices, "a numeric matrix, a sparse matrix or a ",
    "DelayedMatrix", when, "; it is a ", found, ".", hint
  )
}

## as(x, "dgCMatrix") is deprecated or missing for several classes; the
## extra steps turn symmetric, triangular and logical matrices into dgCMatrix
.asDgCMatrix <- function(x) {
  x <- methods::as(x, "CsparseMatrix")
  methods::as(methods::as(x, "generalMatrix"), "dMatrix")
}

## The initialization's PCA and UMAP need a few cells in each batch
.checkBatchSizes <- function(batch, batchGiven, params) {
  # uwot's 15 neighbors in the legacy UMAP; scrapper's PCA needs 3 cells.
  # Without the embedding, the EM's prior on the contamination can't be
  # estimated from one cell.
  minCells <- if (!params$embedding) 2 else if (params$legacyInit) 15 else 3
  legacy <- if (!params$embedding) {
    " to estimate the contamination prior"
  } else if (params$legacyInit) {
    " with 'legacyInit = TRUE'"
  } else {
    ""
  }
  sizes <- table(factor(batch, levels = unique(batch)))
  small <- sizes[sizes < minCells]
  if (length(small) == 0) {
    return(invisible(NULL))
  }
  if (!batchGiven) {
    stop(
      "'x' needs at least ", minCells, " cells", legacy, "; it has ",
      small[[1]], "."
    )
  }
  found <- .firstFew(paste0("batch '", names(small), "' has ", small))
  stop(
    "Each batch needs at least ", minCells, " cells", legacy, "; ",
    found, "."
  )
}


## Without background, the contamination profile of each cluster comes from
## the other clusters, so one cluster can't be fit
.checkOneCluster <- function(z, batch = NULL) {
  if (length(unique(z)) > 1) {
    return(invisible(NULL))
  }
  where <- if (is.null(batch)) {
    ""
  } else {
    paste0(" in batch '", batch, "'")
  }
  stop(
    "'z' has only one cluster", where, ". decontX needs at least two ",
    "clusters to estimate contamination without 'background'. Supply ",
    "'background' to allow one cluster."
  )
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
  if (is.factor(z)) {
    # An unused level would be an empty cluster, which gives NaN estimates
    z <- droplevels(z)
  } else {
    z <- plyr::mapvalues(z, unique(z), seq_along(unique(z)))
    z <- as.factor(z)
  }
  z
}


## The EM numbers the clusters with the codes of res$z. Name the phi and
## eta columns by each cluster's label in 'labels' (one per cell), ordered
## as the levels of 'labels'. Only the columns move, not the values.
.labelEstimates <- function(res, labels) {
  labels <- droplevels(labels)
  codes <- as.integer(res$z)
  codeLabels <- as.character(labels)[match(seq_len(ncol(res$phi)), codes)]
  columns <- match(levels(labels), codeLabels)
  phi <- res$phi[, columns, drop = FALSE]
  eta <- res$eta[, columns, drop = FALSE]
  colnames(phi) <- levels(labels)
  colnames(eta) <- levels(labels)
  list(z = labels, phi = phi, eta = eta)
}


.decontxInitializeZ <- function(counts,
                                varGenes,
                                dbscanEps,
                                estimateCellTypes,
                                seed,
                                legacyInit,
                                numThreads) {
  if (isTRUE(legacyInit)) {
    init <- .decontxInitializeZLegacy(counts,
      varGenes = varGenes,
      seed = seed,
      numThreads = numThreads
    )
  } else {
    init <- .decontxInitializeZScrapper(counts,
      varGenes = varGenes,
      seed = seed,
      numThreads = numThreads
    )
  }
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
      if (is.null(init$pcs)) {
        ## The legacy initialization clusters all genes
        ## transpose before realizing so only one dense copy is allocated;
        ## Matrix::t() because base t() fails on the sparse log counts
        cl <- stats::kmeans(as.matrix(Matrix::t(init$normed)), 2)
      } else {
        ## Cells x PCs, small and dense, instead of cells x all genes
        cl <- stats::kmeans(init$pcs, 2)
      }
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
.decontxInitializeZScrapper <- function(counts, varGenes, seed, numThreads) {
  libSizes <- .checkLibSizes(counts)
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
    num.threads = numThreads
  )
  hvg <- scrapper::chooseHighlyVariableGenes(geneVar$statistics$residuals,
    top = varGenes
  )

  ## 50 PCs matches the internal default of scater::calculateUMAP,
  ## which this pipeline replaces
  pca <- scrapper::runPca(normed[hvg, , drop = FALSE],
    number = 50,
    num.threads = numThreads
  )

  ## scrapper >= 1.5 split runUmap's 'seed' into 'initialize.seed' and
  ## 'optimize.seed'; support both APIs
  umapArgs <- list(pca$components, num.threads = numThreads)
  if (!is.null(seed)) {
    if ("seed" %in% names(formals(scrapper::runUmap))) {
      umapArgs$seed <- seed
    } else {
      umapArgs$initialize.seed <- seed
      umapArgs$optimize.seed <- seed
    }
  }
  resUmap <- do.call(scrapper::runUmap, umapArgs)

  ## The PCs, cells x PCs, are what the k-means fallback clusters
  list(umap = resUmap, pcs = t(pca$components))
}

## Original scater/scuttle-based initialization, kept for backwards
## compatibility so results from decontX 1.11.0 and earlier can be
## reproduced with legacyInit = TRUE. logNormCounts/normalizeCounts are
## deprecated upstream in favor of 'scrapper'; this path keeps working
## for as long as scater/scuttle continue to export them.
.decontxInitializeZLegacy <- function(counts, varGenes, seed, numThreads) {
  if (!.hasScater()) {
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

  resUmap <- .withSeed(
    seed,
    scater::calculateUMAP(sce,
      ntop = varGenes,
      n_threads = numThreads,
      exprs_values = "logcounts"
    )
  )

  list(normed = SingleCellExperiment::logcounts(sce), umap = resUmap)
}

## The library size of each cell, which must be above zero: the
## normalization divides by it, and the EM can't estimate an empty cell
.checkLibSizes <- function(counts) {
  libSizes <- Matrix::colSums(counts)
  if (any(libSizes == 0)) {
    stop(
      "All cells must have at least one count. Remove empty cells or ",
      "droplets before running decontX."
    )
  }
  libSizes
}

## A function of its own so tests can mock a missing scater
.hasScater <- function() {
  requireNamespace("scater", quietly = TRUE)
}


## process varGenes
.processvarGenes <- function(varGenes) {
  if (is.null(varGenes)) {
    varGenes <- 5000
  } else {
    # Inf uses all genes
    if (!(.isWholeNumber(varGenes, min = 2) || .isInf(varGenes))) {
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
    if (!.isNumber(dbscanEps, min = 0)) {
      stop("Parameter 'dbscanEps' needs to be non-negative.")
    }
  }
  dbscanEps
}

.checkDelta <- function(delta) {
  if (!is.numeric(delta) || length(delta) != 2 ||
        !all(is.finite(delta)) || any(delta <= 0)) {
    stop(
      "'delta' needs to be a numeric vector of length 2",
      " containing positive values."
    )
  }
  delta
}


## The background needs droplets, and bgBatch one value per droplet
.checkBackgroundShape <- function(background, bgBatch) {
  if (is.null(dim(background))) {
    stop("'background' must be a matrix or a SingleCellExperiment.")
  }
  if (ncol(background) == 0) {
    stop("'background' has no droplets.")
  }
  if (!is.null(bgBatch)) {
    if (length(bgBatch) != ncol(background)) {
      stop("'bgBatch' must have one value per column of 'background'.")
    }
    if (anyNA(bgBatch)) {
      stop("'bgBatch' must not contain NA.")
    }
  }
  invisible(NULL)
}

.checkBackground <- function(x, background, bgBatch,
                             logfile = NULL, verbose = FALSE) {
  .checkBackgroundShape(background, bgBatch)

  # Remove background barcodes that have already appeared in x
  # If bgBatch param is supplied, also remove duplicate bgBatch
  if (!is.null(colnames(background))) {
    dupBarcode <- colnames(background) %in% colnames(x)
  } else {
    dupBarcode <- FALSE
    warning(
      "No column names were found for the 'background' matrix. ",
      "No checking was performed between the ids in the 'background' ",
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
      bgBatch <- bgBatch[!(dupBarcode)]
    }
    if (ncol(background) == 0) {
      stop(
        "No droplets are left in 'background' after removing the barcodes ",
        "that are also in 'x'."
      )
    }
  }

  re <- list(
    background = background,
    bgBatch = bgBatch
  )

  re
}

## Match the rows of a background matrix to the rows of the counts: by name
## when both have row names, otherwise by position
.matchBackgroundRows <- function(counts,
                                 background,
                                 countsArg = "x",
                                 bgArg = "background",
                                 feature = "genes",
                                 logfile = NULL,
                                 verbose = TRUE) {
  if (is.null(dim(background))) {
    stop("'", bgArg, "' must be a matrix.")
  }
  countsNames <- rownames(counts)
  bgNames <- rownames(background)
  if (is.null(countsNames) || is.null(bgNames)) {
    if (nrow(counts) != nrow(background)) {
      stop(
        "'", bgArg, "' and '", countsArg, "' have a different number of ",
        feature, " (", nrow(background), " and ", nrow(counts), "). Rows ",
        "are matched by name when both have row names, and by position ",
        "otherwise."
      )
    }
    return(background)
  }
  if (identical(countsNames, bgNames)) {
    return(background)
  }
  .checkBackgroundNames(countsNames, bgNames, countsArg, bgArg, feature)
  extra <- setdiff(bgNames, countsNames)
  if (length(extra) > 0) {
    .logMessages(
      date(),
      " .. Dropped ", feature, " in '", bgArg, "' that aren't in '",
      countsArg, "': ", length(extra), " of ", length(bgNames), " (",
      .firstFew(extra), ")",
      sep = "",
      logfile = logfile,
      append = TRUE,
      verbose = verbose
    )
  }
  background[match(countsNames, bgNames), , drop = FALSE]
}

## Every row name of the counts must be in the background, and the names
## that are matched must be unique
.checkBackgroundNames <- function(countsNames, bgNames,
                                  countsArg, bgArg, feature) {
  # Duplicates among the extra background rows are dropped with them
  if (anyDuplicated(countsNames) ||
        anyDuplicated(bgNames[bgNames %in% countsNames])) {
    stop(
      "'", countsArg, "' or '", bgArg, "' has duplicated row names, so ",
      "their rows can't be matched by name. Give both the same row names ",
      "in the same order, or make the names unique."
    )
  }
  missing <- setdiff(countsNames, bgNames)
  if (length(missing) > 0) {
    stop(
      "'", bgArg, "' lacks ", length(missing), " of the ",
      length(countsNames), " ", feature, " in '", countsArg, "': ",
      .firstFew(missing), ". Rows are matched by name."
    )
  }
  invisible(NULL)
}

## The first few values, for messages
.firstFew <- function(x, n = 5) {
  more <- if (length(x) > n) ", ..." else ""
  paste0(paste(x[seq_len(min(n, length(x)))], collapse = ", "), more)
}
