#' Decontaminate using decontPro
#'
#' @name decontPro
#'
#' @param filtered_counts Count matrix NxM (feature x droplet) with only
#' filtered droplets after cell calling, with at least 2 ADTs and 2
#' droplets. The counts must be non-negative whole numbers with no NA.
#' If this is a \linkS4class[SingleCellExperiment]{SingleCellExperiment} or a
#' \linkS4class[SeuratObject]{Seurat} object,
#' expect \code{counts} in the assay slot. Split \code{counts} layers of a
#' Seurat v5 assay are joined first, so the droplets follow
#' \code{Cells()}.
#' @param cellType Vector of length M giving the cell type of each droplet,
#' with no NA. Any labels can be used, such as integers, character strings
#' or a factor. The rows of \code{r_est} in the returned parameters follow
#' the sorted labels, or a factor's levels with unused levels dropped.
#' Alias: \code{cell_type}.
#' @param deltaSd Standard deviation of the half-normal prior on how much
#' the ambient contamination levels vary around each droplet's mean level.
#' Smaller values keep the levels of a droplet closer together. Default
#' 2e-5. Alias: \code{delta_sd}.
#' @param backgroundSd Standard deviation of the half-normal prior on how
#' much the levels of the model's background component vary around each
#' ADT's mean level. This background is part of the model, not the
#' \code{background} counts. Default 2e-6. Alias: \code{background_sd}.
#' @param background Count matrix NxM (feature x droplet) of empty droplets,
#' used to estimate the ambient ADT profile; without it, the profile comes
#' from \code{filtered_counts}. Similar to \code{filtered_counts}, if it is
#' a wrapper object, expect \code{counts} in the assay slot. When both
#' have row names, its rows are matched to the ADTs of
#' \code{filtered_counts} by name: ADTs that aren't in
#' \code{filtered_counts} are dropped, and every ADT of
#' \code{filtered_counts} must be present. Otherwise the rows are matched
#' by position. The ADTs used must have no NA or negative counts.
#' Default NULL. Alias: \code{ambient_counts}.
#' @param ... Only the aliases \code{cell_type}, \code{delta_sd},
#' \code{background_sd} and \code{ambient_counts}, which can be used in
#' place of \code{cellType}, \code{deltaSd}, \code{backgroundSd} and
#' \code{background}. Giving both names of an argument, or any other
#' argument, stops with an error.
#' @param iter Maximum number of iterations of Stan's variational
#' inference. The fit usually stops earlier, when it converges; each
#' iteration takes longer on larger data. Default 50000.
#' @param seed Seed for Stan's random number generator, a whole number
#' from 0 to \code{.Machine$integer.max}. Default 12345.
#' @param verbose Whether to print messages and Stan's progress. With
#' \code{FALSE}, Stan's output is shown only if the fit fails. Default TRUE.
#'
#' @return A list containing decontaminated counts, and estimated parameters.
#' The row names of \code{parameters$r_est}, the native ADT profile of each
#' cell type, are the \code{cellType} labels.
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
#' background <- matrix(
#'   sample(1:2,
#'     1000,
#'     replace = TRUE
#'   ),
#'   ncol = 10
#' )
#'
#' # Decontamination
#' out <- decontPro(counts, k, 1e-2, 1e-2, background)
#'
#' # Decontaminated counts
#' decontaminated_counts <- out$decontaminated_counts
setGeneric("decontPro", function(filtered_counts,
                                 cellType,
                                 ...) {
  standardGeneric("decontPro")
}, signature = "filtered_counts")


#' @export
#' @rdname decontPro
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
setMethod("decontPro", "SingleCellExperiment", function(filtered_counts,
                                                        cellType,
                                                        deltaSd = 2e-5,
                                                        backgroundSd = 2e-6,
                                                        background = NULL,
                                                        ...,
                                                        iter = 50000,
                                                        seed = 12345,
                                                        verbose = TRUE) {
  args <- .decontProArgs(environment(), ...)
  counts <- SummarizedExperiment::assay(filtered_counts, "counts")
  if (!is.null(args[["background"]])) {
    args[["background"]] <- SummarizedExperiment::assay(
      args[["background"]], "counts"
    )
  }
  .decontPro(counts, args)
})


#' @export
#' @rdname decontPro
#' @importClassesFrom SeuratObject Seurat
setMethod("decontPro", "Seurat", function(filtered_counts,
                                          cellType,
                                          deltaSd = 2e-5,
                                          backgroundSd = 2e-6,
                                          background = NULL,
                                          ...,
                                          iter = 50000,
                                          seed = 12345,
                                          verbose = TRUE) {
  args <- .decontProArgs(environment(), ...)
  counts <- .seuratCounts(filtered_counts)
  if (!is.null(args[["background"]])) {
    args[["background"]] <- .seuratCounts(args[["background"]])
  }
  .decontPro(counts, args)
})


#' @export
#' @rdname decontPro
setMethod("decontPro", "ANY", function(filtered_counts,
                                       cellType,
                                       deltaSd = 2e-5,
                                       backgroundSd = 2e-6,
                                       background = NULL,
                                       ...,
                                       iter = 50000,
                                       seed = 12345,
                                       verbose = TRUE) {
  args <- .decontProArgs(environment(), ...)
  .decontPro(filtered_counts, args)
})


# The permanent snake_case alias of each camelCase argument (ADR-0006)
.decontProAliases <- c(
  cellType = "cell_type",
  deltaSd = "delta_sd",
  backgroundSd = "background_sd",
  background = "ambient_counts"
)


# Collects the arguments of a decontPro() method from its frame `env` and
# the aliases in `...`. `argNames` records the name the user gave for each
# argument, so error messages use it.
.decontProArgs <- function(env, ...) {
  # Unknown arguments are reported before anything in `...` is evaluated
  dotNames <- ...names()
  if (is.null(dotNames)) {
    dotNames <- rep("", ...length())
  }
  .checkDecontProDots(dotNames[!dotNames %in% .decontProAliases])
  args <- list()
  argNames <- character()
  for (arg in names(.decontProAliases)) {
    alias <- .decontProAliases[[arg]]
    aliasCount <- sum(dotNames == alias)
    given <- !eval(call("missing", as.name(arg)), envir = env)
    if (aliasCount > 1) {
      stop("'", alias, "' is given more than once.")
    }
    if (given && aliasCount == 1) {
      stop(
        "Give '", arg, "' or its alias '", alias, "', not both. A value ",
        "given by position fills '", arg, "'."
      )
    }
    if (aliasCount == 1) {
      args[arg] <- list(...elt(which(dotNames == alias)))
      argNames[[arg]] <- alias
    } else if (given || arg != "cellType") {
      args[arg] <- list(get(arg, envir = env))
      argNames[[arg]] <- arg
    }
  }
  if (!"cellType" %in% names(args)) {
    stop("'cellType' is missing. Give the cell type of each droplet.")
  }
  args[c("iter", "seed", "verbose")] <- mget(
    c("iter", "seed", "verbose"),
    envir = env
  )
  args$argNames <- argNames
  args
}


# `background = <number>` used to be a partial match for background_sd, so
# a single number most likely means the prior
.checkBackgroundNotSd <- function(background, arg) {
  if (is.numeric(background) && length(background) == 1) {
    stop(
      "'", arg, "' must be a matrix of empty-droplet counts. To set the ",
      "prior of the background component, use 'backgroundSd' (or its ",
      "alias 'background_sd')."
    )
  }
}


# iter, seed and verbose are passed on to rstan
.checkDecontProFitArgs <- function(args) {
  if (!(.isWholeNumber(args$iter, min = 1) &&
          args$iter <= .Machine$integer.max)) {
    stop(
      "'iter' must be a single whole number from 1 to ",
      ".Machine$integer.max."
    )
  }
  if (!(.isWholeNumber(args$seed, min = 0) &&
          args$seed <= .Machine$integer.max)) {
    stop(
      "'seed' must be a single whole number from 0 to ",
      ".Machine$integer.max."
    )
  }
  .checkFlags(list(verbose = args$verbose))
}


# `args` comes from .decontProArgs(), with the counts of `background`
# already extracted from a wrapper object
.decontPro <- function(counts, args) {
  .checkDecontProFitArgs(args)
  ## Prep data
  counts <- as.matrix(counts)
  .checkDecontProCounts(counts)
  N <- nrow(counts)
  M <- ncol(counts)
  .checkDecontProDims(N, M)
  cellType <- .recodeCellType(args$cellType, M, args$argNames[["cellType"]])

  if (is.null(args[["background"]])) {
    temp <- counts
  } else {
    bgArg <- args$argNames[["background"]]
    .checkBackgroundNotSd(args[["background"]], bgArg)
    temp <- .matchBackgroundRows(
      counts = counts,
      background = args[["background"]],
      countsArg = "filtered_counts",
      bgArg = bgArg,
      feature = "ADTs",
      verbose = args$verbose
    )
    .checkCountsDecon(temp, bgArg)
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
    delta_sd = args$deltaSd,
    background_sd = args$backgroundSd
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
  out <- .callStanVb(
    dat, init,
    iter = args$iter, seed = args$seed, verbose = args$verbose
  )


  ## Process Stan output
  dat$cell_type_labels <- cellType$labels
  re <- .processStanVbOut(out, dat)


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


# rstan passes length-1 data to Stan as a scalar, so one ADT or one
# droplet would fail inside Stan with an unclear message (#110)
.checkDecontProDims <- function(N, M) {
  if (N < 2 || M < 2) {
    stop(
      "'filtered_counts' has ", N, ngettext(N, " ADT", " ADTs"), " and ",
      M, ngettext(M, " droplet", " droplets"),
      "; decontPro needs at least 2 of each."
    )
  }
}


# Stan needs cell types coded 1..K with every code used. Labels are sorted,
# or kept in a factor's level order with unused levels dropped. `arg` is
# the name the user gave the argument.
.recodeCellType <- function(cell_type, M, arg = "cellType") {
  if (length(cell_type) != M) {
    stop(
      "'", arg, "' has ", length(cell_type), " entries but ",
      "'filtered_counts' has ", M, " droplets (columns)."
    )
  }
  if (anyNA(cell_type)) {
    stop("'", arg, "' has NA values. Give every droplet a cell type.")
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


# `argNames` holds the names of the arguments in `...` that aren't
# aliases, "" for an unnamed one
.checkDecontProDots <- function(argNames) {
  if (length(argNames) == 0) {
    return(invisible(NULL))
  }
  argNames <- ifelse(argNames == "", "an unnamed value",
    paste0("'", argNames, "'")
  )
  stop(
    "Unknown arguments to decontPro(): ",
    paste(argNames, collapse = ", "), ". See ?decontPro."
  )
}
