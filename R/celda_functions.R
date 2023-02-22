.logMessages <- function(...,
                         sep = " ",
                         logfile = NULL,
                         append = FALSE,
                         verbose = TRUE) {
  if (isTRUE(verbose)) {
    if (!is.null(logfile)) {
      if (!is.character(logfile) || length(logfile) > 1) {
        stop(
          "The log file parameter needs to be a single character",
          " string."
        )
      }
      cat(paste(..., "\n", sep = sep),
          file = logfile,
          append = append
      )
    } else {
      message(paste(..., sep = sep))
    }
  }
}


#' @title Normalization of count data
#' @description Performs normalization, transformation, and/or scaling of a
#'  counts matrix
#' @param counts Integer, Numeric or Sparse matrix. Rows represent features
#' and columns represent cells.
#' @param normalize Character.
#'  Divides counts by the library sizes for each cell. One of 'proportion',
#'  'cpm', 'median', or 'mean'. 'proportion' uses the total counts for each
#'  cell as the library size. 'cpm' divides the library size of each cell by
#'  one million to produce counts per million. 'median' divides the library
#'  size of each cell by the median library size across all cells. 'mean'
#'  divides the library size of each cell by the mean library size across all
#'  cells.
#' @param scaleFactor Numeric. Sets the scale factor for cell-level
#'  normalization. This scale factor is multiplied to each cell after the
#'  library size of each cell had been adjusted in \code{normalize}. Default
#'  \code{NULL} which means no scale factor is applied.
#' @param transformationFun Function. Applys a transformation such as
#'  \link{sqrt}, \link{log}, \link{log2}, \link{log10}, or \link{log1p}.
#'  If NULL, no transformation will be applied. Occurs after normalization.
#'  Default NULL.
#' @param scaleFun Function. Scales the rows of the normalized and transformed
#'  count matrix. For example, 'scale' can be used to z-score normalize the
#'  rows. Default NULL.
#' @param pseudocountNormalize Numeric. Add a pseudocount to counts before
#'  normalization. Default 0.
#' @param pseudocountTransform Numeric. Add a pseudocount to normalized counts
#'  before applying the transformation function. Adding a pseudocount
#'  can be useful before applying a log transformation. Default  0.
#' @return Numeric Matrix. A normalized matrix.
#' @examples
#' data(celdaCGSim)
#' normalizedCounts <- normalizeCounts(celdaCGSim$counts, "proportion",
#'   pseudocountNormalize = 1)
#' @importFrom Matrix colSums
#' @export
normalizeCounts <- function(counts,
                            normalize = c("proportion", "cpm",
                                          "median", "mean"),
                            scaleFactor = NULL,
                            transformationFun = NULL,
                            scaleFun = NULL,
                            pseudocountNormalize = 0,
                            pseudocountTransform = 0) {

  normalize <- match.arg(normalize)

  if (!is.null(transformationFun) &&
      !is.function(transformationFun)) {
    stop("'transformationFun' needs to be of class 'function'")
  }
  if (!is.null(scaleFun) && !is.function(scaleFun)) {
    stop("'scaleFun' needs to be of class 'function'")
  }
  # Perform normalization
  if (normalize == "proportion" & inherits(counts, "matrix")) {
    norm <- fastNormProp(counts, pseudocountNormalize)
  } else {
    counts <- counts + pseudocountNormalize
    cs <- colSums(counts)
    norm <- switch(
      normalize,
      "proportion" = sweep(counts, 2, cs, "/"),
      "cpm" = sweep(counts, 2, cs / 1e6, "/"),
      "median" = sweep(counts, 2, cs / stats::median(cs), "/"),
      "mean" = sweep(counts, 2, cs / mean(cs), "/")
    )
  }

  if (!is.null(scaleFactor)) {
    norm <- norm * scaleFactor
  }

  if (!is.null(transformationFun)) {
    norm <- do.call(
      transformationFun,
      list(norm + pseudocountTransform)
    )
  }
  if (!is.null(scaleFun)) {
    norm <- t(base::apply(norm, 1, scaleFun))
  }
  colnames(norm) <- colnames(counts)
  rownames(norm) <- rownames(counts)
  return(norm)
}

# Generate n random deviates from the Dirichlet function with shape parameters
# alpha. Adapted from gtools v3.5
.rdirichlet <- function(n, alpha) {
  l <- length(alpha)
  x <- matrix(stats::rgamma(l * n, alpha),
              ncol = l,
              byrow = TRUE
  )
  # Check for case where all sampled entries are zero due to round off
  # One entry will be randomly chosen to be one
  isZero <- rowSums(x) == 0
  assignment <- sample(seq(l), size = sum(isZero), replace = TRUE)
  x[cbind(which(isZero), assignment)] <- 1
  # Normalize
  sm <- x %*% rep(1, l)
  y <- x / as.vector(sm)
  return(y)
}

#' @title Retrieve row index for a set of features
#' @description This will return indices of features among the rownames
#' or rowData of a data.frame, matrix, or a \linkS4class{SummarizedExperiment}
#' object including a \linkS4class{SingleCellExperiment}.
#' Partial matching (i.e. grepping) can be used by setting
#' \code{exactMatch = FALSE}.
#' @param features Character vector of feature names to find in the rows of
#' \code{x}.
#' @param x A data.frame, matrix, or \linkS4class{SingleCellExperiment}
#' object to search.
#' @param by Character. Where to search for features in \code{x}. If set to
#' \code{"rownames"} then the features will be searched for among
#' \code{rownames(x)}. If \code{x} inherits from class
#' \linkS4class{SummarizedExperiment}, then \code{by} can be one of the
#' fields in the row annotation data.frame (i.e. one of
#' \code{colnames(rowData(x))}).
#' @param exactMatch Boolean. Whether to only identify exact matches
#' or to identify partial matches using \code{\link{grep}}.
#' @param removeNA Boolean. If set to \code{FALSE}, features not found in
#' \code{x} will be given \code{NA} and the returned vector will be the same
#' length as \code{features}. If set to \code{TRUE}, then the \code{NA}
#' values will be removed from the returned vector. Default \code{FALSE}.
#' @return A vector of row indices for the matching features in \code{x}.
#' @author Yusuke Koga, Joshua Campbell
#' @seealso '\link[scater]{retrieveFeatureInfo}' from package \code{'scater'}
#' and \code{link{regex}} for how to use regular expressions when
#' \code{exactMatch = FALSE}.
#' @examples
#' data(celdaCGSim)
#' retrieveFeatureIndex(c("Gene_1", "Gene_5"), celdaCGSim$counts)
#' retrieveFeatureIndex(c("Gene_1", "Gene_5"), celdaCGSim$counts,
#'                                             exactMatch = FALSE)
#' @export
retrieveFeatureIndex <- function(features,
                                 x,
                                 by = "rownames",
                                 exactMatch = TRUE,
                                 removeNA = FALSE) {

  # Extract vector to search through
  if (by == "rownames") {
    if (is.null(rownames(x))) {
      stop("'rownames' of 'x' are 'NULL'. Please set 'rownames' or change",
           " 'by' to search a different column in 'x'.")
    }
    search <- rownames(x)
  } else if (length(ncol(x)) > 0) {
    if (inherits(x, "SummarizedExperiment")) {
      if (!(by %in% colnames(SummarizedExperiment::rowData(x)))) {
        stop("'by' is not a column in 'rowData(x)'.")
      }
      search <- SummarizedExperiment::rowData(x)[, by]
    } else {
      if (!(by %in% colnames(x))) {
        stop("'by' is not a column in 'x'.")
      }
      search <- x[, by]
    }
  } else {
    search <- as.character(x)
  }

  # Match each element of 'pattern' in vector 'search'
  if (!isTRUE(exactMatch)) {
    featuresIndices <- rep(NA, length(features))
    for (i in seq_along(features)) {
      g <- grep(features[i], search)
      if (length(g) == 1) {
        featuresIndices[i] <- g
      } else if (length(g) > 1) {
        warning(
          "Feature '", features[i], "' matched multiple items in '",
          by, "': ", paste(search[g], collapse = ","),
          ". Only the first match will be selected."
        )
        featuresIndices[i] <- g[1]
      }
    }
  } else {
    featuresIndices <- match(features, search)
  }

  if (sum(is.na(featuresIndices)) > 0) {
    if (sum(is.na(featuresIndices)) == length(features)) {
      if (isTRUE(exactMatch)) {
        stop(
          "None of the provided features had matching",
          " items in '", by, "' within 'x'. ",
          "Check the spelling or try setting",
          " 'exactMatch = FALSE'."
        )
      } else {
        stop(
          "None of the provided features had matching",
          " items in '", by, "' within 'x'. ",
          "Check the spelling and make sure 'by' is set",
          " to the appropriate place in 'x'."
        )
      }
    }
    warning(
      "The following features were not present in 'x': ",
      paste(features[which(is.na(featuresIndices))],
            collapse = ","
      )
    )
  }

  if (isTRUE(removeNA)) {
    featuresIndices <- featuresIndices[!is.na(featuresIndices)]
  }
  return(featuresIndices)
}
