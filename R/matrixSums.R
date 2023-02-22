.colSumByGroup <- function(counts, group, K) {
  if (inherits(counts, "matrix") & is.integer(counts)) {
    res <- .colSumByGroupInteger(counts, group, K)
  } else if (inherits(counts, "matrix") & is.numeric(counts)) {
    res <- .colSumByGroupNumeric(counts, group, K)
  } else if (inherits(counts, "dgCMatrix")) {
    res <- colSumByGroupSparse(counts, group, K)
  } else {
    stop("'counts' must be an integer, numeric, or dgCMatrix matrix.")
  }
  return(res)
}

#' @useDynLib decontX _colSumByGroup_numeric
.colSumByGroupNumeric <- function(x, group, K) {
  group <- factor(group, levels = seq(K))
  res <- .Call("_colSumByGroup_numeric", x, group)
  return(res)
}

#' @useDynLib decontX _colSumByGroup
.colSumByGroupInteger <- function(x, group, K) {
  group <- factor(group, levels = seq(K))
  res <- .Call("_colSumByGroup", x, group)
  return(res)
}
