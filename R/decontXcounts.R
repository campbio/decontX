#' @title Get or set decontaminated counts matrix
#'
#' @description Gets or sets the decontaminated counts matrix from a
#' a \linkS4class[SingleCellExperiment]{SingleCellExperiment} object.
#' @name decontXcounts
#' @param object A
#' \linkS4class[SingleCellExperiment]{SingleCellExperiment} object.
#' @param value A matrix to save as an assay called \code{decontXcounts}
#' @param ... For the generic, further arguments to pass to each method.
#' @return If getting, the assay from \code{object} with the name
#' \code{decontXcounts} will be returned. If setting, a
#' \linkS4class[SingleCellExperiment]{SingleCellExperiment} object will be
#' returned with
#' \code{decontXcounts} listed in the \code{assay} slot.
#' @seealso \code{\link[SummarizedExperiment]{assay}} and
#' \code{\link[SummarizedExperiment]{assay<-}}
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
setMethod("decontXcounts", "SingleCellExperiment", function(object, ...) {
  SummarizedExperiment::assay(object, i = "decontXcounts", ...)
})

#' @export
#' @rdname decontXcounts
setMethod(
  "decontXcounts<-", c("SingleCellExperiment", "ANY"),
  function(object, ..., value) {
    SummarizedExperiment::assay(object, i = "decontXcounts", ...) <- value
    object
  }
)
