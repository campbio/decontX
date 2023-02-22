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
