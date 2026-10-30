# Compares decontX results between two commits with identical(). Each
# commit is installed from `git archive` into its own temporary library, and
# the paths run in a fresh R process per library. It covers the paths the
# committed references can't pin at 1e-10: no z (initialization and UMAP),
# legacyInit, and decontPro's Stan fits. Run it locally, not in CI;
# installing a commit compiles Stan, which takes several minutes.
#
# Usage, from the package root:
#   Rscript dev/compare-commits.R <old-ref> [<new-ref>]
# <new-ref> defaults to HEAD; uncommitted changes are not included. The
# path definitions come from the working tree's
# tests/testthat/helper-references.R plus the extra paths below, so both
# commits run the same code. Exits with status 1 if any path differs.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1 || length(args) > 2) {
  stop("Usage: Rscript dev/compare-commits.R <old-ref> [<new-ref>]")
}
refs <- c(old = args[1], new = if (length(args) == 2) args[2] else "HEAD")
helper <- normalizePath(file.path("tests", "testthat", "helper-references.R"))

# Install one commit into a new temporary library and return the library
installCommit <- function(ref) {
  src <- tempfile("decontX-src-")
  lib <- tempfile("decontX-lib-")
  dir.create(src)
  dir.create(lib)
  archive <- tempfile(fileext = ".tar")
  if (system2("git", c("archive", "--format=tar", "-o", archive, ref)) != 0) {
    stop("git archive failed for ", ref)
  }
  utils::untar(archive, exdir = src)
  rBin <- file.path(R.home("bin"), "R")
  status <- system2(rBin, c("CMD", "INSTALL", paste0("--library=", lib), src))
  if (status != 0) {
    stop("Installing ", ref, " failed")
  }
  lib
}

# Run every path with the decontX installed in 'lib'. Runs in a fresh R
# process, so it can't use anything defined outside the function.
runPaths <- function(helper) {
  env <- new.env(parent = asNamespace("decontX"))
  sys.source(helper, envir = env)
  local(envir = env, {
    proFixture <- function() {
      withr::with_seed(42, {
        counts <- matrix(sample(1:10, 80, replace = TRUE), nrow = 10, ncol = 8)
        ambient <- matrix(sample(1:100, 50, replace = TRUE), nrow = 10)
      })
      list(counts = counts, ambient = ambient, k = c(1, 1, 2, 2, 3, 3, 4, 4))
    }
    extraPaths <- list(
      noZ = function(d) decontX(d$sim$observedCounts, verbose = FALSE),
      noZSce = function(d) {
        sce <- SingleCellExperiment::SingleCellExperiment(
          list(counts = d$sim$observedCounts)
        )
        decontX(sce, verbose = FALSE)
      },
      legacyInit = function(d) {
        if (!requireNamespace("scater", quietly = TRUE)) {
          return("skipped: scater is not installed")
        }
        suppressWarnings(decontX(d$sim$observedCounts,
          legacyInit = TRUE, verbose = FALSE
        ))
      },
      decontPro = function(d) {
        p <- proFixture()
        decontPro(p$counts, p$k, 1e-2, 1e-2)
      },
      decontProAmbient = function(d) {
        p <- proFixture()
        decontPro(p$counts, p$k, 1e-2, 1e-2, ambient_counts = p$ambient)
      }
    )
  })
  d <- env$refData()
  paths <- c(env$refPaths, env$extraPaths)
  lapply(paths, function(path) {
    tryCatch(path(d), error = function(e) {
      paste("error:", conditionMessage(e))
    })
  })
}

results <- list()
for (side in names(refs)) {
  message("Installing ", side, " (", refs[[side]], ")")
  lib <- installCommit(refs[[side]])
  message("Running the paths for ", side)
  results[[side]] <- callr::r(runPaths,
    args = list(helper = helper),
    libpath = c(lib, .libPaths())
  )
}

same <- TRUE
for (path in union(names(results$old), names(results$new))) {
  old <- results$old[[path]]
  new <- results$new[[path]]
  if (identical(old, new)) {
    message(path, ": identical")
  } else {
    same <- FALSE
    difference <- all.equal(old, new, tolerance = 0)
    message(path, ": DIFFERS\n  ", paste(utils::head(difference, 5),
      collapse = "\n  "
    ))
  }
}
quit(status = if (same) 0 else 1)
