test_that(desc = "decontX does not call deprecated normalization functions", {
  deconSim <- simulateContamination(K = 10, delta = c(1, 5))
  # Regression guard for the Bioc 3.24 scuttle deprecations
  # (logNormCounts/normalizeCounts -> scrapper): decontX initialization
  # must not emit these deprecation warnings, which R CMD check flags
  # and which become errors when upstream turns them defunct. Scoped to
  # the functions decontX calls so unrelated upstream deprecations
  # cannot fail this test.
  expect_no_warning(
    decontX(deconSim$observedCounts, seed = 1, maxIter = 2, verbose = FALSE),
    message = "'(normalizeCounts|logNormCounts)' is deprecated"
  )
})

test_that(desc = "decontX initialization without z returns valid results", {
  deconSim <- simulateContamination(K = 10, delta = c(1, 5))
  res <- decontX(deconSim$observedCounts, seed = 1, verbose = FALSE)
  expect_true(all(res$contamination >= 0 & res$contamination <= 1))
  expect_equal(length(res$z), ncol(deconSim$observedCounts))
  expect_gt(length(unique(res$z)), 1)
  umap <- res$estimates$all_cells$UMAP
  expect_equal(dim(umap), c(ncol(deconSim$observedCounts), 2))

  # Same seed twice gives identical results
  res2 <- decontX(deconSim$observedCounts, seed = 1, verbose = FALSE)
  expect_equal(res$contamination, res2$contamination)
  expect_equal(res$z, res2$z)
})

test_that(desc = "decontX initialization handles very sparse data", {
  # Regression test: scrapper::modelGeneVariances' default abundance
  # filter can drop every gene in low-count data (e.g. shallow
  # sequencing); initialization must rank all genes instead of erroring.
  withr::local_seed(7)
  lowCounts <- matrix(stats::rbinom(400 * 100, size = 1, prob = 0.05),
    nrow = 400, ncol = 100,
    dimnames = list(
      paste0("g", seq_len(400)),
      paste0("c", seq_len(100))
    )
  )
  lowCounts[1, Matrix::colSums(lowCounts) == 0] <- 1
  res <- decontX(lowCounts, seed = 1, verbose = FALSE)
  expect_true(all(res$contamination >= 0 & res$contamination <= 1))
})

test_that(desc = "decontX initialization handles duplicated gene names", {
  deconSim <- simulateContamination(K = 10, delta = c(1, 5))
  # Regression test: real 10x data labeled with gene symbols often has
  # duplicated rownames, which scrapper::modelGeneVariances rejects.
  counts <- deconSim$observedCounts
  rownames(counts) <- c(
    "DUP", "DUP",
    paste0("g", seq_len(nrow(counts) - 2))
  )
  res <- decontX(counts, seed = 1, verbose = FALSE)
  expect_true(all(res$contamination >= 0 & res$contamination <= 1))
})

test_that(desc = "decontX gives an informative error for empty cells", {
  deconSim <- simulateContamination(K = 10, delta = c(1, 5))
  counts <- deconSim$observedCounts
  counts[, 1] <- 0
  expect_error(decontX(counts, seed = 1, maxIter = 2, verbose = FALSE),
    regexp = "at least one count"
  )
})

test_that("legacyInit with z given changes only the UMAP", {
  skip_if_not_installed("scater")
  d <- refData()
  x <- d$sim$observedCounts
  # The legacy path calls upstream functions that are deprecated on
  # Bioconductor devel (by design -- it exists to reproduce old results
  # while they still work), so warnings are suppressed here.
  default <- decontX(x, z = d$sim$z, verbose = FALSE)
  legacy <- suppressWarnings(
    decontX(x, z = d$sim$z, legacyInit = TRUE, verbose = FALSE)
  )
  expect_identical(legacy$contamination, default$contamination)
  expect_identical(legacy$decontXcounts, default$decontXcounts)
  umap <- legacy$estimates$all_cells$UMAP
  expect_identical(dim(umap), c(ncol(x), 2L))
  expect_false(anyNA(umap))
  expect_false(isTRUE(all.equal(umap, default$estimates$all_cells$UMAP)))
})

test_that("legacyInit without z estimates the clusters", {
  skip_if_not_installed("scater")
  deconSim <- simulateContamination(K = 10, delta = c(1, 5))
  res <- suppressWarnings(
    decontX(deconSim$observedCounts,
      seed = 1, legacyInit = TRUE,
      verbose = FALSE
    )
  )
  expect_true(all(res$contamination >= 0 & res$contamination <= 1))
  expect_equal(length(res$z), ncol(deconSim$observedCounts))
  expect_gt(length(unique(res$z)), 1)
})

test_that("numeric batch labels are not used as positional indices", {
  # Different seeds, so a misaligned batch would change the results
  sim1 <- simulateContamination(K = 5, delta = c(1, 5), seed = 1)
  sim2 <- simulateContamination(K = 5, delta = c(1, 5), seed = 2)
  combined <- cbind(sim1$observedCounts, sim2$observedCounts)
  z <- c(sim1$z, sim2$z)
  n1 <- ncol(sim1$observedCounts)
  n2 <- ncol(sim2$observedCounts)

  numBatch <- c(rep(3, n1), rep(5, n2))
  charBatch <- c(rep("3", n1), rep("5", n2))

  resNum <- decontX(combined,
    z = z, batch = numBatch,
    seed = 12345, verbose = FALSE
  )
  resChar <- decontX(combined,
    z = z, batch = charBatch,
    seed = 12345, verbose = FALSE
  )

  expect_equal(resNum$contamination, resChar$contamination)
  expect_equal(resNum$decontXcounts, resChar$decontXcounts)
  expect_true(all(resNum$contamination >= 0 & resNum$contamination <= 1))
  expect_false(is.null(resNum$estimates[["3"]]))
  expect_false(is.null(resNum$estimates[["5"]]))
})

test_that("retrieveFeatureIndex never matches an NA feature", {
  counts <- matrix(1:6, nrow = 3, dimnames = list(c("A", NA, "C"), NULL))
  expect_warning(
    idx <- retrieveFeatureIndex(c("A", NA), counts),
    regexp = "not present in 'x': NA"
  )
  expect_equal(idx, c(1L, NA))
  # Only the "not present" warning, not one about NA matching every row
  warnings <- capture_warnings(
    idx <- retrieveFeatureIndex(c("A", NA), counts, exactMatch = FALSE)
  )
  expect_length(warnings, 1)
  expect_match(warnings, "not present in 'x': NA")
  expect_equal(idx, c(1L, NA))
})

test_that("retrieveFeatureIndex needs 'by' to be a single string", {
  counts <- matrix(1:6, nrow = 3, dimnames = list(c("A", "B", "C"), NULL))
  expect_error(retrieveFeatureIndex("A", counts, by = c("rownames", "x")),
    regexp = "'by' must be a single character string"
  )
  expect_error(retrieveFeatureIndex("A", counts, by = NA),
    regexp = "'by' must be a single character string"
  )
})

test_that("runParams holds seed and iterLogLik", {
  deconSim <- simulateContamination(K = 10, delta = c(1, 5))
  res <- decontX(deconSim$observedCounts,
    z = deconSim$z, maxIter = 2,
    iterLogLik = 5, seed = 3, verbose = FALSE
  )
  expect_equal(res$runParams$seed, 3)
  expect_equal(res$runParams$iterLogLik, 5)
})

test_that("runParams always has batchBackground", {
  d <- refData()
  res <- decontX(d$sim$observedCounts, z = d$sim$z, verbose = FALSE)
  expect_true("batchBackground" %in% names(res$runParams))
  expect_null(res$runParams$batchBackground)

  res <- decontX(d$sim$observedCounts,
    z = d$sim$z, background = d$bg,
    verbose = FALSE
  )
  expect_identical(res$runParams$batchBackground, rep("all_cells", 40))
})

test_that("the matrix result has the batch of each cell", {
  d <- refData()
  res <- decontX(d$sim$observedCounts, z = d$sim$z, verbose = FALSE)
  expect_identical(res$batch, rep("all_cells", 100))

  batch <- factor(d$batch)
  res <- decontX(d$sim$observedCounts,
    z = d$sim$z, batch = batch,
    verbose = FALSE
  )
  expect_identical(res$batch, as.character(batch))
  expect_identical(
    names(res),
    c("runParams", "estimates", "decontXcounts", "contamination", "z", "batch")
  )
})

test_that("interleaved batches give each batch's own results", {
  d <- refData()
  counts <- d$sim$observedCounts
  batch <- rep(c("A", "B"), 50)
  res <- decontX(counts, z = d$sim$z, batch = batch, verbose = FALSE)
  expect_identical(dimnames(res$decontXcounts), dimnames(counts))
  for (b in c("A", "B")) {
    inBatch <- batch == b
    alone <- decontX(counts[, inBatch],
      z = d$sim$z[inBatch],
      verbose = FALSE
    )
    expect_identical(
      res$decontXcounts[, inBatch],
      alone$decontXcounts
    )
    expect_identical(res$contamination[inBatch], alone$contamination)
  }
})

test_that("decontXcounts drops the names of the counts' dimnames list", {
  d <- refData()
  counts <- d$sim$observedCounts
  names(dimnames(counts)) <- c("genes", "cells")
  for (batch in list(NULL, d$batch)) {
    res <- decontX(counts, z = d$sim$z, batch = batch, verbose = FALSE)
    expect_identical(dimnames(res$decontXcounts), unname(dimnames(counts)))
  }
})

test_that(".bindBatches and .inCellOrder put the columns in cell order", {
  counts <- Matrix::sparseMatrix(
    i = c(1, 2, 1, 2, 1), j = c(1, 2, 3, 4, 5), x = c(1, 2, 3, 0, 5),
    dims = c(2, 5)
  )
  bindInCellOrder <- function(batch, dimNames) {
    batchCounts <- lapply(unique(batch), function(b) {
      counts[, batch == b, drop = FALSE]
    })
    .inCellOrder(.bindBatches(batchCounts), batch, dimNames)
  }
  withNames <- function(dimNames) {
    expected <- counts
    dimnames(expected) <- dimNames
    expected
  }

  # Interleaved batches, one of them labelled with an empty string
  dimNames <- list(c("g1", "g2"), paste0("c", seq_len(5)))
  for (batch in list(c("A", "B", "A", "B", "A"), c("", "B", "", "B", ""))) {
    expect_identical(bindInCellOrder(batch, dimNames), withNames(dimNames))
  }

  # Cells already in batch order, and no column names
  dimNames <- list(c("g1", "g2"), NULL)
  expect_identical(
    bindInCellOrder(c("A", "A", "A", "B", "B"), dimNames),
    withNames(dimNames)
  )

  # One batch comes back as it is
  expect_identical(.bindBatches(list(counts)), counts)

  # One gene stays a matrix
  oneGene <- counts[1, , drop = FALSE]
  batchCounts <- list(
    oneGene[, c(1, 3, 5), drop = FALSE],
    oneGene[, c(2, 4), drop = FALSE]
  )
  expected <- oneGene
  dimnames(expected) <- list("g1", NULL)
  expect_identical(
    .inCellOrder(.bindBatches(batchCounts), c("A", "B", "A", "B", "A"),
      list("g1", NULL)
    ),
    expected
  )
})

test_that("decontX keeps the SCE results out of the metadata", {
  d <- refData()
  sce <- SingleCellExperiment::SingleCellExperiment(
    list(counts = d$sim$observedCounts)
  )
  sce <- decontX(sce, z = d$sim$z, batch = d$batch, verbose = FALSE)
  expect_identical(
    names(S4Vectors::metadata(sce)$decontX),
    c("runParams", "estimates", "contamination", "z")
  )
})

test_that("batched SCE runs store one padded UMAP per batch", {
  d <- refData()
  sce <- SingleCellExperiment::SingleCellExperiment(
    list(counts = d$sim$observedCounts)
  )
  sce <- decontX(sce, z = d$sim$z, batch = d$batch, verbose = FALSE)
  expect_identical(
    SingleCellExperiment::reducedDimNames(sce),
    c("decontX_A_UMAP", "decontX_B_UMAP")
  )
  estimates <- S4Vectors::metadata(sce)$decontX$estimates
  for (b in c("A", "B")) {
    name <- paste0("decontX_", b, "_UMAP")
    umap <- SingleCellExperiment::reducedDim(sce, name)
    inBatch <- d$batch == b
    expect_identical(dimnames(umap), list(colnames(sce), c("UMAP_1", "UMAP_2")))
    expect_true(all(is.na(umap[!inBatch, ])))
    expect_identical(unname(umap[inBatch, ]), unname(estimates[[b]]$UMAP))
  }
})

test_that("the decontX core and its params constructor have no defaults", {
  noDefault <- function(fn) {
    fmls <- formals(fn)
    vapply(names(fmls), function(n) identical(fmls[[n]], quote(expr = )),
      logical(1)
    )
  }
  expect_true(all(noDefault(.decontX)))
  expect_true(all(noDefault(.decontXParams)))
})

test_that("the decontX methods have the same arguments and defaults", {
  sceArgs <- decontXMethodFormals("SingleCellExperiment")
  anyArgs <- decontXMethodFormals("ANY")
  containerArgs <- c("x", "assayName", "bgAssayName")
  expect_identical(
    sceArgs[setdiff(names(sceArgs), containerArgs)],
    anyArgs[setdiff(names(anyArgs), containerArgs)]
  )
  dataArgs <- c("z", "batch", "background", "bgBatch")
  expect_identical(
    setdiff(names(anyArgs), c(containerArgs, dataArgs)),
    names(formals(.decontXParams))
  )
  # Every parameter reaches the list that runParams records
  expect_identical(names(defaultParams()), names(formals(.decontXParams)))
})

test_that(".decontXUMAPs pads each batch's UMAP to every cell", {
  umapA <- matrix(1:4, 2, dimnames = list(c("c1", "c3"), c("u1", "u2")))
  umapB <- matrix(5:6, 1, dimnames = list("c2", c("u1", "u2")))
  single <- list(estimates = list(all_cells = list(UMAP = umapA)))
  expect_identical(
    .decontXUMAPs(single, c("c1", "c3")),
    list(decontX_UMAP = umapA)
  )

  batched <- list(
    estimates = list(A = list(UMAP = umapA), B = list(UMAP = umapB)),
    batch = c("A", "B", "A")
  )
  cells <- c("c1", "c2", "c3")
  padded <- function(rows, values) {
    m <- matrix(NA_integer_, 3, 2,
      dimnames = list(cells, c("UMAP_1", "UMAP_2"))
    )
    m[rows, ] <- values
    m
  }
  expect_identical(
    .decontXUMAPs(batched, cells),
    list(
      decontX_A_UMAP = padded(c(1, 3), 1:4),
      decontX_B_UMAP = padded(2, 5:6)
    )
  )
})

test_that(desc = "decontX recovers simulated contamination (oracle)", {
  # Behavioral oracle: without this, a no-op decontX (returning its input
  # unchanged) would pass the suite. Uses the simulator's ground truth to
  # assert decontX actually estimates and removes contamination. Clusters
  # (z) are supplied so this isolates the EM/decontamination from the
  # stochastic clustering step.
  sim <- simulate_oracle()
  # iterLogLik = 1 records the log-likelihood every iteration so the
  # monotonicity check below sees the full EM sequence.
  res <- decontX(sim$observedCounts,
    z = sim$z, seed = 12345, iterLogLik = 1,
    verbose = FALSE
  )

  # (a) Estimated per-cell contamination tracks the true fraction.
  expect_gt(stats::cor(res$contamination, sim$contamination), 0.7)

  # (b) Decontaminated counts are closer to the native (true) counts than
  #     the observed counts are, and never exceed the observed counts.
  decont <- as.matrix(res$decontXcounts)
  observed <- as.matrix(sim$observedCounts)
  native <- as.matrix(sim$nativeCounts)
  expect_lt(sum(abs(decont - native)), sum(abs(observed - native)))
  expect_true(all(decont <= observed + 1e-6))

  # (c) EM log-likelihood is non-decreasing across iterations.
  ll <- res$estimates$all_cells$logLikelihood
  expect_gt(length(ll), 1)
  expect_true(all(diff(ll) >= -1e-6))
})

## DecontX
test_that("decontX on a matrix gives results the plots show", {
  s <- simulateContamination()
  res <- decontX(s$observedCounts, verbose = FALSE)
  expect_equal(dim(res$decontXcounts), dim(s$observedCounts))
  expect_true(all(res$contamination >= 0 & res$contamination <= 1))
  expect_equal(length(res$z), ncol(s$observedCounts))

  # The contamination plot colours each cell by its contamination
  b <- ggplot2::ggplot_build(plotDecontXContamination(res))
  expect_equal(b$plot$data$Contamination, unname(res$contamination))

  # Decontamination only removes counts, so no cluster has more cells
  # expressing a marker afterwards
  percent <- function(counts) {
    p <- plotDecontXMarkerPercentage(counts, z = res$z, markers = s$markers)
    ggplot2::ggplot_build(p)$plot$data$percent
  }
  before <- percent(s$observedCounts)
  after <- percent(res$decontXcounts)
  expect_true(all(after <= before))
  expect_true(any(after < before))

  b <- ggplot2::ggplot_build(
    plotDecontXMarkerExpression(s$observedCounts, s$markers[[1]], z = s$z)
  )
  expect_setequal(as.character(b$plot$data$Marker), s$markers[[1]])
  expect_equal(
    sort(b$plot$data$Expression),
    sort(as.vector(s$observedCounts[s$markers[[1]], ]))
  )

  # With background, eta is the background's gene profile
  bg <- s$observedCounts[, 1:5]
  colnames(bg) <- paste(colnames(bg), "_", sep = "")
  resBg <- decontX(s$observedCounts, background = bg, verbose = FALSE)
  expect_equal(dim(resBg$decontXcounts), dim(s$observedCounts))
  expect_true(all(resBg$contamination >= 0 & resBg$contamination <= 1))
  eta <- resBg$estimates$all_cells$eta
  expect_equal(eta, matrix(bgProfile(bg), nrow(eta), ncol(eta)),
    ignore_attr = TRUE
  )
})

test_that("marker plots build for matrix and SCE input without warnings", {
  s <- simulateContamination(seed = 12345)
  sce <- SingleCellExperiment::SingleCellExperiment(
    list(counts = s$observedCounts)
  )
  sce <- decontX(sce, z = s$z, seed = 12345, verbose = FALSE)
  groups <- list(A = 1, B = 2)

  # ggplot evaluates aesthetics lazily, so build each plot to exercise them
  expect_no_warning(ggplot2::ggplot_build(plotDecontXMarkerPercentage(
    s$observedCounts, s$markers, z = s$z, labelBars = TRUE
  )))
  expect_no_warning(ggplot2::ggplot_build(plotDecontXMarkerExpression(
    s$observedCounts, s$markers[[1]], groups, z = s$z
  )))
  sparse <- methods::as(s$observedCounts, "CsparseMatrix")
  expect_no_warning(ggplot2::ggplot_build(plotDecontXMarkerPercentage(
    sparse, s$markers, z = s$z
  )))
  expect_no_warning(ggplot2::ggplot_build(plotDecontXMarkerExpression(
    sparse, s$markers[[1]], groups, z = s$z
  )))

  p <- plotDecontXMarkerPercentage(sce, s$markers, groups,
    assayName = c("counts", "decontXcounts"), labelBars = TRUE
  )
  expect_no_warning(b <- ggplot2::ggplot_build(p))
  expect_length(b$data, 2)
  expect_length(unique(b$data[[1]]$fill), 2)

  p <- plotDecontXMarkerExpression(sce, s$markers[[1]], groups)
  expect_no_warning(b <- ggplot2::ggplot_build(p))
  expect_length(unique(b$data[[1]]$fill), 2)
  expect_equal(levels(b$layout$layout$Cell_Type), names(groups))
  expect_no_warning(ggplot2::ggplot_build(plotDecontXContamination(sce)))
})

test_that("decontX on an SCE stores results the plots read", {
  s <- simulateContamination()
  sce <- SingleCellExperiment::SingleCellExperiment(
    list(counts = s$observedCounts)
  )
  sce <- decontX(sce, verbose = FALSE)
  expect_true("decontXcounts" %in% SummarizedExperiment::assayNames(sce))
  expect_equal(dim(decontXcounts(sce)), dim(s$observedCounts))
  contam <- sce$decontX_contamination
  expect_equal(length(contam), ncol(sce))
  expect_true(all(contam >= 0 & contam <= 1))
  expect_equal(length(sce$decontX_clusters), ncol(sce))
  expect_equal(
    dim(SingleCellExperiment::reducedDim(sce, "decontX_UMAP")),
    c(ncol(sce), 2)
  )

  b <- ggplot2::ggplot_build(plotDecontXContamination(sce))
  expect_equal(b$plot$data$Contamination, unname(contam))

  # z as a colData column gives the same percentages as the labels
  percent <- function(z) {
    p <- plotDecontXMarkerPercentage(sce,
      z = z, markers = s$markers,
      assayName = "decontXcounts"
    )
    ggplot2::ggplot_build(p)$plot$data$percent
  }
  # With K = 3, "X1".."X3" sort like 1..3, so the rows line up
  sce$newz2 <- paste0("X", s$z)
  expect_identical(percent("newz2"), percent(s$z))

  # The marker plot defaults to the decontX clusters
  b <- ggplot2::ggplot_build(plotDecontXMarkerExpression(sce, s$markers[[1]]))
  expect_setequal(
    as.character(b$layout$layout$Cell_Type),
    levels(sce$decontX_clusters)
  )

  # With estimateDelta = FALSE, delta stays at its starting value
  sce <- decontX(sce, estimateDelta = FALSE, verbose = FALSE)
  delta <- S4Vectors::metadata(sce)$decontX$estimates$all_cells$delta
  expect_equal(delta, c(10, 10), ignore_attr = TRUE)

  # With background, eta is the background's gene profile
  bg <- sce[, 1:5]
  colnames(bg) <- paste(colnames(bg), "_", sep = "")
  sce <- decontX(sce, background = bg, verbose = FALSE)
  expect_true("decontXcounts" %in% SummarizedExperiment::assayNames(sce))
  profile <- bgProfile(SummarizedExperiment::assay(bg, "counts"))
  eta <- S4Vectors::metadata(sce)$decontX$estimates$all_cells$eta
  expect_equal(eta, matrix(profile, nrow(eta), ncol(eta)), ignore_attr = TRUE)
})


test_that(desc = "decontX checks delta, z and counts", {
  deconSim <- simulateContamination(K = 10, delta = c(1, 5))
  expect_error(
    decontX(
      x = deconSim$observedCounts,
      z = deconSim$z, delta = c(1, -1), verbose = FALSE
    ),
    regexp = "'delta' needs to be a numeric vector of length 2"
  )
  expect_error(
    decontX(
      x = deconSim$observedCounts,
      z = deconSim$z, delta = c(1, 1, 1), verbose = FALSE
    ),
    regexp = "'delta' needs to be a numeric vector of length 2"
  )
  expect_error(
    decontX(
      x = deconSim$observedCounts,
      z = c(deconSim$z, 1), verbose = FALSE
    ),
    paste0(
      "'z' must be of the same length as the number of cells in the",
      " 'counts' matrix."
    )
  )
  expect_error(
    decontX(
      x = deconSim$observedCounts,
      z = rep(1, ncol(
        deconSim$observedCounts
      )),
      verbose = FALSE
    ),
    "'z' has only one cluster. decontX needs at least two clusters"
  )
  countsNA <- deconSim$observedCounts
  countsNA[1, 1] <- NA
  expect_error(
    decontX(x = countsNA, z = deconSim$z, verbose = FALSE),
    "Missing value in 'counts' matrix."
  )
})

test_that("decontX stops when batch and bgBatch don't go together", {
  d <- refData()
  expect_error(
    decontX(d$sim$observedCounts,
      z = d$sim$z, background = d$bg,
      bgBatch = d$bgBatch, verbose = FALSE
    ),
    regexp = "'bgBatch' was given without 'batch'"
  )
  expect_error(
    decontX(d$sim$observedCounts,
      z = d$sim$z, batch = d$batch,
      background = d$bg, verbose = FALSE
    ),
    regexp = "Please also supply background batch"
  )
})

test_that("decontX stops when bgBatch is given without background", {
  d <- refData()
  expect_error(
    decontX(d$sim$observedCounts,
      z = d$sim$z, batch = d$batch,
      bgBatch = d$bgBatch, verbose = FALSE
    ),
    regexp = "'bgBatch' was given without 'background'; give both or neither"
  )
  sce <- SingleCellExperiment::SingleCellExperiment(
    list(counts = d$sim$observedCounts)
  )
  expect_error(
    decontX(sce,
      z = d$sim$z, batch = d$batch,
      bgBatch = d$bgBatch, verbose = FALSE
    ),
    regexp = "'bgBatch' was given without 'background'; give both or neither"
  )
})

test_that("decontX writes the same log text", {
  # Pins every line the core logs, trailing spaces included, so the log
  # stays the same when the logging code changes. convergence = 0 runs
  # every iteration, so the number of lines is fixed
  d <- refData()
  logfile <- withr::local_tempfile(fileext = ".log")
  res <- decontX(d$sim$observedCounts,
    z = d$sim$z, maxIter = 3, iterLogLik = 1, convergence = 0,
    logfile = logfile
  )
  expect_snapshot(writeLines(maskLog(readLines(logfile))))

  # A batched DelayedMatrix run, with a background barcode that is also in
  # the counts and a background gene that isn't
  bg <- rbind(d$bg, extra = 1)
  colnames(bg)[1] <- colnames(d$sim$observedCounts)[1]
  logfile <- withr::local_tempfile(fileext = ".log")
  res <- decontX(DelayedArray::DelayedArray(d$sim$observedCounts),
    z = d$sim$z, batch = d$batch, background = bg, bgBatch = d$bgBatch,
    maxIter = 3, iterLogLik = 1, convergence = 0, logfile = logfile
  )
  expect_snapshot(writeLines(maskLog(readLines(logfile))))

  # Without z, decontX also estimates the cell types
  logfile <- withr::local_tempfile(fileext = ".log")
  res <- decontX(d$sim$observedCounts,
    maxIter = 3, iterLogLik = 1, convergence = 0, logfile = logfile
  )
  expect_snapshot(writeLines(maskLog(readLines(logfile))))
})

test_that("decontX writes the same messages to the console", {
  d <- refData()
  messages <- testthat::capture_messages(
    decontX(d$sim$observedCounts,
      z = d$sim$z, maxIter = 3, iterLogLik = 1, convergence = 0
    )
  )
  expect_snapshot(writeLines(maskLog(sub("\n$", "", messages))))
})

test_that("decontX checks varGenes and dbscanEps even when z is given", {
  d <- refData()
  expect_error(
    decontX(d$sim$observedCounts, z = d$sim$z, varGenes = 1, verbose = FALSE),
    regexp = "'varGenes' must be an integer larger than 1"
  )
  expect_error(
    decontX(d$sim$observedCounts,
      z = d$sim$z, dbscanEps = -1,
      verbose = FALSE
    ),
    regexp = "'dbscanEps' needs to be non-negative"
  )
})

test_that("the bgBatch error and the cluster warning read correctly", {
  # The snapshots show the wording #74 fixed
  d <- refData()
  expect_snapshot(
    decontX(d$sim$observedCounts,
      z = d$sim$z, background = d$bg,
      bgBatch = d$bgBatch, verbose = FALSE
    ),
    error = TRUE
  )
  expect_snapshot(sim <- simulateContamination(K = 101, C = 10))
})

test_that("sparse and DelayedArray input give the dense-input results", {
  # The EM is deterministic, so two iterations are enough to compare
  d <- refData()
  dense <- decontX(d$sim$observedCounts,
    z = d$sim$z, maxIter = 2,
    verbose = FALSE
  )
  # Dense input returns a dgCMatrix, which #71 kept (A1-14)
  expect_s4_class(dense$decontXcounts, "dgCMatrix")

  sparseCounts <- methods::as(d$sim$observedCounts, "CsparseMatrix")
  sparse <- decontX(sparseCounts, z = d$sim$z, maxIter = 2, verbose = FALSE)
  expect_identical(sparse$contamination, dense$contamination)
  expect_identical(sparse$decontXcounts, dense$decontXcounts)

  delayedCounts <- DelayedArray::DelayedArray(sparseCounts)
  delayed <- decontX(delayedCounts,
    z = d$sim$z, maxIter = 2,
    verbose = FALSE
  )
  expect_s4_class(delayed$decontXcounts, "DelayedMatrix")
  expect_identical(delayed$contamination, dense$contamination)
  expect_identical(
    methods::as(delayed$decontXcounts, "dgCMatrix"),
    dense$decontXcounts
  )
})

test_that("other sparse classes give the dense-input results as a dgCMatrix", {
  d <- refData()
  dense <- decontX(d$sim$observedCounts,
    z = d$sim$z, maxIter = 2,
    verbose = FALSE
  )
  for (cls in c("RsparseMatrix", "TsparseMatrix", "SVT_SparseMatrix")) {
    counts <- methods::as(d$sim$observedCounts, cls)
    expect_no_warning(
      res <- decontX(counts, z = d$sim$z, maxIter = 2, verbose = FALSE)
    )
    expect_identical(res$contamination, dense$contamination)
    expect_identical(res$decontXcounts, dense$decontXcounts)
  }
})

test_that("a DelayedMatrix over a dense matrix gives a DelayedMatrix", {
  d <- refData()
  dense <- decontX(d$sim$observedCounts,
    z = d$sim$z, maxIter = 2,
    verbose = FALSE
  )
  delayed <- decontX(DelayedArray::DelayedArray(d$sim$observedCounts),
    z = d$sim$z, maxIter = 2,
    verbose = FALSE
  )
  expect_s4_class(delayed$decontXcounts, "DelayedMatrix")
  expect_s4_class(DelayedArray::seed(delayed$decontXcounts), "dgCMatrix")
  expect_identical(delayed$contamination, dense$contamination)
  expect_identical(
    methods::as(delayed$decontXcounts, "dgCMatrix"),
    dense$decontXcounts
  )
})

test_that("decontX names unsupported input classes", {
  d <- refData()
  x <- d$sim$observedCounts
  expect_error(
    decontX(as.data.frame(x), z = d$sim$z, verbose = FALSE),
    regexp = "'x' must be a SingleCellExperiment, .*; it is a data\\.frame\\."
  )
  se <- SummarizedExperiment::SummarizedExperiment(list(counts = x))
  expect_error(
    decontX(se, z = d$sim$z, verbose = FALSE),
    regexp = paste0(
      "it is a SummarizedExperiment\\. .*",
      "as\\(x, \"SingleCellExperiment\"\\)"
    )
  )
  sceBg <- SingleCellExperiment::SingleCellExperiment(list(counts = d$bg))
  expect_error(
    decontX(x, z = d$sim$z, background = sceBg, verbose = FALSE),
    regexp = paste0(
      "'background' must be .* when 'x' is a matrix; it is a ",
      "SingleCellExperiment\\."
    )
  )
  # Counts stored as text, with the dimnames kept
  textX <- x
  storage.mode(textX) <- "character"
  textBg <- d$bg
  storage.mode(textBg) <- "character"
  expect_error(
    decontX(textX, z = d$sim$z, verbose = FALSE),
    regexp = "'x' must be .*; it is a character matrix\\."
  )
  expect_error(
    decontX(DelayedArray::DelayedArray(textX), z = d$sim$z, verbose = FALSE),
    regexp = "'x' must be .*; it is a character DelayedMatrix\\."
  )
  # Presence/absence isn't counts, whatever the storage
  expect_error(
    decontX(methods::as(x > 0, "CsparseMatrix"),
      z = d$sim$z, verbose = FALSE
    ),
    regexp = "'x' must be .*; it is a lgCMatrix\\."
  )
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = x))
  textSce <- SingleCellExperiment::SingleCellExperiment(list(counts = textX))
  expect_error(
    decontX(textSce, z = d$sim$z, verbose = FALSE),
    regexp = "The 'x' assay must be .*; it is a character matrix\\."
  )
  textBgSce <- SingleCellExperiment::SingleCellExperiment(
    list(counts = textBg)
  )
  expect_error(
    decontX(sce, z = d$sim$z, background = textBgSce, verbose = FALSE),
    regexp = "The 'background' assay must be .*; it is a character matrix\\."
  )
  expect_error(
    decontX(sce, z = d$sim$z, background = d$bg, verbose = FALSE),
    regexp = paste0(
      "'background' must be a SingleCellExperiment when 'x' is one; ",
      "it is a matrix\\."
    )
  )
})

test_that("decontX stops on batches too small to initialize", {
  d <- refData()
  x <- d$sim$observedCounts
  # Without z, nothing else catches a one-cell batch
  expect_error(
    decontX(x, batch = c(rep("A", 99), "B"), verbose = FALSE),
    regexp = "Each batch needs at least 3 cells; batch 'B' has 1\\."
  )
  # With background, one cluster is allowed, so the size check stops it
  expect_error(
    decontX(x,
      z = d$sim$z, batch = c(rep("A", 99), "B"), background = d$bg,
      bgBatch = c(rep("A", 39), "B"), verbose = FALSE
    ),
    regexp = "Each batch needs at least 3 cells; batch 'B' has 1\\."
  )
  expect_error(
    decontX(x, batch = c(rep("A", 98), "B", "B"), verbose = FALSE),
    regexp = "batch 'B' has 2\\."
  )
  expect_error(
    decontX(x[, 1:2], z = c(1, 2), verbose = FALSE),
    regexp = "'x' needs at least 3 cells; it has 2\\."
  )

  # scrapper caps its neighbor count for tiny batches, with a warning
  res <- suppressWarnings(
    decontX(x,
      z = c(d$sim$z[1:97], 1, 2, 2), batch = c(rep("A", 97), rep("B", 3)),
      maxIter = 2, verbose = FALSE
    )
  )
  expect_true(all(is.finite(res$contamination)))
  res <- suppressWarnings(
    decontX(x[, 1:3], z = c(1, 2, 2), maxIter = 2, verbose = FALSE)
  )
  expect_true(all(is.finite(res$contamination)))
})

test_that("legacyInit needs 15 cells in each batch", {
  skip_if_not_installed("scater")
  d <- refData()
  x <- d$sim$observedCounts
  z <- rep(1:2, length.out = 15)
  # The legacy path's upstream deprecation warnings are suppressed
  expect_error(
    suppressWarnings(
      decontX(x[, 1:14], z = z[1:14], legacyInit = TRUE, verbose = FALSE)
    ),
    regexp = paste0(
      "'x' needs at least 15 cells with 'legacyInit = TRUE'; ",
      "it has 14\\."
    )
  )
  res <- suppressWarnings(
    decontX(x[, 1:15],
      z = z, legacyInit = TRUE, maxIter = 2,
      verbose = FALSE
    )
  )
  expect_true(all(is.finite(res$contamination)))
})

test_that("decontXInitialize stops on cluster labels below 1", {
  # Not symmetric, so it becomes a dgCMatrix
  counts <- methods::as(matrix(as.numeric(1:9), 3, 3), "CsparseMatrix")
  for (z in list(c(0L, 1L, 2L), c(NA, 1L, 2L), c(-1L, 1L, 2L))) {
    expect_error(
      decontXInitialize(counts, theta = rep(0.5, 3), z = z, pseudocount = 0),
      regexp = "The entries in 'z' need to be"
    )
  }
})

test_that("decontX checks that z and batch have one value per cell", {
  d <- refData()
  x <- d$sim$observedCounts
  zNA <- d$sim$z
  zNA[1:5] <- NA
  batchNA <- d$batch
  batchNA[1] <- NA
  lengthMessage <- "must be of the same length as the number of cells"
  expect_error(decontX(x, z = d$sim$z[1:90], verbose = FALSE),
    regexp = paste("'z'", lengthMessage)
  )
  expect_error(decontX(x, z = "cluster", verbose = FALSE),
    regexp = paste("'z'", lengthMessage)
  )
  expect_error(decontX(x, z = zNA, verbose = FALSE),
    regexp = "'z' must not contain NA"
  )
  expect_error(
    decontX(x, z = d$sim$z, batch = d$batch[1:90], verbose = FALSE),
    regexp = paste("'batch'", lengthMessage)
  )
  expect_error(decontX(x, z = d$sim$z, batch = batchNA, verbose = FALSE),
    regexp = "'batch' must not contain NA"
  )
})

test_that("decontX stops on an empty batch name (#102)", {
  d <- refData()
  x <- d$sim$observedCounts
  # [[""]] never matches, so the batch's estimates couldn't be read back
  batchEmpty <- replace(d$batch, d$batch == "A", "")
  for (batch in list(batchEmpty, factor(batchEmpty))) {
    expect_error(decontX(x, z = d$sim$z, batch = batch, verbose = FALSE),
      regexp = "'batch' must not contain empty strings"
    )
    expect_error(decontX(x, batch = batch, verbose = FALSE),
      regexp = "'batch' must not contain empty strings"
    )
  }
})

test_that("decontX checks its scalar arguments", {
  d <- refData()
  x <- d$sim$observedCounts
  z <- d$sim$z
  checks <- list(
    list(args = list(maxIter = 0), regexp = "'maxIter' must be a single"),
    list(args = list(maxIter = "a"), regexp = "'maxIter' must be a single"),
    list(args = list(maxIter = 2.5), regexp = "'maxIter' must be a single"),
    list(args = list(iterLogLik = 0), regexp = "'iterLogLik' must be a single"),
    list(
      args = list(convergence = "x"),
      regexp = "'convergence' must be a single non-negative number"
    ),
    list(
      args = list(convergence = -1),
      regexp = "'convergence' must be a single non-negative number"
    ),
    list(
      args = list(estimateDelta = "yes"),
      regexp = "'estimateDelta' must be TRUE or FALSE"
    ),
    list(
      args = list(estimateDelta = NA),
      regexp = "'estimateDelta' must be TRUE or FALSE"
    ),
    list(
      args = list(legacyInit = "yes"),
      regexp = "'legacyInit' must be TRUE or FALSE"
    ),
    list(
      args = list(verbose = "no"),
      regexp = "'verbose' must be TRUE or FALSE"
    ),
    list(args = list(seed = "a"), regexp = "'seed' must be NULL or a single"),
    list(
      args = list(seed = c(1, 2)),
      regexp = "'seed' must be NULL or a single"
    ),
    list(
      args = list(varGenes = NA),
      regexp = "'varGenes' must be an integer larger than 1"
    ),
    list(
      args = list(varGenes = "a"),
      regexp = "'varGenes' must be an integer larger than 1"
    ),
    list(
      args = list(dbscanEps = NA),
      regexp = "'dbscanEps' needs to be non-negative"
    ),
    list(
      args = list(dbscanEps = "a"),
      regexp = "'dbscanEps' needs to be non-negative"
    ),
    list(
      args = list(logfile = 1),
      regexp = "'logfile' must be NULL or a single character string"
    ),
    list(
      args = list(embedding = "yes"),
      regexp = "'embedding' must be TRUE or FALSE"
    ),
    list(
      args = list(embedding = NA),
      regexp = "'embedding' must be TRUE or FALSE"
    )
  )
  for (bad in list(0, 1.5, NA, "2", c(1, 2), Inf, -1)) {
    checks <- c(checks, list(list(
      args = list(numThreads = bad),
      regexp = "'numThreads' must be a single whole number of at least 1"
    )))
  }
  for (check in checks) {
    call <- utils::modifyList(list(x = x, z = z, verbose = FALSE), check$args)
    expect_error(do.call(decontX, call), regexp = check$regexp)
  }
})

test_that("decontX requires two finite, positive delta values", {
  d <- refData()
  deltaMessage <- "'delta' needs to be a numeric vector of length 2"
  for (delta in list(c(0, 0), c(10, NA), c(10, Inf))) {
    expect_error(
      decontX(d$sim$observedCounts,
        z = d$sim$z, delta = delta,
        verbose = FALSE
      ),
      regexp = deltaMessage
    )
  }
})

test_that("decontX stops on negative counts", {
  d <- refData()
  counts <- d$sim$observedCounts
  counts[1, 1:3] <- -5L
  expect_error(decontX(counts, z = d$sim$z, verbose = FALSE),
    regexp = "Negative value in 'counts' matrix"
  )
})

test_that("decontX checks assay names before running", {
  d <- refData()
  sce <- SingleCellExperiment::SingleCellExperiment(
    list(counts = d$sim$observedCounts)
  )
  bg <- SingleCellExperiment::SingleCellExperiment(list(counts = d$bg))
  expect_error(decontX(sce, assayName = "raw", z = d$sim$z, verbose = FALSE),
    regexp = "'assayName' \\(\"raw\"\\) is not an assay of 'x'"
  )
  expect_error(
    decontX(sce,
      z = d$sim$z, background = bg, bgAssayName = "raw",
      verbose = FALSE
    ),
    regexp = "'bgAssayName' \\(\"raw\"\\) is not an assay of 'background'"
  )
})

test_that("the zero-count error doesn't suggest supplying z", {
  # Initialization runs by default even when z is supplied
  d <- refData()
  counts <- d$sim$observedCounts
  counts[, 1] <- 0
  err <- expect_error(decontX(counts, z = d$sim$z, verbose = FALSE),
    regexp = "at least one count"
  )
  expect_false(grepl("'z'", conditionMessage(err)))
})

## embedding and numThreads (#79)

# A decontX result without the parts the embedding and the new arguments
# change, to compare runs with and without the embedding
withoutEmbedding <- function(res) {
  res$runParams$embedding <- NULL
  res$runParams$numThreads <- NULL
  res$estimates <- lapply(res$estimates, function(e) {
    e$UMAP <- NULL
    e
  })
  res
}

test_that("embedding = FALSE gives the same decontamination without a UMAP", {
  d <- refData()
  x <- d$sim$observedCounts
  for (batch in list(NULL, d$batch)) {
    with <- decontX(x, z = d$sim$z, batch = batch, verbose = FALSE)
    without <- decontX(x,
      z = d$sim$z, batch = batch, embedding = FALSE,
      verbose = FALSE
    )
    for (e in without$estimates) {
      expect_true("UMAP" %in% names(e))
      expect_null(e$UMAP)
    }
    expect_false(without$runParams$embedding)
    expect_identical(withoutEmbedding(without), withoutEmbedding(with))
  }
})

test_that("embedding = FALSE with legacyInit gives the same decontamination", {
  skip_if_not_installed("scater")
  d <- refData()
  x <- d$sim$observedCounts
  with <- suppressWarnings(
    decontX(x, z = d$sim$z, legacyInit = TRUE, verbose = FALSE)
  )
  without <- decontX(x,
    z = d$sim$z, legacyInit = TRUE, embedding = FALSE,
    verbose = FALSE
  )
  expect_null(without$estimates$all_cells$UMAP)
  expect_identical(withoutEmbedding(without), withoutEmbedding(with))
})

test_that("embedding = FALSE writes no UMAP and keeps an earlier one", {
  d <- refData()
  sce <- SingleCellExperiment::SingleCellExperiment(
    list(counts = d$sim$observedCounts)
  )
  res <- decontX(sce, z = d$sim$z, embedding = FALSE, verbose = FALSE)
  expect_identical(SingleCellExperiment::reducedDimNames(res), character())
  expect_identical(
    res$decontX_contamination,
    decontX(sce, z = d$sim$z, verbose = FALSE)$decontX_contamination
  )

  # A UMAP from an earlier run stays, and the batched run adds none
  earlier <- decontX(sce, z = d$sim$z, verbose = FALSE)
  res <- decontX(earlier,
    z = d$sim$z, batch = d$batch, embedding = FALSE,
    verbose = FALSE
  )
  expect_identical(SingleCellExperiment::reducedDimNames(res), "decontX_UMAP")
  expect_identical(
    SingleCellExperiment::reducedDim(res, "decontX_UMAP"),
    SingleCellExperiment::reducedDim(earlier, "decontX_UMAP")
  )
})

test_that("embedding = FALSE stores the same SCE results without a UMAP", {
  d <- refData()
  sce <- SingleCellExperiment::SingleCellExperiment(
    list(counts = d$sim$observedCounts)
  )
  for (batch in list(NULL, d$batch)) {
    with <- decontX(sce, z = d$sim$z, batch = batch, verbose = FALSE)
    without <- decontX(sce,
      z = d$sim$z, batch = batch, embedding = FALSE,
      verbose = FALSE
    )
    expect_identical(decontXcounts(without), decontXcounts(with))
    expect_identical(
      SummarizedExperiment::colData(without),
      SummarizedExperiment::colData(with)
    )
    expect_identical(
      withoutEmbedding(S4Vectors::metadata(without)$decontX),
      withoutEmbedding(S4Vectors::metadata(with)$decontX)
    )
    expect_length(SingleCellExperiment::reducedDimNames(without), 0)
  }
})

test_that("embedding = FALSE needs z", {
  d <- refData()
  expect_error(
    decontX(d$sim$observedCounts, embedding = FALSE, verbose = FALSE),
    regexp = "'embedding = FALSE' needs 'z'"
  )
})

test_that("embedding = FALSE still stops on cells without counts", {
  d <- refData()
  counts <- d$sim$observedCounts
  counts[, 1] <- 0
  expect_error(
    decontX(counts, z = d$sim$z, embedding = FALSE, verbose = FALSE),
    regexp = "All cells must have at least one count"
  )
})

test_that("embedding = FALSE needs 2 cells in each batch", {
  d <- refData()
  x <- d$sim$observedCounts
  # The 3-cell (15 with legacyInit) minimum is for the PCA and UMAP only
  res <- decontX(x[, 1:2],
    z = c(1, 2), embedding = FALSE, maxIter = 2, verbose = FALSE
  )
  expect_true(all(is.finite(res$contamination)))
  res <- decontX(x,
    z = c(d$sim$z[1:98], 1, 2), batch = c(rep("A", 98), "B", "B"),
    embedding = FALSE, maxIter = 2, verbose = FALSE
  )
  expect_true(all(is.finite(res$contamination)))
  # One cell can't estimate the contamination prior, even with background,
  # where one cluster is allowed
  expect_error(
    decontX(x,
      z = d$sim$z, batch = c(rep("A", 99), "B"), background = d$bg,
      bgBatch = c(rep("A", 39), "B"), embedding = FALSE, verbose = FALSE
    ),
    regexp = paste0(
      "Each batch needs at least 2 cells to estimate the contamination ",
      "prior; batch 'B' has 1\\."
    )
  )
  expect_error(
    decontX(x[, 1, drop = FALSE],
      z = 1, background = d$bg, embedding = FALSE, verbose = FALSE
    ),
    regexp = "'x' needs at least 2 cells to estimate the contamination prior"
  )
})

test_that("embedding = FALSE has no 15-cell minimum with legacyInit", {
  skip_if_not_installed("scater")
  d <- refData()
  res <- decontX(d$sim$observedCounts[, 1:10],
    z = rep(1:2, 5), legacyInit = TRUE, embedding = FALSE, maxIter = 2,
    verbose = FALSE
  )
  expect_true(all(is.finite(res$contamination)))
})

test_that("embedding = FALSE leaves the UMAP out of the log", {
  d <- refData()
  logfile <- withr::local_tempfile(fileext = ".log")
  res <- decontX(d$sim$observedCounts,
    z = d$sim$z, embedding = FALSE, maxIter = 3, iterLogLik = 1,
    convergence = 0, logfile = logfile
  )
  expect_snapshot(writeLines(maskLog(readLines(logfile))))
})

test_that("runParams records embedding and numThreads", {
  d <- refData()
  res <- decontX(d$sim$observedCounts, z = d$sim$z, verbose = FALSE)
  expect_true(res$runParams$embedding)
  expect_identical(res$runParams$numThreads, 1)
  res <- decontX(d$sim$observedCounts,
    z = d$sim$z, embedding = FALSE, numThreads = 2L, verbose = FALSE
  )
  expect_false(res$runParams$embedding)
  expect_identical(res$runParams$numThreads, 2L)
})

test_that("numThreads reaches every threaded initialization call", {
  d <- refData()
  seen <- list()
  record <- function(name, value) {
    seen[[name]] <<- value
  }
  # Each mock is the scrapper function with a first line that records
  # num.threads, so it keeps the function's formals: decontX reads them to
  # choose runUmap's seed arguments
  recording <- function(fn, name) {
    body(fn) <- bquote({
      .(record)(.(name), num.threads)
      .(body(fn))
    })
    fn
  }
  local_mocked_bindings(
    modelGeneVariances = recording(
      scrapper::modelGeneVariances, "modelGeneVariances"
    ),
    runPca = recording(scrapper::runPca, "runPca"),
    runUmap = recording(scrapper::runUmap, "runUmap"),
    .package = "scrapper"
  )
  res <- decontX(d$sim$observedCounts,
    z = d$sim$z, numThreads = 3, maxIter = 2, verbose = FALSE
  )
  expect_identical(
    seen,
    list(modelGeneVariances = 3, runPca = 3, runUmap = 3)
  )
})

test_that("numThreads reaches the legacy UMAP", {
  skip_if_not_installed("scater")
  d <- refData()
  threads <- NULL
  local_mocked_bindings(
    calculateUMAP = function(x, ..., n_threads) {
      threads <<- n_threads
      matrix(0, ncol(x), 2)
    },
    .package = "scater"
  )
  res <- suppressWarnings(decontX(d$sim$observedCounts,
    z = d$sim$z, legacyInit = TRUE, numThreads = 3, maxIter = 2,
    verbose = FALSE
  ))
  expect_identical(threads, 3)
})

test_that("numThreads > 1 gives the same legacy results with z", {
  skip_if_not_installed("scater")
  d <- refData()
  x <- d$sim$observedCounts
  one <- suppressWarnings(
    decontX(x, z = d$sim$z, legacyInit = TRUE, verbose = FALSE)
  )
  two <- suppressWarnings(decontX(x,
    z = d$sim$z, legacyInit = TRUE, numThreads = 2,
    verbose = FALSE
  ))
  expect_identical(two$decontXcounts, one$decontXcounts)
  expect_identical(two$contamination, one$contamination)
  expect_identical(
    two$estimates$all_cells$UMAP,
    one$estimates$all_cells$UMAP
  )
})

test_that("numThreads > 1 gives the same decontamination with z", {
  d <- refData()
  x <- d$sim$observedCounts
  one <- decontX(x, z = d$sim$z, verbose = FALSE)
  two <- decontX(x, z = d$sim$z, numThreads = 2, verbose = FALSE)
  expect_identical(two$decontXcounts, one$decontXcounts)
  expect_identical(two$contamination, one$contamination)
  expect_identical(dim(two$estimates$all_cells$UMAP), c(ncol(x), 2L))
})

test_that("the contamination plot uses a reducedDim for every cell", {
  d <- refData()
  sce <- SingleCellExperiment::SingleCellExperiment(
    list(counts = d$sim$observedCounts)
  )
  coords <- withr::with_seed(1, matrix(rnorm(2 * ncol(sce)), ncol = 2))
  colnames(coords) <- c("PC1", "PC2")
  SingleCellExperiment::reducedDim(sce, "PCA") <- coords
  sce <- decontX(sce,
    z = d$sim$z, batch = d$batch, embedding = FALSE,
    verbose = FALSE
  )

  # Every cell, with the reducedDim's column names as the axis titles
  b <- ggplot2::ggplot_build(
    plotDecontXContamination(sce, reducedDimName = "PCA")
  )
  expect_equal(b$data[[1]]$x, unname(coords[, 1]))
  expect_equal(b$data[[1]]$y, unname(coords[, 2]))
  expect_equal(b$plot$data$Contamination, sce$decontX_contamination)
  expect_identical(b$plot$labels$x, "PC1")
  expect_identical(b$plot$labels$y, "PC2")

  # One batch's cells
  b <- ggplot2::ggplot_build(
    plotDecontXContamination(sce, reducedDimName = "PCA", batch = "B")
  )
  inB <- d$batch == "B"
  expect_equal(b$data[[1]]$x, unname(coords[inB, 1]))
  expect_equal(b$plot$data$Contamination, sce$decontX_contamination[inB])

  # Unnamed columns get generic axis titles
  SingleCellExperiment::reducedDim(sce, "PCA") <- unname(coords)
  b <- ggplot2::ggplot_build(
    plotDecontXContamination(sce, reducedDimName = "PCA")
  )
  expect_identical(b$plot$labels$x, "Dim_1")
  expect_identical(b$plot$labels$y, "Dim_2")

  expect_error(
    plotDecontXContamination(sce, reducedDimName = "TSNE"),
    regexp = "'TSNE' is not a reducedDim of 'x'\\. Available: 'PCA'"
  )
  for (bad in list(1, NA_character_, c("PCA", "PCA"))) {
    expect_error(
      plotDecontXContamination(sce, reducedDimName = bad),
      regexp = "'reducedDimName' must be a single reducedDim name"
    )
  }
  expect_error(
    plotDecontXContamination(sce, batch = c("A", "B")),
    regexp = "'batch' must be NULL or a single batch name"
  )

  # A subset of the cells: every cell works, a batch needs the run's cells
  b <- ggplot2::ggplot_build(
    plotDecontXContamination(sce[, 1:10], reducedDimName = "PCA")
  )
  expect_equal(b$data[[1]]$x, unname(coords[1:10, 1]))
  expect_error(
    plotDecontXContamination(sce[, 1:10], reducedDimName = "PCA", batch = "A"),
    regexp = "'x' has 10 cells and the run had 100"
  )

  # The decontX UMAP is missing from the SCE's estimates too
  expect_error(plotDecontXContamination(sce),
    regexp = "Batch 'A' has no decontX UMAP"
  )

  noContamination <- sce
  noContamination$decontX_contamination <- NULL
  expect_error(
    plotDecontXContamination(noContamination, reducedDimName = "PCA"),
    regexp = "'x' has no 'decontX_contamination' column"
  )

  SingleCellExperiment::reducedDim(sce, "PCA") <- coords[, 1, drop = FALSE]
  expect_error(
    plotDecontXContamination(sce, reducedDimName = "PCA"),
    regexp = "reducedDim 'PCA' needs at least two columns"
  )
  SingleCellExperiment::reducedDims(sce) <- list()
  expect_error(
    plotDecontXContamination(sce, reducedDimName = "PCA"),
    regexp = "'x' has no reducedDims"
  )
  expect_error(
    plotDecontXContamination(sce, reducedDimName = "PCA", batch = "C"),
    regexp = "'C' is not one of the batches in 'x'"
  )
})

test_that("the contamination plot explains a missing decontX UMAP", {
  d <- refData()
  res <- decontX(d$sim$observedCounts,
    z = d$sim$z, embedding = FALSE, verbose = FALSE
  )
  expect_error(plotDecontXContamination(res),
    regexp = "no decontX UMAP.*embedding = FALSE.*reducedDimName"
  )
  expect_error(plotDecontXContamination(res, reducedDimName = "PCA"),
    regexp = "'reducedDimName' needs 'x' to be a SingleCellExperiment"
  )
})

test_that("decontX still accepts integer, factor and Inf arguments", {
  d <- refData()
  counts <- d$sim$observedCounts
  expect_no_error(.checkCellValues(factor(d$sim$z), "z", ncol(counts)))
  expect_no_error(.checkCellValues(factor(d$batch), "batch", ncol(counts)))
  expect_no_error(.decontXParams(
    maxIter = 500L,
    delta = c(10L, 10L),
    estimateDelta = TRUE,
    convergence = 0L,
    iterLogLik = 10L,
    varGenes = 5000L,
    dbscanEps = 0,
    seed = 1L,
    legacyInit = FALSE,
    logfile = NULL,
    verbose = FALSE,
    embedding = TRUE,
    numThreads = 2L
  ))
  expect_no_error(.checkCellValues(as.character(d$sim$z), "z", ncol(counts)))
  expect_no_error(.decontXParams(
    maxIter = Inf,
    delta = c(10, 10),
    estimateDelta = FALSE,
    convergence = 0.001,
    iterLogLik = Inf,
    varGenes = Inf,
    dbscanEps = NULL,
    seed = NULL,
    legacyInit = TRUE,
    logfile = tempfile(),
    verbose = TRUE,
    embedding = FALSE,
    numThreads = 1
  ))
})

test_that("decontX rejects settings that never stop or can't seed", {
  d <- refData()
  x <- d$sim$observedCounts
  expect_error(
    decontX(x, z = d$sim$z, maxIter = Inf, convergence = 0, verbose = FALSE),
    regexp = "'convergence' must be above 0 when 'maxIter' is Inf"
  )
  expect_error(decontX(x, z = d$sim$z, seed = 1e10, verbose = FALSE),
    regexp = "'seed' must be NULL or a single number in R's integer range"
  )
})

test_that("decontX checks background counts for NA and negative values", {
  d <- refData()
  bgNegative <- d$bg
  bgNegative[1, 1] <- -1L
  bgMissing <- d$bg
  bgMissing[1, 1] <- NA
  expect_error(
    decontX(d$sim$observedCounts,
      z = d$sim$z, background = bgNegative,
      verbose = FALSE
    ),
    regexp = "Negative value in 'background' matrix"
  )
  sce <- SingleCellExperiment::SingleCellExperiment(
    list(counts = d$sim$observedCounts)
  )
  bgSce <- SingleCellExperiment::SingleCellExperiment(
    list(counts = bgMissing)
  )
  expect_error(
    decontX(sce, z = d$sim$z, background = bgSce, verbose = FALSE),
    regexp = "Missing value in 'background' matrix"
  )
})

test_that("decontX takes assays by index and names the assay argument", {
  d <- refData()
  sce <- SingleCellExperiment::SingleCellExperiment(
    list(counts = d$sim$observedCounts)
  )
  byName <- decontX(sce, z = d$sim$z, maxIter = 2, verbose = FALSE)
  byIndex <- decontX(sce,
    assayName = 1, z = d$sim$z, maxIter = 2,
    verbose = FALSE
  )
  expect_identical(byIndex$decontX_contamination, byName$decontX_contamination)
  expect_error(decontX(sce, assayName = 2, z = d$sim$z, verbose = FALSE),
    regexp = "'assayName' \\(2\\) is not an assay of 'x'"
  )

  # bgAssayName defaults to assayName, so the message names assayName
  sce <- SingleCellExperiment::SingleCellExperiment(
    list(raw = d$sim$observedCounts)
  )
  bg <- SingleCellExperiment::SingleCellExperiment(list(counts = d$bg))
  expect_error(
    decontX(sce,
      assayName = "raw", z = d$sim$z, background = bg,
      verbose = FALSE
    ),
    regexp = "'assayName' \\(\"raw\"\\) is not an assay of 'background'"
  )
})

test_that("decontX matches background rows to the counts by gene name", {
  d <- refData()
  x <- d$sim$observedCounts
  run <- function(bg, counts = x, verbose = FALSE) {
    decontX(counts,
      z = d$sim$z, background = bg, maxIter = 5,
      verbose = verbose
    )
  }
  matched <- run(d$bg)$contamination

  reversed <- d$bg[rev(seq_len(nrow(d$bg))), ]
  expect_identical(run(reversed)$contamination, matched)
  extra <- rbind(d$bg, Extra_1 = 1L, Extra_2 = 2L)
  expect_identical(run(extra)$contamination, matched)
  messages <- capture_messages(run(extra, verbose = TRUE))
  expect_match(messages,
    regexp = paste0(
      "Dropped genes in 'background' that aren't in 'x': 2 of 62 ",
      "\\(Extra_1, Extra_2\\)"
    ),
    all = FALSE
  )

  # Without row names on either, rows are matched by position
  unnamed <- d$bg
  rownames(unnamed) <- NULL
  expect_identical(run(unnamed)$contamination, matched)

  # Duplicated names that are identical in both are used as they are
  dupX <- x
  rownames(dupX)[2] <- rownames(dupX)[1]
  dupBg <- d$bg
  rownames(dupBg) <- rownames(dupX)
  expect_identical(run(dupBg, counts = dupX)$contamination, matched)

  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = x))
  for (bg in list(reversed, extra)) {
    bgSce <- SingleCellExperiment::SingleCellExperiment(list(counts = bg))
    res <- decontX(sce,
      z = d$sim$z, background = bgSce, maxIter = 5,
      verbose = FALSE
    )
    expect_identical(res$decontX_contamination, matched)
  }
})

test_that("decontX stops when background rows can't be matched", {
  d <- refData()
  x <- d$sim$observedCounts
  run <- function(bg, counts = x) {
    decontX(counts, z = d$sim$z, background = bg, verbose = FALSE)
  }
  expect_error(run(d$bg[-c(1, 2), ]),
    regexp = paste0(
      "'background' lacks 2 of the 60 genes in 'x': ",
      rownames(x)[1], ", ", rownames(x)[2]
    )
  )
  noOverlap <- d$bg
  rownames(noOverlap) <- paste0("other_", seq_len(nrow(noOverlap)))
  expect_error(run(noOverlap),
    regexp = "'background' lacks 60 of the 60 genes in 'x'"
  )
  unnamed <- d$bg[-1, ]
  rownames(unnamed) <- NULL
  expect_error(run(unnamed),
    regexp = "'background' and 'x' have a different number of genes"
  )
  dupX <- x
  rownames(dupX)[2] <- rownames(dupX)[1]
  expect_error(run(d$bg, counts = dupX),
    regexp = "duplicated row names, so their rows can't be matched by name"
  )
})

test_that("decontX checks bgBatch against batch and the background", {
  d <- refData()
  x <- d$sim$observedCounts
  run <- function(bgBatch, bg = d$bg) {
    decontX(x,
      z = d$sim$z, batch = d$batch, background = bg,
      bgBatch = bgBatch, verbose = FALSE
    )
  }
  expect_error(run(c("A", "B")),
    regexp = "'bgBatch' must have one value per column of 'background'"
  )
  expect_error(run(replace(d$bgBatch, 1, NA)),
    regexp = "'bgBatch' must not contain NA"
  )
  expect_error(run(replace(d$bgBatch, 1, "")),
    regexp = "'bgBatch' must not contain empty strings"
  )
  expect_error(run(rep(c("C", "A", "B"), length.out = 40)),
    regexp = "'bgBatch' has values that aren't in 'batch': C"
  )
  expect_error(run(rep("A", 40)),
    regexp = "Every batch needs droplets in 'background'; none for: B"
  )

  # A batch can also lose its droplets to barcodes that are in x
  bg <- d$bg
  colnames(bg)[21:40] <- colnames(x)[1:20]
  expect_error(run(d$bgBatch, bg = bg),
    regexp = "Every batch needs droplets in 'background'; none for: B"
  )
})

test_that("decontX stops when no background droplets are left", {
  d <- refData()
  x <- d$sim$observedCounts
  bg <- d$bg
  colnames(bg) <- colnames(x)[seq_len(ncol(bg))]
  expect_error(
    decontX(x, z = d$sim$z, background = bg, verbose = FALSE),
    regexp = "No droplets are left in 'background' after removing"
  )
  expect_error(
    decontX(x, z = d$sim$z, background = d$bg[, 0], verbose = FALSE),
    regexp = "'background' has no droplets"
  )
  expect_error(
    decontX(x, z = d$sim$z, background = d$bg[, 1], verbose = FALSE),
    regexp = "'background' must be a matrix"
  )
})

test_that("decontX handles one droplet per batch and duplicated extra genes", {
  d <- refData()
  x <- d$sim$observedCounts
  expect_no_error(
    decontX(x,
      z = d$sim$z, batch = d$batch, background = d$bg,
      bgBatch = c(rep("A", 39), "B"), maxIter = 5, verbose = FALSE
    )
  )

  # Duplicated names only matter among the genes that are matched to x
  run <- function(bg) {
    decontX(x, z = d$sim$z, background = bg, maxIter = 5, verbose = FALSE)
  }
  extra <- rbind(d$bg, Extra = 1L, Extra = 2L)
  expect_identical(run(extra)$contamination, run(d$bg)$contamination)
})

test_that("decontX allows one cluster when the batch has background", {
  d <- refData()
  x <- d$sim$observedCounts
  res <- decontX(x,
    z = rep(1, ncol(x)), background = d$bg, maxIter = 5,
    verbose = FALSE
  )
  expect_length(unique(res$z), 1)
  expect_true(all(is.finite(res$contamination)))
  expect_true(all(res$contamination >= 0 & res$contamination <= 1))

  # One batch down to one cluster, the other with several
  z <- ifelse(d$batch == "A", 1, d$sim$z)
  res <- decontX(x,
    z = z, batch = d$batch, background = d$bg, bgBatch = d$bgBatch,
    maxIter = 5, verbose = FALSE
  )
  expect_length(unique(res$z[d$batch == "A"]), 1)
  expect_true(all(is.finite(res$contamination)))
})

test_that("decontX names the batch with one cluster and no background", {
  d <- refData()
  x <- d$sim$observedCounts
  expect_error(
    decontX(x, z = rep(1, ncol(x)), verbose = FALSE),
    regexp = paste0(
      "'z' has only one cluster\\. decontX needs at least two clusters ",
      "to estimate contamination without 'background'"
    )
  )
  z <- ifelse(d$batch == "B", 1, d$sim$z)
  expect_error(
    decontX(x, z = z, batch = d$batch, verbose = FALSE),
    regexp = "'z' has only one cluster in batch 'B'"
  )
  # A batch the user named "all_cells" is named too
  batch <- ifelse(d$batch == "A", "all_cells", "S2")
  z <- ifelse(batch == "all_cells", 1, d$sim$z)
  expect_error(
    decontX(x, z = z, batch = batch, verbose = FALSE),
    regexp = "'z' has only one cluster in batch 'all_cells'"
  )
})

test_that("the k-means fallback clusters the PCA components", {
  d <- refData()
  x <- d$sim$observedCounts
  received <- NULL
  kmeans <- stats::kmeans
  local_mocked_bindings(
    kmeans = function(x, centers, ...) {
      received <<- x
      kmeans(x, centers, ...)
    },
    .package = "stats"
  )
  # dbscanEps = 0 skips dbscan, so the fallback runs
  decontX(x, dbscanEps = 0, maxIter = 5, verbose = FALSE)
  # One row per cell, and at most 50 PCs rather than all 60 genes
  expect_identical(nrow(received), ncol(x))
  expect_lte(ncol(received), 50L)
  # The PCs of the default initialization (varGenes = 5000, seed = 12345);
  # decontX runs it on a dgCMatrix, so the last digits can differ
  pcs <- .decontxInitializeZScrapper(x,
    varGenes = 5000, seed = 12345,
    numThreads = 1
  )$pcs
  expect_equal(received, pcs)
})

test_that("the legacy k-means fallback clusters all genes", {
  skip_if_not_installed("scater")
  d <- refData()
  x <- d$sim$observedCounts
  received <- NULL
  kmeans <- stats::kmeans
  local_mocked_bindings(
    kmeans = function(x, centers, ...) {
      received <<- x
      kmeans(x, centers, ...)
    },
    .package = "stats"
  )
  # dbscanEps = 0 skips dbscan, so the fallback runs. The legacy path
  # calls deprecated upstream functions, so warnings are suppressed.
  res <- suppressWarnings(
    decontX(x,
      dbscanEps = 0, maxIter = 5, legacyInit = TRUE,
      verbose = FALSE
    )
  )
  # A dense cells x genes matrix
  expect_true(is.matrix(received))
  expect_identical(dim(received), rev(dim(x)))
  expect_identical(length(unique(res$z)), 2L)
})

test_that("decontX returns the user's cluster labels as a factor", {
  d <- refData()
  x <- d$sim$observedCounts
  labels <- c("T", "B", "Mono")[d$sim$z]
  res <- decontX(x, z = labels, maxIter = 5, verbose = FALSE)
  expect_identical(res$z, factor(labels))

  # Labels found in only one batch each come back without a prefix
  own <- paste0(labels, ifelse(d$batch == "A", "", "2"))
  res <- decontX(x, z = own, batch = d$batch, maxIter = 5, verbose = FALSE)
  expect_identical(res$z, factor(own))

  # Labels found in more than one batch get the batch prefix, ordered by
  # batch and then by label
  res <- decontX(x,
    z = labels, batch = d$batch, maxIter = 5, verbose = FALSE
  )
  expected <- c("A-B", "A-Mono", "A-T", "B-B", "B-Mono", "B-T")
  expect_identical(
    res$z,
    factor(paste0(d$batch, "-", labels), levels = expected)
  )
  # The estimates of each batch keep the labels without the prefix
  expect_identical(colnames(res$estimates$B$phi), c("B", "Mono", "T"))

  # A factor keeps its level order, without unused levels
  z <- factor(labels, levels = c("Mono", "unused", "T", "B"))
  res <- decontX(x, z = z, maxIter = 5, verbose = FALSE)
  expect_identical(res$z, droplevels(z))
  res <- decontX(x, z = z, batch = d$batch, maxIter = 5, verbose = FALSE)
  expected <- c("A-Mono", "A-T", "A-B", "B-Mono", "B-T", "B-B")
  expect_identical(
    res$z,
    factor(paste0(d$batch, "-", z), levels = expected)
  )

  # More than nine character labels, which used to be coded as text
  z12 <- rep(paste0("c", 1:12), length.out = ncol(x))
  res <- decontX(x, z = z12, maxIter = 5, verbose = FALSE)
  expect_identical(as.character(res$z), z12)
  # Each cell's column of phi and eta is the one named by its label
  est <- res$estimates$all_cells
  expect_identical(colnames(est$phi), levels(est$z))
  expect_identical(colnames(est$phi)[as.integer(est$z)], z12)
  native <- .nativeCounts(methods::as(x, "dgCMatrix"), est)
  expect_equal(unname(as.matrix(native)), unname(as.matrix(res$decontXcounts)))

  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = x))
  sceOne <- decontX(sce, z = labels, maxIter = 5, verbose = FALSE)
  expect_identical(sceOne$decontX_clusters, factor(labels))
  sceBatched <- decontX(sce,
    z = own, batch = d$batch, maxIter = 5, verbose = FALSE
  )
  expect_identical(sceBatched$decontX_clusters, factor(own))
  sceBatched <- decontX(sce,
    z = labels, batch = d$batch, maxIter = 5, verbose = FALSE
  )
  expect_identical(
    as.character(sceBatched$decontX_clusters),
    paste0(d$batch, "-", labels)
  )
})

test_that("phi and eta columns are named and ordered by the labels in z", {
  d <- refData()
  x <- d$sim$observedCounts
  res <- decontX(x,
    z = d$sim$z, batch = d$batch, maxIter = 5, verbose = FALSE
  )
  for (est in res$estimates) {
    expect_true(is.factor(est$z))
    expect_identical(colnames(est$phi), levels(est$z))
    expect_identical(colnames(est$eta), levels(est$z))
  }

  # Relabelling the clusters reorders the columns but not the values
  relabelled <- decontX(x,
    z = c(30, 10, 20)[d$sim$z], maxIter = 5, verbose = FALSE
  )
  res <- decontX(x, z = d$sim$z, maxIter = 5, verbose = FALSE)
  est <- relabelled$estimates$all_cells
  expect_identical(colnames(est$phi), c("10", "20", "30"))
  expect_identical(
    unname(est$phi[, c("30", "10", "20")]),
    unname(res$estimates$all_cells$phi[, c("1", "2", "3")])
  )
  expect_identical(relabelled$contamination, res$contamination)
  expect_identical(relabelled$decontXcounts, res$decontXcounts)
})

test_that("estimated clusters come back as a factor", {
  d <- refData()
  x <- d$sim$observedCounts
  res <- decontX(x, maxIter = 5, verbose = FALSE)
  expect_true(is.factor(res$z))
  expect_identical(levels(res$z), as.character(seq_len(nlevels(res$z))))

  res <- decontX(x, batch = d$batch, maxIter = 5, verbose = FALSE)
  expect_true(is.factor(res$z))
  expect_true(all(startsWith(as.character(res$z), paste0(d$batch, "-"))))
  expect_identical(levels(res$z), c(
    paste0("A-", levels(res$estimates$A$z)),
    paste0("B-", levels(res$estimates$B$z))
  ))
})

test_that("unused cluster levels don't give NaN estimates", {
  d <- refData()
  x <- d$sim$observedCounts
  # An unused level before used ones
  z <- factor(d$sim$z, levels = c(1, 4, 2, 3))
  res <- decontX(x, z = z, maxIter = 5, verbose = FALSE)
  expect_identical(levels(res$estimates$all_cells$z), c("1", "2", "3"))
  expect_false(anyNA(res$estimates$all_cells$phi))

  # A level that batch B doesn't have
  z <- d$sim$z
  z[d$batch == "B" & z == 2] <- 1
  res <- decontX(x,
    z = factor(z), batch = d$batch, maxIter = 5, verbose = FALSE
  )
  expect_identical(levels(res$estimates$B$z), c("1", "3"))
  expect_false(anyNA(res$estimates$B$phi))
})

test_that("marker percentages are right for any cluster labels", {
  d <- refData()
  x <- d$sim$observedCounts
  markers <- d$sim$markers
  ref <- plotDecontXMarkerPercentage(x, markers, z = d$sim$z)$data
  for (codes in list(10:12, c(2, 5, 9), 0:2)) {
    p <- plotDecontXMarkerPercentage(x, markers, z = codes[d$sim$z])$data
    expect_identical(p$percent, ref$percent)
    expect_identical(
      as.character(p$cellTypeLabels),
      as.character(codes[ref$cellType])
    )
  }

  # An unused level is left out
  z <- factor(d$sim$z, levels = 0:3)
  p <- plotDecontXMarkerPercentage(x, markers, z = z)$data
  expect_identical(p$percent, ref$percent)
  expect_identical(levels(p$cellTypeLabels), c("1", "2", "3"))

  # One marker gene
  one <- list(One = markers[[1]][1])
  p <- plotDecontXMarkerPercentage(x, one, z = d$sim$z)$data
  expect_identical(nrow(p), 3L)

  # groupClusters that keeps one cell
  z <- c(1, rep(2, ncol(x) - 1))
  p <- plotDecontXMarkerPercentage(x, markers, list(A = 1), z = z)$data
  expect_identical(as.character(unique(p$cellTypeLabels)), "A")
  df <- plotDecontXMarkerExpression(x, markers[[1]], list(A = 1), z = z)$data
  expect_identical(nrow(df), length(markers[[1]]))
})

test_that("the marker expression plot gives each cell its own cluster", {
  d <- refData()
  order <- c(51:100, 1:50)
  x <- d$sim$observedCounts[, order]
  z <- d$sim$z[order]
  # Numeric-like cell names, permuted or not 1..n
  for (cells in list(as.character(order), as.character(1000 + order))) {
    colnames(x) <- cells
    df <- plotDecontXMarkerExpression(x, d$sim$markers[[1]], z = z)$data
    expected <- z[match(as.character(df$Cell), cells)]
    expect_identical(as.character(df$Cluster), as.character(expected))
  }

  # Two assays of a SingleCellExperiment, with and without groupClusters
  sce <- SingleCellExperiment::SingleCellExperiment(
    list(counts = x, decontXcounts = x)
  )
  sce$cluster <- z
  df <- plotDecontXMarkerExpression(sce, d$sim$markers[[1]],
    z = "cluster"
  )$data
  expected <- z[match(as.character(df$Cell), cells)]
  expect_identical(as.character(df$Cluster), as.character(expected))
  groups <- list(A = 1, B = 2)
  df <- plotDecontXMarkerExpression(sce, d$sim$markers[[1]], groups,
    z = "cluster"
  )$data
  expected <- z[match(as.character(df$Cell), cells)]
  expect_identical(
    as.character(df$Cell_Type),
    c("A", "B")[expected]
  )
})

test_that("the groupClusters error lists the labels in z", {
  d <- refData()
  x <- d$sim$observedCounts
  markers <- d$sim$markers
  # Batched decontX runs prefix the labels with the batch
  z <- paste0("b1-", d$sim$z)
  groups <- list(T = c(1, 2))
  expect_error(
    plotDecontXMarkerPercentage(x, markers, groups, z = z),
    regexp = "not found in 'z': 1,2. Labels in 'z': b1-1, b1-2, b1-3"
  )
  expect_error(
    plotDecontXMarkerExpression(x, markers[[1]], groups, z = z),
    regexp = "Labels in 'z': b1-1, b1-2, b1-3",
    fixed = TRUE
  )
  # Only the first 10 labels are listed
  z <- paste0("b", rep(1:4, each = 25), "-", d$sim$z)
  expect_error(
    plotDecontXMarkerPercentage(x, markers, groups, z = z),
    regexp = "b4-1, \\.\\.\\.$"
  )
})

test_that("the contamination plot reads the UMAP columns by position", {
  umap <- withr::with_seed(1, matrix(rnorm(10), ncol = 2))
  res <- list(estimates = list(all_cells = list(
    contamination = seq(0.1, 0.5, by = 0.1),
    UMAP = umap
  )))
  # Unnamed columns, as after an h5ad round trip
  b <- ggplot2::ggplot_build(plotDecontXContamination(res))
  expect_equal(b$data[[1]]$x, umap[, 1])
  expect_equal(b$data[[1]]$y, umap[, 2])
  expect_identical(b$plot$labels$x, "UMAP_1")
  expect_identical(b$plot$labels$y, "UMAP_2")

  # So are empty and NA names
  for (names in list(c("", ""), c(NA, NA))) {
    colnames(res$estimates[[1]]$UMAP) <- names
    b <- ggplot2::ggplot_build(plotDecontXContamination(res))
    expect_identical(b$plot$labels$x, "UMAP_1")
    expect_identical(b$plot$labels$y, "UMAP_2")
  }

  # Named columns are the axis titles
  colnames(res$estimates[[1]]$UMAP) <- c("DecontX_UMAP_1", "DecontX_UMAP_2")
  b <- ggplot2::ggplot_build(plotDecontXContamination(res))
  expect_identical(b$plot$labels$x, "DecontX_UMAP_1")
  expect_identical(b$plot$labels$y, "DecontX_UMAP_2")
})

test_that("legacyInit asks for scater when it isn't installed", {
  d <- refData()
  local_mocked_bindings(.hasScater = function() FALSE)
  expect_error(
    decontX(d$sim$observedCounts,
      z = d$sim$z, legacyInit = TRUE,
      verbose = FALSE
    ),
    regexp = "'legacyInit = TRUE' requires the 'scater' package"
  )
})

test_that("seed = NULL takes the EM's random start from the session", {
  d <- refData()
  x <- d$sim$observedCounts
  seeded <- decontX(x, z = d$sim$z, seed = 7, verbose = FALSE)
  unseeded <- withr::with_seed(7,
    decontX(x, z = d$sim$z, seed = NULL, verbose = FALSE)
  )
  expect_identical(unseeded$contamination, seeded$contamination)
  expect_identical(unseeded$decontXcounts, seeded$decontXcounts)

  expect_identical(
    withr::with_seed(3, simulateContamination(C = 50, G = 40, seed = NULL)),
    simulateContamination(C = 50, G = 40, seed = 3)
  )
})

test_that("seed = NULL with legacyInit uses the session's random numbers", {
  skip_if_not_installed("scater")
  d <- refData()
  x <- d$sim$observedCounts
  # The legacy path's upstream deprecation warnings are suppressed
  seeded <- suppressWarnings(
    decontX(x, z = d$sim$z, seed = 7, legacyInit = TRUE, verbose = FALSE)
  )
  unseeded <- suppressWarnings(withr::with_seed(7,
    decontX(x, z = d$sim$z, seed = NULL, legacyInit = TRUE, verbose = FALSE)
  ))
  expect_equal(
    unseeded$estimates$all_cells$UMAP,
    seeded$estimates$all_cells$UMAP
  )
})

test_that(".withSeed seeds the code only when given a seed", {
  withr::local_seed(1)
  expect_identical(
    .withSeed(10, stats::runif(2)),
    withr::with_seed(10, stats::runif(2))
  )
  # A seeded run leaves the caller's stream where it was
  withr::with_seed(5, {
    .withSeed(10, stats::runif(3))
    afterSeeded <- stats::runif(1)
  })
  expect_identical(afterSeeded, withr::with_seed(5, stats::runif(1)))
  # NULL draws from the caller's stream
  withr::with_seed(5, unseeded <- .withSeed(NULL, stats::runif(2)))
  expect_identical(unseeded, withr::with_seed(5, stats::runif(2)))
})

test_that("decontX logs to standard output when logfile is empty", {
  # As cat(file = "") did
  d <- refData()
  expect_output(
    decontX(d$sim$observedCounts, z = d$sim$z, maxIter = 1, logfile = ""),
    "Starting DecontX"
  )
})

test_that("decontX writes its log to logfile, not the console", {
  d <- refData()
  logfile <- withr::local_tempfile(fileext = ".log")
  expect_silent(decontX(d$sim$observedCounts, z = d$sim$z, logfile = logfile))
  log <- readLines(logfile)
  expect_match(log, "Starting DecontX", all = FALSE)
  expect_match(log, "Estimating contamination", all = FALSE)
  expect_match(log, "Completed DecontX", all = FALSE)
})

test_that("decontX warns when the background has no column names", {
  d <- refData()
  bg <- d$bg
  colnames(bg) <- NULL
  expect_warning(
    decontX(d$sim$observedCounts,
      z = d$sim$z, background = bg,
      verbose = FALSE
    ),
    regexp = "No column names were found for the 'background' matrix"
  )
})

test_that("decontX says when an SCE has no named assays", {
  d <- refData()
  sce <- SingleCellExperiment::SingleCellExperiment(
    list(d$sim$observedCounts)
  )
  expect_error(decontX(sce, z = d$sim$z, verbose = FALSE),
    regexp = "is not an assay of 'x'\\. It has no named assays\\."
  )
})

test_that("varGenes = NULL uses 5000 genes", {
  d <- refData()
  expect_no_error(decontX(d$sim$observedCounts,
    z = d$sim$z, varGenes = NULL, maxIter = 2,
    verbose = FALSE
  ))
  # The fixture has fewer genes than that, so check the number itself
  expect_identical(.processvarGenes(NULL), 5000)
})

test_that("retrieveFeatureIndex searches rowData and data.frame columns", {
  counts <- matrix(1:12, nrow = 4, dimnames = list(paste0("g", 1:4), NULL))
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
  SummarizedExperiment::rowData(sce)$symbol <- c("A", "B", "C", "D")
  expect_identical(
    retrieveFeatureIndex(c("C", "A"), sce, by = "symbol"),
    c(3L, 1L)
  )
  expect_error(retrieveFeatureIndex("A", sce, by = "id"),
    regexp = "'by' is not a column in 'rowData\\(x\\)'\\."
  )

  df <- data.frame(symbol = c("A", "B", "C", "D"))
  expect_identical(retrieveFeatureIndex(c("D", "B"), df, by = "symbol"),
    c(4L, 2L)
  )
  expect_error(retrieveFeatureIndex("A", df, by = "id"),
    regexp = "'by' is not a column in 'x'\\."
  )

  # A vector is searched itself
  expect_identical(
    retrieveFeatureIndex(c("D", "B"), c("A", "B", "C", "D"), by = "symbol"),
    c(4L, 2L)
  )
})

test_that("retrieveFeatureIndex reports features it can't place", {
  counts <- matrix(1:16,
    nrow = 4,
    dimnames = list(c("CD3D", "CD3E", "CD4", "CD8A"), NULL)
  )
  expect_error(retrieveFeatureIndex("CD4", unname(counts)),
    regexp = "'rownames' of 'x' are 'NULL'"
  )

  # A partial match with several rows keeps the first
  expect_warning(
    res <- retrieveFeatureIndex(c("CD3", "CD4"), counts, exactMatch = FALSE),
    regexp = "Feature 'CD3' matched multiple items in 'rownames': CD3D,CD3E"
  )
  expect_identical(res, c(1L, 3L))

  expect_error(retrieveFeatureIndex(c("X", "Y"), counts),
    regexp = "None of the provided features .* try setting 'exactMatch"
  )
  expect_error(retrieveFeatureIndex(c("X", "Y"), counts, exactMatch = FALSE),
    regexp = "None of the provided features .* make sure 'by' is set"
  )

  # removeNA drops the features that weren't found
  expect_warning(
    res <- retrieveFeatureIndex(c("CD4", "X"), counts, removeNA = TRUE),
    regexp = "not present in 'x': X"
  )
  expect_identical(res, 3L)
})

test_that("plotDecontXContamination plots the batch asked for", {
  d <- refData()
  res <- decontX(d$sim$observedCounts,
    z = d$sim$z, batch = d$batch,
    verbose = FALSE
  )
  b <- ggplot2::ggplot_build(plotDecontXContamination(res, batch = "B"))
  expect_equal(b$data[[1]]$x, unname(res$estimates$B$UMAP[, 1]))
  expect_equal(
    b$plot$data$Contamination,
    unname(res$estimates$B$contamination)
  )

  expect_error(plotDecontXContamination(res, batch = "C"),
    regexp = "'C' is not one of the batches in 'x'\\. .* plotting: 'A,B'"
  )
  expect_error(plotDecontXContamination(list()),
    regexp = "decontX estimates not found"
  )
})

test_that("the marker plots check markers, groupClusters and z", {
  d <- refData()
  x <- d$sim$observedCounts
  z <- d$sim$z
  markers <- d$sim$markers
  expect_error(plotDecontXMarkerPercentage(x, unlist(markers), z = z),
    regexp = "'markers' needs to be a named list\\."
  )
  expect_error(
    plotDecontXMarkerPercentage(x, markers, list(1, 2), z = z),
    regexp = "'groupClusters' needs to be a named list\\."
  )
  expect_error(
    plotDecontXMarkerPercentage(x, markers, list(A = 1:2, B = 2), z = z),
    regexp = "clusters: 2\\. Clusters need be assigned to a single group\\."
  )
  expect_error(plotDecontXMarkerPercentage(x, markers, z = z[-1]),
    regexp = "the length of 'z' needs to\\s+be the same as the number"
  )

  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = x))
  expect_error(plotDecontXMarkerExpression(sce, markers[[1]]),
    regexp = "'decontX_clusters' not found in 'colData\\(x\\)'"
  )
  expect_error(plotDecontXMarkerExpression(sce, markers[[1]], z = "cluster"),
    regexp = "'cluster' not found in 'colData\\(x\\)'\\."
  )
})

test_that("the marker expression plot can log1p and add dots", {
  d <- refData()
  x <- d$sim$observedCounts
  markers <- d$sim$markers[[1]]
  p <- plotDecontXMarkerExpression(x, markers,
    z = d$sim$z,
    log1p = TRUE, plotDots = TRUE
  )
  b <- ggplot2::ggplot_build(p)
  expect_equal(sort(b$plot$data$Expression), sort(log1p(x[markers, ])))
  expect_identical(b$plot$labels$y, "Expression (log1p)")
  expect_s3_class(b$plot$layers[[2]]$position, "PositionJitter")

  # Without groups, each cluster is its own group
  expect_setequal(
    as.character(b$layout$layout$Cell_Type),
    as.character(unique(d$sim$z))
  )
})

test_that("internal checks stop on input no exported function passes", {
  # The exported functions check these first; the internal checks guard
  # other callers
  expect_error(.logMessages("text", logfile = 1),
    regexp = "The log file parameter needs to be a single character string"
  )
  expect_error(.checkCountsDecon(1:10),
    regexp = "At least 2 genes need to have non-zero expressions"
  )
  expect_error(.processCellLabels(1:3, numCells = 4),
    regexp = "'z' must be of the same length as the number of cells"
  )
})

test_that("decontXcounts are the counts times the native probability", {
  sim <- simulateContamination(C = 30, G = 40, K = 3, seed = 1)
  counts <- methods::as(sim$observedCounts, "dgCMatrix")
  # Dense R form of the native counts, with the EM's pseudocount
  expectedNative <- function(counts, est) {
    k <- as.integer(est$z)
    theta <- rep(est$theta, each = nrow(counts))
    native <- (est$phi[, k] + 1e-20) * (theta + 1e-20)
    contaminant <- (est$eta[, k] + 1e-20) * (1 - theta + 1e-20)
    unname(as.matrix(counts) * native / (native + contaminant))
  }

  res <- decontX(counts, z = sim$z, embedding = FALSE, verbose = FALSE)
  expect_equal(unname(as.matrix(res$decontXcounts)),
    expectedNative(counts, res$estimates$all_cells),
    tolerance = 1e-12
  )
  # Same sparsity pattern as the counts
  expect_identical(res$decontXcounts@i, counts@i)
  expect_identical(res$decontXcounts@p, counts@p)

  # Interleaved batches: each batch's counts go back to their own cells
  batch <- rep(c("a", "b"), length.out = ncol(counts))
  res <- decontX(counts,
    z = sim$z, batch = batch, embedding = FALSE, verbose = FALSE
  )
  expected <- matrix(0, nrow(counts), ncol(counts))
  for (b in c("a", "b")) {
    cells <- batch == b
    expected[, cells] <- expectedNative(
      counts[, cells, drop = FALSE], res$estimates[[b]]
    )
  }
  expect_equal(unname(as.matrix(res$decontXcounts)), expected,
    tolerance = 1e-12
  )

  # A subclass of dgCMatrix gives a plain dgCMatrix
  classEnv <- new.env()
  countsClass <- methods::setClass("decontXTestCounts",
    contains = "dgCMatrix", slots = c(note = "character"), where = classEnv
  )
  subCounts <- countsClass(counts, note = "extra")
  res <- decontX(subCounts, z = sim$z, embedding = FALSE, verbose = FALSE)
  expect_identical(as.character(class(res$decontXcounts)), "dgCMatrix")
})

test_that("the C++ EM functions check the shapes of their arguments", {
  counts <- methods::as(matrix(as.numeric(1:12), 4, 3), "CsparseMatrix")
  valid <- list(
    counts = counts,
    counts_colsums = Matrix::colSums(counts),
    theta = rep(0.5, 3),
    estimate_eta = TRUE,
    eta = matrix(0.25, 4, 2),
    phi = matrix(0.25, 4, 2),
    z = c(1L, 2L, 1L),
    estimate_delta = FALSE,
    delta = c(10, 10),
    pseudocount = 1e-20
  )
  # Each change breaks one check; the checks run in this order
  changes <- list(
    theta = list(theta = rep(0.5, 2), regexp = "Length of 'theta'"),
    z = list(z = c(1L, 2L), regexp = "Length of 'z'"),
    colsums = list(
      counts_colsums = c(1, 2),
      regexp = "Length of 'counts_colsums'"
    ),
    phiRows = list(phi = matrix(0.25, 3, 2), regexp = "rows in 'phi'"),
    etaRows = list(eta = matrix(0.25, 3, 2), regexp = "rows in 'eta'"),
    etaCols = list(
      eta = matrix(0.25, 4, 3),
      regexp = "columns in 'eta' must be equal to the number of columns"
    ),
    zRange = list(
      z = c(1L, 3L, 1L),
      regexp = "between 1 and the number of columns in eta and phi"
    ),
    delta = list(delta = c(10, -1), regexp = "'delta' must be a numeric")
  )
  checks <- list(
    decontXEM = names(changes),
    decontXLogLik = c("theta", "z", "phiRows", "etaRows", "etaCols", "zRange"),
    calculateNativeCounts = c(
      "theta", "z", "phiRows", "etaRows", "etaCols", "zRange"
    ),
    decontXInitialize = c("theta", "z")
  )
  for (fun in names(checks)) {
    f <- get(fun)
    # A renamed argument would otherwise be passed as NULL
    expect_true(all(names(formals(f)) %in% names(valid)), info = fun)
    for (check in checks[[fun]]) {
      change <- changes[[check]]
      args <- utils::modifyList(valid, change[names(change) != "regexp"])
      expect_error(do.call(f, args[names(formals(f))]),
        regexp = change$regexp, info = paste(fun, check)
      )
    }
    # The valid arguments pass every check
    expect_no_error(do.call(f, valid[names(formals(f))]))
  }
})

test_that("cluster labels are numbered in the order the EM has used", {
  # Labels are numbered by first appearance. Character labels stay
  # character, so with 10 or more of them the codes sort as text
  z <- paste0("c", c(1:12, 1:12))
  labels <- .processCellLabels(z, numCells = 24)
  expect_identical(
    levels(labels),
    c("1", "10", "11", "12", "2", "3", "4", "5", "6", "7", "8", "9")
  )
  expect_identical(as.character(labels), as.character(c(1:12, 1:12)))
  expect_identical(
    as.integer(labels)[1:12],
    c(1L, 5:12, 2:4)
  )

  # Numbers and logicals are numbered as numbers
  labels <- .processCellLabels(c(12, 3, 3, 0, 12), numCells = 5)
  expect_identical(labels, factor(c(1, 2, 2, 3, 1)))
  labels <- .processCellLabels(c(TRUE, FALSE, TRUE), numCells = 3)
  expect_identical(labels, factor(c(1, 2, 1)))

  # Factors keep their labels, without unused levels
  z <- factor(c("b", "a", "b"), levels = c("c", "b", "a"))
  expect_identical(
    .processCellLabels(z, numCells = 3),
    factor(c("b", "a", "b"), levels = c("b", "a"))
  )
})

test_that("the marker plots' data have one row per marker and cell", {
  x <- matrix(c(0, 3, 5, 1, 2, 0),
    nrow = 2,
    dimnames = list(c("CD3", "CD19"), c("10", "2", "3"))
  )
  markers <- list(T = "CD3", B = "CD19")
  # Cell names that look like numbers become numbers
  df <- plotDecontXMarkerExpression(x, c("CD3", "CD19"), z = c(1, 2, 1))$data
  expect_identical(df[1:3], data.frame(
    Marker = factor(rep(c("CD3", "CD19"), 3), levels = c("CD3", "CD19")),
    Cell = rep(c(10L, 2L, 3L), each = 2),
    Expression = c(0, 3, 5, 1, 2, 0)
  ))
  expect_identical(as.character(df$Cluster), rep(c("1", "2", "1"), each = 2))

  # Other names become factors in their order
  colnames(x) <- c("c", "a", "b")
  df <- plotDecontXMarkerExpression(x, "CD19", z = c(1, 2, 1))$data
  expect_identical(df$Cell, factor(c("c", "a", "b"), levels = c("c", "a", "b")))
  expect_identical(df$Expression, c(3, 1, 0))

  # Without cell names, the cells are numbered
  sce <- SingleCellExperiment::SingleCellExperiment(
    list(counts = unname(x) + 1, decontXcounts = unname(x))
  )
  rownames(sce) <- rownames(x)
  df <- plotDecontXMarkerExpression(sce, "CD3", z = c(1, 2, 1),
    assayName = c("counts", "decontXcounts")
  )$data
  expect_identical(df$Cell, rep(1:3, 2))
  expect_identical(df$Expression, c(1, 6, 3, 0, 5, 2))
  expect_identical(df$assay, factor(rep(c("counts", "decontXcounts"),
    each = 3
  ), levels = c("counts", "decontXcounts")))

  # Markers and clusters are numbered
  df <- plotDecontXMarkerPercentage(x, markers, z = c(1, 2, 1))$data
  expect_identical(df$markers, c(1L, 2L, 1L, 2L))
  expect_identical(df$cellType, c(1L, 1L, 2L, 2L))
  expect_identical(df$percent, c(50, 50, 100, 100))
})
