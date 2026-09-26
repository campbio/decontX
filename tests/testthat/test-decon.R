library(decontX)

deconSim <- simulateContamination(K = 10, delta = c(1, 5))

test_that(desc = "decontX does not call deprecated normalization functions", {
  # Regression guard for the Bioc 3.24 scuttle deprecations
  # (logNormCounts/normalizeCounts -> scrapper): decontX initialization
  # must not emit these deprecation warnings, which R CMD check flags
  # and which become errors when upstream turns them defunct. Scoped to
  # the functions decontX calls so unrelated upstream deprecations
  # cannot fail this test.
  expect_no_warning(
    decontX(deconSim$observedCounts, seed = 1, maxIter = 2),
    message = "'(normalizeCounts|logNormCounts)' is deprecated"
  )
})

test_that(desc = "decontX initialization without z returns valid results", {
  res <- decontX(deconSim$observedCounts, seed = 1, maxIter = 2)
  expect_true(all(res$contamination >= 0 & res$contamination <= 1))
  expect_equal(length(res$z), ncol(deconSim$observedCounts))
  expect_gt(length(unique(res$z)), 1)
  umap <- res$estimates$all_cells$UMAP
  expect_equal(dim(umap), c(ncol(deconSim$observedCounts), 2))

  # Same seed twice gives identical results
  res2 <- decontX(deconSim$observedCounts, seed = 1, maxIter = 2)
  expect_equal(res$contamination, res2$contamination)
  expect_equal(res$z, res2$z)
})

test_that(desc = "decontX initialization handles very sparse data", {
  # Regression test: scrapper::modelGeneVariances' default abundance
  # filter can drop every gene in low-count data (e.g. shallow
  # sequencing); initialization must rank all genes instead of erroring.
  set.seed(7)
  lowCounts <- matrix(stats::rbinom(400 * 100, size = 1, prob = 0.05),
    nrow = 400, ncol = 100,
    dimnames = list(
      paste0("g", seq_len(400)),
      paste0("c", seq_len(100))
    )
  )
  lowCounts[1, Matrix::colSums(lowCounts) == 0] <- 1
  res <- decontX(lowCounts, seed = 1, maxIter = 2)
  expect_true(all(res$contamination >= 0 & res$contamination <= 1))
})

test_that(desc = "decontX initialization handles duplicated gene names", {
  # Regression test: real 10x data labeled with gene symbols often has
  # duplicated rownames, which scrapper::modelGeneVariances rejects.
  counts <- deconSim$observedCounts
  rownames(counts) <- c(
    "DUP", "DUP",
    paste0("g", seq_len(nrow(counts) - 2))
  )
  res <- decontX(counts, seed = 1, maxIter = 2)
  expect_true(all(res$contamination >= 0 & res$contamination <= 1))
})

test_that(desc = "decontX gives an informative error for empty cells", {
  counts <- deconSim$observedCounts
  counts[, 1] <- 0
  expect_error(decontX(counts, seed = 1, maxIter = 2),
    regexp = "at least one count"
  )
})

test_that(desc = "legacyInit reproduces the scater-based initialization", {
  skip_if_not_installed("scater")
  # The legacy path calls upstream functions that are deprecated on
  # Bioconductor devel (by design -- it exists to reproduce old results
  # while they still work), so warnings are suppressed here.
  res <- suppressWarnings(
    decontX(deconSim$observedCounts,
      seed = 1, maxIter = 2,
      legacyInit = TRUE
    )
  )
  expect_true(all(res$contamination >= 0 & res$contamination <= 1))
  expect_equal(length(res$z), ncol(deconSim$observedCounts))
})

test_that("numeric batch labels are not used as positional indices", {
  sim1 <- simulateContamination(K = 5, delta = c(1, 5))
  sim2 <- simulateContamination(K = 5, delta = c(1, 5))
  combined <- cbind(sim1$observedCounts, sim2$observedCounts)
  z <- c(sim1$z, sim2$z)
  n1 <- ncol(sim1$observedCounts)
  n2 <- ncol(sim2$observedCounts)

  numBatch <- c(rep(3, n1), rep(5, n2))
  charBatch <- c(rep("3", n1), rep("5", n2))

  resNum <- decontX(combined,
    z = z, batch = numBatch,
    maxIter = 2, seed = 12345
  )
  resChar <- decontX(combined,
    z = z, batch = charBatch,
    maxIter = 2, seed = 12345
  )

  expect_equal(resNum$contamination, resChar$contamination)
  expect_equal(resNum$decontXcounts, resChar$decontXcounts)
  expect_true(all(resNum$contamination >= 0 & resNum$contamination <= 1))
  expect_false(is.null(resNum$estimates[["3"]]))
  expect_false(is.null(resNum$estimates[["5"]]))
})

test_that(desc = "Testing simulateContamination", {
  expect_equivalent(
    object = colSums(deconSim$observedCounts),
    expected = deconSim$NByC
  )
  expect_equal(
    object = dim(deconSim$phi),
    expected = dim(deconSim$eta)
  )
  expect_equal(typeof(deconSim$observedCounts), "integer")
  expect_warning(simulateContamination(K = 101, C = 10))
  expect_error(simulateContamination(K = 3, G = 2, numMarkers = 10))

  # ADR-0003 Stage 1: the contamination distribution eta is now normalized
  # by the package's own fastNormProp instead of celda::normalizeCounts.
  # Guard that (a) each cell-type column is a proper proportion and (b) the
  # simulator no longer emits the celda/scuttle deprecation warning that
  # R CMD check flags in examples.
  expect_equal(colSums(deconSim$eta), rep(1, ncol(deconSim$eta)))
  expect_no_warning(
    simulateContamination(K = 3),
    message = "'normalizeCounts' is deprecated"
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
  res <- decontX(sim$observedCounts, z = sim$z, seed = 12345, iterLogLik = 1)

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
test_that(desc = "Testing DecontX on counts matrix", {
  s <- simulateContamination()
  res <- decontX(s$observedCounts)
  expect_equal(dim(res$decontXcounts), dim(s$observedCounts))
  expect_true(all(res$contamination >= 0 & res$contamination <= 1))
  expect_equal(length(res$z), ncol(s$observedCounts))

  p <- plotDecontXMarkerPercentage(s$observedCounts,
    z = res$z,
    markers = s$markers
  )
  expect_s3_class(p, "ggplot")
  p <- plotDecontXMarkerPercentage(res$decontXcounts,
    z = res$z,
    markers = s$markers
  )
  expect_s3_class(p, "ggplot")
  p <- plotDecontXMarkerExpression(s$observedCounts,
    s$markers[[1]],
    z = s$z
  )
  expect_s3_class(p, "ggplot")
  p <- plotDecontXContamination(res)
  expect_s3_class(p, "ggplot")

  # test with background input
  b <- s$observedCounts[, 1:5]
  colnames(b) <- paste(colnames(b), "_", sep = "")
  resBg <- decontX(s$observedCounts,
    background = b
  )
  expect_equal(dim(resBg$decontXcounts), dim(s$observedCounts))
  expect_true(all(resBg$contamination >= 0 & resBg$contamination <= 1))
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

test_that(desc = "Testing DecontX on SCE", {
  s <- simulateContamination()
  sce <- SingleCellExperiment::SingleCellExperiment(
    list(counts = s$observedCounts)
  )
  sce <- decontX(sce)
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

  p <- plotDecontXContamination(sce)
  expect_s3_class(p, "ggplot")
  p <- plotDecontXMarkerPercentage(sce,
    z = s$z,
    markers = s$markers,
    assayName = "decontXcounts"
  )
  expect_s3_class(p, "ggplot")
  p <- plotDecontXMarkerExpression(sce, s$markers[[1]])
  expect_s3_class(p, "ggplot")
  newz <- paste0("X", s$z)
  sce$newz2 <- newz
  p <- plotDecontXMarkerPercentage(sce,
    z = "newz2",
    markers = s$markers,
    assayName = "decontXcounts"
  )
  expect_s3_class(p, "ggplot")
  sce <- decontX(sce, estimateDelta = FALSE)

  # test with background input
  bg <- sce[, 1:5]
  colnames(bg) <- paste(colnames(bg), "_", sep = "")
  sce <- decontX(sce, background = bg)
  expect_true("decontXcounts" %in% SummarizedExperiment::assayNames(sce))
})


## .decontXoneBatch
test_that(desc = "Testing .decontXoneBatch", {
  expect_error(decontX(
    x = deconSim$observedCounts,
    z = deconSim$z, delta = c(1, -1)
  ))
  expect_error(decontX(
    x = deconSim$observedCounts,
    z = deconSim$z, delta = c(1, 1, 1)
  ))
  expect_error(
    decontX(
      x = deconSim$observedCounts,
      z = c(deconSim$z, 1)
    ),
    paste0(
      "'z' must be of the same length as the number of cells in the",
      " 'counts' matrix."
    )
  )
  expect_error(
    .decontXoneBatch(
      counts = deconSim$observedCounts,
      z = rep(1, ncol(
        deconSim$observedCounts
      ))
    ),
    "No need to decontaminate when only one cluster is in the dataset."
  )
  countsNA <- deconSim$observedCounts
  countsNA[1, 1] <- NA
  expect_error(
    .decontXoneBatch(counts = countsNA, z = deconSim$z),
    "Missing value in 'counts' matrix."
  )
})
