library(decontX)
context("Testing DecontX functions")

deconSim <- simulateContamination(K = 10, delta = c(1, 5))
modelDecontXoneBatch <- decontX(deconSim$observedCounts,
        z = deconSim$z,
        maxIter = 2)

deconSim2 <- simulateContamination(K = 10, delta = c(1, 5))
batchDecontX <- decontX(cbind(deconSim$observedCounts,
    deconSim2$observedCounts),
        z = c(deconSim$z, deconSim2$z),
        batch = rep(seq(2), each = ncol(deconSim$observedCounts)),
        maxIter = 2)

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
                      dimnames = list(paste0("g", seq_len(400)),
                                      paste0("c", seq_len(100))))
  lowCounts[1, Matrix::colSums(lowCounts) == 0] <- 1
  res <- decontX(lowCounts, seed = 1, maxIter = 2)
  expect_true(all(res$contamination >= 0 & res$contamination <= 1))
})

test_that(desc = "decontX initialization handles duplicated gene names", {
  # Regression test: real 10x data labeled with gene symbols often has
  # duplicated rownames, which scrapper::modelGeneVariances rejects.
  counts <- deconSim$observedCounts
  rownames(counts) <- c("DUP", "DUP",
                        paste0("g", seq_len(nrow(counts) - 2)))
  res <- decontX(counts, seed = 1, maxIter = 2)
  expect_true(all(res$contamination >= 0 & res$contamination <= 1))
})

test_that(desc = "decontX gives an informative error for empty cells", {
  counts <- deconSim$observedCounts
  counts[, 1] <- 0
  expect_error(decontX(counts, seed = 1, maxIter = 2),
               regexp = "at least one count")
})

test_that(desc = "legacyInit reproduces the scater-based initialization", {
  skip_if_not_installed("scater")
  # The legacy path calls upstream functions that are deprecated on
  # Bioconductor devel (by design -- it exists to reproduce old results
  # while they still work), so warnings are suppressed here.
  res <- suppressWarnings(
    decontX(deconSim$observedCounts, seed = 1, maxIter = 2,
            legacyInit = TRUE)
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

  resNum <- decontX(combined, z = z, batch = numBatch,
                    maxIter = 2, seed = 12345)
  resChar <- decontX(combined, z = z, batch = charBatch,
                     maxIter = 2, seed = 12345)

  expect_equal(resNum$contamination, resChar$contamination)
  expect_equal(resNum$decontXcounts, resChar$decontXcounts)
  expect_true(all(resNum$contamination >= 0 &
                    resNum$contamination <= 1))
  expect_false(is.null(resNum$estimates[["3"]]))
  expect_false(is.null(resNum$estimates[["5"]]))
})

test_that(desc = "Testing simulateContamination", {
    expect_equivalent(object = colSums(deconSim$observedCounts),
        expected = deconSim$NByC)
    expect_equal(object = dim(deconSim$phi),
        expected = dim(deconSim$eta))
    expect_equal(typeof(deconSim$observedCounts), "integer")
    expect_warning(simulateContamination(K = 101, C = 10))
    expect_error(simulateContamination(K = 3, G = 2, numMarkers = 10))
})

## DecontX
test_that(desc = "Testing DecontX on counts matrix", {
  s <- simulateContamination()
  res <- decontX(s$observedCounts)
  p <- plotDecontXMarkerPercentage(s$observedCounts,
                                   z = res$z,
                                   markers = s$markers)
  p <- plotDecontXMarkerPercentage(res$decontXcounts,
                                   z = res$z,
                                   markers = s$markers)
  p <- plotDecontXMarkerExpression(s$observedCounts,
                                   s$markers[[1]],
                                   z = s$z)
  p <- plotDecontXContamination(res)

  # test with background input
  b <- s$observedCounts[, 1:5]
  colnames(b) <- paste(colnames(b), "_", sep = "")
  res <- decontX(s$observedCounts,
                 background = b)
})

test_that(desc = "Testing DecontX on SCE", {
  s <- simulateContamination()
  sce <- SingleCellExperiment::SingleCellExperiment(
                               list(counts = s$observedCounts))
  sce <- decontX(sce)
  p <- plotDecontXContamination(sce)
  p <- plotDecontXMarkerPercentage(sce,
                                   z = s$z,
                                   markers = s$markers,
                                   assayName = "decontXcounts")
  p <- plotDecontXMarkerExpression(sce, s$markers[[1]])
  newz <- paste0("X", s$z)
  sce$newz2 <- newz
  p <- plotDecontXMarkerPercentage(sce,
                                   z = "newz2",
                                   markers = s$markers,
                                   assayName = "decontXcounts")
  sce <- decontX(sce, estimateDelta = FALSE)

  # test with background input
  bg <- sce[, 1:5]
  colnames(bg) <- paste(colnames(bg), "_", sep = "")
  sce <- decontX(sce, background = bg)
})


## .decontXoneBatch
test_that(desc = "Testing .decontXoneBatch", {
    expect_error(decontX(x = deconSim$observedCounts,
        z = deconSim$z, delta = c(1, -1)))
    expect_error(decontX(x = deconSim$observedCounts,
        z = deconSim$z, delta = c(1, 1, 1)))
    expect_error(decontX(x = deconSim$observedCounts,
        z = c(deconSim$z, 1)),
        paste0("'z' must be of the same length as the number of cells in the",
            " 'counts' matrix."))
    expect_error(.decontXoneBatch(counts = deconSim$observedCounts,
        z = rep(1, ncol(
            deconSim$observedCounts))),
        "No need to decontaminate when only one cluster is in the dataset.")
    countsNA <- deconSim$observedCounts
    countsNA[1, 1] <- NA
    expect_error(.decontXoneBatch(counts = countsNA, z = deconSim$z),
        "Missing value in 'counts' matrix.")
})
