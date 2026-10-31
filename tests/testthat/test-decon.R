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
  expect_equal(
    object = colSums(deconSim$observedCounts),
    expected = deconSim$NByC,
    ignore_attr = TRUE
  )
  expect_equal(
    object = dim(deconSim$phi),
    expected = dim(deconSim$eta)
  )
  expect_equal(typeof(deconSim$observedCounts), "integer")
  expect_warning(simulateContamination(K = 101, C = 10),
    regexp = "Only 10 clusters are simulated"
  )
  expect_error(simulateContamination(K = 3, G = 2, numMarkers = 10),
    regexp = "cannot be greater than the number of genes"
  )

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
  expect_error(
    decontX(
      x = deconSim$observedCounts,
      z = deconSim$z, delta = c(1, -1)
    ),
    regexp = "'delta' needs to be a numeric vector of length 2"
  )
  expect_error(
    decontX(
      x = deconSim$observedCounts,
      z = deconSim$z, delta = c(1, 1, 1)
    ),
    regexp = "'delta' needs to be a numeric vector of length 2"
  )
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

test_that("decontX stops when batch and bgBatch don't go together", {
  d <- refData()
  expect_error(
    decontX(d$sim$observedCounts,
      z = d$sim$z, background = d$bg,
      bgBatch = d$bgBatch, verbose = FALSE
    ),
    regexp = "background should also default to no batch"
  )
  expect_error(
    decontX(d$sim$observedCounts,
      z = d$sim$z, batch = d$batch,
      background = d$bg, verbose = FALSE
    ),
    regexp = "Please also supply background batch"
  )
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

test_that("misspelled messages keep their current wording", {
  # #74 fixes the spelling; the snapshots show that change
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
  # Today dense input returns a dgCMatrix (A1-14; #71 changes this)
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
    )
  )
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
  # Initialization runs even when z is supplied (#79 adds an opt-out)
  d <- refData()
  counts <- d$sim$observedCounts
  counts[, 1] <- 0
  err <- expect_error(decontX(counts, z = d$sim$z, verbose = FALSE),
    regexp = "at least one count"
  )
  expect_false(grepl("'z'", conditionMessage(err)))
})

test_that("decontX still accepts integer, factor and Inf arguments", {
  d <- refData()
  counts <- d$sim$observedCounts
  expect_no_error(.checkDecontXArgs(
    counts = counts,
    z = factor(d$sim$z),
    batch = factor(d$batch),
    maxIter = 500L,
    convergence = 0L,
    iterLogLik = 10L,
    delta = c(10L, 10L),
    estimateDelta = TRUE,
    varGenes = 5000L,
    dbscanEps = 0,
    seed = 1L,
    legacyInit = FALSE,
    logfile = NULL,
    verbose = FALSE
  ))
  expect_no_error(.checkDecontXArgs(
    counts = counts,
    z = as.character(d$sim$z),
    batch = NULL,
    maxIter = Inf,
    convergence = 0.001,
    iterLogLik = Inf,
    delta = c(10, 10),
    estimateDelta = FALSE,
    varGenes = Inf,
    dbscanEps = NULL,
    seed = NULL,
    legacyInit = TRUE,
    logfile = tempfile(),
    verbose = TRUE
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
