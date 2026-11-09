test_that(desc = "Testing simulateContamination", {
  deconSim <- simulateContamination(K = 10, delta = c(1, 5))
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

test_that("simulateContamination uses NRange when its ends are equal", {
  sim <- simulateContamination(C = 20, G = 30, K = 2, NRange = c(500, 500))
  expect_equal(sim$NByC, rep(500L, 20))
  expect_equal(colSums(sim$observedCounts), rep(500, 20), ignore_attr = TRUE)
})

test_that("simulateContamination needs at least two populations", {
  expect_error(simulateContamination(K = 1),
    regexp = "'K' must be at least 2"
  )
  expect_error(simulateContamination(K = 2, C = 1),
    regexp = "Only 1 cell population was simulated"
  )
})

test_that("simulateContamination uses one delta value for both shapes", {
  expect_identical(
    simulateContamination(C = 50, G = 40, delta = 2),
    simulateContamination(C = 50, G = 40, delta = c(2, 2))
  )
})

test_that("simulateContamination numbers fewer clusters than K from 1", {
  expect_warning(
    sim <- simulateContamination(C = 10, G = 40, K = 50),
    regexp = "clusters are simulated"
  )
  expect_type(sim$z, "integer")
  expect_identical(unique(sim$z), seq_len(max(sim$z)))
})
