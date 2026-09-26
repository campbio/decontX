library(decontX)

test_that(desc = "Testing decontPro on count matrix", {
  # Real Stan VB end-to-end run. Kept intentionally small (10 features x
  # 8 droplets) so the MCMC path stays fast; the decomposition-sums-to-input
  # invariant it checks is algebraic and holds regardless of convergence.
  set.seed(42)
  counts <- matrix(sample(1:10, 80, replace = TRUE), nrow = 10, ncol = 8)

  k <- c(1, 1, 2, 2, 3, 3, 4, 4)

  out <- decontPro(counts, k, 1e-2, 1e-2)

  # Sum decomposed matrices
  matsum <- out$decontaminated_counts +
    out$ambient_counts +
    out$background_counts

  expect_equal(matsum, counts)
})

test_that("decontPro does not keep posterior draws that scale with data", {
  # Issue #46: rstan's vb defaults to 1000 draws of every parameter, and
  # delta/background have one value per ADT x droplet, so 215 x 52222 data
  # produced a ~215 GB temp CSV. decontPro only uses the approximation mean.
  set.seed(42)
  counts <- matrix(sample(1:10, 80, replace = TRUE), nrow = 10, ncol = 8)
  k <- c(1, 1, 2, 2, 3, 3, 4, 4)

  fit <- NULL
  process <- .process_stan_vb_out
  local_mocked_bindings(.process_stan_vb_out = function(stan_vb_output, dat) {
    fit <<- stan_vb_output
    process(stan_vb_output, dat)
  })
  suppressWarnings(decontPro(counts, k, 1e-2, 1e-2))

  # Fewer retained values than there are ADT x droplet entries (80)
  expect_lt(length(unlist(fit@sim$samples)), nrow(counts) * ncol(counts))
  # The mean estimate for every parameter is still available
  expect_equal(dim(fit@sim$est$delta), c(8, 10))
  expect_equal(dim(fit@sim$est$background), c(10, 8))
})

test_that("decontPro results match a Stan fit with default output draws", {
  set.seed(42)
  counts <- matrix(sample(1:10, 80, replace = TRUE), nrow = 10, ncol = 8)
  k <- c(1, 1, 2, 2, 3, 3, 4, 4)

  out <- suppressWarnings(decontPro(counts, k, 1e-2, 1e-2))
  local_mocked_bindings(.call_stan_vb = function(data, initial_condition) {
    rstan::vb(stanmodels$shrinkage,
      init = initial_condition, data = data,
      seed = 12345, iter = 50000
    )
  })
  reference <- suppressWarnings(decontPro(counts, k, 1e-2, 1e-2))

  expect_identical(out, reference)
})

test_that("decontPro does not surface rstan's draw-based Pareto k warning", {
  set.seed(42)
  counts <- matrix(sample(1:10, 80, replace = TRUE), nrow = 10, ncol = 8)
  k <- c(1, 1, 2, 2, 3, 3, 4, 4)

  expect_no_warning(decontPro(counts, k, 1e-2, 1e-2), message = "Pareto k")
})

test_that("decontPro errors when ambient_counts feature count differs", {
  # Guard in .decontPro: this stops before any Stan call, so the test is fast.
  counts <- matrix(sample(1:10, 80, replace = TRUE), nrow = 10, ncol = 8)
  k <- c(1, 1, 2, 2, 3, 3, 4, 4)
  ambient <- matrix(sample(1:10, 40, replace = TRUE), nrow = 5, ncol = 8)
  expect_error(decontPro(counts, k, 1e-2, 1e-2, ambient_counts = ambient),
               regexp = "different number of ADTs")
})

test_that("decontPro errors on droplets with zero total counts", {
  # Guard in .decontPro: stops before any Stan call, so the test is fast.
  counts <- matrix(sample(1:10, 80, replace = TRUE), nrow = 10, ncol = 8)
  counts[, 1] <- 0L
  k <- c(1, 1, 2, 2, 3, 3, 4, 4)
  expect_error(decontPro(counts, k, 1e-2, 1e-2),
               regexp = "0 total counts")
})

test_that("decontPro SingleCellExperiment method extracts the counts assay", {
  # Stan is mocked out so this exercises only the SCE method dispatch and
  # assay extraction (fast, no MCMC).
  mat <- matrix(sample(1:10, 80, replace = TRUE), nrow = 10, ncol = 8,
                dimnames = list(paste0("g", seq_len(10)),
                                paste0("c", seq_len(8))))
  k <- c(1, 1, 2, 2, 3, 3, 4, 4)
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = mat))

  received <- NULL
  local_mocked_bindings(
    .decontPro = function(counts, cell_type, ...) {
      received <<- counts
      list(decontaminated_counts = as.matrix(counts))
    }
  )
  decontPro(sce, k, 1e-2, 1e-2)
  expect_equal(as.matrix(received), mat)
})

test_that("decontPro Seurat method uses LayerData not GetAssayData", {
  skip_if_not_installed("SeuratObject")

  set.seed(42)
  mat <- matrix(sample(1:10, 1000, replace = TRUE), ncol = 10,
                dimnames = list(paste0("g", seq_len(100)),
                                paste0("c", seq_len(10))))
  k <- c(1, 1, 2, 2, 2, 3, 3, 4, 4, 4)
  obj <- SeuratObject::CreateSeuratObject(counts = mat)

  received <- NULL
  local_mocked_bindings(
    .decontPro = function(counts, cell_type, ...) {
      received <<- counts
      nc <- ncol(counts)
      nr <- nrow(counts)
      list(
        decontaminated_counts = as.matrix(counts),
        ambient_counts = matrix(0, nr, nc),
        background_counts = matrix(0, nr, nc)
      )
    }
  )

  res <- decontPro(obj, k, 1e-2, 1e-2)
  expect_equal(dim(received), dim(mat))
  expect_equal(as.matrix(received), mat)
})
