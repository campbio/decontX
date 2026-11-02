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
  decontPro(counts, k, 1e-2, 1e-2)

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

  out <- decontPro(counts, k, 1e-2, 1e-2)
  # Reference: rstan's default 1000 draws. Keep seed/iter in sync with
  # .call_stan_vb() so only the draw settings differ.
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

test_that("decontPro takes the ambient profile from ambient_counts", {
  # Stan is mocked out: this pins the ambient profile passed to Stan (row
  # sums of ambient_counts, matched by position) for the matrix and SCE
  # methods.
  p <- refProData()
  mat <- p$counts
  ambient <- p$ambient
  k <- p$k
  expected <- rowSums(ambient) / sum(ambient)

  received <- NULL
  local_mocked_bindings(
    .call_stan_vb = function(data, initial_condition) {
      received <<- data
      NULL
    },
    .process_stan_vb_out = function(stan_vb_output, dat) list()
  )

  decontPro(mat, k, 1e-2, 1e-2, ambient_counts = ambient)
  expect_equal(received$p, expected)
  expect_equal(received$counts, mat)

  received <- NULL
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = mat))
  ambientSce <- SingleCellExperiment::SingleCellExperiment(
    list(counts = ambient)
  )
  decontPro(sce, k, 1e-2, 1e-2, ambient_counts = ambientSce)
  expect_equal(received$p, expected)
})

test_that("decontPro matches ambient_counts rows to the counts by name", {
  p <- refProData()
  mat <- p$counts
  ambient <- p$ambient
  k <- p$k
  expected <- rowSums(ambient) / sum(ambient)

  received <- NULL
  local_mocked_bindings(
    .call_stan_vb = function(data, initial_condition) {
      received <<- data
      NULL
    },
    .process_stan_vb_out = function(stan_vb_output, dat) list()
  )

  decontPro(mat, k, 1e-2, 1e-2, ambient_counts = ambient[10:1, ])
  expect_equal(received$p, expected)

  received <- NULL
  extra <- rbind(ambient, extra_ADT = 50L)
  expect_message(
    decontPro(mat, k, 1e-2, 1e-2, ambient_counts = extra),
    regexp = paste0(
      "Dropped ADTs in 'ambient_counts' that aren't in 'filtered_counts': ",
      "1 of 11 \\(extra_ADT\\)"
    )
  )
  expect_equal(received$p, expected)

  expect_error(
    decontPro(mat, k, 1e-2, 1e-2, ambient_counts = ambient[-3, ]),
    regexp = "'ambient_counts' lacks 1 of the 10 ADTs in 'filtered_counts': g3"
  )
})

test_that("decontPro names ambient_counts when it isn't a matrix", {
  p <- refProData()
  expect_error(
    decontPro(p$counts, p$k, 1e-2, 1e-2, ambient_counts = rowSums(p$ambient)),
    regexp = "'ambient_counts' must be a matrix"
  )
})

test_that("decontPro joins split Seurat layers in Cells() order", {
  skip_if_not_installed("SeuratObject")
  p <- refProData()
  sample <- rep(c("a", "b"), 4)
  obj <- SeuratObject::CreateSeuratObject(
    counts = as(p$counts, "CsparseMatrix")
  )
  obj[["RNA"]] <- split(obj[["RNA"]], f = sample)
  ambientObj <- SeuratObject::CreateSeuratObject(
    counts = as(p$ambient, "CsparseMatrix")
  )
  ambientSample <- c("a", "b", "a", "b", "a")
  ambientObj[["RNA"]] <- split(ambientObj[["RNA"]], f = ambientSample)

  received <- NULL
  local_mocked_bindings(
    .decontPro = function(counts, cell_type, delta_sd, background_sd,
                          ambient_counts) {
      received <<- list(counts = counts, ambient = ambient_counts)
      list()
    }
  )
  expect_no_warning(
    decontPro(obj, p$k, 1e-2, 1e-2, ambient_counts = ambientObj)
  )
  expect_identical(colnames(received$counts), SeuratObject::Cells(obj))
  expect_equal(as.matrix(received$counts), p$counts)
  expect_equal(as.matrix(received$ambient), p$ambient)
  # The user's object keeps its layers
  expect_identical(SeuratObject::Layers(obj), c("counts.a", "counts.b"))
})

test_that("decontPro reads a Seurat v3 assay", {
  skip_if_not_installed("SeuratObject")
  withr::local_options(Seurat.object.assay.version = "v3")
  p <- refProData()
  obj <- SeuratObject::CreateSeuratObject(
    counts = as(p$counts, "CsparseMatrix")
  )
  expect_false(methods::is(obj[["RNA"]], "Assay5"))

  received <- NULL
  local_mocked_bindings(
    .decontPro = function(counts, cell_type, ...) {
      received <<- counts
      list()
    }
  )
  decontPro(obj, p$k, 1e-2, 1e-2)
  expect_equal(as.matrix(received), p$counts)
})

test_that("decontPro checks the length of cell_type and NA", {
  p <- refProData()
  expect_error(
    decontPro(p$counts, p$k[-1], 1e-2, 1e-2),
    regexp = "'cell_type' has 7 entries but 'filtered_counts' has 8 droplets"
  )
  expect_error(
    decontPro(p$counts, c(p$k, 1), 1e-2, 1e-2),
    regexp = "'cell_type' has 9 entries but 'filtered_counts' has 8 droplets"
  )
  expect_error(
    decontPro(p$counts, c(p$k[-8], NA), 1e-2, 1e-2),
    regexp = "'cell_type' has NA values"
  )
})

test_that("decontPro recodes cell_type labels to 1..K for Stan", {
  p <- refProData()
  received <- NULL
  local_mocked_bindings(
    .call_stan_vb = function(data, initial_condition) {
      received <<- data
      NULL
    },
    .process_stan_vb_out = function(stan_vb_output, dat) list()
  )
  codes <- c(1, 1, 2, 2, 3, 3, 4, 4)
  labelSets <- list(
    gaps = c(1, 1, 3, 3, 5, 5, 7, 7),
    zeroBased = c(0, 0, 1, 1, 2, 2, 3, 3),
    character = c("B", "B", "M", "M", "NK", "NK", "T", "T"),
    factor = factor(c("T", "T", "B", "B", "NK", "NK", "M", "M"),
      levels = c("T", "B", "unused", "NK", "M")
    ),
    oneRowMatrix = matrix(c(2, 2, 4, 4, 6, 6, 8, 8), nrow = 1)
  )
  for (labels in labelSets) {
    decontPro(p$counts, labels, 1e-2, 1e-2)
    expect_equal(received$cell_type, codes)
    expect_identical(received$K, 4L)
  }
})

test_that("decontPro names the r_est rows by cell_type label", {
  # One real Stan fit per labelling (10 x 8, about a second each)
  p <- refProData()
  out <- decontPro(p$counts, p$k, 1e-2, 1e-2)
  gapped <- decontPro(p$counts, p$k * 2 + 1, 1e-2, 1e-2)

  expect_identical(rownames(out$parameters$r_est), c("1", "2", "3", "4"))
  expect_identical(rownames(gapped$parameters$r_est), c("3", "5", "7", "9"))
  rownames(gapped$parameters$r_est) <- rownames(out$parameters$r_est)
  expect_identical(gapped, out)
})

test_that("decontPro requires whole-number counts", {
  p <- refProData()
  expect_error(
    decontPro(p$counts + 0.5, p$k, 1e-2, 1e-2),
    regexp = "'filtered_counts' must hold whole-number counts"
  )
  negative <- p$counts
  negative[1, 1] <- -1
  expect_error(
    decontPro(negative, p$k, 1e-2, 1e-2),
    regexp = "'filtered_counts' has negative values"
  )
  infinite <- p$counts + 0
  infinite[1, 1] <- Inf
  expect_error(
    decontPro(infinite, p$k, 1e-2, 1e-2),
    regexp = "'filtered_counts' must hold whole-number counts"
  )
  missing <- p$counts
  missing[1, 1] <- NA
  expect_error(
    decontPro(missing, p$k, 1e-2, 1e-2),
    regexp = "'filtered_counts' has NA values"
  )
  expect_error(
    decontPro(p$counts > 5, p$k, 1e-2, 1e-2),
    regexp = "'filtered_counts' must be numeric"
  )
})

test_that("decontPro reports a failed Stan fit", {
  # When rstan::vb() fails it returns a list or an empty stanfit (mode 2)
  p <- refProData()
  local_mocked_bindings(vb = function(...) list(), .package = "rstan")
  expect_error(
    decontPro(p$counts, p$k, 1e-2, 1e-2),
    regexp = "Stan's variational inference failed"
  )
  local_mocked_bindings(
    vb = function(...) methods::new("stanfit", mode = 2L),
    .package = "rstan"
  )
  expect_error(
    decontPro(p$counts, p$k, 1e-2, 1e-2),
    regexp = "Stan's variational inference failed"
  )
})

test_that("decontPro rejects unknown arguments", {
  p <- refProData()
  expect_error(
    decontPro(p$counts, p$k, seed = 1),
    regexp = "Unknown arguments to decontPro\\(\\): 'seed'"
  )
  expect_error(
    decontPro(p$counts, p$k, delta_SD = 1e-2, batch = 1),
    regexp = "Unknown arguments to decontPro\\(\\): 'delta_SD', 'batch'"
  )
  expect_error(
    decontPro(p$counts, p$k, 1e-2, 1e-2, NULL, 5),
    regexp = "Unknown arguments to decontPro\\(\\): an unnamed value"
  )

  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = p$counts))
  expect_error(
    decontPro(sce, p$k, altExpName = "ADT"),
    regexp = "Unknown arguments to decontPro\\(\\): 'altExpName'"
  )
  expect_error(
    decontPro(sce, p$k, 1e-2, 1e-2, NULL, 5),
    regexp = "Unknown arguments to decontPro\\(\\): an unnamed value"
  )

  skip_if_not_installed("SeuratObject")
  obj <- SeuratObject::CreateSeuratObject(
    counts = as(p$counts, "CsparseMatrix")
  )
  expect_error(
    decontPro(obj, p$k, assay = "ADT"),
    regexp = "Unknown arguments to decontPro\\(\\): 'assay'"
  )
  expect_error(
    decontPro(obj, p$k, 1e-2, 1e-2, NULL, 5),
    regexp = "Unknown arguments to decontPro\\(\\): an unnamed value"
  )
})

test_that("decontPro rejects NA and negative values in ambient_counts", {
  # Checked before any Stan call, so no fit runs
  p <- refProData()
  missing <- p$ambient
  missing[2, 3] <- NA
  expect_error(
    decontPro(p$counts, p$k, 1e-2, 1e-2, ambient_counts = missing),
    regexp = "Missing value in 'ambient_counts' matrix"
  )
  negative <- p$ambient
  negative[2, 3] <- -1L
  expect_error(
    decontPro(p$counts, p$k, 1e-2, 1e-2, ambient_counts = negative),
    regexp = "Negative value in 'ambient_counts' matrix"
  )
  expect_error(
    decontPro(p$counts, p$k, 1e-2, 1e-2,
      ambient_counts = as(negative, "CsparseMatrix")
    ),
    regexp = "Negative value in 'ambient_counts' matrix"
  )
})
