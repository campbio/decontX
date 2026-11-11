test_that("decontPro splits the counts from one compact Stan fit", {
  # One real Stan VB run on the 10 x 8 fixture (about a second)
  p <- refProData()
  fit <- NULL
  process <- .processStanVbOut
  local_mocked_bindings(.processStanVbOut = function(stan_vb_output, dat) {
    fit <<- stan_vb_output
    process(stan_vb_output, dat)
  })
  # rstan's draw-based Pareto k warning isn't surfaced
  expect_no_warning(
    out <- quietStan(decontPro(p$counts, p$k, 1e-2, 1e-2)),
    message = "Pareto k"
  )

  # The three parts add up to the counts
  matsum <- out$decontaminated_counts +
    out$ambient_counts +
    out$background_counts
  expect_equal(matsum, p$counts)

  # Issue #46: rstan's vb defaults to 1000 draws of every parameter, and
  # delta/background have one value per ADT x droplet, so 215 x 52222 data
  # produced a ~215 GB temp CSV. decontPro only uses the approximation
  # mean, so fewer values are kept than there are ADT x droplet entries
  expect_lt(length(unlist(fit@sim$samples)), length(p$counts))
  # The mean estimate for every parameter is still available
  expect_equal(dim(fit@sim$est$delta), c(8, 10))
  expect_equal(dim(fit@sim$est$background), c(10, 8))
})

test_that("decontPro splits each count by the estimated rates", {
  # A hand-built fit for 2 ADTs x 3 droplets, without Stan. Droplet 1 has
  # no ambient (delta 0), droplet 2 is all ambient (delta 1) and droplet 3
  # is half and half; ADT 2 in droplet 1 has a background of 0.5
  counts <- matrix(c(6, 12, 8, 4, 12, 16), nrow = 2)
  dat <- list(
    N = 2, M = 3, counts = counts, OC = c(18, 12, 28),
    p = c(0.25, 0.75), cell_type = c(1, 2, 2),
    cell_type_labels = c("B", "T"), delta_sd = 0.1, background_sd = 0.2
  )
  est <- list(
    r = rbind(c(0.5, 0.5), c(0.75, 0.25)),
    # Stan gives delta as droplets x ADTs
    delta = rbind(c(0, 0), c(1, 1), c(0.5, 0.5)),
    delta_mean = 0.3,
    background = matrix(c(0, 0.5, 0, 0, 0, 0), nrow = 2),
    background_mean = 0.1
  )
  fit <- methods::new("stanfit", sim = list(est = est))
  out <- .processStanVbOut(fit, dat)

  # Droplet 1: ADT 1 all native; ADT 2 native 0.25 vs background 0.5.
  # Droplet 3: ambient p vs native r for cell type 2
  expect_equal(out$decontaminated_counts, matrix(c(6, 4, 0, 0, 9, 4), 2))
  expect_equal(out$ambient_counts, matrix(c(0, 0, 8, 4, 3, 12), 2))
  expect_equal(out$background_counts, matrix(c(0, 8, 0, 0, 0, 0), 2))

  params <- out$parameters
  expect_named(params, c(
    "delta_sd", "background_sd", "r_est", "p_est", "delta_est",
    "delta_mean_est", "background_est", "background_mean_est"
  ))
  expect_identical(rownames(params$r_est), c("B", "T"))
  expect_identical(params$delta_est, t(est$delta))
  expect_identical(params$p_est, dat$p)
  expect_identical(params$delta_sd, 0.1)
  expect_identical(params$background_sd, 0.2)
  expect_identical(params$background_est, est$background)
  expect_identical(params$delta_mean_est, 0.3)
  expect_identical(params$background_mean_est, 0.1)
})

test_that("decontPro results match a Stan fit with default output draws", {
  p <- refProData()
  out <- quietStan(decontPro(p$counts, p$k, 1e-2, 1e-2))
  # Reference: rstan's default 1000 draws. Keep seed/iter in sync with
  # .callStanVb() so only the draw settings differ.
  local_mocked_bindings(.callStanVb = function(data, initial_condition) {
    rstan::vb(stanmodels$shrinkage,
      init = initial_condition, data = data,
      seed = 12345, iter = 50000
    )
  })
  reference <- suppressWarnings(
    quietStan(decontPro(p$counts, p$k, 1e-2, 1e-2))
  )

  expect_identical(out, reference)
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
  withr::local_seed(42)
  mat <- matrix(sample(1:10, 1000, replace = TRUE), ncol = 10,
                dimnames = list(paste0("g", seq_len(100)),
                                paste0("c", seq_len(10))))
  k <- c(1, 1, 2, 2, 2, 3, 3, 4, 4, 4)
  obj <- SeuratObject::CreateSeuratObject(
    counts = as(mat, "CsparseMatrix")
  )

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

  decontPro(obj, k, 1e-2, 1e-2)
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
    .callStanVb = function(data, initial_condition) {
      received <<- data
      NULL
    },
    .processStanVbOut = function(stan_vb_output, dat) list()
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
    .callStanVb = function(data, initial_condition) {
      received <<- data
      NULL
    },
    .processStanVbOut = function(stan_vb_output, dat) list()
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

test_that("decontPro needs at least 2 ADTs and 2 droplets", {
  # rstan passes length-1 data to Stan as a scalar, so a single ADT or
  # droplet used to fail inside Stan with an unclear message (#110)
  local_mocked_bindings(.callStanVb = function(data, initial_condition) {
    stop("reached Stan")
  })
  oneAdt <- matrix(20, 1, 8, dimnames = list("a1", paste0("d", 1:8)))
  expect_error(
    decontPro(oneAdt, rep(1:2, 4), 1e-2, 1e-2),
    regexp = paste(
      "'filtered_counts' has 1 ADT and 8 droplets;",
      "decontPro needs at least 2 of each"
    )
  )
  oneDroplet <- matrix(20, 10, 1, dimnames = list(paste0("a", 1:10), "d1"))
  expect_error(
    decontPro(oneDroplet, 1, 1e-2, 1e-2),
    regexp = "'filtered_counts' has 10 ADTs and 1 droplet;"
  )
  twoByTwo <- matrix(20, 2, 2, dimnames = list(c("a1", "a2"), c("d1", "d2")))
  expect_error(decontPro(twoByTwo, 1:2, 1e-2, 1e-2), regexp = "reached Stan")
})

test_that("decontPro recodes cell_type labels to 1..K for Stan", {
  p <- refProData()
  received <- NULL
  local_mocked_bindings(
    .callStanVb = function(data, initial_condition) {
      received <<- data
      NULL
    },
    .processStanVbOut = function(stan_vb_output, dat) list()
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
  out <- quietStan(decontPro(p$counts, p$k, 1e-2, 1e-2))
  gapped <- quietStan(decontPro(p$counts, p$k * 2 + 1, 1e-2, 1e-2))

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

test_that("plotDensity smooths both curves of a panel with one bandwidth", {
  counts <- withr::with_seed(1, matrix(rnbinom(200, mu = 50, size = 2),
    nrow = 2, dimnames = list(c("A", "B"), paste0("c", seq_len(100)))
  ))
  # Most droplets drop to zero, so the two curves differ in spread
  decontaminated <- counts
  decontaminated[, 1:60] <- 0
  p <- plotDensity(counts, decontaminated, c("A", "B"))
  for (i in 1:2) {
    values <- c(counts[i, ], decontaminated[i, ])
    # The bandwidth is computed on the scale the density is drawn on
    scale <- ggplot2::layer_scales(p[[i]])$x$get_transformation()
    expect_equal(scale$transform(values), asinh(values / 2))
    expect_equal(
      p[[i]]$layers[[1]]$stat_params$bw,
      stats::bw.nrd0(asinh(values / 2))
    )
  }

  # NA and Inf are left out, as geom_density() drops them; with fewer
  # than 2 values left, the default bandwidth is used. ggplot2 warns about
  # the rows it drops
  decontaminated[1, 1] <- NA
  decontaminated[1, 2] <- Inf
  p <- suppressWarnings(plotDensity(counts, decontaminated, "A"))
  values <- c(counts[1, ], decontaminated[1, -(1:2)])
  expect_equal(
    p[[1]]$layers[[1]]$stat_params$bw,
    stats::bw.nrd0(asinh(values / 2))
  )
  p <- suppressWarnings(plotDensity(counts[, 1, drop = FALSE],
    decontaminated[, 1, drop = FALSE], "A"
  ))
  expect_identical(p[[1]]$layers[[1]]$stat_params$bw, "nrd0")
})

test_that("plotDensity and plotBoxByCluster save to the file given", {
  counts <- matrix(c(60, 72, 52, 49, 89, 112),
    nrow = 2, dimnames = list(c("CD3", "CD4"), paste0("c", 1:3))
  )
  decontaminated <- counts - 1
  plots <- list(
    density = function(file) {
      plotDensity(counts, decontaminated, "CD3", file = file)
    },
    box = function(file) {
      plotBoxByCluster(counts, decontaminated, c(1, 2, 1), "CD3", file = file)
    }
  )
  dir <- withr::local_tempdir()
  for (name in names(plots)) {
    # A name ending in ".pdf", in any case, is used as given
    for (ext in c(".pdf", ".PDF")) {
      path <- file.path(dir, paste0(name, ext))
      res <- withVisible(plots[[name]](path))
      expect_false(res$visible)
      expect_identical(res$value, path)
      expect_true(file.exists(path))
      expect_false(file.exists(paste0(path, ".pdf")))
    }
    # Otherwise ".pdf" is added
    path <- file.path(dir, paste0(name, "_plot"))
    res <- withVisible(plots[[name]](path))
    expect_false(res$visible)
    expect_identical(res$value, paste0(path, ".pdf"))
    expect_true(file.exists(paste0(path, ".pdf")))
  }
})

test_that("plotBoxByCluster returns a box plot per feature", {
  counts <- matrix(c(60, 72, 52, 49, 89, 112),
    nrow = 2, dimnames = list(c("CD3", "CD4"), paste0("c", 1:3))
  )
  decontaminated <- counts - 1
  p <- plotBoxByCluster(counts, decontaminated, c(1, 2, 1), c("CD4", "CD3"))
  expect_s3_class(p, "patchwork")
  expect_length(p, 2)
  cd3 <- p[[2]]
  expect_identical(cd3$labels$title, "CD3 counts by cell clusters")
  expect_equal(
    cd3$data$value,
    unname(c(counts["CD3", ], decontaminated["CD3", ]))
  )
  expect_identical(
    as.character(cd3$data$variable),
    rep(c("con", "decon"), each = 3)
  )
  expect_identical(cd3$data$cell_type, factor(c(1, 2, 1, 1, 2, 1)))
})

test_that("plotDensity and plotBoxByCluster share scales and layout", {
  # A is concentrated, so its density peaks above the cap of 5; B and C
  # are spread out
  counts <- matrix(
    c(
      10, 10, 10, 10, 10, 12,
      0, 5, 20, 100, 300, 1000,
      1, 2, 3, 50, 60, 70
    ),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("A", "B", "C"), paste0("c", 1:6))
  )
  decontaminated <- counts
  decontaminated["A", ] <- 8
  decontaminated[c("B", "C"), ] <- counts[c("B", "C"), ] / 2
  cellType <- c(1, 1, 2, 2, 3, 3)
  breaks <- c(1, 5, 10, 100, 1000, 10000)
  plots <- list(
    density = plotDensity(counts, decontaminated, c("A", "B", "C")),
    box = plotBoxByCluster(counts, decontaminated, cellType, c("A", "B", "C"))
  )
  axis <- c(density = "x", box = "y")
  for (name in names(plots)) {
    p <- plots[[name]]
    # Three panels in round(sqrt(3)) = 2 columns, one legend for all
    expect_s3_class(p, "patchwork")
    expect_length(p, 3)
    expect_identical(p$patches$layout$ncol, 2)
    expect_identical(p$patches$layout$guides, "collect")
    for (i in 1:3) {
      expect_identical(p[[i]]$theme$legend.position, "bottom")
      expect_identical(
        p[[i]]$theme$plot.margin,
        ggplot2::unit(c(3, 3, 2, 1), "pt")
      )
      # Pseudo-log scale on the counts axis
      scale <- ggplot2::layer_scales(p[[i]])[[axis[[name]]]]
      expect_identical(scale$breaks, breaks)
      expect_equal(
        scale$get_transformation()$transform(breaks),
        asinh(breaks / 2)
      )
      # One fill colour and label for each of the two counts
      fill <- p[[i]]$scales$get_scales("fill")
      expect_identical(fill$labels, c("Original", "Decontaminated"))
      built <- ggplot2::layer_data(p[[i]])
      expect_identical(
        unique(built$fill),
        c("#E64B35B2", "#4DBBD5B2")
      )
    }
  }

  # The density's y axis stops at 5, or at the top of the curves if lower
  density <- plots$density
  expect_identical(density[[1]]$coordinates$limits$y[2], 5)
  for (i in 2:3) {
    top <- max(ggplot2::layer_data(density[[i]])$y)
    expect_lt(top, 5)
    expect_equal(density[[i]]$coordinates$limits$y, c(0, top))
  }

  # Box statistics for B, on the pseudo-log scale
  box <- ggplot2::layer_data(plots$box[[2]])
  values <- asinh(c(counts["B", ], decontaminated["B", ]) / 2)
  type <- rep(cellType, 2)
  variable <- rep(c("con", "decon"), each = 6)
  expected <- vapply(split(values, list(variable, type)), stats::median, 1)
  expect_equal(box$middle, unname(expected))
})

test_that("plotDensity and plotBoxByCluster stack both counts", {
  counts <- matrix(c(60L, 72L, 52L, 49L, 89L, 112L),
    nrow = 2, dimnames = list(c("CD3", "CD4"), paste0("c", 1:3))
  )
  decontaminated <- counts - 0.5
  variable <- factor(rep(c("con", "decon"), each = 3),
    levels = c("con", "decon")
  )
  p <- plotDensity(counts, decontaminated, "CD4")
  expect_identical(p[[1]]$data, data.frame(
    variable = variable,
    value = c(72, 49, 112, 71.5, 48.5, 111.5)
  ))
  p <- plotBoxByCluster(counts, decontaminated, c("b", "a", "b"), "CD3")
  expect_identical(p[[1]]$data, data.frame(
    cell_type = factor(rep(c("b", "a", "b"), 2)),
    variable = variable,
    value = c(60, 52, 89, 59.5, 51.5, 88.5)
  ))
})

test_that("plotDensity and plotBoxByCluster plot a factor feature by name", {
  counts <- matrix(c(60, 72, 52, 49, 89, 112),
    nrow = 2, dimnames = list(c("CD3", "CD4"), paste0("c", 1:3))
  )
  decontaminated <- counts - 1
  # factor("CD4") has the integer code 1, which is the row of CD3
  feature <- factor("CD4")
  expected <- unname(c(counts["CD4", ], decontaminated["CD4", ]))
  p <- plotDensity(counts, decontaminated, feature)
  expect_equal(p[[1]]$data$value, expected)
  expect_identical(p[[1]]$labels$title, "CD4")
  p <- plotBoxByCluster(counts, decontaminated, c(1, 2, 1), feature)
  expect_equal(p[[1]]$data$value, expected)
  expect_identical(p[[1]]$labels$title, "CD4 counts by cell clusters")
})

test_that("plotBoxByCluster stops when cell_type isn't one per droplet", {
  counts <- matrix(seq_len(12),
    nrow = 2, dimnames = list(c("CD3", "CD4"), paste0("c", 1:6))
  )
  decontaminated <- counts - 1
  # Length 3 divides the 12 stacked rows, so it used to be recycled
  # silently; length 4 doesn't
  for (cellType in list(c(1, 2, 3), c(1, 2, 3, 4))) {
    expect_error(
      plotBoxByCluster(counts, decontaminated, cellType, "CD3"),
      regexp = paste0(
        "'cell_type' must have one entry per droplet \\(6\\), ",
        "not ", length(cellType)
      )
    )
  }
})

test_that("plotDensity and plotBoxByCluster stop for empty or unknown names", {
  counts <- matrix(c(60, 72, 52, 49, 89, 112),
    nrow = 2, dimnames = list(c("CD3", "CD4"), paste0("c", 1:3))
  )
  decontaminated <- counts - 1
  plots <- list(
    density = function(features) {
      plotDensity(counts, decontaminated, features)
    },
    box = function(features) {
      plotBoxByCluster(counts, decontaminated, c(1, 2, 1), features)
    }
  )
  for (plot in plots) {
    expect_error(
      plot(character(0)),
      regexp = "'features' must name at least one feature"
    )
    expect_error(
      plot(c("CD3", "typo", "CD8")),
      regexp = paste0(
        "'features' not found in the rownames of 'counts' or ",
        "'decontaminated_counts': typo, CD8"
      )
    )
  }
  # Features are matched by name, so row numbers aren't accepted
  expect_error(
    plotDensity(counts, decontaminated, 2),
    regexp = "'decontaminated_counts': 2$"
  )
  # A feature missing from only one of the two matrices is reported too
  expect_error(
    plotDensity(counts, decontaminated[1, , drop = FALSE], "CD4"),
    regexp = "'decontaminated_counts': CD4"
  )
})
