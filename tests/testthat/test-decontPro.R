library(decontX)
context("Testing decontPro function")

test_that(desc = "Testing decontPro on count matrix", {
  set.seed(42)
  counts <- matrix(sample(1:10,
                          1000,
                          replace = TRUE),
                   ncol = 10)

  k <- c(1, 1, 2, 2, 2, 3, 3, 4, 4, 4)

  out <- decontPro(counts, k, 1e-2, 1e-2)

  # Sum decomposed matrices
  matsum <- out$decontaminated_counts +
    out$ambient_counts +
    out$background_counts

  expect_equal(matsum, counts)
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
