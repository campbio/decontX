library(decontX)
context("Testing decontPro function")

test_that(desc = "Testing decontPro on count matrix", {
  set.seed(42)
  counts <- matrix(sample(1:10,
                          1000,
                          replace = TRUE),
                   ncol = 10)

  k <- c(1,1,2,2,2,3,3,4,4,4)

  out <- decontPro(counts, k, 1e-2, 1e-2)

  # Sum decomposed matrices
  matsum = out$decontaminated_counts +
    out$ambient_counts +
    out$background_counts

  expect_equal(matsum, counts)
})
