# Compares today's results with the stored references in
# fixtures/decontx-refs.rds (see helper-references.R). The tolerance allows
# for libm differences between platforms; integer parts match exactly.
refs <- readRDS(test_path("fixtures", "decontx-refs.rds"))
d <- refData()

test_that("the stored references cover every reference path", {
  expect_setequal(names(refs), names(refPaths))
})

for (path in names(refPaths)) {
  test_that(paste("results match the stored reference:", path), {
    res <- refPaths[[path]](d)
    expect_equal(res, refs[[path]], tolerance = 1e-10)

    # Cluster labels, iteration counts and simulated counts are integers
    if (!is.null(res$estimates)) {
      expect_identical(res$z, refs[[path]]$z)
      expect_identical(
        lapply(res$estimates, `[[`, "iteration"),
        lapply(refs[[path]]$estimates, `[[`, "iteration")
      )
    }
    if (!is.null(res$observedCounts)) {
      expect_identical(res$observedCounts, refs[[path]]$observedCounts)
    }
  })
}
