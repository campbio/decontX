# Compares today's results with the stored references in
# fixtures/decontx-refs.rds (see helper-references.R). The tolerance allows
# for libm differences between platforms; integer parts match exactly.
test_that("the stored references cover every reference path", {
  refs <- readRDS(test_path("fixtures", "decontx-refs.rds"))
  expect_setequal(names(refs), names(refPaths))
})

for (path in names(refPaths)) {
  test_that(paste("results match the stored reference:", path), {
    refs <- readRDS(test_path("fixtures", "decontx-refs.rds"))
    ref <- refs[[path]]
    res <- refPaths[[path]](refData())
    expect_equal(res, ref, tolerance = 1e-10)

    # Integer parts (cluster labels, iteration counts, simulated counts)
    # must match exactly, type included
    if (is.list(ref) && !is.data.frame(ref)) {
      isInteger <- vapply(ref, is.integer, logical(1))
      expect_identical(res[isInteger], ref[isInteger])
    }
    if (!is.null(ref$estimates)) {
      expect_identical(
        lapply(res$estimates, `[[`, "iteration"),
        lapply(ref$estimates, `[[`, "iteration")
      )
    }
  })
}
