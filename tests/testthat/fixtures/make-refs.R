# Regenerates decontx-refs.rds, the reference outputs that
# test-references.R compares against. Run from the package root:
#   Rscript tests/testthat/fixtures/make-refs.R
# Rerun it only when a change alters results on purpose, and say in the PR
# which references changed and why.
pkgload::load_all(".", quiet = TRUE)
source(file.path("tests", "testthat", "helper-references.R"))

d <- refData()
refs <- lapply(refPaths, function(path) path(d))
saveRDS(refs, file.path("tests", "testthat", "fixtures", "decontx-refs.rds"),
  compress = "xz"
)
