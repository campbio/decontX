# Reference paths pinned by test-references.R. fixtures/make-refs.R stores
# their output in fixtures/decontx-refs.rds, and dev/compare-commits.R runs
# them on two installed commits. Each path takes the shared data from
# refData() and returns only the parts of the result that are pinned.

# Simulated counts, a second simulation used as the background, and the
# batch labels shared by every path
refData <- function() {
  sim <- simulateContamination(C = 100, G = 60, K = 3, seed = 1)
  bg <- simulateContamination(C = 40, G = 60, K = 3, seed = 2)$observedCounts
  rownames(bg) <- rownames(sim$observedCounts)
  colnames(bg) <- paste0("bg_", seq_len(ncol(bg)))
  list(
    sim = sim,
    bg = bg,
    batch = rep(c("A", "B"), each = 50),
    bgBatch = rep(c("A", "B"), each = 20)
  )
}

# The 10 x 8 decontPro fixture, with an ambient matrix for the same ADTs
refProData <- function() {
  withr::with_seed(42, {
    counts <- matrix(sample(1:10, 80, replace = TRUE), nrow = 10, ncol = 8,
      dimnames = list(paste0("g", seq_len(10)), paste0("c", seq_len(8)))
    )
    ambient <- matrix(sample(1:100, 50, replace = TRUE), nrow = 10, ncol = 5,
      dimnames = list(paste0("g", seq_len(10)), paste0("e", seq_len(5)))
    )
  })
  list(counts = counts, ambient = ambient, k = c(1, 1, 2, 2, 3, 3, 4, 4))
}

# The UMAP is left out: with z supplied it doesn't feed the EM, and it
# depends on the platform
refKeep <- function(res) {
  estimates <- lapply(res$estimates, function(e) {
    e[c("theta", "phi", "eta", "delta", "logLikelihood", "iteration")]
  })
  list(
    contamination = res$contamination,
    z = res$z,
    decontXcounts = res$decontXcounts,
    estimates = estimates
  )
}

refRun <- function(d, ...) {
  res <- decontX(d$sim$observedCounts, z = d$sim$z, verbose = FALSE, ...)
  refKeep(res)
}

refPaths <- list(
  single = function(d) refRun(d),
  batched = function(d) refRun(d, batch = d$batch),
  background = function(d) refRun(d, background = d$bg),
  bgBatched = function(d) {
    refRun(d, batch = d$batch, background = d$bg, bgBatch = d$bgBatch)
  },
  fixedDelta = function(d) refRun(d, estimateDelta = FALSE),
  simulateDefault = function(d) simulateContamination(),
  simulateSeed1 = function(d) d$sim,
  markerPercent = function(d) {
    p <- plotDecontXMarkerPercentage(d$sim$observedCounts,
      markers = d$sim$markers, z = d$sim$z
    )
    p$data
  }
)
