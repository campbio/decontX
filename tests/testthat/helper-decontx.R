# Shared fixtures/helpers for decontX tests.
# testthat automatically sources helper-*.R before running the test files.

# A deterministic simulated dataset with known ground truth, used for the
# behavioral (oracle) assertions in test-decon.R: recovery of the true
# per-cell contamination fraction, the decontaminated <= observed invariant,
# and EM log-likelihood monotonicity. Strong contamination (delta = c(1, 10))
# and well-separated clusters give a clear signal to recover.
simulate_oracle <- function(seed = 12345,
                            C = 300,
                            G = 100,
                            K = 3,
                            delta = c(1, 10)) {
  simulateContamination(C = C, G = G, K = K, delta = delta, seed = seed)
}

# Runs a real decontPro Stan fit without printing rstan's progress. If the
# fit fails, Stan's output is printed, since it holds the reason
quietStan <- function(expr) {
  output <- withr::local_tempfile()
  withr::with_output_sink(output, {
    res <- tryCatch(expr, error = function(e) e)
  })
  if (inherits(res, "error")) {
    writeLines(readLines(output))
    stop(res)
  }
  res
}

# The arguments of a decontX() method, with their defaults, including those
# the methods package moves into an inner .local function
decontXMethodFormals <- function(signature) {
  method <- methods::getMethod("decontX", signature)
  formals(methods::unRematchDefinition(method))
}

# The decontX parameters with the ANY method's defaults, for tests that
# call the internal functions directly
defaultParams <- function() {
  defaults <- decontXMethodFormals("ANY")[names(formals(.decontXParams))]
  do.call(.decontXParams, lapply(defaults, eval))
}

# The gene profile decontX uses as eta when given background counts
bgProfile <- function(background) {
  profile <- Matrix::rowSums(background) + 1e-20
  profile / sum(profile)
}
