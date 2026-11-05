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

# decontX log lines with the parts that change between runs masked: the
# date, the total time and the per-iteration numbers. A "|" marks the end
# of each line, so the snapshots show trailing spaces
maskLog <- function(lines) {
  lines <- sub(
    "^[A-Z][a-z]{2} [A-Z][a-z]{2} [ 0-9][0-9] [0-9:]{8} [0-9]{4}",
    "<date>", lines
  )
  lines <- sub("Total time: [0-9.]+ [a-z]+", "Total time: <time>", lines)
  lines <- sub(
    "iteration: [0-9]+ \\| converge: [0-9.e-]+",
    "iteration: <n> | converge: <x>", lines
  )
  paste0(lines, "|")
}

# The gene profile decontX uses as eta when given background counts
bgProfile <- function(background) {
  profile <- Matrix::rowSums(background) + 1e-20
  profile / sum(profile)
}
