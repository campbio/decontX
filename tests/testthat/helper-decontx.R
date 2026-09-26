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
