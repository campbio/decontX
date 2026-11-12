#' Call Stan variational bayes for inference
#'
#' @param data A list of input data for Stan.
#' @param initial_condition Initial values for Stan params.
#' @param iter Maximum number of iterations.
#' @param seed Seed for Stan's random number generator.
#' @param verbose Whether to print Stan's output. If FALSE, it is captured
#' and shown with message() only when the fit fails.
#'
#' @return Stan output
#' @noRd
.callStanVb <- function(data, initial_condition, iter, seed, verbose) {
  # rstan writes the draws of every parameter to this CSV, about 100 bytes
  # per ADT x droplet, and would keep it until R exits (A5-08)
  sampleFile <- tempfile(fileext = ".csv")
  on.exit(unlink(sampleFile), add = TRUE)
  fitStan <- function() {
    # decontPro only uses the mean of the variational approximation
    # (`@sim$est`), which does not depend on the posterior draws. delta and
    # background have one value per ADT x droplet, so rstan's default of
    # 1000 draws can take hundreds of GB on large data (issue #46). Keep
    # the minimum rstan allows (2) and retain draws for a single scalar
    # parameter only. rstan's Pareto k diagnostic is estimated from the
    # draws; with 2 draws it is always Inf and carries no information, so
    # its warning is dropped.
    withCallingHandlers(
      rstan::vb(
        object = stanmodels$shrinkage,
        init = initial_condition,
        data = data,
        seed = seed,
        iter = iter,
        output_samples = 2,
        pars = "tau_a",
        sample_file = sampleFile
      ),
      warning = function(w) {
        if (startsWith(conditionMessage(w), "Pareto k diagnostic value")) {
          invokeRestart("muffleWarning")
        }
      }
    )
  }
  if (verbose) {
    fit <- fitStan()
  } else {
    stanOutput <- utils::capture.output(
      fit <- tryCatch(fitStan(), error = function(e) e)
    )
  }
  failed <- inherits(fit, "error") ||
    !methods::is(fit, "stanfit") || fit@mode != 0L
  # Stan's output holds the reason for a failure. It goes to message()
  # because an error message is cut at 1000 characters.
  if (failed && !verbose) {
    message(paste(stanOutput, collapse = "\n"))
  }
  if (inherits(fit, "error")) {
    stop(fit)
  }
  # On failure rstan prints the reason and returns a list or an empty
  # stanfit (mode 2). rstan exports no accessor for the mode.
  if (failed) {
    stop(
      "Stan's variational inference failed; Stan printed the reason ",
      if (verbose) "above." else "in the message above."
    )
  }
  fit
}


#' Process Stan output.
#'
#' @param stan_vb_output Stan variational bayes output
#' @param dat List of data input to stan vb
#'
#' @return Decomposed counts based on Stan estimate.
#' @noRd
.processStanVbOut <- function(stan_vb_output, dat) {
  val <- stan_vb_output@sim$est

  r_est <- val$r
  p_est <- dat$p
  delta_est <- t(val$delta)
  delta_mean_est <- val$delta_mean
  background_est <- val$background
  background_mean_est <- val$background_mean


  processed_stan <- list()


  # Save raw estimations
  parameters <- list()
  parameters[["delta_sd"]] <- dat$delta_sd
  parameters[["background_sd"]] <- dat$background_sd
  parameters[["r_est"]] <- r_est
  rownames(parameters[["r_est"]]) <- dat$cell_type_labels
  parameters[["p_est"]] <- p_est
  parameters[["delta_est"]] <- delta_est
  parameters[["delta_mean_est"]] <- delta_mean_est
  parameters[["background_est"]] <- background_est
  parameters[["background_mean_est"]] <- background_mean_est

  processed_stan[["parameters"]] <- parameters


  ## Split each count by its share of the total rate. The three rates of
  ## a droplet all carry its total count OC as a factor, which cancels in
  ## the shares, so it is left out (#111). The profile and rate matrices,
  ## each a full ADT x droplet matrix, are removed once they are used.
  not_background <- 1 - background_est

  #  Ambient: p_est recycles down each droplet's column
  ambient_rate_est <- p_est * delta_est * not_background

  #  Native
  native_profile <- t(r_est[dat$cell_type, , drop = FALSE])
  cell_rate_est <- native_profile * (1 - delta_est) * not_background
  rm(native_profile, not_background)

  #  Background: its rate is background_est itself. counts_per_rate comes
  #  first in each split so the results keep the counts' dimnames.
  counts_per_rate <- dat$counts /
    (ambient_rate_est + cell_rate_est + background_est)
  decontaminated_counts <- counts_per_rate * cell_rate_est
  rm(cell_rate_est)
  ambient_counts <- counts_per_rate * ambient_rate_est
  rm(ambient_rate_est)
  background_counts <- counts_per_rate * background_est


  processed_stan[["decontaminated_counts"]] <- decontaminated_counts
  processed_stan[["ambient_counts"]] <- ambient_counts
  processed_stan[["background_counts"]] <- background_counts

  processed_stan
}
