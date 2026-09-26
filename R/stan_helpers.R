#' Call Stan variational bayes for inference
#'
#' @param data A list of input data for Stan.
#' @param initial_condition Initial values for Stan params.
#'
#' @return Stan output
.call_stan_vb <- function(data, initial_condition) {
  # decontPro only uses the mean of the variational approximation
  # (`@sim$est`), which does not depend on the posterior draws. delta and
  # background have one value per ADT x droplet, so rstan's default of 1000
  # draws can take hundreds of GB on large data (issue #46). Keep the minimum
  # rstan allows (2) and retain draws for a single scalar parameter only.
  # rstan's Pareto k diagnostic is estimated from the draws; with 2 draws it
  # is always Inf and carries no information, so its warning is dropped.
  withCallingHandlers(
    rstan::vb(
      object = stanmodels$shrinkage,
      init = initial_condition,
      data = data,
      seed = 12345,
      iter = 50000,
      output_samples = 2,
      pars = "tau_a"
    ),
    warning = function(w) {
      if (startsWith(conditionMessage(w), "Pareto k diagnostic value")) {
        invokeRestart("muffleWarning")
      }
    }
  )
}


#' Process Stan output.
#'
#' @param stan_vb_output Stan variational bayes output
#' @param dat List of data input to stan vb
#'
#' @return Decomposed counts based on Stan estimate.
.process_stan_vb_out <- function(stan_vb_output, dat) {
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
  parameters[["p_est"]] <- p_est
  parameters[["delta_est"]] <- delta_est
  parameters[["delta_mean_est"]] <- delta_mean_est
  parameters[["background_est"]] <- background_est
  parameters[["background_mean_est"]] <- background_mean_est

  processed_stan[["parameters"]] <- parameters


  ## Decontaminated matrix
  #  Ambient
  unscaled_rates <- matrix(dat$p,
    nrow = dat$N,
    ncol = dat$M
  )

  scaling_factor <- matrix(dat$OC,
    nrow = dat$N,
    ncol = dat$M,
    byrow = TRUE
  )

  ambient_rate_est <- scaling_factor *
    unscaled_rates *
    delta_est *
    (1 - background_est)

  #  Native
  unscaled_rates <- t(r_est[dat$cell_type, ])

  cell_rate_est <- scaling_factor *
    unscaled_rates *
    (1 - delta_est) *
    (1 - background_est)

  #  Background
  background_est <- val$background

  background_rate_est <- scaling_factor *
    background_est


  # Decontaminated counts
  counts <- dat$counts
  decontaminated_counts <- counts * cell_rate_est /
    (ambient_rate_est + cell_rate_est + background_rate_est)

  ambient_counts <- counts * ambient_rate_est /
    (ambient_rate_est + cell_rate_est + background_rate_est)

  background_counts <- counts * background_rate_est /
    (ambient_rate_est + cell_rate_est + background_rate_est)


  processed_stan[["decontaminated_counts"]] <- decontaminated_counts
  processed_stan[["ambient_counts"]] <- ambient_counts
  processed_stan[["background_counts"]] <- background_counts

  processed_stan
}
