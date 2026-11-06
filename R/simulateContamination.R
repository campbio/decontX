#' @title Simulate contaminated count matrix
#' @description This function generates a list containing two count matrices --
#'  one for real expression, the other one for contamination, as well as other
#'  parameters used in the simulation which can be useful for running
#'  decontamination.
#' @param C Integer. Number of cells to be simulated. Default \code{300}.
#' @param G Integer. Number of genes to be simulated. Default \code{100}.
#' @param K Integer. Number of cell populations to be simulated, at least 2,
#' since the contamination comes from the other populations. With few cells,
#' fewer populations may be drawn; at least 2 are still needed.
#' Default \code{3}.
#' @param NRange Integer vector. A vector of length 2 that specifies the lower
#'  and upper bounds of the number of counts generated for each cell. Default
#'  \code{c(500, 1000)}.
#' @param beta Numeric. Concentration parameter for Phi. Default \code{0.1}.
#' @param delta Numeric or Numeric vector. Concentration parameter for Theta.
#'  If input as a single numeric value, symmetric values for beta
#'  distribution are specified; if input as a vector of length 2, the two
#'  values will be the shape1 and shape2 parameters of the beta distribution
#'  respectively. Default \code{c(1, 5)}.
#' @param numMarkers Integer. Number of markers for each cell population.
#' Default \code{3}.
#' @param seed Integer. Passed to \code{\link[withr]{with_seed}}.
#' For reproducibility, a default value of 12345 is used. If NULL, no calls to
#'  \code{\link[withr]{with_seed}} are made.
#' @return A list containing the \code{nativeCounts} (real expression),
#' \code{observedCounts} (real expression + contamination), as well as other
#' parameters used in the simulation.
#' @author Shiyi Yang, Yuan Yin, Joshua Campbell
#' @examples
#' contaminationSim <- simulateContamination(K = 3, delta = c(1, 10))
#' @export
simulateContamination <- function(C = 300,
                                  G = 100,
                                  K = 3,
                                  NRange = c(500, 1000),
                                  beta = 0.1,
                                  delta = c(1, 10),
                                  numMarkers = 3,
                                  seed = 12345) {
  .withSeed(
    seed,
    .simulateContaminatedMatrix(
      C = C,
      G = G,
      K = K,
      NRange = NRange,
      beta = beta,
      delta = delta,
      numMarkers = numMarkers
    )
  )
}


.simulateContaminatedMatrix <- function(C = 300,
                                        G = 100,
                                        K = 3,
                                        NRange = c(500, 1000),
                                        beta = 0.5,
                                        delta = c(1, 2),
                                        numMarkers = 3) {
  # The contamination of each population comes from the other populations
  if (K < 2) {
    stop(
      "'K' must be at least 2, since contamination comes from the other ",
      "cell populations."
    )
  }
  if (length(delta) == 1) {
    cpByC <- stats::rbeta(
      n = C,
      shape1 = delta,
      shape2 = delta
    )
  } else {
    cpByC <- stats::rbeta(
      n = C,
      shape1 = delta[1],
      shape2 = delta[2]
    )
  }

  z <- sample(seq(K), size = C, replace = TRUE)
  if (length(unique(z)) == 1) {
    stop(
      "Only 1 cell population was simulated. Increase the number of ",
      "cells 'C'."
    )
  }
  if (length(unique(z)) < K) {
    warning(
      "Only ",
      length(unique(z)),
      " clusters are simulated. Try to increase number of cells 'C' if",
      " more clusters are needed"
    )
    K <- length(unique(z))
    z <- plyr::mapvalues(z, unique(z), seq_along(unique(z)))
  }

  # sample() on a single value would draw from 1 to that value instead
  nValues <- seq(min(NRange), max(NRange))
  NbyC <- nValues[sample.int(length(nValues), size = C, replace = TRUE)]
  cNbyC <- vapply(seq(C), function(i) {
    stats::rbinom(
      n = 1,
      size = NbyC[i],
      p = cpByC[i]
    )
  }, integer(1))
  rNbyC <- NbyC - cNbyC

  phi <- .rdirichlet(K, rep(beta, G))

  ## Select random genes to be markers in each cell population
  ## by setting their values to zero.
  if (K * numMarkers > G) {
    stop(
      "The number of markers ('numMarkers') times the number of cell",
      " populations ('K') cannot be greater than the number of",
      " genes ('G')."
    )
  }
  markerKIndex <- rep(seq(K), each = numMarkers)
  markerRowIndex <- sample(seq(G), numMarkers * K)
  for (i in seq(K)) {
    ix <- markerRowIndex[markerKIndex == i]
    phi[i, ix] <- max(phi[i, ])
    for (j in setdiff(seq(K), i)) {
      phi[j, ix] <- 0
    }
  }
  phi <- prop.table(phi, margin = 1)

  ## sample real expressed count matrix
  cellRmat <- vapply(seq(C), function(i) {
    stats::rmultinom(1, size = rNbyC[i], prob = phi[z[i], ])
  }, integer(G))

  rownames(cellRmat) <- paste0("Gene_", seq(G))
  colnames(cellRmat) <- paste0("Cell_", seq(C))

  ## Get list of marker names
  markerNames <- list()
  for (i in seq(K)) {
    markerNames[[i]] <- rownames(cellRmat)[markerRowIndex[markerKIndex == i]]
  }
  names(markerNames) <- paste0("CellType_", seq(K), "_Markers")

  ## sample contamination count matrix
  nGByK <-
    rowSums(cellRmat) - .colSumByGroup(cellRmat, group = z, K = K)
  ## Column-wise proportion normalization. This reproduces
  ## celda::normalizeCounts(normalize = "proportion") exactly (verified to
  ## the bit on the simulation path) via the package's own fastNormProp
  ## C++ routine with a zero pseudocount, removing the celda dependency
  ## (ADR-0003 Stage 1).
  eta <- fastNormProp(nGByK, 0)

  cellCmat <- vapply(seq(C), function(i) {
    stats::rmultinom(1, size = cNbyC[i], prob = eta[, z[i]])
  }, integer(G))
  cellOmat <- cellRmat + cellCmat
  contamination <- colSums(cellCmat) / colSums(cellOmat)

  rownames(cellOmat) <- paste0("Gene_", seq(G))
  colnames(cellOmat) <- paste0("Cell_", seq(C))

  list(
    "nativeCounts" = cellRmat,
    "observedCounts" = cellOmat,
    "NByC" = NbyC,
    "z" = z,
    "eta" = eta,
    "phi" = t(phi),
    "markers" = markerNames,
    "numMarkers" = numMarkers,
    "contamination" = contamination
  )
}
