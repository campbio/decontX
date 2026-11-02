#' Density of each ADT, raw counts overlapped with decontaminated counts
#'
#' Both curves of an ADT are smoothed with the same bandwidth, computed from
#' the raw and decontaminated counts together, so their heights can be
#' compared.
#'
#' @param counts original count matrix of nADT x nDroplet.
#' @param decontaminated_counts decontaminated count matrix.
#' @param features names of ADT to plot
#' @param file file name to save the plot into as a pdf. ".pdf" is added
#' unless \code{file} already ends in ".pdf". If omitted, return a
#' \code{ggplot} object.
#'
#' @return A \code{ggplot} object, or, if \code{file} is given, the path of
#' the saved pdf, invisibly.
#'
#' @examples
#' # Simulate a dataset with 3 cells and 2 ADTs
#' counts <- matrix(c(60, 72, 52, 49, 89, 112),
#'   nrow = 2,
#'   dimnames = list(
#'     c("CD3", "CD4"),
#'     c(
#'       "CTGTTTACACCGCTAG",
#'       "CTCTACGGTGTGGCTC",
#'       "AGCAGCCAGGCTCATT"
#'     )
#'   )
#' )
#'
#' decontaminated_counts <- matrix(c(58, 36, 26, 45, 88, 110),
#'   nrow = 2,
#'   dimnames = list(
#'     c("CD3", "CD4"),
#'     c(
#'       "CTGTTTACACCGCTAG",
#'       "CTCTACGGTGTGGCTC",
#'       "AGCAGCCAGGCTCATT"
#'     )
#'   )
#' )
#'
#' plotDensity(
#'   counts,
#'   decontaminated_counts,
#'   c("CD3", "CD4")
#' )
#' @export
plotDensity <- function(counts,
                        decontaminated_counts,
                        features,
                        file = NULL) {
  p <- list()

  for (i in seq_along(features)) {
    feature <- features[i]


    df <- data.frame(
      con = counts[feature, ],
      decon = decontaminated_counts[feature, ]
    )
    df.m <- reshape2::melt(df, measure.var = c("con", "decon"))

    # One bandwidth for both curves, on the pseudo-log scale the density is
    # computed on: asinh(x / 2) is "pseudo_log" with its default sigma and
    # base. NAs are left out, as geom_density() drops them
    bw <- stats::bw.nrd0(asinh(stats::na.omit(df.m$value) / 2))

    # Plot
    p1 <- ggplot2::ggplot(
      df.m,
      ggplot2::aes(.data$value, fill = .data$variable)
    ) +
      ggplot2::geom_density(alpha = 0.7, bw = bw) +
      ggplot2::scale_x_continuous(
        transform = "pseudo_log",
        breaks = c(1, 5, 10^seq(1, 4))
      ) +
      ggplot2::scale_fill_manual(
        values = c("#E64B35B2", "#4DBBD5B2"),
        labels = c("Original", "Decontaminated")
      ) +
      ggplot2::ggtitle(feature) +
      ggplot2::labs(x = "", fill = "") +
      ggplot2::theme_classic() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          angle = 45,
          vjust = 1,
          hjust = 0.9
        ),
        legend.margin = ggplot2::margin(t = -10)
      )

    ylimit <- ggplot2::layer_scales(p1)$y$get_limits()
    ylimit[2] <- min(ylimit[2], 5)
    p1 <- p1 + ggplot2::coord_cartesian(ylim = ylimit)

    p[[i]] <- p1
  }

  # Combine plots
  p_wrap <- patchwork::wrap_plots(p,
    ncol = round(sqrt(length(features))),
    guides = "collect"
  ) +
    patchwork::plot_layout() & ggplot2::theme(
    legend.position = "bottom",
    plot.margin = ggplot2::unit(c(3, 3, 2, 1), "pt")
  )


  # Output
  if (is.null(file)) {
    p_wrap
  } else {
    .savePlotPdf(p_wrap, file)
  }
}


#' Boxplot of features grouped by cell type
#'
#' @param counts original count matrix of nADT x nDroplet.
#' @param decontaminated_counts decontaminated count matrix.
#' @param cell_type 1xnDroplet vector of cell_type.
#' @param features names of ADT to plot
#' @param file file name to save the plot into as a pdf. ".pdf" is added
#' unless \code{file} already ends in ".pdf". If omitted, return a
#' \code{ggplot} object.
#'
#' @return A \code{ggplot} object, or, if \code{file} is given, the path of
#' the saved pdf, invisibly.
#'
#' @examples
#' # Simulate a dataset with 3 cells and 2 ADTs
#' counts <- matrix(c(60, 72, 52, 49, 89, 112),
#'   nrow = 2,
#'   dimnames = list(
#'     c("CD3", "CD4"),
#'     c(
#'       "CTGTTTACACCGCTAG",
#'       "CTCTACGGTGTGGCTC",
#'       "AGCAGCCAGGCTCATT"
#'     )
#'   )
#' )
#'
#' decontaminated_counts <- matrix(c(58, 36, 26, 45, 88, 110),
#'   nrow = 2,
#'   dimnames = list(
#'     c("CD3", "CD4"),
#'     c(
#'       "CTGTTTACACCGCTAG",
#'       "CTCTACGGTGTGGCTC",
#'       "AGCAGCCAGGCTCATT"
#'     )
#'   )
#' )
#'
#' plotBoxByCluster(
#'   counts,
#'   decontaminated_counts,
#'   c(1, 2, 1),
#'   c("CD3", "CD4")
#' )
#' @export
plotBoxByCluster <- function(counts,
                             decontaminated_counts,
                             cell_type,
                             features,
                             file = NULL) {
  p <- list()

  for (i in seq_along(features)) {
    feature <- features[i]


    df <- data.frame(
      con = counts[feature, ],
      decon = decontaminated_counts[feature, ],
      cell_type = as.factor(cell_type)
    )

    df.m <- reshape2::melt(df, measure.var = c("con", "decon"))


    # Plot
    p1 <- ggplot2::ggplot(
      df.m,
      ggplot2::aes(
        x = .data$cell_type,
        y = .data$value,
        fill = .data$variable
      )
    ) +
      ggplot2::geom_boxplot(
        lwd = 0.2,
        outlier.size = 0.5,
        alpha = 0.7
      ) +
      ggplot2::scale_y_continuous(
        transform = "pseudo_log",
        breaks = c(1, 5, 10^seq(1, 4))
      ) +
      ggplot2::scale_fill_manual(
        values = c("#E64B35B2", "#4DBBD5B2"),
        labels = c("Original", "Decontaminated")
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        legend.position = "bottom",
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(
          angle = 45,
          vjust = 1,
          hjust = 0.9
        )
      ) +
      ggplot2::labs(x = "", y = "", fill = "") +
      ggplot2::ggtitle(paste0(feature, " counts by cell clusters"))


    p[[i]] <- p1
  }

  # Combine plots
  p_wrap <- patchwork::wrap_plots(p,
    ncol = round(sqrt(length(features))),
    guides = "collect"
  ) +
    patchwork::plot_layout() & ggplot2::theme(
    legend.position = "bottom",
    plot.margin =
      ggplot2::unit(
        c(3, 3, 2, 1),
        "pt"
      )
  )

  if (is.null(file)) {
    p_wrap
  } else {
    .savePlotPdf(p_wrap, file)
  }
}


## Save a plot as a letter-size pdf, adding ".pdf" unless 'file' has it,
## and return the path invisibly
.savePlotPdf <- function(p, file) {
  if (!grepl("\\.pdf$", file, ignore.case = TRUE)) {
    file <- paste0(file, ".pdf")
  }
  ggplot2::ggsave(file,
    p,
    width = 8.5,
    height = 11,
    units = "in"
  )
  invisible(file)
}
