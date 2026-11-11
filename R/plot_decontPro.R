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
#' unless \code{file} already ends in ".pdf". If omitted, return the plot.
#'
#' @return A \code{patchwork} object with one \code{ggplot} per feature, or,
#' if \code{file} is given, the path of the saved pdf, invisibly.
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
  p <- .plotEachFeature(
    counts,
    decontaminated_counts,
    features,
    function(df.m, feature) {
      xScale <- .pseudoLogScale("x")

      # One bandwidth for both curves, on the scale the density is computed
      # on. Values that geom_density() drops are left out; with fewer than
      # 2 left, it keeps its default
      values <- xScale$get_transformation()$transform(df.m$value)
      values <- values[is.finite(values)]
      bw <- if (length(values) >= 2) stats::bw.nrd0(values) else "nrd0"

      # Plot
      p1 <- ggplot2::ggplot(
        df.m,
        ggplot2::aes(.data$value, fill = .data$variable)
      ) +
        ggplot2::geom_density(alpha = 0.7, bw = bw) +
        xScale +
        .conDeconFillScale() +
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
      p1 + ggplot2::coord_cartesian(ylim = ylimit)
    }
  )

  .wrapFeaturePlots(p, file)
}


#' Boxplot of features grouped by cell type
#'
#' @param counts original count matrix of nADT x nDroplet.
#' @param decontaminated_counts decontaminated count matrix.
#' @param cell_type 1xnDroplet vector of cell_type.
#' @param features names of ADT to plot
#' @param file file name to save the plot into as a pdf. ".pdf" is added
#' unless \code{file} already ends in ".pdf". If omitted, return the plot.
#'
#' @return A \code{patchwork} object with one \code{ggplot} per feature, or,
#' if \code{file} is given, the path of the saved pdf, invisibly.
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
  p <- .plotEachFeature(
    counts,
    decontaminated_counts,
    features,
    function(df.m, feature) {
      ggplot2::ggplot(
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
        .pseudoLogScale("y") +
        .conDeconFillScale() +
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
    },
    as.factor(cell_type)
  )

  .wrapFeaturePlots(p, file)
}


## One plot per feature: the feature's counts before and after
## decontamination are stacked by .meltConDecon(), with the cell types when
## given, and passed to 'plotOne' with the feature name
.plotEachFeature <- function(counts,
                             decontaminated_counts,
                             features,
                             plotOne,
                             cellType = NULL) {
  lapply(features, function(feature) {
    df.m <- .meltConDecon(
      counts[feature, ],
      decontaminated_counts[feature, ],
      cellType
    )
    plotOne(df.m, feature)
  })
}


## Pseudo-log scale for the counts axis ("x" or "y")
.pseudoLogScale <- function(axis = c("x", "y")) {
  axis <- match.arg(axis)
  breaks <- c(1, 5, 10^seq(1, 4))
  if (axis == "x") {
    ggplot2::scale_x_continuous(transform = "pseudo_log", breaks = breaks)
  } else {
    ggplot2::scale_y_continuous(transform = "pseudo_log", breaks = breaks)
  }
}


## Fill colours and labels of the counts before and after decontamination
.conDeconFillScale <- function() {
  ggplot2::scale_fill_manual(
    values = c("#E64B35B2", "#4DBBD5B2"),
    labels = c("Original", "Decontaminated")
  )
}


## Combine the per-feature plots in a square-ish grid with one legend at the
## bottom, then return the plot or, if 'file' is given, save it as a pdf
.wrapFeaturePlots <- function(plots, file = NULL) {
  p_wrap <- patchwork::wrap_plots(plots,
    ncol = round(sqrt(length(plots))),
    guides = "collect"
  ) +
    patchwork::plot_layout() & ggplot2::theme(
    legend.position = "bottom",
    plot.margin = ggplot2::unit(c(3, 3, 2, 1), "pt")
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


## The counts before ("con") and after ("decon") decontamination stacked in
## one column, as reshape2::melt() gave them, with the cell types repeated
## when given
.meltConDecon <- function(con, decon, cellType = NULL) {
  df <- data.frame(
    variable = factor(rep(c("con", "decon"), each = length(con)),
      levels = c("con", "decon")
    ),
    value = unname(c(con, decon))
  )
  if (!is.null(cellType)) {
    df <- data.frame(cell_type = rep(cellType, 2), df)
  }
  df
}
