# Builds the two overview figures that have simulated data layers:
#
# - man/figures/decontX-overview.svg (README and site home page) from
#   decontX-overview.tmpl.svg, with UMAP points before and after
#   decontamination;
# - vignettes/decontPro-model.svg ("Get started") from
#   decontPro-model.tmpl.svg, with ADT density curves.
#
# The data are simulated for illustration only; they are not decontX
# output. vignettes/decontX-model.svg has no data layer, so it is edited
# directly. Run from the repository root:
#
#   Rscript dev/figures/make-figures.R

set.seed(84)
here <- "dev/figures"

# UMAP: three cell types; after decontamination they form tight clusters.
# Before, each cell drifts toward the middle in proportion to its
# contamination, and non-T cells pick up the T cell marker.
clusters <- list(
  T = list(cx = 90, cy = 80, sx = 20, sy = 15, n = 90),
  B = list(cx = 250, cy = 70, sx = 16, sy = 13, n = 60),
  M = list(cx = 175, cy = 180, sx = 22, sy = 12, n = 70)
)
cells <- do.call(rbind, lapply(names(clusters), function(k) {
  z <- clusters[[k]]
  data.frame(
    type = k, cx = z$cx, cy = z$cy,
    xa = stats::rnorm(z$n, z$cx, z$sx), ya = stats::rnorm(z$n, z$cy, z$sy)
  )
}))
n <- nrow(cells)
cells$contam <- stats::rbeta(n, 1.5, 6)
mid <- c(
  mean(vapply(clusters, `[[`, numeric(1), "cx")),
  mean(vapply(clusters, `[[`, numeric(1), "cy"))
)
pull <- 1.1 * cells$contam
cells$xb <- cells$xa + pull * (mid[1] - cells$cx) + stats::rnorm(n, 0, 6)
cells$yb <- cells$ya + pull * (mid[2] - cells$cy) + stats::rnorm(n, 0, 6)
is_t <- cells$type == "T"
cells$after <- ifelse(is_t, stats::runif(n, 0.7, 1), stats::runif(n, 0, 0.04))
cells$before <- ifelse(
  is_t, cells$after * (1 - 0.4 * cells$contam), pmin(1, 1.4 * cells$contam)
)

ramp <- grDevices::colorRamp(c("#E3E3E3", "#B9A3D3", "#54278F"))
to_hex <- function(v) {
  m <- ramp(pmin(pmax(v, 0), 1))
  sprintf("#%02X%02X%02X", round(m[, 1]), round(m[, 2]), round(m[, 3]))
}
umap_points <- function(x, y, v) {
  o <- order(v) # draw high expression on top
  paste0(
    "    ",
    sprintf(
      '<circle cx="%.1f" cy="%.1f" r="3.6" fill="%s"/>',
      x[o], y[o], to_hex(v[o])
    ),
    collapse = "\n"
  )
}

# Density of one T cell ADT across all cells, on a 260 x 120 box.
xs <- seq(0, 1, length.out = 120)
mixture <- function(w, m, s) {
  rowSums(vapply(
    seq_along(w), function(i) w[i] * stats::dnorm(xs, m[i], s[i]),
    numeric(length(xs))
  ))
}
density_path <- function(d, ymax) {
  py <- 120 - d / ymax * 105
  paste0(
    "M0,120 ", paste(sprintf("L%.1f,%.1f", xs * 260, py), collapse = " "),
    " L260,120 Z"
  )
}
before <- mixture(c(0.62, 0.38), c(0.32, 0.76), c(0.07, 0.07))
after <- mixture(c(0.62, 0.38), c(0.07, 0.73), c(0.03, 0.07))
ymax <- max(before, after)

fill_template <- function(template, out, values) {
  svg <- paste(readLines(file.path(here, template)), collapse = "\n")
  for (key in names(values)) {
    svg <- sub(paste0("{{", key, "}}"), values[[key]], svg, fixed = TRUE)
  }
  writeLines(svg, out)
}
fill_template(
  "decontX-overview.tmpl.svg", "man/figures/decontX-overview.svg",
  list(
    UMAP_BEFORE = umap_points(cells$xb, cells$yb, cells$before),
    UMAP_AFTER = umap_points(cells$xa, cells$ya, cells$after)
  )
)
fill_template(
  "decontPro-model.tmpl.svg", "vignettes/decontPro-model.svg",
  list(
    DENS_BEFORE = density_path(before, ymax),
    DENS_AFTER = density_path(after, ymax)
  )
)
