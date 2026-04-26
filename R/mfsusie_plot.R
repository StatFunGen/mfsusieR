# Plot helpers and `mfsusie_plot()` for fits returned by mfsusie() /
# fsusie(). Base R only (graphics + grDevices). No ggplot2 / cowplot.
#
# Style notes drawn from `susieR::susie_plot()` and
# `mvsusieR::mvsusie_plot()`:
#   * susieR uses a bold high-saturation palette and exposes
#     `pos` as either a vector or a list(attr=, start=, end=).
#     We adopt the vector form; the list form is rarely used in
#     practice and easy to add later.
#   * susieR draws `abline(h = 0.95, lty = 3)` on PIP plots; we do
#     the same for visual continuity.
#   * mvsusieR tiles outcomes with cowplot::plot_grid; we use
#     base-R `par(mfrow = ...)` to avoid the ggplot stack.
#
# Color scheme: Okabe-Ito (colorblind-friendly, modern). Recycles
# beyond 8 CSes.

.cs_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                 "#0072B2", "#D55E00", "#CC79A7", "#999999")

#' Per-CS color vector
#' @keywords internal
#' @noRd
mf_cs_colors <- function(n_cs) {
  if (n_cs <= 0L) return(character(0))
  rep(.cs_palette, length.out = n_cs)
}

# Per-variable color for the PIP panel: CS color if in a CS, grey
# otherwise.
.pip_colors <- function(fit) {
  p   <- length(fit$pip)
  cs  <- fit$sets$cs %||% list()
  col <- rep("grey60", p)
  pal <- mf_cs_colors(length(cs))
  for (i in seq_along(cs)) col[cs[[i]]] <- pal[i]
  col
}

# Position-space curve for effect l, outcome m. Reads
# `fit$effect_curves[[m]][[l]]` if present (post-processed). Else
# inverts via the existing coef() path.
.effect_curve <- function(fit, l, m) {
  if (!is.null(fit$effect_curves) &&
      !is.null(fit$effect_curves[[m]]) &&
      !is.null(fit$effect_curves[[m]][[l]])) {
    return(fit$effect_curves[[m]][[l]])
  }
  coef_lm <- coef(fit)[[m]]   # L x T_basis
  coef_lm[l, ]
}

# Optional credible band: T x 2 matrix or NULL.
.credible_band <- function(fit, l, m) {
  if (is.null(fit$credible_bands) ||
      is.null(fit$credible_bands[[m]]) ||
      is.null(fit$credible_bands[[m]][[l]])) {
    return(NULL)
  }
  fit$credible_bands[[m]][[l]]
}

# Internal: PIP panel.
.draw_pip <- function(fit, pos = NULL, main = "PIP",
                      xlab = "variable", ylab = "PIP",
                      cex = 0.9, add_legend = TRUE) {
  pip <- fit$pip
  if (is.null(pos)) pos <- seq_along(pip)
  col <- .pip_colors(fit)
  pch <- ifelse(col == "grey60", 1L, 19L)
  plot(pos, pip, type = "p", pch = pch, col = col, cex = cex,
       xlab = xlab, ylab = ylab, ylim = c(0, 1), main = main, las = 1)
  abline(h = 0.95, lty = 3, col = "grey50")
  cs <- fit$sets$cs %||% list()
  if (add_legend && length(cs) > 0L) {
    pal <- mf_cs_colors(length(cs))
    legend("topleft", legend = paste0("CS", seq_along(cs)),
           col = pal, pch = 19L, bty = "n", cex = 0.75)
  }
}

# Internal: effect-curve panel for one outcome.
.draw_effect <- function(fit, m, pos = NULL, main = NULL,
                         show_grid_dots = FALSE, lwd = 1.5,
                         add_legend = TRUE) {
  T_basis <- fit$dwt_meta$T_basis[m]
  if (is.null(pos)) pos <- fit$dwt_meta$pos[[m]]
  if (is.null(main)) {
    main <- if (length(fit$dwt_meta$T_basis) == 1L) {
      sprintf("Effect curves (T = %d)", T_basis)
    } else {
      sprintf("Outcome %d (T = %d)", m, T_basis)
    }
  }
  cs   <- fit$sets$cs %||% list()
  if (length(cs) == 0L) {
    plot.new()
    title(main = paste(main, "\n(no credible sets)"))
    return(invisible(NULL))
  }
  cs_l <- fit$sets$cs_index %||% seq_along(cs)
  pal  <- mf_cs_colors(length(cs))

  # Scalar outcome: dot plot of per-effect mean.
  if (T_basis == 1L) {
    eff <- vapply(cs_l, function(l) .effect_curve(fit, l, m), numeric(1))
    plot(seq_along(cs_l), eff, type = "p", pch = 19L,
         col = pal, cex = 1.4,
         xlab = "credible set", ylab = "effect",
         xaxt = "n", main = main, las = 1)
    axis(1, at = seq_along(cs_l), labels = paste0("CS", seq_along(cs_l)))
    abline(h = 0, lty = 2, col = "grey60")
    return(invisible(NULL))
  }

  curves <- lapply(cs_l, function(l) .effect_curve(fit, l, m))
  bands  <- lapply(cs_l, function(l) .credible_band(fit, l, m))
  yrange <- range(unlist(curves), unlist(lapply(bands, c)), 0,
                  na.rm = TRUE)
  plot(NA, xlim = range(pos), ylim = yrange,
       xlab = "position", ylab = "effect", main = main, las = 1)
  abline(h = 0, lty = 2, col = "grey60")

  for (i in seq_along(cs_l)) {
    band <- bands[[i]]
    if (is.null(band)) next
    polygon(c(pos, rev(pos)),
            c(band[, 1L], rev(band[, 2L])),
            col = adjustcolor(pal[i], alpha.f = 0.25), border = NA)
  }
  for (i in seq_along(cs_l)) {
    lines(pos, curves[[i]], col = pal[i], lwd = lwd)
    if (show_grid_dots) {
      points(pos, curves[[i]], col = pal[i], pch = 21L,
             bg = "white", cex = 0.6)
    }
  }
  if (add_legend) {
    legend("topright", legend = paste0("CS", seq_along(cs_l)),
           col = pal, lwd = lwd, bty = "n", cex = 0.75)
  }
}

#' Plot an mfsusie() / fsusie() fit
#'
#' One unified plot for both single-outcome (`fsusie()`) and
#' multi-outcome (`mfsusie()`) fits. When `M = 1` the layout is a
#' single 2-panel column (PIP on top, effect curves on bottom). When
#' `M > 1` the layout tiles per-outcome panels on a dense grid with
#' the PIP plot in the top-left slot.
#'
#' Effect curves come from `fit$effect_curves` when post-processed
#' smoothing has been applied (see `mf_post_smooth()`); otherwise
#' from the wavelet-domain inverse via `coef()`. Pointwise credible
#' bands are drawn as transparent ribbons when
#' `fit$credible_bands` is populated.
#'
#' For the PIP plot only, see also `susie_plot(fit, y = "PIP")`
#' (re-exported from susieR), which is the standard SuSiE-family
#' PIP plotter.
#'
#' @param fit a fit returned by `mfsusie()` or `fsusie()`.
#' @param m optional integer index. When supplied (and `M > 1`),
#'   plot only that outcome's effect panel (no PIP). Default
#'   `NULL`: full layout.
#' @param pos optional length-p vector for the PIP x-axis.
#' @param show_grid_dots logical, draw circles at each post-remap
#'   grid point on top of the curves. Default `FALSE`.
#' @param lwd numeric, curve line width. Default `1.5`.
#' @param add_legend logical, show CS legend on each panel.
#'   Default `TRUE`.
#' @param ... reserved.
#' @return Called for side effect; returns `invisible(NULL)`.
#' @export
mfsusie_plot <- function(fit, m = NULL, pos = NULL,
                         show_grid_dots = FALSE, lwd = 1.5,
                         add_legend = TRUE, ...) {
  if (!inherits(fit, "mfsusie")) {
    stop("`fit` must be an `mfsusie` (or `fsusie`) fit object.")
  }
  M <- length(fit$dwt_meta$T_basis)

  # Single-outcome focus: just one effect panel.
  if (!is.null(m)) {
    if (m < 1L || m > M) stop(sprintf("`m` must be in 1..%d.", M))
    op <- par(mar = c(4, 4, 2.5, 1))
    on.exit(par(op), add = TRUE)
    .draw_effect(fit, m = m, show_grid_dots = show_grid_dots,
                 lwd = lwd, add_legend = add_legend)
    return(invisible(NULL))
  }

  # M = 1: simple 2-panel column.
  if (M == 1L) {
    op <- par(mfrow = c(2L, 1L), mar = c(4, 4, 2.5, 1))
    on.exit(par(op), add = TRUE)
    .draw_pip(fit, pos = pos, add_legend = add_legend)
    .draw_effect(fit, m = 1L, show_grid_dots = show_grid_dots,
                 lwd = lwd, add_legend = add_legend)
    return(invisible(NULL))
  }

  # M > 1: dense grid; PIP top-left, then per-outcome effect panels.
  n_panels <- M + 1L
  cols     <- ceiling(sqrt(n_panels))
  rows     <- ceiling(n_panels / cols)
  op <- par(mfrow = c(rows, cols), mar = c(4, 4, 2.5, 1))
  on.exit(par(op), add = TRUE)

  .draw_pip(fit, pos = pos, add_legend = add_legend)
  for (mi in seq_len(M)) {
    .draw_effect(fit, m = mi, show_grid_dots = show_grid_dots,
                 lwd = lwd, add_legend = add_legend)
  }
  remaining <- rows * cols - n_panels
  if (remaining > 0L) {
    for (i in seq_len(remaining)) plot.new()
  }
  invisible(NULL)
}

#' @method plot mfsusie
#' @export
plot.mfsusie <- function(x, ...) {
  mfsusie_plot(x, ...)
}
