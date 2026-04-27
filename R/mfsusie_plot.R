# Plot helpers and `mfsusie_plot()` for fits returned by mfsusie() /
# fsusie(). Base R only (graphics + grDevices). No ggplot2 / cowplot.
#
# Style notes:
#   * Color scheme: Okabe-Ito (colorblind-friendly).
#   * Layout uses `layout()` so stack-mode (one row per CS) and
#     dense-tile mode (M > 1) compose cleanly.
#   * Effect panel supports two styles:
#       - `effect_style = "band"`  -- mean curve + ribbon
#       - `effect_style = "errorbar"` -- per-position dot + bar
#   * Facet of CSes:
#       - `facet_cs = "overlay"` -- all CSes on one panel
#       - `facet_cs = "stack"`   -- one row per CS
#       - `facet_cs = "auto"`    -- stack when length(cs) >= 3 or
#                                   affected-region masks disjoint;
#                                   overlay otherwise

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

# Position-space curve for effect l, outcome m. Reads from a
# pre-resolved `smoothed` payload (the chosen entry of
# `fit$smoothed`) if non-NULL, else inverts via the raw coef()
# path.
.effect_curve <- function(fit, l, m, smoothed = NULL) {
  if (!is.null(smoothed) &&
      !is.null(smoothed$effect_curves) &&
      !is.null(smoothed$effect_curves[[m]]) &&
      !is.null(smoothed$effect_curves[[m]][[l]])) {
    return(smoothed$effect_curves[[m]][[l]])
  }
  coef_lm <- coef(fit)[[m]]   # L x T_basis
  coef_lm[l, ]
}

# Optional credible band: T x 2 matrix or NULL.
.credible_band <- function(fit, l, m, smoothed = NULL) {
  if (is.null(smoothed) ||
      is.null(smoothed$credible_bands) ||
      is.null(smoothed$credible_bands[[m]]) ||
      is.null(smoothed$credible_bands[[m]][[l]])) {
    return(NULL)
  }
  smoothed$credible_bands[[m]][[l]]
}

# Optional lfsr curve: numeric vector of length T or NULL. Populated
# only by the HMM smoother.
.lfsr_curve <- function(fit, l, m, smoothed = NULL) {
  if (is.null(smoothed) ||
      is.null(smoothed$lfsr_curves) ||
      is.null(smoothed$lfsr_curves[[m]]) ||
      is.null(smoothed$lfsr_curves[[m]][[l]])) {
    return(NULL)
  }
  smoothed$lfsr_curves[[m]][[l]]
}

# Contiguous (start, end) runs where the credible band excludes
# zero. Returns a list of integer index pairs into `pos`.
.affected_runs <- function(band) {
  if (is.null(band)) return(list())
  flag <- band[, 1L] > 0 | band[, 2L] < 0
  if (!any(flag)) return(list())
  rle_f <- rle(flag)
  ends   <- cumsum(rle_f$lengths)
  starts <- ends - rle_f$lengths + 1L
  on <- which(rle_f$values)
  lapply(on, function(i) c(starts[i], ends[i]))
}

# Affected-region mask (logical T) for a single CS. Returns a
# zero-length vector when band is NULL.
.affected_mask <- function(band) {
  if (is.null(band)) return(logical(0))
  band[, 1L] > 0 | band[, 2L] < 0
}

# Resolve facet_cs: "auto" -> "stack" or "overlay". Stack when
# length(cs) >= 3 OR when affected-region masks are pairwise
# disjoint (and length(cs) >= 2). Overlay otherwise.
.resolve_facet <- function(facet_cs, fit, m, smoothed = NULL) {
  if (facet_cs %in% c("stack", "overlay")) return(facet_cs)
  cs <- fit$sets$cs %||% list()
  K  <- length(cs)
  if (K < 2L) return("overlay")
  if (K >= 3L) return("stack")
  cs_l <- fit$sets$cs_index %||% seq_along(cs)
  bands <- lapply(cs_l, function(l) .credible_band(fit, l, m, smoothed))
  if (any(vapply(bands, is.null, logical(1)))) return("overlay")
  masks <- lapply(bands, .affected_mask)
  if (any(vapply(masks, length, integer(1)) == 0L)) return("overlay")
  if (!any(masks[[1L]] & masks[[2L]])) return("stack")
  "overlay"
}

# Build the PIP-panel title from `fit$sets`. Default: "Posterior
# inclusion probability and credible sets" with the requested
# coverage and the minimum CS purity (smallest min-abs-corr
# across the displayed CSes).
.pip_title <- function(fit) {
  base <- "Non-zero effects"
  cov  <- fit$sets$requested_coverage
  pur  <- fit$sets$purity
  parts <- character()
  if (!is.null(cov) && length(cov) == 1L && is.finite(cov))
    parts <- c(parts, sprintf("coverage %g%%", 100 * cov))
  pur_min <- if (is.data.frame(pur) && "min.abs.corr" %in% colnames(pur))
    suppressWarnings(min(pur$min.abs.corr, na.rm = TRUE)) else NA_real_
  if (is.finite(pur_min))
    parts <- c(parts, sprintf("min |r| %.2f", pur_min))
  if (length(parts) > 0L)
    sprintf("%s\n(%s)", base, paste(parts, collapse = ", "))
  else
    base
}

# Internal: PIP panel.
.draw_pip <- function(fit, pos = NULL, main = NULL,
                      xlab = "variable", ylab = "PIP",
                      cex = 1.2, add_legend = TRUE) {
  pip <- fit$pip
  if (is.null(pos)) pos <- seq_along(pip)
  if (is.null(main)) main <- .pip_title(fit)
  col <- .pip_colors(fit)
  pch <- ifelse(col == "grey60", 1L, 19L)
  plot(pos, pip, type = "p", pch = pch, col = col, cex = cex,
       lwd = 1.6,
       xlab = xlab, ylab = ylab, ylim = c(0, 1), main = main,
       cex.main = 1.05, font.main = 2L, las = 1)
  abline(h = 0.95, lty = 3, col = "grey50")
  cs <- fit$sets$cs %||% list()
  if (add_legend && length(cs) > 0L) {
    pal <- mf_cs_colors(length(cs))
    legend("topleft", legend = paste0("CS", seq_along(cs)),
           col = pal, pch = 19L, bty = "n", cex = 0.75)
  }
}

# =====================================================================
# Effect-panel cell drawer
# =====================================================================
#
# A single `(outcome, cs_subset)` panel is rendered into the
# current device cell by `.draw_effect_panel`. It does not own
# the layout: the caller decides whether the cell is one of M*K
# stacked rows, the lower half of a 2-row PIP+effect column, or
# a tile in an M-outcome grid.
#
# Two effect styles share the cell layout:
#   - "band":     ribbon polygon + line for each CS posterior mean
#   - "errorbar": vertical capped bar + filled point per position
# Everything else (axes, h=0 reference, optional affected-region
# bars, optional lfsr right-axis overlay, optional truth dots,
# legend) is identical, factored into the helpers below.

# y-axis range that fits all curves, bands, optional truth dots,
# and 0.
.effect_yrange <- function(curves, bands, truth_per_cs = NULL) {
  has_truth <- !is.null(truth_per_cs) &&
    any(!vapply(truth_per_cs, is.null, logical(1)))
  range(unlist(curves), unlist(lapply(bands, c)),
        if (has_truth) unlist(truth_per_cs) else NULL,
        0, na.rm = TRUE)
}

# Per-CS truth values as dark-grey dots.
.overlay_truth_dots <- function(pos, truth_per_cs) {
  if (is.null(truth_per_cs)) return(invisible(NULL))
  for (i in seq_along(truth_per_cs)) {
    tr <- truth_per_cs[[i]]
    if (is.null(tr)) next
    pos_i <- pos[seq_along(tr)]
    points(pos_i, tr, col = "grey30", pch = 20L, cex = 0.7)
  }
  invisible(NULL)
}

# Bottom-of-panel affected-region bars (band style).
.overlay_affected_bars <- function(pos, bands, pal, yrange) {
  bar_y    <- yrange[1L]
  bar_step <- 0.03 * diff(yrange)
  for (i in seq_along(bands)) {
    runs <- .affected_runs(bands[[i]])
    if (length(runs) == 0L) next
    y_i <- bar_y + (i - 1L) * bar_step
    for (rg in runs) {
      segments(pos[rg[1L]], y_i, pos[rg[2L]], y_i,
               col = pal[i], lwd = 4.5, lend = "butt")
    }
  }
}

# Affected-region dots on the zero line (errorbar style).
.overlay_affected_dots <- function(pos, bands) {
  for (i in seq_along(bands)) {
    mask <- .affected_mask(bands[[i]])
    if (length(mask) == 0L || !any(mask)) next
    points(pos[mask], rep(0, sum(mask)),
           col = "black", pch = 20L, cex = 1.2)
  }
}

# Right-axis lfsr overlay (band style only).
.overlay_lfsr_axis <- function(pos, lfsrs, pal, lfsr_threshold) {
  op2 <- par(new = TRUE); on.exit(par(op2), add = TRUE)
  plot(NA, xlim = range(pos), ylim = c(0, 1),
       axes = FALSE, xlab = "", ylab = "")
  axis(4, at = c(0, lfsr_threshold, 0.5, 1),
       labels = c("0", sprintf("%.2g", lfsr_threshold), "0.5", "1"),
       las = 1)
  mtext("lfsr", side = 4, line = 2.4)
  abline(h = lfsr_threshold, lty = 2, col = "firebrick", lwd = 1.2)
  label_x <- pos[1L] + 0.02 * diff(range(pos))
  text(label_x, lfsr_threshold,
       labels = sprintf("lfsr = %.2g", lfsr_threshold),
       col = "firebrick", cex = 0.75, pos = 3)
  for (i in seq_along(lfsrs)) {
    if (is.null(lfsrs[[i]])) next
    lines(pos, lfsrs[[i]], col = pal[i], lwd = 1.4, lty = 2)
  }
}

# CS legend with optional "truth" header row.
.effect_cs_legend <- function(cs_subset, pal, has_truth, lwd, style) {
  K <- length(cs_subset)
  lab <- paste0("CS", seq_len(K))
  cols <- pal
  if (style == "band") {
    lwds <- rep(lwd, K); ltys <- rep(1L, K); pchs <- rep(NA_integer_, K)
  } else {
    lwds <- rep(NA_real_, K); ltys <- rep(NA_integer_, K); pchs <- rep(16L, K)
  }
  if (has_truth) {
    lab  <- c("truth", lab)
    cols <- c("grey30", cols)
    lwds <- c(NA_real_, lwds)
    ltys <- c(NA_integer_, ltys)
    pchs <- c(20L, pchs)
  }
  legend("topright", legend = lab, col = cols, lwd = lwds,
         lty = ltys, pch = pchs, bty = "n", cex = 0.75)
}

# Body draw for the band style: ribbon polygons + mean curves.
.draw_effect_body_band <- function(pos, curves, bands, pal, lwd,
                                   show_grid_dots) {
  for (i in seq_along(curves)) {
    band <- bands[[i]]
    if (is.null(band)) next
    polygon(c(pos, rev(pos)),
            c(band[, 1L], rev(band[, 2L])),
            col = adjustcolor(pal[i], alpha.f = 0.25), border = NA)
  }
  for (i in seq_along(curves)) {
    lines(pos, curves[[i]], col = pal[i], lwd = lwd)
    if (show_grid_dots) {
      points(pos, curves[[i]], col = pal[i], pch = 21L,
             bg = "white", cex = 0.85, lwd = 1.4)
    }
  }
}

# Body draw for the errorbar style: capped vertical bars +
# filled points per position.
.draw_effect_body_errorbar <- function(pos, curves, bands, pal) {
  cap_len <- min(0.04, 0.3 / max(length(pos), 1L))
  for (i in seq_along(curves)) {
    band    <- bands[[i]]
    if (!is.null(band)) {
      suppressWarnings(
        arrows(pos, band[, 1L], pos, band[, 2L],
               col = pal[i], lwd = 1.4,
               length = cap_len, angle = 90, code = 3L))
    }
    points(pos, curves[[i]], col = pal[i], pch = 16L, cex = 1.0)
  }
}

# Single-cell effect panel. `effect_style` selects the body
# drawer; everything else (frame, lfsr, truth, legend) is shared.
.draw_effect_panel <- function(fit, m, cs_subset, effect_style,
                                pos, pal, lwd, smoothed,
                                show_grid_dots, show_affected_region,
                                show_lfsr_curve, lfsr_threshold,
                                add_legend, main,
                                xlab = "outcome position", xaxt = "s",
                                truth_per_cs = NULL) {
  curves <- lapply(cs_subset, function(l) .effect_curve(fit, l, m, smoothed))
  bands  <- lapply(cs_subset, function(l) .credible_band(fit, l, m, smoothed))
  lfsrs  <- lapply(cs_subset, function(l) .lfsr_curve(fit, l, m, smoothed))

  if (effect_style == "errorbar" &&
      all(vapply(bands, is.null, logical(1)))) {
    stop("`effect_style = \"errorbar\"` requires post-smoothed ",
         "credible bands. Call `mf_post_smooth()` first.")
  }

  has_lfsr  <- effect_style == "band" && show_lfsr_curve &&
    any(!vapply(lfsrs, is.null, logical(1)))
  has_truth <- !is.null(truth_per_cs) &&
    any(!vapply(truth_per_cs, is.null, logical(1)))

  yrange <- .effect_yrange(curves, bands, truth_per_cs)
  plot(NA, xlim = range(pos), ylim = yrange,
       xlab = xlab, ylab = "effect", main = main, las = 1,
       xaxt = xaxt, cex.main = 1.05, font.main = 2L)
  abline(h = 0, lty = 2, col = "grey60")

  if (effect_style == "band") {
    .draw_effect_body_band(pos, curves, bands, pal, lwd, show_grid_dots)
    if (show_affected_region) .overlay_affected_bars(pos, bands, pal, yrange)
  } else {
    .draw_effect_body_errorbar(pos, curves, bands, pal)
    if (show_affected_region) .overlay_affected_dots(pos, bands)
  }

  # Truth dots use effect-axis coordinates and must be drawn
  # BEFORE the lfsr overlay, which switches the panel's
  # plot.window to (0, 1) y-coordinates for the right axis.
  # Anything plotted after `.overlay_lfsr_axis` would inherit
  # those (0, 1) coords and be misplaced.
  .overlay_truth_dots(pos, truth_per_cs)
  if (has_lfsr) .overlay_lfsr_axis(pos, lfsrs, pal, lfsr_threshold)

  if (add_legend) {
    if (length(cs_subset) > 1L || has_truth) {
      .effect_cs_legend(cs_subset, pal, has_truth,
                        lwd = if (effect_style == "band") lwd else NA_real_,
                        style = effect_style)
    }
    if (has_lfsr) {
      legend("topleft",
             legend = c("effect (left axis)", "lfsr (right axis)"),
             lty = c(1, 2), lwd = c(lwd, 1.4),
             col = "grey30", bty = "n", cex = 0.7,
             seg.len = 2.0)
    }
  }
}

# Internal: default outcome label for the effect panel title.
.outcome_main <- function(fit, m) {
  T_basis <- fit$dwt_meta$T_basis[m]
  nm      <- fit$dwt_meta$outcome_names[m]
  label   <- if (!is.null(nm) && nzchar(nm)) nm else sprintf("Outcome %d", m)
  if (length(fit$dwt_meta$T_basis) == 1L)
    sprintf("Effect sizes (T = %d)", T_basis)
  else
    sprintf("%s (T = %d)", label, T_basis)
}

# Internal: scalar (T = 1) outcome dot plot.
.draw_scalar_effect <- function(fit, m, smoothed, main) {
  cs   <- fit$sets$cs %||% list()
  cs_l <- fit$sets$cs_index %||% seq_along(cs)
  pal  <- mf_cs_colors(length(cs))
  eff  <- vapply(cs_l, function(l) .effect_curve(fit, l, m, smoothed),
                 numeric(1))
  plot(seq_along(cs_l), eff, type = "p", pch = 19L,
       col = pal, cex = 1.4,
       xlab = "credible set", ylab = "effect",
       xaxt = "n", main = main, las = 1)
  axis(1, at = seq_along(cs_l), labels = paste0("CS", seq_along(cs_l)))
  abline(h = 0, lty = 2, col = "grey60")
}

# Internal: single-cell overlay effect panel (all CSes on one cell).
# Always draws into the current device cell. Caller owns the layout.
.draw_effect_in_cell <- function(fit, m, smoothed,
                                  effect_style, pos, lwd,
                                  show_grid_dots, show_affected_region,
                                  show_lfsr_curve, lfsr_threshold,
                                  add_legend, main = NULL,
                                  truth_per_cs = NULL) {
  T_basis <- fit$dwt_meta$T_basis[m]
  if (is.null(pos))  pos  <- fit$dwt_meta$pos[[m]]
  if (is.null(main)) main <- .outcome_main(fit, m)
  cs <- fit$sets$cs %||% list()
  if (length(cs) == 0L) {
    plot.new(); title(main = paste(main, "\n(no credible sets)"))
    return(invisible(NULL))
  }
  if (T_basis == 1L) {
    .draw_scalar_effect(fit, m, smoothed, main)
    return(invisible(NULL))
  }
  cs_l <- fit$sets$cs_index %||% seq_along(cs)
  pal  <- mf_cs_colors(length(cs))
  .draw_effect_panel(fit, m, cs_subset = cs_l,
                     effect_style = effect_style,
                     pos = pos, pal = pal, lwd = lwd,
                     smoothed = smoothed,
                     show_grid_dots = show_grid_dots,
                     show_affected_region = show_affected_region,
                     show_lfsr_curve = show_lfsr_curve,
                     lfsr_threshold = lfsr_threshold,
                     add_legend = add_legend,
                     main = main,
                     truth_per_cs = truth_per_cs)
}

# Internal: K stacked per-CS effect cells. Caller has already
# allocated K cells via `layout()`. `mar_cs` lets the M = 1 path
# tighten margins for the stacked sub-panels.
.draw_effect_stack_cells <- function(fit, m, smoothed, K,
                                     effect_style, pos, lwd,
                                     show_grid_dots, show_affected_region,
                                     show_lfsr_curve, lfsr_threshold,
                                     add_legend,
                                     mar_cs = NULL,
                                     truth_per_cs = NULL,
                                     main_prefix = NULL,
                                     last_cell = TRUE) {
  cs   <- fit$sets$cs %||% list()
  cs_l <- fit$sets$cs_index %||% seq_along(cs)
  pal  <- mf_cs_colors(length(cs))
  if (is.null(pos)) pos <- fit$dwt_meta$pos[[m]]
  if (!is.null(mar_cs)) {
    op <- par(mar = mar_cs); on.exit(par(op), add = TRUE)
  }
  for (i in seq_len(K)) {
    is_last <- last_cell && (i == K)
    title_i <- if (is.null(main_prefix)) sprintf("CS%d", i)
               else sprintf("%s — CS%d", main_prefix, i)
    .draw_effect_panel(fit, m, cs_subset = cs_l[i],
                       effect_style = effect_style,
                       pos = pos, pal = pal[i], lwd = lwd,
                       smoothed = smoothed,
                       show_grid_dots = show_grid_dots,
                       show_affected_region = show_affected_region,
                       show_lfsr_curve = show_lfsr_curve,
                       lfsr_threshold = lfsr_threshold,
                       add_legend = (i == 1L) && add_legend,
                       main = title_i,
                       xaxt = if (is_last) "s" else "n",
                       xlab = if (is_last) "outcome position" else "",
                       truth_per_cs = if (!is.null(truth_per_cs))
                         truth_per_cs[i] else NULL)
  }
}

# Internal: normalize the user-facing `truth` argument into a
# canonical list[M] of list[K] of length-T_m vectors (or NULL each).
# Accepted shapes:
#
#  M = 1, K = K:
#    NULL                    -> no truth
#    numeric vector length-T -> single truth, replicated for every CS
#    list of K vectors       -> per-CS truth
#    list(vec)               -> single truth, replicated for every CS
#
#  M > 1:
#    NULL                    -> no truth
#    list of length M, each entry NULL / vector / list-of-K-vectors
.normalize_truth <- function(truth, M, K) {
  empty <- replicate(M, vector("list", K), simplify = FALSE)
  if (is.null(truth)) return(empty)
  per_cs <- function(x) {
    if (is.null(x)) return(vector("list", K))
    if (is.numeric(x) && !is.list(x)) return(rep(list(x), K))
    if (is.list(x) && length(x) == K) return(x)
    if (is.list(x) && length(x) == 1L && is.numeric(x[[1L]]))
      return(rep(list(x[[1L]]), K))
    stop("`truth` entry must be NULL, a length-T_m vector, or a ",
         "length-K list of length-T_m vectors.")
  }
  if (M == 1L) {
    if (is.numeric(truth) && !is.list(truth)) return(list(rep(list(truth), K)))
    if (is.list(truth) && length(truth) == K) return(list(truth))
    if (is.list(truth) && length(truth) == 1L)  return(list(per_cs(truth[[1L]])))
    stop("`truth` for M = 1 must be NULL, a length-T_1 vector, ",
         "or a length-K list of length-T_1 vectors.")
  }
  if (!is.list(truth) || length(truth) != M)
    stop("`truth` for M > 1 must be a length-M list.")
  lapply(truth, per_cs)
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
#' from the wavelet-domain inverse via `coef()`. With
#' `effect_style = "band"` (default) credible bands are drawn as
#' transparent ribbons; with `effect_style = "errorbar"` each
#' position gets a dot at the mean and a vertical bar to the
#' credible band (errorbar mode requires post-smoothed credible
#' bands and errors otherwise).
#'
#' Multiple credible sets can be drawn overlaid (`facet_cs =
#' "overlay"`) or in stacked rows (`facet_cs = "stack"`); the
#' default `"auto"` picks stack when there are at least 3 CSes or
#' when their affected-region masks are pairwise disjoint.
#'
#' After `mf_post_smooth(method = "HMM")` the per-position lfsr
#' curves are overlaid on a secondary axis (band mode) with a
#' threshold line at `lfsr_threshold` (default `0.01`). For the
#' per-CS lfsr bubble grid see `mfsusie_plot_lfsr()`.
#'
#' @param fit a fit returned by `mfsusie()` or `fsusie()`.
#' @param m optional integer index. When supplied, plot only that
#'   outcome's effect panel (no PIP). Default `NULL`: full layout.
#' @param pos optional length-p vector for the PIP x-axis.
#' @param effect_style `"band"` (default) or `"errorbar"`. Style
#'   of the effect panel.
#' @param facet_cs `"auto"` (default), `"stack"`, or `"overlay"`.
#'   Layout of multiple credible sets within an effect panel.
#' @param show_grid_dots logical, draw circles at each post-remap
#'   grid point on top of the curves. Default `FALSE`.
#' @param show_lfsr_curve logical, overlay HMM lfsr curves on a
#'   secondary axis when `fit$lfsr_curves` is populated. Only
#'   honored in `effect_style = "band"`. Default `TRUE`. The
#'   solid line is the effect (left axis), the dashed line is
#'   the lfsr (right axis, 0 to 1). For a dedicated per-CS
#'   lfsr view see `mfsusie_plot_lfsr()`.
#' @param show_affected_region logical, mark contiguous segments
#'   where each CS's credible band excludes zero with thick bars
#'   at the bottom of the effect panel (band mode) or with black
#'   dots on the zero line (errorbar mode). Default `TRUE`.
#' @param lfsr_threshold numeric, dashed reference line on the lfsr
#'   secondary axis. Default `0.01`.
#' @param lwd numeric, curve line width. Default `1.5`.
#' @param add_legend logical, show CS legend on each panel.
#'   Default `TRUE`.
#' @param truth optional ground-truth overlay. For `M = 1`: a
#'   length-`T_1` numeric vector (replicated across all CS panels)
#'   or a length-`K` list of length-`T_1` vectors (per-CS truth).
#'   For `M > 1`: a length-`M` list, each entry one of those two
#'   shapes (or `NULL` to skip an outcome). Truth values are
#'   overlaid as dark-grey dots and added to the legend. Default
#'   `NULL` (no overlay).
#' @param save optional file path. When non-NULL the plot is
#'   written to the file at the dimensions returned by
#'   `mfsusie_plot_dimensions()`. The graphics device is selected
#'   from the file extension: `.pdf`, `.png`, `.jpg`/`.jpeg`, or
#'   `.svg`. Default `NULL` (draw to the current device).
#' @param ... reserved.
#' @return Called for side effect; returns `invisible(NULL)`.
#' @export
mfsusie_plot <- function(fit, m = NULL, pos = NULL,
                         effect_style = c("band", "errorbar"),
                         facet_cs = c("auto", "stack", "overlay"),
                         show_grid_dots = FALSE,
                         show_lfsr_curve = TRUE,
                         show_affected_region = TRUE,
                         lfsr_threshold = 0.01,
                         lwd = 2.0,
                         add_legend = TRUE,
                         truth = NULL,
                         smooth_method = NULL,
                         save = NULL, ...) {
  if (!inherits(fit, "mfsusie")) {
    stop("`fit` must be an `mfsusie` (or `fsusie`) fit object.")
  }
  effect_style <- match.arg(effect_style)
  facet_cs     <- match.arg(facet_cs)

  picked   <- .pick_smooth_method(fit, smooth_method)
  smoothed <- if (!is.null(picked)) fit$smoothed[[picked]] else NULL

  M <- length(fit$dwt_meta$T_basis)
  K <- length(fit$sets$cs %||% list())

  # Optional file output: open a sized device and close it on
  # exit so the figure renders at the recommended dimensions.
  if (!is.null(save)) {
    dims <- mfsusie_plot_dimensions(fit, m = m, facet_cs = facet_cs,
                                    smooth_method = smooth_method)
    .open_save_device(save, dims)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  truth_norm <- .normalize_truth(truth, M, max(K, 1L))
  truth_for  <- function(mi) {
    if (K < 1L) return(NULL)
    truth_norm[[mi]]
  }

  # Resolve facet for the outcome we are about to render. Stack
  # only fires when there are at least 2 credible sets and the
  # caller did not force overlay.
  resolve_facet_for <- function(mi) {
    if (K < 2L) return("overlay")
    .resolve_facet(facet_cs, fit, mi, smoothed)
  }

  draw_overlay_cell <- function(mi, main = NULL) {
    .draw_effect_in_cell(fit, m = mi, smoothed = smoothed,
                         effect_style = effect_style,
                         pos = NULL, lwd = lwd,
                         show_grid_dots = show_grid_dots,
                         show_affected_region = show_affected_region,
                         show_lfsr_curve = show_lfsr_curve,
                         lfsr_threshold = lfsr_threshold,
                         add_legend = add_legend, main = main,
                         truth_per_cs = truth_for(mi))
  }
  draw_stack_cells <- function(mi, main_outer, last_cell = TRUE,
                                main_prefix = NULL) {
    .draw_effect_stack_cells(fit, m = mi, smoothed = smoothed, K = K,
                             effect_style = effect_style,
                             pos = NULL, lwd = lwd,
                             show_grid_dots = show_grid_dots,
                             show_affected_region = show_affected_region,
                             show_lfsr_curve = show_lfsr_curve,
                             lfsr_threshold = lfsr_threshold,
                             add_legend = add_legend,
                             mar_cs = c(2.5, 4, 1.8, 4),
                             truth_per_cs = truth_for(mi),
                             main_prefix = main_prefix,
                             last_cell = last_cell)
    if (!is.null(main_outer))
      mtext(main_outer, side = 3, outer = TRUE, line = 0.2,
            cex = 1.05, font = 2L)
  }

  # Single-outcome focus: one effect region. May expand to K
  # stacked sub-panels when the resolver picks stack.
  if (!is.null(m)) {
    if (m < 1L || m > M) stop(sprintf("`m` must be in 1..%d.", M))
    facet_res <- resolve_facet_for(m)
    if (facet_res == "stack" && K >= 2L) {
      layout(matrix(seq_len(K), ncol = 1L), heights = rep(1, K))
      on.exit(layout(1L), add = TRUE)
      op <- par(mar = c(2.5, 4, 1.8, 4),
                oma = c(2.5, 0, 1.8, 0))
      on.exit(par(op), add = TRUE)
      draw_stack_cells(m, main_outer = .outcome_main(fit, m))
    } else {
      op <- par(mar = c(4, 4, 2.5, 4))
      on.exit(par(op), add = TRUE)
      draw_overlay_cell(m)
    }
    return(invisible(NULL))
  }

  # Spacing policy: in any layout that combines a PIP cell with
  # effect cells, the PIP cell is allocated 1.5x the height of a
  # single effect cell. The PIP cell carries more axis content
  # (the variable axis, the CS legend, the threshold line) than
  # an effect cell and looks visibly squeezed at equal weight,
  # especially when the effect region is split into K stacked
  # CS panels.
  pip_weight <- 1.5

  # M = 1: PIP on top, effect region below. Effect region is
  # one overlay cell or K stacked cells per the facet resolver.
  if (M == 1L) {
    facet_res       <- resolve_facet_for(1L)
    n_effect_panels <- if (facet_res == "stack" && K >= 2L) K else 1L
    layout(matrix(seq_len(1L + n_effect_panels), ncol = 1L),
           heights = c(pip_weight, rep(1.0, n_effect_panels)))
    on.exit(layout(1L), add = TRUE)
    op <- par(mar = c(4, 4, 2.5, 4))
    on.exit(par(op), add = TRUE)
    .draw_pip(fit, pos = pos, add_legend = add_legend)
    if (n_effect_panels == 1L) {
      draw_overlay_cell(1L)
    } else {
      draw_stack_cells(1L, main_outer = NULL)
    }
    return(invisible(NULL))
  }

  # M > 1: two layouts.
  #   facet_cs = "stack" with K >= 2: vertical stack of 1 PIP cell
  #     and M*K per-(outcome, CS) cells. Use this when each CS
  #     deserves its own panel per outcome (e.g., DNAm-CS1, ...,
  #     RNA-CS3 stacked vertically).
  #   else: tiled grid, PIP in the top-left slot, M overlay cells.
  if (facet_cs == "stack" && K >= 2L) {
    n_cells <- 1L + M * K
    layout(matrix(seq_len(n_cells), ncol = 1L),
           heights = c(pip_weight, rep(1.0, M * K)))
    on.exit(layout(1L), add = TRUE)
    op <- par(mar = c(4, 4, 2.5, 4))
    on.exit(par(op), add = TRUE)
    .draw_pip(fit, pos = pos, add_legend = add_legend)
    for (mi in seq_len(M)) {
      is_last_outcome <- (mi == M)
      draw_stack_cells(mi, main_outer = NULL,
                       last_cell = is_last_outcome,
                       main_prefix = .outcome_main(fit, mi))
    }
    return(invisible(NULL))
  }

  n_panels <- M + 1L
  cols     <- ceiling(sqrt(n_panels))
  rows     <- ceiling(n_panels / cols)
  op <- par(mfrow = c(rows, cols), mar = c(4, 4, 2.5, 4))
  on.exit(par(op), add = TRUE)

  .draw_pip(fit, pos = pos, add_legend = add_legend)
  for (mi in seq_len(M)) draw_overlay_cell(mi)
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

# Internal: open a graphics device sized to `dims` if `save` is
# non-NULL. Returns NULL when no device is opened. The caller is
# responsible for `dev.off()` (typically via `on.exit`).
.open_save_device <- function(save, dims) {
  if (is.null(save)) return(invisible(NULL))
  ext <- tolower(tools::file_ext(save))
  if (ext == "pdf") {
    grDevices::pdf(save, width = dims$width, height = dims$height)
  } else if (ext == "png") {
    grDevices::png(save, width = dims$width, height = dims$height,
                   units = "in", res = 150)
  } else if (ext %in% c("jpg", "jpeg")) {
    grDevices::jpeg(save, width = dims$width, height = dims$height,
                    units = "in", res = 150, quality = 95)
  } else if (ext == "svg") {
    if (!requireNamespace("svglite", quietly = TRUE)) {
      grDevices::svg(save, width = dims$width, height = dims$height)
    } else {
      svglite::svglite(save, width = dims$width, height = dims$height)
    }
  } else {
    stop("`save` must end in .pdf, .png, .jpg/.jpeg, or .svg; got '",
         ext, "'.")
  }
  invisible(save)
}

#' Recommended figure dimensions for `mfsusie_plot()`
#'
#' Returns a width and height (in inches) sized to the number of
#' panels `mfsusie_plot()` will draw on the supplied fit. The
#' recommendation accounts for the PIP cell, the per-outcome
#' effect cells, and the per-CS sub-cells when `facet_cs =
#' "stack"` is used with `K >= 2` credible sets. Useful for
#' setting `fig.width` / `fig.height` in a knitr chunk so each
#' panel has enough vertical space to render legibly.
#'
#' Heuristics: `min_pip_height` (default 1.6 in) for the PIP
#' cell, `min_cell_height` (default 1.5 in) for each effect
#' cell. The PIP cell carries 1.5x the height weight of an
#' effect cell. The recommended figure width defaults to
#' `width = 8` inches.
#'
#' @param fit a fit returned by `mfsusie()` or `fsusie()`.
#' @param m optional integer index. When supplied, recommend
#'   dimensions for the focused single-outcome view.
#' @param facet_cs `"auto"`, `"stack"`, or `"overlay"`. Default
#'   `"auto"` matches `mfsusie_plot()`'s default.
#' @param smooth_method optional smoother name to resolve
#'   `facet_cs = "auto"` against the smoothed credible bands.
#' @param min_pip_height minimum height (in) for the PIP cell.
#' @param min_cell_height minimum height (in) for each effect
#'   cell.
#' @param width figure width in inches. Default `8`.
#' @return a list with components `width`, `height` (numeric,
#'   inches) and `n_cells` (integer, total panel count).
#' @examples
#' \donttest{
#' set.seed(1L)
#' n <- 100; p <- 20
#' X <- matrix(rnorm(n * p), n)
#' Y <- list(matrix(rnorm(n * 32), n), matrix(rnorm(n * 32), n))
#' fit <- mfsusie(X, Y, L = 4, verbose = FALSE)
#' fit_s <- mf_post_smooth(fit)
#' mfsusie_plot_dimensions(fit_s, facet_cs = "stack")
#' }
#' @export
mfsusie_plot_dimensions <- function(fit, m = NULL,
                                    facet_cs = c("auto", "stack",
                                                 "overlay"),
                                    smooth_method = NULL,
                                    min_pip_height = 2.2,
                                    min_cell_height = 2.0,
                                    width = 8) {
  if (!inherits(fit, "mfsusie")) {
    stop("`fit` must be an `mfsusie` (or `fsusie`) fit object.")
  }
  facet_cs <- match.arg(facet_cs)
  M <- length(fit$dwt_meta$T_basis)
  K <- length(fit$sets$cs %||% list())

  # Resolve facet_cs the same way mfsusie_plot does.
  picked   <- .pick_smooth_method(fit, smooth_method)
  smoothed <- if (!is.null(picked)) fit$smoothed[[picked]] else NULL
  facet_for <- function(mi) {
    if (K < 2L) return("overlay")
    .resolve_facet(facet_cs, fit, mi, smoothed)
  }

  # Cell counts per layout branch in mfsusie_plot.
  if (!is.null(m)) {
    if (m < 1L || m > M)
      stop(sprintf("`m` must be in 1..%d.", M))
    facet_res <- facet_for(m)
    n_cells <- if (facet_res == "stack" && K >= 2L) K else 1L
    has_pip <- FALSE
  } else if (M == 1L) {
    facet_res <- facet_for(1L)
    n_eff   <- if (facet_res == "stack" && K >= 2L) K else 1L
    n_cells <- 1L + n_eff
    has_pip <- TRUE
  } else if (facet_cs == "stack" && K >= 2L) {
    n_cells <- 1L + M * K
    has_pip <- TRUE
  } else {
    # Tiled grid layout: ceiling(sqrt(M+1)) x ceiling((M+1)/cols)
    n_panels <- M + 1L
    cols     <- ceiling(sqrt(n_panels))
    rows     <- ceiling(n_panels / cols)
    return(list(
      width   = max(width, cols * min_cell_height),
      height  = max(min_cell_height, rows * min_cell_height),
      n_cells = n_panels))
  }

  effect_cells <- n_cells - as.integer(has_pip)
  height <- (if (has_pip) min_pip_height else 0) +
            effect_cells * min_cell_height
  list(width = width, height = height, n_cells = n_cells)
}

#' Recommended figure dimensions for `mfsusie_plot_lfsr()`
#'
#' Returns a `(width, height)` (in inches) sized to the per-CS
#' bubble grid so each bubble row has enough vertical space and
#' each outcome panel has enough horizontal space. Mirrors
#' `mfsusie_plot_dimensions()` for the lfsr companion plot. Use
#' the result for `fig.width` / `fig.height` in a knitr chunk so
#' the legend and bubbles do not collide.
#'
#' @param fit an `mfsusie` / `fsusie()` fit object that has been
#'   post-processed with `mf_post_smooth(method = "HMM")`.
#' @param add_legend logical, whether the lfsr legend strip is
#'   shown (matches the corresponding `mfsusie_plot_lfsr()` arg).
#'   Default `TRUE`.
#' @return a list with components `width`, `height` (numeric, in
#'   inches), `cols`, `rows` (integer, the M-grid layout used).
#' @export
mfsusie_plot_lfsr_dimensions <- function(fit, add_legend = TRUE) {
  if (!inherits(fit, "mfsusie")) {
    stop("`fit` must be an `mfsusie` (or `fsusie`) fit object.")
  }
  M    <- length(fit$dwt_meta$T_basis)
  K    <- length(fit$sets$cs %||% list())
  cols <- ceiling(sqrt(M))
  rows <- ceiling(M / cols)
  width  <- max(6, 4 * cols)
  # Extra header strip when legend is on: title (~1 line) +
  # legend row (~1 line) + spacing -> ~2.5 in.
  height <- max(3, 0.6 * max(K, 1L) * rows +
                  (if (isTRUE(add_legend)) 2.5 else 1))
  list(width = width, height = height, cols = cols, rows = rows)
}

# =====================================================================
# Per-CS lfsr bubble grid
# =====================================================================
#
# Layout: rows are CSes, columns are positions. Dot size encodes
# `-log10(lfsr)` clamped to `[0, cex_max]`. Color marks lfsr at or
# below `lfsr_threshold` (or, when `truth` is supplied, the
# ground-truth affected mask). Multi-outcome fits tile one panel
# per outcome; the function takes no `m` argument and chooses the
# layout from `length(fit$dwt_meta$T_basis)`.

# Internal: draw a single bubble panel for one outcome.
.draw_lfsr_bubble <- function(fit, m, lfsr_threshold, truth_mask,
                              cex_max, main = NULL,
                              smoothed = NULL,
                              add_legend = TRUE) {
  cs   <- fit$sets$cs %||% list()
  cs_l <- fit$sets$cs_index %||% seq_along(cs)
  K    <- length(cs_l)
  if (K == 0L) {
    plot.new()
    title(main = paste(main %||% "", "\n(no credible sets)"))
    return(invisible(NULL))
  }
  pos <- fit$dwt_meta$pos[[m]]
  T_m <- length(pos)

  lfsrs <- lapply(cs_l, function(l) .lfsr_curve(fit, l, m, smoothed))
  if (all(vapply(lfsrs, is.null, logical(1)))) {
    stop("`mfsusie_plot_lfsr()` requires HMM-smoothed lfsr ",
         "curves. Call `mf_post_smooth(method = \"HMM\")` first.")
  }

  pal <- mf_cs_colors(K)
  if (is.null(main)) {
    nm <- fit$dwt_meta$outcome_names[m]
    label <- if (!is.null(nm) && nzchar(nm)) nm else
      sprintf("Outcome %d", m)
    main <- if (length(fit$dwt_meta$T_basis) == 1L) {
      "Per-CS local false sign rate"
    } else {
      sprintf("%s: per-CS lfsr", label)
    }
  }

  plot(NA, xlim = range(pos),
       ylim = c(0.5, K + 0.5),
       xlab = "outcome position", ylab = "credible set",
       yaxt = "n", main = "", las = 1)
  # Title sits above the legend strip so neither overlaps the
  # plot box. `mtext` uses margin-line coordinates: line = 4 is
  # the top of the reserved top margin, line = 1.2 is the strip
  # between the panel and the title used by the legend rows.
  mtext(main, side = 3L, line = 4.0, cex = 1.05, font = 2L)
  axis(2, at = seq_len(K), labels = paste0("CS", seq_len(K)),
       las = 1)
  abline(h = seq_len(K), lty = 3, col = "grey85")

  for (i in seq_len(K)) {
    lfsr_i <- lfsrs[[i]]
    if (is.null(lfsr_i)) next

    sig <- lfsr_i <= lfsr_threshold

    # Color rule: if `truth` is supplied, color by truth; else by
    # `lfsr <= threshold`.
    if (!is.null(truth_mask) && length(truth_mask) >= i &&
        !is.null(truth_mask[[i]])) {
      mask_i <- as.logical(truth_mask[[i]])
      if (length(mask_i) != T_m) {
        stop(sprintf(
          "`truth[[%d]]` length %d does not match T_m = %d.",
          i, length(mask_i), T_m))
      }
      col_i <- ifelse(mask_i, pal[i], "grey60")
    } else {
      col_i <- ifelse(sig, pal[i], "grey60")
    }

    raw <- -log10(pmax(lfsr_i, .Machine$double.eps))
    points(pos, rep(i, T_m), pch = 1L,
           cex = .lfsr_cex(raw, cex_max), col = col_i, lwd = 2.0)
  }

  if (add_legend) {
    # Both legends sit horizontally just below the title and
    # ABOVE the plot box, in the reserved top margin. `xpd = NA`
    # lets `legend()` draw beyond the panel; positioning uses
    # `par("cxy")` (character size in user coordinates) so the
    # offset is measured in margin lines, not panel fractions.
    op_xpd <- par(xpd = NA); on.exit(par(op_xpd), add = TRUE)
    usr <- par("usr")
    cxy <- par("cxy")
    x_left  <- usr[1L]
    x_mid   <- usr[1L] + 0.55 * (usr[2L] - usr[1L])
    # Bottom of the legend rests ~1.4 char-heights above the
    # top of the panel; the title sits ~2.6 char-heights above
    # that (at line = 4.0).
    y_above <- usr[4L] + 1.4 * cxy[2L]

    size_breaks <- c(1.3, round(cex_max / 2, 1), cex_max)
    legend(x = x_left, y = y_above,
           legend = sprintf("-log10(lfsr) = %g", size_breaks),
           pch = 1L, pt.cex = .lfsr_cex(size_breaks, cex_max),
           pt.lwd = 2.0, col = "black", bty = "n", cex = 0.75,
           horiz = TRUE)

    color_label <- if (!is.null(truth_mask) &&
                       any(!vapply(truth_mask, is.null, logical(1L))))
      c("truly affected", "not affected")
    else
      c(sprintf("lfsr <= %g", lfsr_threshold),
        sprintf("lfsr > %g",  lfsr_threshold))
    legend(x = x_mid, y = y_above,
           legend = color_label,
           pch = 1L, pt.cex = 1.4, pt.lwd = 2.0,
           col = c(pal[1L], "grey60"), bty = "n", cex = 0.75,
           horiz = TRUE)
  }
}

# Map -log10(lfsr) values to a tame visual cex range. The raw
# `-log10(lfsr)` runs from 0 to >> 1, and using it directly as
# `cex` produces dots that dominate the panel at high
# significance and disappear near the threshold. Linearly map
# `[0, cex_max]` to `[0.7, 2.5]` so the legend and panel both
# stay readable.
.lfsr_cex <- function(raw, cex_max,
                      cex_range = c(0.7, 2.5)) {
  raw_c <- pmin(pmax(raw, 0), cex_max)
  cex_range[1L] +
    (cex_range[2L] - cex_range[1L]) * raw_c / cex_max
}

#' Per-CS local false sign rate bubble grid
#'
#' Visualize the per-position local false sign rate (lfsr)
#' produced by `mf_post_smooth(method = "HMM")` as a bubble
#' grid: rows are credible sets, columns are positions, dot
#' size is `-log10(lfsr)` clamped to `[0, cex_max]`. Color
#' marks positions where lfsr is at or below `lfsr_threshold`.
#' For simulated data where the truly affected positions are
#' known, pass `truth` and dots are recolored by the
#' ground-truth mask.
#'
#' Companion to `mfsusie_plot()` for HMM-smoothed fits. The
#' function handles single-outcome (`M = 1`) and multi-outcome
#' (`M > 1`) fits transparently and takes no `m` argument; for
#' `M > 1` it tiles one bubble panel per outcome using the same
#' layout policy as `mfsusie_plot()`.
#'
#' @param fit an `mfsusie` / `fsusie()` fit object that has
#'   been post-processed with `mf_post_smooth(method = "HMM")`
#'   (so `fit$lfsr_curves` is populated).
#' @param lfsr_threshold numeric, color cutoff. Positions with
#'   lfsr at or below `lfsr_threshold` use the CS color; the
#'   rest are grey. Default `0.01`.
#' @param truth optional ground-truth mask for simulated data.
#'   For `M = 1`: a length-`T_1` boolean vector OR a length-K
#'   list of length-`T_1` boolean vectors (one per CS). For
#'   `M > 1`: a length-M list, each entry a length-`T_m`
#'   boolean vector or list-of-CS-vectors.
#' @param cex_max numeric, upper clamp on `-log10(lfsr)` for
#'   the dot size. Default `6`.
#' @param save optional file path. When non-NULL the bubble grid
#'   is written to the file at a size proportional to the number
#'   of outcomes and credible sets. The graphics device is
#'   selected from the file extension: `.pdf`, `.png`,
#'   `.jpg`/`.jpeg`, or `.svg`.
#' @param ... reserved.
#' @return Called for side effect; returns `invisible(NULL)`.
#' @export
mfsusie_plot_lfsr <- function(fit,
                              lfsr_threshold = 0.01,
                              truth          = NULL,
                              cex_max        = 6,
                              add_legend     = TRUE,
                              smooth_method  = NULL,
                              save           = NULL, ...) {
  if (!inherits(fit, "mfsusie")) {
    stop("`fit` must be an `mfsusie` (or `fsusie`) fit object.")
  }
  # `mfsusie_plot_lfsr` only makes sense for HMM-smoothed fits.
  # When the user does not name a smoothing, default to HMM if
  # present; otherwise emit a clear error.
  if (is.null(smooth_method)) {
    if (!is.null(fit$smoothed) && !is.null(fit$smoothed[["HMM"]])) {
      smooth_method <- "HMM"
    }
  }
  if (is.null(smooth_method) ||
      is.null(fit$smoothed[[smooth_method]]$lfsr_curves)) {
    stop("`mfsusie_plot_lfsr()` requires HMM-smoothed lfsr ",
         "curves. Call `mf_post_smooth(method = \"HMM\")` first.")
  }
  smoothed <- fit$smoothed[[smooth_method]]

  M <- length(fit$dwt_meta$T_basis)
  K <- length(fit$sets$cs %||% list())

  # Optional file output: open a sized device and close it on
  # exit. Sizing computed by `mfsusie_plot_lfsr_dimensions()`.
  if (!is.null(save)) {
    .open_save_device(save,
                      mfsusie_plot_lfsr_dimensions(fit, add_legend))
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  # Normalize `truth` into a length-M list of length-K lists of
  # length-T_m boolean vectors (or NULL). Acceptable inputs:
  #   M = 1: bool vec, list of bool vecs (per CS), or NULL.
  #   M > 1: list of (bool vec or list-of-bool-vecs) per outcome,
  #          or NULL.
  normalize_truth <- function(tr, M) {
    if (is.null(tr)) return(replicate(M, NULL, simplify = FALSE))
    if (M == 1L) {
      if (is.list(tr)) return(list(tr))
      return(list(list(tr)))
    }
    if (!is.list(tr) || length(tr) != M) {
      stop("`truth` must be a length-M list when M > 1.")
    }
    lapply(tr, function(x) if (is.list(x)) x else list(x))
  }
  per_outcome_truth <- normalize_truth(truth, M)

  # Top-margin strip holds the size + colour legend (under the
  # title). The panel mar reserves an extra ~3 lines on top
  # when the legend is enabled.
  # Top-margin lines: when the legend strip is shown, reserve
  # enough room for (a) the title at line 4.0, (b) the legend
  # row at ~1.4 char-heights above the panel, plus a small pad.
  # When no legend, the title sits at line 1 and only ~2.5 lines
  # of margin are needed.
  top_mar <- if (isTRUE(add_legend)) 6.5 else 2.5

  if (M == 1L) {
    op <- par(mar = c(4, 5, top_mar, 2))
    on.exit(par(op), add = TRUE)
    .draw_lfsr_bubble(fit, m = 1L,
                      lfsr_threshold = lfsr_threshold,
                      truth_mask = per_outcome_truth[[1L]],
                      cex_max = cex_max,
                      smoothed = smoothed,
                      add_legend = add_legend)
    return(invisible(NULL))
  }

  cols <- ceiling(sqrt(M))
  rows <- ceiling(M / cols)
  op <- par(mfrow = c(rows, cols), mar = c(4, 5, top_mar, 2))
  on.exit(par(op), add = TRUE)
  for (mi in seq_len(M)) {
    .draw_lfsr_bubble(fit, m = mi,
                      lfsr_threshold = lfsr_threshold,
                      truth_mask = per_outcome_truth[[mi]],
                      cex_max = cex_max,
                      smoothed = smoothed,
                      add_legend = add_legend && (mi == 1L))
  }
  remaining <- rows * cols - M
  if (remaining > 0L) {
    for (i in seq_len(remaining)) plot.new()
  }
  invisible(NULL)
}
