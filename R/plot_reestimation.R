# =====================================================================
# plot_reestimation.R
# ---------------------------------------------------------------------
# WHAT THIS FILE DOES
#   Draws, for each declared transition, the rate the user supplied beside
#   the rate recovered from the simulated population: one figure per
#   transition, two panels sharing the x and the y axis.
#
#     [ supplied rate ] [ recovered: bands, 95% CIs, smooth, reference ]
#     [               legend                                           ]
#
# DESIGN
#   Uses base R graphics with Okabe-Ito group colours and a light grid.
#   Figures use the same display rules across transitions:
#   * one y axis for both panels, chosen from the supplied rate, the
#     smooth and the well-supported bands (so one band with tiny exposure
#     cannot squash the plot); bands above it are marked at the top edge;
#   * a logarithmic axis when the values span more than 2.5 orders of
#     magnitude (mortality), bands with no events then marked at the
#     bottom edge;
#   * forced transitions (rate = Inf) shaded where infinite;
#   * empty risk sets, no events and failed fits said in words.
#   * the size of the device is checked before any graphics call;
#     only the graphics parameters that are changed are restored.
#   * colours come from plot_style(), PLOT_STYLE or the user's style.
# =====================================================================


#' Plot supplied against recovered transition rates
#'
#' @description
#' Draws one figure per declared transition. The left panel shows the rate
#' the model was given; the right panel shows the rate recovered from the
#' simulated population ([reestimate()]): the occurrence/exposure rate of
#' each band as a point whose size reflects its person-years, with its
#' exact 95% Poisson interval, the smooth Poisson-spline fit with its 95%
#' band, and the supplied rate as a dashed reference. Both panels share the
#' x and the y axis.
#'
#' @details
#' If the rate depends on a second time scale (age and calendar time, for
#' example), each period of that scale has its own colour in both panels,
#' and the supplied rate is drawn at the middle of each period.
#'
#' The dashed reference in the right panel is the supplied rate itself
#' when it depends only on the time scale on the x axis. Otherwise it is
#' the supplied rate averaged over the exposure in each band, which is
#' exactly what the band estimate estimates.
#'
#' The y axis is logarithmic when the values span more than 2.5 orders of
#' magnitude, unless `log` says otherwise. Points beyond the axis are
#' marked by a triangle at the edge. Shaded regions are where the rate is
#' infinite (a forced transition); there is no hazard to smooth there, and
#' the band estimates reflect the timing rule rather than a rate.
#'
#' Several figures are drawn in turn; on an interactive device the plot
#' waits for confirmation between them (`ask`). To keep them, draw to a
#' file device, e.g. `pdf("rates.pdf", width = 11, height = 5)`, or select
#' one with `which`. A size of about 11 by 5 inches suits them.
#'
#' Before drawing, the size of the device is checked. If the full layout
#' does not fit, a compact one without legend, subtitle and notes is used.
#' If that does not fit either, an error says so and nothing is drawn.
#'
#' Use `panels` and `layout = FALSE` to draw into a grid of your own,
#' e.g. `par(mfrow = c(2, 2))`; `plot(result, which = 1:4, panels = "recovered",
#' layout = FALSE)`.
#'
#' @param x A `micsim_reestimation`, a `micsim_result`, or one element of a
#'   reestimation.
#' @param which Transitions to draw: positions or labels (the names of the
#'   reestimation). All by default.
#' @param log `TRUE` or `FALSE` to force a logarithmic or linear y axis;
#'   `NULL` decides per transition.
#' @param ask Wait between figures? Defaults to `TRUE` on an interactive
#'   device when there is more than one.
#' @param panels "both" (default), "supplied" or "recovered": which panels
#' to draw.
#' @param legend Draw the legend? NULL (default) means TRUE in the full
#' layout and FALSE with layout = FALSE.
#' @param layout TRUE (default) arranges each figure on its own page.
#' FALSE draws the panels into the figure regions already set up,
#' e.g. with par(mfrow = ).
#' @param style Named list of colours overriding the defaults: groups
#' (up to three), truth, points, smooth, forced, grid, axis, muted.
#' @param ... Graphical arguments main and sub (one per figure, recycled),
#' xlab, ylab, xlim, ylim (in rate units, also on a log axis) and cex.
#' Others are ignored with a warning.
#' @return `x`, invisibly.
#'
#' @seealso [reestimate()], [simulate.micsim_model()]
#'
#' @examples
#' mort <- function(age) 0.0005 * exp(0.07 * age)
#' model <- micsim_model(
#'   state_space(sex = c("m", "f")),
#'   transitions(transition(NA, "dead", mort)),
#'   data.frame(ID = 1:500, birthDate = "19400101", initState = "f"),
#'   horizon = c(20000101, 20291231)
#' )
#' result <- simulate(model, seed = 1)
#' plot(result)
#'
#' @export
plot.micsim_reestimation <- function(x, which = NULL, log = NULL, ask = NULL,
                                     panels = c("both", "supplied", "recovered"),
                                     legend = NULL, layout = TRUE,
                                     style = NULL, ...) {
  entries <- unclass(x)
  if (!is.null(which)) {
    if (is.character(which)) {
      unknown <- setdiff(which, names(entries))
      if (length(unknown) > 0L) {
        micsim_abort("micsim_reestimation_error", paste0(
          "No transition labelled ", toString(paste0("\"", unknown, "\"")),
          ".\n  Labels: ", toString(paste0("\"", names(entries), "\"")), "."))
      }
    } else if (!is.numeric(which) || any(which < 1 | which > length(entries))) {
      micsim_abort("micsim_reestimation_error", paste0(
        "`which` must be positions from 1 to ", length(entries),
        " or transition labels."))
    }
    entries <- entries[which]
  }
  if (length(entries) == 0L) return(invisible(x))
  opts <- plot_options(panels, legend, layout, style, log, list(...))
  if (is.null(ask)) {
    ask <- isTRUE(opts$layout) && length(entries) > 1L &&
      grDevices::dev.interactive()
  }
  if (isTRUE(ask)) {
    old_ask <- grDevices::devAskNewPage(TRUE)
    on.exit(grDevices::devAskNewPage(old_ask), add = TRUE)
  }
  # `main` and `sub` may give one title per figure.
  drawn <- character(length(entries))
  for (k in seq_along(entries)) {
    drawn[k] <- plot_rate_check(entries[[k]], log = log,
                                opts = figure_options(opts, k))
  }
  note_compact(sum(drawn == "compact"), length(entries))
  invisible(x)
}


#' @rdname plot.micsim_reestimation
#' @export
plot.micsim_result <- function(x, which = NULL, log = NULL, ask = NULL,
                               panels = c("both", "supplied", "recovered"),
                               legend = NULL, layout = TRUE, style = NULL,
                               ...) {
  r <- x$reestimation
  if (is.null(r)) r <- reestimate(x)
  plot(r, which = which, log = log, ask = ask, panels = panels,
       legend = legend, layout = layout, style = style, ...)
  invisible(x)
}


#' @rdname plot.micsim_reestimation
#' @export
plot.micsim_rate_check <- function(x, log = NULL,
                                   panels = c("both", "supplied", "recovered"),
                                   legend = NULL, layout = TRUE,
                                   style = NULL, ...) {
  opts <- plot_options(panels, legend, layout, style, log, list(...))
  drawn <- plot_rate_check(x, log = log, opts = figure_options(opts, 1L))
  note_compact(sum(drawn == "compact"), 1L)
  invisible(x)
}


#' Print one transition's reestimation
#' @param x A `micsim_rate_check`.
#' @param ... Not used.
#' @return `x`, invisibly.
#' @export
print.micsim_rate_check <- function(x, ...) {
  cat("<micsim_rate_check> ", x$label, "\n", sep = "")
  cat("  rate       : ", rate_signature(x), "\n", sep = "")
  cat("  shown along: ", axis_label(x$axis),
      if (!is.na(x$group_by)) paste0(", by ", scale_name(x$group_by)), "\n", sep = "")
  cat("  events     : ", format(x$events, big.mark = ","), " in ",
      format(round(x$exposure), big.mark = ","), " person-years\n", sep = "")
  if (!is.na(x$ratio)) {
    cat("  expected   : ", format(round(x$expected, 1), big.mark = ","),
        " (observed/expected ", sprintf("%.3f", x$ratio), ", 95% CI ",
        sprintf("%.3f", x$ratio_ci[1L]), " to ", sprintf("%.3f", x$ratio_ci[2L]),
        ")\n", sep = "")
  }
  for (n in x$notes) cat("  note       : ", n, "\n", sep = "")
  invisible(x)
}


# ---------------------------------------------------------------------
# The figure
# ---------------------------------------------------------------------

#' Colours and constants of the figure
#' @keywords internal
#' @noRd
PLOT_STYLE <- list(
  groups  = c("#0072B2", "#D55E00", "#009E73"),
  truth   = "#1A1A1A",
  points  = "#0072B2",
  smooth  = "#D55E00",
  forced  = "#E69F00",
  grid    = "#E6E6E6",
  axis    = "#4D4D4D",
  muted   = "#6B6B6B"
)


#' Draw one transition's figure
#'
#' @param e A `micsim_rate_check`.
#' @param log `TRUE`, `FALSE` or `NULL`.
#' @param opts Options from figure_options(plot_options(...))
#'
#' @return "full", "compact" or "into regions", invisibly.
#'
#' @keywords internal
#' @noRd
plot_rate_check <- function(e, log = NULL,
                            opts = figure_options(plot_options(), 1L)) {
  # The style for this figure, read by the drawing helpers through
  # plot_style(), and put back afterwards.
  old_style <- set_plot_style(opts$style)
  on.exit(set_plot_style(old_style), add = TRUE)
  st <- plot_style()

  n_groups <- nrow(e$groups)
  grouped <- n_groups > 1L
  xlim <- if (is.null(opts$labs$xlim)) x_window(e) else opts$labs$xlim
  scale <- y_scale(e, log, xlim = xlim)
  if (!is.null(opts$labs$ylim)) scale <- user_ylim(scale, opts$labs$ylim)
  labs <- figure_labels(e, scale, opts$labs)
  n_panels <- if (opts$panels == "both") 2L else 1L

  if (!isTRUE(opts$layout)) {
    # Draw into the figure regions the caller has set up, e.g. with
    # par(mfrow = ): only the axis style is set here, and restored.
    check_room(n_panels, compact = TRUE, margins_in = graphics::par("mai"),
               region_in = next_figure_size())
    op <- graphics::par(c("mgp", "tcl", "las"))
    on.exit(graphics::par(op), add = TRUE)
    graphics::par(mgp = c(2.6, 0.6, 0), tcl = -0.3, las = 1)
    draw_panels(e, opts, xlim, scale, labs, grouped, n_groups,
                header_line = 0.3, first_title = labs$main)
    if (isTRUE(opts$legend)) {
      inner_legend(legend_items(e, grouped, scale))
    }
    return(invisible("into regions"))
  }

  # Plan the figure before any graphics call, so that a device too small
  # for it is reported, or a compact layout used, without leaving a half-
  # drawn figure in the device's display list.
  cex <- if (is.null(opts$labs$cex)) 1 else opts$labs$cex
  key <- if (isTRUE(opts$legend)) legend_plan(legend_items(e, grouped, scale), cex = 0.8 * cex) else NULL
  notes <- wrapped_notes(e$notes, cex)
  line_in <- graphics::par("cin")[2L] * cex
  key_in <- if (is.null(key)) 0 else line_in * (1.5 + key$rows * 0.95 + if (grouped) 2.5 else 0)
  top_lines <- 4.2 + 1.1 * length(notes)
  full_fits <- fits_device(n_panels, line_in,
                           width_lines = 5.2 + 1.2, height_lines = top_lines + 4.2 + 2.4,
                           extra_height_in = key_in)
  compact <- !full_fits
  if (compact) {
    check_room(n_panels, compact = TRUE,
               margins_in = line_in * c(4, 4.6, 2.6, 0.8),
               region_in = graphics::par("din") / c(n_panels, 1))
  }

  # Restore only what is changed here; restoring `mfrow` also ends the
  # layout, as filled.contour() does.
  op <- graphics::par(c("mar", "oma", "mgp", "tcl", "las", "family", "cex",
                        "mfrow"))
  on.exit(graphics::par(op), add = TRUE)

  if (compact) {
    graphics::par(mfrow = c(1L, n_panels), oma = c(0, 0, 1.4, 0),
                  mar = c(4, 4.6, 2.6, 0.8), mgp = c(2.6, 0.6, 0),
                  tcl = -0.3, las = 1, family = "sans", cex = cex)
    draw_panels(e, opts, xlim, scale, labs, grouped, n_groups,
                header_line = 0.5)
    graphics::mtext(labs$main, side = 3, outer = TRUE, line = 0.2,
                    cex = 1.1, adj = 0.02)
    return(invisible("compact"))
  }

  if (is.null(key)) {
    graphics::layout(matrix(seq_len(n_panels), nrow = 1L))
  } else {
    graphics::layout(
      matrix(c(seq_len(n_panels), rep(n_panels + 1L, n_panels)),
             nrow = 2L, byrow = TRUE),
      # The legend's share of the device height, as a fraction rather than
      # an absolute height: if the device is later redrawn smaller (an
      # IDE plot pane being resized), the strip shrinks with it instead
      # of leaving the panels no room.
      heights = c(1 - key_in / graphics::par("din")[2L],
                  key_in / graphics::par("din")[2L]))
  }
  graphics::par(oma = c(0, 0, top_lines, 0), mgp = c(2.6, 0.6, 0),
                tcl = -0.3, las = 1, family = "sans", cex = cex,
                mar = c(4.2, 5.2, 2.4, 1.2))
  draw_panels(e, opts, xlim, scale, labs, grouped, n_groups,
              header_line = 0.7)

  # -- title, subtitle and notes ----------------------------------------
  top <- 1.1 * length(notes)
  graphics::mtext(labs$main, side = 3, outer = TRUE, line = 2.2 + top,
                  cex = 1.3, adj = 0.02)
  graphics::mtext(labs$sub, side = 3, outer = TRUE, line = 0.7 + top,
                  cex = 0.85, adj = 0.02, col = st$muted)
  for (i in seq_along(notes)) {
    graphics::mtext(notes[i], side = 3, outer = TRUE,
                    line = 0.5 + 1.1 * (length(notes) - i), cex = 0.8,
                    adj = 0.02, col = "#9A3D00", font = 3)
  }

  # -- legend -----------------------------------------------------------
  if (!is.null(key)) {
    graphics::par(mar = c(0, 0, 0, 0))
    graphics::plot.new()
    draw_legend(e, grouped, key)
  }
  invisible("full")
}


#' The panels of one figure, in order: supplied rate, recovered rate.
#'
#' @description
#' Each panel starts with plot.new(), so it goes into the next figure
#' region of whatever layout is in effect. header_line places the panel
#' headers; first_title, if given, is drawn above the first panel (used
#' when the figure has no outer title, with layout = FALSE).
#'
#' @param e The transition-rate comparison.
#' @param opts The options to use.
#' @param xlim The x range to use.
#' @param scale The scale of the axis.
#' @param labs The labels to use.
#' @param grouped Whether the figure is grouped.
#' @param n_groups The number of groups.
#' @param header_line The line to place the panel headers.
#' @param first_title The title to draw above the first panel.
#'
#' @return Invisible(NULL).
#'
#' @keywords internal
#' @noRd
draw_panels <- function(e, opts, xlim, scale, labs, grouped, n_groups,
                        header_line, first_title = NULL) {
  st <- plot_style()
  col_of <- function(g) if (grouped) st$groups[g] else NA

  # -- left: supplied rate ----------------------------------------------
  if (opts$panels %in% c("both", "supplied")) {
    frame_panel(xlim, scale, e$axis, xlab = labs$xlab, ylab = labs$ylab)
    put_title(first_title)
    first_title <- NULL
    shade_forced(e$truth, scale)
    if (shows_range(e, scale)) {
      tr <- e$truth
      draw_ribbon(tr$x, tr$lower, tr$upper, scale,
                  col = with_alpha(st$truth, 0.12), floor_missing = FALSE)
    }
    for (g in seq_len(n_groups)) {
      tr <- e$truth[e$truth$group == g, , drop = FALSE]
      draw_curve(tr$x, tr$rate, scale,
                 col = if (grouped) col_of(g) else st$truth, lwd = 2.2)
    }
    left_note <- if (all(is.na(e$truth$rate))) "could not be evaluated" else
      if (length(e$varying) > 0L) {
        paste0("averaged over the person-years at risk",
               if (grouped) " in each period" else "")
      } else NULL
    panel_header("Supplied rate", left_note, line = header_line)
    if (all(is.na(e$truth$rate))) {
      centre_note("The supplied rate could not be evaluated.")
    }
  }

  # -- right: recovered rate --------------------------------------------
  if (opts$panels %in% c("both", "recovered")) {
    frame_panel(xlim, scale, e$axis, xlab = labs$xlab, ylab = labs$ylab)
    put_title(first_title)
    shade_forced(e$truth, scale)
    b <- e$bands
    max_exp <- max(c(b$exposure, 0), na.rm = TRUE)
    for (g in seq_len(n_groups)) {
      bg <- b[b$group == g & b$exposure > 0, , drop = FALSE]
      pc <- if (grouped) col_of(g) else st$points
      sc <- if (grouped) col_of(g) else st$smooth
      rc <- if (grouped) col_of(g) else st$truth
      # smooth band and curve below the points
      sm <- e$smooth[e$smooth$group == g, , drop = FALSE]
      if (nrow(sm) > 0L) {
        draw_ribbon(sm$x, sm$lower, sm$upper, scale, col = with_alpha(sc, 0.16))
      }
      # interval whiskers, then points
      if (nrow(bg) > 0L) {
        # A band with no events where the supplied rate is zero says nothing
        # (a structural zero); it gets neither a whisker nor a marker.
        ref <- band_reference(e, bg)
        bg <- bg[bg$events > 0L | is.na(ref) | ref > 0, , drop = FALSE]
      }
      if (nrow(bg) > 0L) {
        # Intervals only for bands with events: for a band without any, the
        # upper limit is 3.7 / person-years, which says only that the band
        # was small, and drawn for every such band it buries the figure.
        # Fainter for bands with fewer person-years, for the same reason.
        wk <- bg[bg$events > 0L, , drop = FALSE]
        if (nrow(wk) > 0L) {
          lo <- to_axis(wk$ci_lower, scale); hi <- to_axis(wk$ci_upper, scale)
          lo[is.na(lo)] <- scale$lim[1L]
          a <- 0.12 + 0.28 * sqrt(wk$exposure / max_exp)
          graphics::segments(wk$mid, lo, wk$mid, hi,
                             col = vapply(a, function(z) with_alpha(pc, z), character(1L)),
                             lwd = 1)
        }
        draw_points(bg, scale, pc, max_exp)
      }
      if (nrow(sm) > 0L) draw_curve(sm$x, sm$rate, scale, col = sc, lwd = 2.4)
      # reference: the supplied rate
      reference_line(e, g, scale, rc)
    }
    panel_header("Recovered from the simulation", NULL, line = header_line)
    if (any(grepl("^Nobody", e$notes))) {
      centre_note("Nobody was at risk of this transition.")
    } else if (any(grepl("^No events", e$notes))) {
      centre_note("No events were observed.")
    }
  }
  invisible(NULL)
}


# ---------------------------------------------------------------------
# Options, labels and sizes
# ---------------------------------------------------------------------

#' A figure title above the current panel, left-aligned; nothing for NULL.
#'
#' @param title The title to plot.
#'
#' @return Invisible(NULL).
#'
#' @keywords internal
#' @noRd
put_title <- function(title) {
  if (!is.null(title)) {
    graphics::title(main = title, line = 1.6, cex.main = 1, adj = 0)
  }
  invisible(NULL)
}


#' Graphical arguments `plot()` accepts through `...`.
#' 
#' @keywords internal
#' @noRd
PLOT_DOTS <- c("main", "sub", "xlab", "ylab", "xlim", "ylim", "cex")

#' Check and collect the plot options shared by all figures of a call.
#'
#' @param panels The panels to plot.
#' @param legend Whether to plot a legend. legend = NULL means TRUE with layout = TRUE, FALSE otherwise.
#' @param layout Whether to layout the figure.
#' @param style The style to use.
#' @param log Whether to plot on a logarithmic scale.
#' @param dots The dots to use. Unknown dots are ignored.
#'
#' @return A list with the plot options.
#'
#' @keywords internal
#' @noRd
plot_options <- function(panels = c("both", "supplied", "recovered"),
                         legend = NULL, layout = TRUE, style = NULL,
                         log = NULL, dots = list()) {
  abort_plot <- function(msg) micsim_abort("micsim_plot_error", msg)
  # By default the legend is drawn in the full layout, where it has a strip
  # of its own, and left out when drawing into the caller's regions, where
  # it would cover the data.
  if (is.null(legend)) legend <- isTRUE(layout)
  panels <- tryCatch(match.arg(panels),
                     error = function(e) abort_plot(
                       "`panels` must be \"both\", \"supplied\" or \"recovered\"."))
  for (flag in c("legend", "layout")) {
    v <- get(flag)
    if (!isTRUE(v) && !isFALSE(v)) abort_plot(paste0("`", flag, "` must be TRUE or FALSE."))
  }
  if (!is.null(log) && !isTRUE(log) && !isFALSE(log)) {
    abort_plot("`log` must be TRUE, FALSE or NULL.")
  }
  nms <- names(dots)
  if (length(dots) > 0L && (is.null(nms) || any(!nzchar(nms)))) {
    abort_plot("Arguments passed through `...` must be named, e.g. main = \"Mortality\".")
  }
  ignored <- setdiff(nms, PLOT_DOTS)
  if (length(ignored) > 0L) {
    warning(micsim_warning("micsim_plot_warning", paste0(
      "Ignoring argument(s) ", toString(paste0("`", ignored, "`")), ".\n",
      "  plot() takes ", toString(paste0("`", PLOT_DOTS, "`")),
      " through `...`; colours\n  are set with `style`.")))
  }
  labs <- dots[intersect(nms, PLOT_DOTS)]
  for (lim in c("xlim", "ylim")) {
    v <- labs[[lim]]
    if (!is.null(v) && !(is.numeric(v) && length(v) == 2L && all(is.finite(v)) && v[1L] < v[2L])) {
      abort_plot(paste0("`", lim, "` must be two increasing finite numbers."))
    }
  }
  if (!is.null(labs$cex) && !(is.numeric(labs$cex) && length(labs$cex) == 1L &&
                                is.finite(labs$cex) && labs$cex > 0)) {
    abort_plot("`cex` must be a single positive number.")
  }
  list(panels = panels, legend = legend, layout = layout,
       style = check_plot_style(style), labs = labs)
}

#' The options of figure k: `main` and `sub` may hold one title per figure
#' and are recycled.
#'
#' @param opts The options to use.
#' @param k The figure number.
#'
#' @return The options with the titles and subtitles recycled.
#'
#' @keywords internal
#' @noRd
figure_options <- function(opts, k) {
  for (f in c("main", "sub")) {
    v <- opts$labs[[f]]
    if (!is.null(v) && !is.language(v) && length(v) > 1L) {
      opts$labs[[f]] <- v[[(k - 1L) %% length(v) + 1L]]
    }
  }
  opts
}

#' Titles and axis labels of a figure, the user's where given.
#'
#' @param e The transition-rate comparison.
#' @param scale The scale of the axis.
#' @param labs The labels to use.
#'
#' @return A list with the main title, subtitle, x-axis label and y-axis label.
#'
#' @keywords internal
#' @noRd
figure_labels <- function(e, scale, labs) {
  pick <- function(f, default) if (is.null(labs[[f]])) default else labs[[f]]
  list(
    main = pick("main", plot_title(e)),
    sub  = pick("sub", subtitle(e)),
    xlab = pick("xlab", axis_label(e$axis)),
    ylab = pick("ylab", if (scale$log) "Rate per person-year (log scale)" else
      "Rate per person-year")
  )
}

#' A user-given y range, on the axis's own scale.
#'
#' @param scale The scale of the axis.
#' @param ylim The y range to set.
#'
#' @return The scale with the y range set.
#'
#' @keywords internal
#' @noRd
user_ylim <- function(scale, ylim) {
  if (scale$log) {
    if (any(ylim <= 0)) {
      micsim_abort("micsim_plot_error",
                   "On a logarithmic axis, `ylim` must be positive. Use log = FALSE for a linear axis.")
    }
    scale$lim <- log10(ylim)
  } else {
    scale$lim <- ylim
  }
  scale
}

#' The notes of a figure, wrapped to the device width.
#'
#' @param notes The notes to wrap.
#' @param cex The character expansion factor.
#'
#' @return The wrapped notes.
#'
#' @keywords internal
#' @noRd
wrapped_notes <- function(notes, cex = 1) {
  if (length(notes) == 0L) return(character(0))
  char_in <- graphics::strwidth("n", units = "inches", cex = 0.8 * cex)
  width <- max(20L, floor(0.95 * graphics::par("din")[1L] / char_in))
  strwrap(paste(notes, collapse = "  "), width = width)
}

#' Does the full figure fit on the current device?
#'
#' @description
#' n_panels side by side, each needing width_lines of margin and 1.0
#' inch of plotting width; height_lines of margins, extra_height_in
#' (the legend strip) and 0.8 inch of plotting height.
#' 
#' @param n_panels The number of panels in the figure.
#' @param line_in The height of one line of margin text in inches.
#' @param width_lines The margin lines one panel needs horizontally.
#' @param height_lines The margin lines the figure needs vertically.
#' @param extra_height_in The extra height in inches.
#'
#' @return `TRUE` if the figure fits on the device, otherwise `FALSE`.
#'
#' @keywords internal
#' @noRd
fits_device <- function(n_panels, line_in, width_lines, height_lines,
                        extra_height_in = 0) {
  din <- graphics::par("din")
  need_w <- n_panels * (width_lines * line_in + 1.0)
  need_h <- height_lines * line_in + extra_height_in + 0.8
  din[1L] >= need_w && din[2L] >= need_h
}

#' Stop, before anything is drawn, if a panel's region cannot hold its
#' margins and a minimal plotting area.
#'
#' @param n_panels The number of panels in the figure.
#' @param compact Whether the figure is drawn in the compact layout.
#' @param margins_in c(bottom, left, top, right) in inches.
#' @param region_in width and height of the region one panel gets, in inches.
#'
#' @return `TRUE` invisibly or raises `micsim_plot_error` if the figure does
#' not fit on the device.
#'
#' @keywords internal
#' @noRd
check_room <- function(n_panels, compact, margins_in, region_in) {
  need <- c(margins_in[2L] + margins_in[4L], margins_in[1L] + margins_in[3L]) + 0.6
  if (any(region_in < need)) {
    micsim_abort("micsim_plot_error", paste0(
      "The plot device is too small to draw this figure (",
      n_panels, " panel(s) of ", sprintf("%.1f x %.1f", region_in[1L], region_in[2L]),
      " in; at least ", sprintf("%.1f x %.1f", need[1L], need[2L]),
      " in each are needed).\n",
      "  Enlarge the plot window, draw a single panel with panels = \"recovered\",\n",
      "  or draw to a file, e.g. pdf(\"rates.pdf\", width = 11, height = 5)."))
  }
  invisible(TRUE)
}

#' The size, in inches, of the next figure region of an mfrow/mfcol grid.
#' For a layout() set up by the caller this is an approximation.
#'
#' @keywords internal
#' @noRd
next_figure_size <- function() {
  mf <- graphics::par("mfrow")
  graphics::par("din") / rev(mf)
}

#' Say, once per plot() call, that figures were drawn in the compact
#' layout because the device was too small for the full one.
#'
#' @keywords internal
#' @noRd
note_compact <- function(n_compact, n_figures) {
  if (n_compact == 0L) return(invisible(NULL))
  message(structure(
    class = c("micsim_plot_message", "message", "condition"),
    list(message = paste0(
      if (n_figures == 1L) "The figure was" else
        paste0(n_compact, " of ", n_figures, " figures were"),
      " drawn without the legend, subtitle and notes,\n",
      "  because the plot device is too small for the full layout. Enlarge\n",
      "  the device, or draw to a file, e.g.\n",
      "  pdf(\"rates.pdf\", width = 11, height = 5).\n"),
      call = NULL)))
  invisible(NULL)
}


#' A legend inside the current panel, for layout = FALSE.
#'
#' @keywords internal
#' @noRd
inner_legend <- function(items) {
  get <- function(f) vapply(items, function(it) as.character(it[[f]]), character(1L))
  num <- function(f) as.numeric(get(f))
  graphics::legend(
    "topright", legend = get("label"), col = get("col"), lty = num("lty"),
    lwd = num("lwd"), pch = num("pch"), pt.bg = get("bg"),
    pt.cex = 0.8 * num("cex"), bty = "n", cex = 0.6, seg.len = 1.8,
    inset = 0.01)
}

.plot_style <- new.env(parent = emptyenv())

#' The style in effect: PLOT_STYLE, or the user's style merged into it.
#' 
#' @keywords internal
#' @noRd
plot_style <- function() {
  s <- .plot_style$current
  if (is.null(s)) PLOT_STYLE else s
}
#' Set the plot style.
#'
#' @keywords internal
#' @noRd
set_plot_style <- function(style) {
  old <- .plot_style$current
  .plot_style$current <- style
  old
}

#' Check a user's `style` and merge it into PLOT_STYLE. NULL means the
#' default style.
#'
#' @keywords internal
#' @noRd
check_plot_style <- function(style) {
  if (is.null(style)) return(NULL)
  if (!is.list(style) || is.null(names(style)) || any(!nzchar(names(style)))) {
    micsim_abort("micsim_plot_error", paste0(
      "`style` must be a named list of colours, e.g. list(smooth = \"darkred\").\n",
      "  Names: ", toString(names(PLOT_STYLE)), "."))
  }
  unknown <- setdiff(names(style), names(PLOT_STYLE))
  if (length(unknown) > 0L) {
    micsim_abort("micsim_plot_error", paste0(
      "Unknown `style` element(s) ", toString(unknown), ".\n",
      "  Names: ", toString(names(PLOT_STYLE)), "."))
  }
  for (n in names(style)) {
    ok <- tryCatch({ grDevices::col2rgb(style[[n]]); TRUE }, error = function(e) FALSE)
    if (!ok || length(style[[n]]) == 0L) {
      micsim_abort("micsim_plot_error", paste0("`style$", n, "` is not a colour."))
    }
  }
  if (!is.null(style$groups) && length(style$groups) < 3L) {
    style$groups <- rep_len(style$groups, 3L)
  }
  out <- PLOT_STYLE
  out[names(style)] <- style
  out
}


#' The x range worth showing
#'
#' Where the supplied rate is positive (or infinite) or an event occurred,
#' padded by a band on each side and kept within the range observed. A
#' fertility rate, zero outside ages 15 to 45 in a risk set aged 0 to 100,
#' is then shown over the ages where it acts.
#' 
#' @param e A transition-rate comparison.
#' @return Length-two numeric axis range.
#' @keywords internal
#' @noRd
x_window <- function(e) {
  full <- range(c(e$bands$lower, e$bands$upper))
  t <- e$truth
  active <- c(t$x[!is.na(t$rate) & t$rate > 0],
              e$bands$lower[e$bands$events > 0], e$bands$upper[e$bands$events > 0])
  if (length(active) == 0L) return(full)
  w <- max(e$bands$upper - e$bands$lower)
  win <- range(active) + c(-1, 1) * max(w, 0.03 * diff(range(active)))
  win <- c(max(full[1L], win[1L]), min(full[2L], win[2L]))
  if (diff(win) <= 0) full else win
}


#' Decide the y axis: linear or log, and its limits
#' @param e A transition-rate comparison.
#' @param log `TRUE`, `FALSE`, or `NULL` for automatic choice.
#' @param xlim The visible x range; NULL uses x_window(e)
#' @return List with `log` and `lim`; limits are log10 values for a log axis.
#' @keywords internal
#' @noRd
y_scale <- function(e, log = NULL, xlim = NULL) {
  # The axis is chosen from what is visible: the user's x range if given.
  win <- if (is.null(xlim)) x_window(e) else xlim
  b <- e$bands[e$bands$upper > win[1L] & e$bands$lower < win[2L], , drop = FALSE]
  inside <- e$truth$x >= win[1L] & e$truth$x <= win[2L]
  truth <- e$truth$rate[inside & is.finite(e$truth$rate)]
  if (has_range(e)) {
    truth <- c(truth, e$truth$upper[inside & is.finite(e$truth$upper)])
  }
  # Use fitted values to set the axis only when their upper/lower
  # interval ratio is below four. Wide spline intervals should not
  # compress the band estimates and reference curves.
  firm <- is.finite(e$smooth$upper) & is.finite(e$smooth$lower) &
    e$smooth$lower > 0 & e$smooth$upper / e$smooth$lower < 4
  sm <- e$smooth$rate[firm]
  ref <- reference_values(e, b)
  supported <- b$exposure > 0 &
    (b$events >= 3L | b$exposure >= 0.25 * stats::median(b$exposure[b$exposure > 0]))
  dots <- b$rate[supported & !is.na(b$rate)]
  if (e$forced) dots <- b$rate[!is.na(b$rate)]
  other <- c(truth, ref, dots)
  other <- other[is.finite(other)]
  # and even then at most half again above everything else
  if (length(other) > 0L && max(other) > 0) sm <- pmin(sm, 1.5 * max(other))
  vals <- c(other, sm)
  vals <- vals[is.finite(vals)]
  pos <- vals[vals > 0]
  hi <- if (length(vals) == 0L) 0 else max(vals)
  if (!is.finite(hi) || hi <= 0) hi <- if (length(pos)) max(pos) else 1
  use_log <- if (!is.null(log)) isTRUE(log) else
    length(pos) >= 2L && max(pos) / min(pos) > 10^2.5
  if (use_log && length(pos) == 0L) use_log <- FALSE
  if (use_log) {
    lo <- max(min(pos), hi * 1e-6)
    list(log = TRUE, lim = c(log10(lo) - 0.08 * log10(hi / lo) - 0.02,
                             log10(hi) + 0.06 * log10(hi / lo) + 0.02))
  } else {
    list(log = FALSE, lim = c(0, hi * 1.07))
  }
}


#' Values on the axis's own scale (log10 when logarithmic)
#' @param v Numeric rate values.
#' @param scale Axis settings from `y_scale()`.
#' @return Values unchanged on a linear axis, or log10 values with
#'   non-positive and missing inputs represented as `NA`.
#' @keywords internal
#' @noRd
to_axis <- function(v, scale) {
  if (!scale$log) return(v)
  out <- rep(NA_real_, length(v))
  ok <- !is.na(v) & v > 0
  out[ok] <- log10(v[ok])
  out
}


#' An empty panel with grid and axes
#' @param xlim Length-two x-axis limits.
#' @param scale Y-axis settings from `y_scale()`.
#' @param axis Clock name for the x-axis label.
#' @param xlab,ylab Axis titles; NULL ylab chooses one from the axis.
#' @return Invisibly returns the result of drawing the axis title.
#' @keywords internal
#' @noRd
frame_panel <- function(xlim, scale, axis, xlab = axis_label(axis),
                        ylab = NULL) {
  st <- plot_style()
  graphics::plot.new()
  graphics::plot.window(xlim = xlim, ylim = scale$lim, xaxs = "i", yaxs = "i")
  xat <- pretty(xlim, n = 7)
  xat <- xat[xat >= xlim[1L] & xat <= xlim[2L]]
  yt <- y_ticks(scale)
  graphics::abline(v = xat, h = yt$at, col = st$grid, lwd = 0.8)
  graphics::axis(1, at = xat, col = NA, col.ticks = st$axis, col.axis = st$axis,
                 cex.axis = 0.85)
  graphics::axis(2, at = yt$at, labels = yt$labels, col = NA,
                 col.ticks = st$axis, col.axis = st$axis, cex.axis = 0.85)
  graphics::box(bty = "l", col = st$axis, lwd = 0.8)
  if (is.null(ylab)) {
    ylab <- if (scale$log) "Rate per person-year (log scale)" else
      "Rate per person-year"
  }
  graphics::title(xlab = xlab, col.lab = st$axis, cex.lab = 0.95)
  # The y label sits inside the left margin, whatever its width.
  graphics::title(ylab = ylab, col.lab = st$axis, cex.lab = 0.95,
                  line = min(3.9, graphics::par("mar")[2L] - 1))
}


#' Tick positions and labels of the y axis
#' @param scale Y-axis settings.
#' @return List with numeric positions `at` and formatted `labels`.
#' @keywords internal
#' @noRd
y_ticks <- function(scale) {
  if (!scale$log) {
    at <- pretty(scale$lim, n = 6)
    at <- at[at >= scale$lim[1L] & at <= scale$lim[2L]]
    return(list(at = at, labels = format_rate(at)))
  }
  e <- seq(floor(scale$lim[1L]), ceiling(scale$lim[2L]))
  decades <- length(e)
  mult <- if (decades <= 3L) c(1, 2, 5) else 1
  at <- log10(as.vector(outer(mult, 10^e)))
  at <- sort(at[at >= scale$lim[1L] & at <= scale$lim[2L]])
  if (length(at) < 2L) {
    at <- pretty(scale$lim, n = 4)
    at <- at[at >= scale$lim[1L] & at <= scale$lim[2L]]
  }
  list(at = at, labels = format_rate(10^at))
}


#' Format rates for an axis: plain decimals, scientific when tiny
#' @param v Numeric rate values.
#' @return Character vector of labels.
#' @keywords internal
#' @noRd
format_rate <- function(v) {
  vapply(v, function(z) {
    if (z == 0) return("0")
    if (abs(z) < 1e-4) return(formatC(z, format = "e", digits = 0))
    format(signif(z, 6), scientific = FALSE, drop0trailing = TRUE, trim = TRUE)
  }, character(1L))
}


#' Axis title of a time scale
#' @param axis Internal clock name.
#' @return Axis-label string for a recognised clock.
#' @keywords internal
#' @noRd
axis_label <- function(axis) {
  switch(axis,
         age = "Age (years)",
         calTime = "Calendar time (year)",
         duration = "Time since entering the state (years)")
}


#' Draw a line, breaking it at missing or infinite values
#' @param x,y Curve coordinates in original units.
#' @param scale Y-axis settings.
#' @param col,lwd,lty Line colour, width and type.
#' @return `NULL`, invisibly, after drawing.
#' @keywords internal
#' @noRd
draw_curve <- function(x, y, scale, col, lwd = 2, lty = 1) {
  yy <- to_axis(y, scale)
  yy[is.infinite(y)] <- NA_real_
  graphics::lines(x, yy, col = col, lwd = lwd, lty = lty, lend = "round",
                  ljoin = "round")
}


#' A translucent band between two curves
#' @param x Axis values.
#' @param lo,hi Lower and upper bounds in rate units.
#' @param scale Y-axis settings.
#' @param col Fill colour, including transparency if required.
#' @param floor_missing Place missing lower bounds at the panel floor;
#'   otherwise draw separate runs with complete bounds.
#' @return `NULL`, invisibly, after drawing.
#' @keywords internal
#' @noRd
draw_ribbon <- function(x, lo, hi, scale, col, floor_missing = TRUE) {
  ylo <- to_axis(lo, scale); yhi <- to_axis(hi, scale)
  # On a log axis a lower bound of zero has no position: rest it on the
  # floor of the panel. Elsewhere a missing bound means no ribbon there.
  zero_lo <- !is.na(lo) & lo <= 0
  ylo[zero_lo | (floor_missing & is.na(ylo))] <- scale$lim[1L]
  yhi[!is.na(yhi)] <- pmin(yhi[!is.na(yhi)], scale$lim[2L])
  if (!floor_missing) {
    ok <- !is.na(ylo) & !is.na(yhi)
    runs <- rle(ok); ends <- cumsum(runs$lengths); starts <- ends - runs$lengths + 1L
    for (k in which(runs$values)) {
      idx <- starts[k]:ends[k]
      if (length(idx) < 2L) next
      graphics::polygon(c(x[idx], rev(x[idx])), c(ylo[idx], rev(yhi[idx])),
                        col = col, border = NA)
    }
    return(invisible(NULL))
  }
  yhi <- pmin(yhi, scale$lim[2L])
  ok <- !is.na(yhi)
  if (sum(ok) < 2L) return(invisible(NULL))
  graphics::polygon(c(x[ok], rev(x[ok])), c(ylo[ok], rev(yhi[ok])),
                    col = col, border = NA)
}


#' The band estimates as points sized by exposure
#'
#' Points beyond the top of the axis become triangles at the top edge; on
#' a log axis, bands with no events become open triangles at the bottom.
#' @param bg Bands for one group, restricted to positive exposure.
#' @param scale Y-axis settings.
#' @param col Point colour.
#' @param max_exp Largest band exposure across groups, for shared sizing.
#' @return Called for its plotting side effects.
#' @keywords internal
#' @noRd
draw_points <- function(bg, scale, col, max_exp) {
  size <- 0.45 + 1.35 * sqrt(bg$exposure / max_exp)
  y <- to_axis(bg$rate, scale)
  top <- !is.na(y) & y > scale$lim[2L]
  zero <- bg$events == 0L
  normal <- !is.na(y) & !top & !(scale$log & zero)
  graphics::points(bg$mid[normal], y[normal], pch = 21, cex = size[normal],
                   bg = with_alpha(col, 0.45), col = with_alpha(col, 0.75), lwd = 0.6)
  edge <- diff(scale$lim) * 0.012
  if (any(top)) {
    graphics::points(bg$mid[top], rep(scale$lim[2L] - edge, sum(top)),
                     pch = 24, cex = 0.8, bg = with_alpha(col, 0.6), col = col,
                     lwd = 0.6)
  }
  if (scale$log && any(zero)) {
    graphics::points(bg$mid[zero], rep(scale$lim[1L] + edge, sum(zero)),
                     pch = 25, cex = 0.55, col = with_alpha(plot_style()$muted, 0.6),
                     lwd = 0.7)
  }
}


#' Is the range of the supplied rate over the other time scales drawn?
#'
#' Only with a single group (with several, the colours already show the
#' dependence) and only where the range is not negligible.
#' @param e A transition-rate comparison.
#' @return Logical scalar; visibility is checked by `shows_range()`.
#' @keywords internal
#' @noRd
has_range <- function(e) {
  if (length(e$varying) == 0L || nrow(e$groups) > 1L) return(FALSE)
  t <- e$truth
  !is.null(t$lower) && any(is.finite(t$lower) & is.finite(t$upper))
}

#' @rdname has_range
#' @param scale The y axis, from `y_scale()`: the range is drawn only where
#'   it would be visible, at least 1.5% of the axis wide somewhere.
#' @keywords internal
#' @noRd
shows_range <- function(e, scale) {
  if (!has_range(e)) return(FALSE)
  t <- e$truth
  lo <- to_axis(t$lower, scale); hi <- to_axis(t$upper, scale)
  if (scale$log) lo[!is.na(t$lower) & t$lower <= 0] <- scale$lim[1L]
  ok <- is.finite(lo) & is.finite(hi)
  any(ok) && max(hi[ok] - lo[ok]) > 0.015 * diff(scale$lim)
}


#' Values held fixed, as text (Inf durations of absorbing moves left out)
#' @param held Named list of reference clock values.
#' @return One string, or `NULL` when no finite values are present.
#' @keywords internal
#' @noRd
held_text <- function(held) {
  held <- held[vapply(held, function(v) is.finite(v), logical(1L))]
  if (length(held) == 0L) return(NULL)
  paste0(vapply(names(held), function(n) paste0(n, " held at ",
                                                format(round(held[[n]], 1), nsmall = 1)),
                character(1L)), collapse = ", ")
}


#' The figure title: the transition exactly, with its condition
#'
#' @param e A transition-rate comparison.
#'
#' @return Plotmath expression, including the origin condition if present.
#'
#' @keywords internal
#' @noRd
plot_title <- function(e) {
  if (is.na(e$where)) return(e$title)
  cond <- paste0("  where ", e$where)
  bquote(.(e$title) ~ plain(.(cond)))
}


#' The values of the reference line, for the axis limits
#'
#' @param e A transition-rate comparison.
#' @param b Bands used for the axis calculation.
#'
#' @return Finite expected/exposure values, or an empty vector when the
#'   supplied comparison curve is used instead.
#'
#' @keywords internal
#' @noRd
reference_values <- function(e, b = e$bands) {
  if (uses_reference_curve(e)) return(numeric(0))
  r <- b$expected / b$exposure
  r[is.finite(r)]
}


#' Is the reference the supplied rate itself (TRUE) or its band average?
#'
#' @param e A transition-rate comparison.
#' @return Logical scalar.
#'
#' @keywords internal
#' @noRd
uses_reference_curve <- function(e) {
  any(is.finite(e$truth$rate)) || e$forced || all(is.na(e$bands$expected))
}


#' Obtain a supplied-rate reference for each band
#'
#' Uses expected events / exposure where finite. Otherwise uses the
#' maximum non-missing supplied-grid value inside the band, if available.
#'
#' @param e A transition-rate comparison.
#' @param bg Bands to match by group and interval.
#'
#' @return Numeric reference per band, or `NA` where unavailable.
#'
#' @keywords internal
#' @noRd
band_reference <- function(e, bg) {
  r <- bg$expected / bg$exposure
  miss <- !is.finite(r)
  if (any(miss)) {
    for (i in which(miss)) {
      t <- e$truth[e$truth$group == bg$group[i] &
                     e$truth$x >= bg$lower[i] & e$truth$x <= bg$upper[i], ]
      r[i] <- if (nrow(t) == 0L || all(is.na(t$rate))) NA_real_ else max(t$rate, na.rm = TRUE)
    }
  }
  r
}


#' The supplied rate as a dashed reference in the right panel
#' @param e A transition-rate comparison.
#' @param g Group index.
#' @param scale Y-axis settings.
#' @param col Line colour.
#' @return Called for its plotting side effects.
#' @keywords internal
#' @noRd
reference_line <- function(e, g, scale, col) {
  if (uses_reference_curve(e)) {
    tr <- e$truth[e$truth$group == g, , drop = FALSE]
    draw_curve(tr$x, tr$rate, scale, col = with_alpha(col, 0.85), lwd = 1.5,
               lty = "22")
    return(invisible(NULL))
  }
  b <- e$bands[e$bands$group == g & e$bands$exposure > 0, , drop = FALSE]
  if (nrow(b) == 0L) return(invisible(NULL))
  r <- b$expected / b$exposure
  y <- to_axis(r, scale)
  graphics::segments(b$lower, y, b$upper, y, col = with_alpha(col, 0.85),
                     lwd = 1.5, lty = "22")
}


#' Shade where the supplied rate is infinite
#' @param truth Supplied-grid data with `group`, `x` and `rate`.
#' @param scale Y-axis settings.
#' @return Called for its plotting side effects. Isolated infinite
#'   points receive a minimum visible width, rather than an exact duration.
#' @keywords internal
#' @noRd
shade_forced <- function(truth, scale) {
  inf <- is.infinite(truth$rate)
  if (!any(inf)) return(invisible(NULL))
  span <- diff(range(truth$x))
  min_w <- span * 0.006
  for (g in unique(truth$group[inf])) {
    t <- truth[truth$group == g, , drop = FALSE]
    r <- rle(is.infinite(t$rate))
    ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1L
    for (k in which(r$values)) {
      x0 <- t$x[starts[k]]; x1 <- t$x[ends[k]]
      if (x1 - x0 < min_w) { m <- (x0 + x1) / 2; x0 <- m - min_w / 2; x1 <- m + min_w / 2 }
      graphics::rect(x0, scale$lim[1L], x1, scale$lim[2L],
                     col = with_alpha(plot_style()$forced, 0.22), border = NA)
      graphics::text((x0 + x1) / 2, scale$lim[2L] - diff(scale$lim) * 0.035,
                     "rate = Inf", pos = 4, offset = 0.3, cex = 0.72,
                     col = "#9A6700")
    }
  }
}


#' Panel title, left; an optional note, right
#' @param title Panel title.
#' @param note Optional note aligned to the right.
#' @param line Margin line of the header.
#' @return Called for its plotting side effects.
#' @keywords internal
#' @noRd
panel_header <- function(title, note = NULL, line = 0.7) {
  graphics::mtext(title, side = 3, line = line, adj = 0, font = 2, cex = 0.9,
                  col = "#262626")
  if (!is.null(note)) {
    graphics::mtext(note, side = 3, line = line, adj = 1, cex = 0.72,
                    col = plot_style()$muted, font = 3)
  }
}


#' A message in the middle of a panel
#' @param text Message to draw at the panel centre.
#' @return Called for its plotting side effects.
#' @keywords internal
#' @noRd
centre_note <- function(text) {
  u <- graphics::par("usr")
  graphics::text(mean(u[1:2]), mean(u[3:4]), text, col = plot_style()$muted,
                 cex = 0.95, font = 3)
}


#' The rate's call signature, e.g. "mortRates(age, calTime)"
#' @param e A transition-rate comparison.
#' @return Function name and argument names as one string.
#' @keywords internal
#' @noRd
rate_signature <- function(e) {
  name <- if (is.na(e$rate_name)) "rate" else e$rate_name
  paste0(name, "(", paste(e$arguments, collapse = ", "), ")")
}


#' The line under the title
#' @param e A transition-rate comparison.
#' @return One string describing the rate, events, exposure and ratio.
#' @keywords internal
#' @noRd
subtitle <- function(e) {
  big <- function(n) format(n, big.mark = ",", scientific = FALSE, trim = TRUE)
  parts <- c(
    rate_signature(e),
    if (e$n_edges > 1L) paste0(e$n_edges, " origin states"),
    paste0(big(e$events), " event", if (e$events == 1L) "" else "s", " in ",
           big(round(e$exposure)), " person-years"),
    if (!is.na(e$ratio) && e$expected >= 1) sprintf("observed/expected %.3f (95%% CI %.3f to %.3f)",
                                 e$ratio, e$ratio_ci[1L], e$ratio_ci[2L]),
    if (length(e$varying) == 0L) held_text(e$held)
  )
  paste(parts, collapse = "   |   ")
}


#' The entries of the legend
#' @param e A transition-rate comparison.
#' @param grouped Whether several groups have separate colours.
#' @param scale Y-axis settings.
#' @return List of labels and drawing attributes for the visible layers.
#' @keywords internal
#' @noRd
legend_items <- function(e, grouped, scale) {
  st <- plot_style()
  ink <- "#4D4D4D"
  pick <- function(own) if (grouped) ink else own
  ref_text <- if (!uses_reference_curve(e)) "supplied rate, averaged over each band (dashed)" else
    if (length(e$varying) > 0L) "supplied rate, averaged over those at risk (dashed)" else
      "supplied rate (dashed)"
  items <- list(
    list(label = "supplied rate", col = pick(st$truth), lty = 1, lwd = 2.2,
         pch = NA, bg = NA, cex = NA),
    if (shows_range(e, scale))
      list(label = paste0("10th to 90th percentile over ",
                          paste(vapply(e$varying, scale_name, character(1L)), collapse = " and ")),
           col = with_alpha(st$truth, 0.18), lty = NA, lwd = NA, pch = 15, bg = NA, cex = 2.2),
    if (any(e$bands$exposure > 0))
      list(label = "band estimate, size by person-years", col = with_alpha(pick(st$points), 0.75),
           lty = NA, lwd = NA, pch = 21, bg = with_alpha(pick(st$points), 0.45), cex = 1.4),
    if (any(e$bands$events > 0L))
      list(label = "exact 95% interval", col = with_alpha(pick(st$points), 0.5),
           lty = 1, lwd = 1, pch = NA, bg = NA, cex = NA),
    if (nrow(e$smooth) > 0L)
      list(label = "Poisson spline fit with 95% band", col = pick(st$smooth),
           lty = 1, lwd = 2.4, pch = NA, bg = NA, cex = NA),
    list(label = ref_text, col = pick(st$truth), lty = 2, lwd = 1.5,
         pch = NA, bg = NA, cex = NA),
    if (e$forced)
      list(label = "rate = Inf (forced)", col = with_alpha(st$forced, 0.45),
           lty = NA, lwd = NA, pch = 15, bg = NA, cex = 2.2),
    if (scale$log && any(e$bands$events == 0L & e$bands$exposure > 0))
      list(label = "band with no events", col = with_alpha(st$muted, 0.6),
           lty = NA, lwd = NA, pch = 25, bg = NA, cex = 0.8)
  )
  Filter(Negate(is.null), items)
}


#' How the legend is laid out on the current device
#'
#' Uniform columns, as many as fit across the device (at most 4), measured
#' in inches so that the plan can be made before the layout.
#' @param items Legend entries from `legend_items()`.
#' @param cex Legend text size.
#' @return List with entries, column and row counts, and text size.
#' @keywords internal
#' @noRd
legend_plan <- function(items, cex = 0.8) {
  labels <- vapply(items, function(it) it$label, character(1L))
  tw <- max(graphics::strwidth(labels, units = "inches", cex = cex))
  key <- graphics::strwidth("MMMMMMMMM", units = "inches", cex = cex)
  avail <- 0.94 * graphics::par("din")[1L]
  ncol <- 1L
  for (k in rev(seq_len(min(4L, length(labels))))) {
    if (k * (tw + key) <= avail) { ncol <- k; break }
  }
  list(items = items, ncol = ncol, rows = ceiling(length(labels) / ncol),
       cex = cex)
}


#' The legend strip
#' @param e A transition-rate comparison.
#' @param grouped Whether to draw a second legend for group colours.
#' @param key Layout plan from `legend_plan()`.
#' @return Called for its plotting side effects.
#' @keywords internal
#' @noRd
draw_legend <- function(e, grouped, key) {
  items <- key$items
  get <- function(f) vapply(items, function(it) as.character(it[[f]]), character(1L))
  num <- function(f) as.numeric(get(f))
  labels <- get("label")
  # legend() does not leave room between a column's text and the next
  # column's symbol; the padding does.
  tw <- max(graphics::strwidth(labels, cex = key$cex)) +
    graphics::strwidth("MMMM", cex = key$cex)
  graphics::legend(
    x = "top", legend = labels, col = get("col"), lty = num("lty"),
    lwd = num("lwd"), pch = num("pch"), pt.bg = get("bg"), pt.cex = num("cex"),
    ncol = key$ncol, bty = "n", cex = key$cex, text.col = "#333333",
    seg.len = 2.2, x.intersp = 0.7, text.width = tw, inset = c(0, 0.02),
    xpd = NA)
  if (grouped) {
    g <- e$groups
    graphics::legend(
      x = "bottom", legend = g$label,
      col = plot_style()$groups[seq_len(nrow(g))], lwd = 4, horiz = TRUE,
      bty = "n", cex = key$cex, seg.len = 1.6, text.col = "#333333",
      inset = c(0, 0.02), xpd = NA,
      title = paste0("Colour: periods of ", scale_name(e$group_by),
                     " with equal person-years"),
      title.col = plot_style()$muted)
  }
}


#' A colour with transparency
#' @param col R colour specification(s).
#' @param a Opacity between zero and one.
#' @return Colour strings including an alpha channel.
#' @keywords internal
#' @noRd
with_alpha <- function(col, a) {
  rgb <- grDevices::col2rgb(col)
  grDevices::rgb(rgb[1L, ], rgb[2L, ], rgb[3L, ], alpha = round(a * 255),
                 maxColorValue = 255)
}
