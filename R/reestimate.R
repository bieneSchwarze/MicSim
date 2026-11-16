# =====================================================================
# reestimate.R
# ---------------------------------------------------------------------
# WHAT THIS FILE DOES
#   Estimates every transition rate of a model from the population it
#   simulated and compares each estimate with the rate the user supplied.
#   This checks that the simulation reproduces the supplied rates.
#
# PLACE IN THE PIPELINE
#   simulate() --> micsim_result --> reestimate() --> result$reestimation
#                                                       |
#                                                       v
#                                                     plot()
#
# THE METHOD  (classical occurrence/exposure on a Lexis diagram)
#   1. The life courses in `pop` are cut into SPELLS. A spell is an
#      uninterrupted stretch of time in one state. It starts at entry (start 
#      of the simulation, arrival or birth) or at the last transition, and 
#      ends at the next transition or the end of observation (the end of the 
#      horizon, or reaching max_age, as in the engine).
#   2. For each transition as DECLARED (one transition() call, however
#      many origin states it expands into), the risk set is every spell in
#      one of its origin states. A spell ending in the declared
#      destination is an event; any other ending is censoring.
#   3. The rate is shown along its MAIN TIME SCALE: age if it uses age, else 
#      duration, else calendar time (the `axis` argument can choose another).
#      This scale is cut into BANDS. They are 1 year wide unless the observed
#      range is short, and the scale is also the x axis of the plot. Wherever
#      the supplied rate is evaluated, it is given every time scale it takes. 
#      Exposure in each band is summed exactly and each event is counted in
#      the band the person was at risk in just before it. Occurrence / exposure
#      is the maximum-likelihood estimate of a rate that is constant within the
#      band, and it comes with an exact Poisson confidence interval.
#   4. A smooth curve is fitted to the band counts by Poisson regression with
#      log(exposure) offset and a natural cubic spline in the time scale.
#      Knots sit at quantiles of the event times, and their number is chosen
#      by AIC. Fitting to band totals gives the same estimates as fitting to
#      individual records split at the band boundaries, as long as the spline is
#      evaluated at one point per band (midpoint). The curve is fitted only
#      where events occurred, never extrapolated, and left out when there are
#      too few events.
#   5. If the rate also depends on another time scale (e.g. age and calendar
#      time), the estimate is compared with the supplied rate averaged, at each
#      value of the main time scale, over the people at risk there 
#      (conditional_truth()), as that average is what a band estimates.
#      If the rate varies by more than 25% along the second scale, 
#      the risk set is also split into up to three periods of it, 
#      and every step is done per period.
#   6. The expected number of events in each band is the supplied rate
#      integrated over the exposure in that band, at each person's age, 
#      calendar time and duration. Observed / expected near 1 means the 
#      simulation reproduces the supplied rate.
#
# =====================================================================


#' Reestimate a model's transition rates from its simulated population
#'
#' @description
#' Estimates every transition rate of the model from the life courses it
#' produced and compares each estimate with the rate that was supplied.
#' [simulate()][simulate.micsim_model] does this automatically and stores
#' the result as `result$reestimation`. You can call `reestimate()` to redo it
#' with other settings. `plot()` supplies comparison plots.
#'
#' @details
#' For each transition declared with [transition()]:
#'
#' * the **risk set** is the time spent in any of a transition's origin 
#'   states, and an **event** is a move into its destination;
#' * the rate is shown along its **main time scale**: age, or
#'   duration in the state if the rate does not use age, or calendar time if
#'   it uses neither (`axis` chooses another);
#' * the main time scale is cut into **bands** (1 year unless the observed 
#'   range is short), and in each band the **occurrence/exposure rate**
#'   (events divided by person-years) is computed, with an exact 95%
#'   Poisson confidence interval. Although the bands lie along one scale,
#'   the supplied rate is always evaluated at every time scale it takes;
#' * a **smooth curve** is fitted to the bands by Poisson regression with a
#'   natural cubic spline, knots at quantiles of the event times and their
#'   number chosen by AIC, over the range in which events occurred, with a
#'   pointwise 95% confidence interval;
#' * the **expected** number of events in each band is each person's
#'   exposure there times the supplied rate at their own age, calendar time
#'   and duration, summed over everyone at risk. `observed / expected` near 1
#'   means the simulated population reproduces the supplied rate.
#'
#' **Rates with more than one time scale.** Suppose a mortality rate
#' depends on age and calendar time, and age is its main time scale. The
#' band estimate at age 60 then reflects everyone who was 60 during the
#' run, in whatever calendar years that happened. So it is compared with
#' the supplied rate averaged over those people, not with the rate in one
#' particular year. The plot also shows the 10th to 90th percentile of the
#' supplied rate among them, so you can see how much it varies with
#' calendar time.
#'
#' If that variation is large (more than 25%), the whole comparison is
#' also done separately for up to three periods of calendar time, each
#' with the same number of person-years, so the plot shows how the rate
#' changes over time.
#'
#' Rates that are `Inf` somewhere (forced transitions, such as school
#' enrolment at a set age) have no finite rate to estimate. They are reported
#' with their bands, but without a smooth curve or an expected count.
#'
#' The band estimates treat the rate as constant within each band. They
#' therefore differ from the supplied rate by Monte Carlo error, which the
#' confidence intervals show, and by however much the rate changes within
#' a band.
#'
#' @param x A `micsim_result`.
#' @param width Width of the bands, in years. `NULL` (the default) uses 1
#'   year, or a finer width when the observed range is shorter than ten
#'   years.
#' @param groups Maximum number of periods of a second time scale, from 1 to 3.
#'   Fewer are used when there are few events (each period needs at least 150
#'   events) or when the rate changes by less than 25% along the second scale.
#' @param axis The main time scale for rates that depend on more
#'   than one: `"age"`, `"calTime"` or `"duration"`. It sets both the scale the
#'   bands are laid along and the x axis of the plot. A single value applies
#'   to every transition whose rate uses that scale. A named vector can set it
#'   per transition, with names as in the result. By default the first of age,
#'   duration, and calendar time that the rate uses.
#' @param ... Not used.
#'
#' @return An object of class `micsim_reestimation`. This is a named list with
#'   one element per declared transition. Each element holds
#'   \describe{
#'     \item{`label`, `title`}{the transition as text and as a plot title.}
#'     \item{`axis`, `group_by`, `groups`}{the main time scale, the second
#'       time scale the bands are grouped by (or `NA`), and a data frame
#'       with one row per period of the second scale: its limits (`lo`,
#'       `hi`), its person-year-weighted average (`at`)
#'       and the label used in the plot. 
#'       Without a second scale there is a single row.}
#'     \item{`bands`}{one row per band and group: `group`, `lower`,
#'       `upper`, `mid`, `events`, `exposure`, `rate`, `ci_lower`,
#'       `ci_upper`, `expected`.}
#'     \item{`smooth`}{the fitted curve: `group`, `x`, `rate`, `lower`,
#'       `upper`.}
#'     \item{`truth`}{the supplied rate on a fine grid: `group`, `x`,
#'       `rate`, and, for a rate that also depends on another time scale,
#'       `lower` and `upper`, the 10th and 90th percentile over the people at
#'       risk (`rate` is then their mean).}
#'     \item{`events`, `exposure`, `expected`, `ratio`, `ratio_ci`}{totals,
#'       observed / expected, and its exact 95% Poisson interval. If the
#'       simulation reproduces the supplied rate, the interval covers 1 in
#'       at least 95% of runs.}
#'     \item{`forced`, `notes`}{whether the rate is `Inf` somewhere, and
#'       anything the plot should say.}
#'   }
#'   `as.data.frame()` stacks the bands of all transitions.
#'
#' @seealso [plot.micsim_reestimation()], [simulate.micsim_model()]
#'
#' @examples
#' mort <- function(age) 0.0005 * exp(0.07 * age)
#' model <- micsim_model(
#'   state_space(sex = c("m", "f")),
#'   transitions(transition(NA, "dead", mort)),
#'   data.frame(ID = 1:300, birthDate = "19500101", initState = "f"),
#'   horizon = c(20000101, 20291231)
#' )
#' result <- simulate(model, seed = 1)
#' result$reestimation
#' head(as.data.frame(result$reestimation))
#'
#' @export
reestimate <- function(x, ...) {
  UseMethod("reestimate")
}


#' @rdname reestimate
#' @export
reestimate.micsim_result <- function(x, width = NULL, groups = 3L,
                                     axis = NULL, ...) {
  if (!is.null(width) &&
        !(is.numeric(width) && length(width) == 1L && is.finite(width) && width > 0)) {
    micsim_abort("micsim_reestimation_error",
                 "`width` must be a single positive number of years, or NULL.")
  }
  groups <- suppressWarnings(as.integer(groups))
  if (length(groups) != 1L || is.na(groups) || groups < 1L || groups > 3L) {
    micsim_abort("micsim_reestimation_error",
                 "`groups` must be 1, 2 or 3.")
  }
  if (!is.null(axis)) {
    bad <- setdiff(axis, c("age", "calTime", "duration"))
    if (!is.character(axis) || length(bad) > 0L ||
          (length(axis) > 1L && is.null(names(axis)))) {
      micsim_abort("micsim_reestimation_error", paste0(
        "`axis` must be \"age\", \"calTime\" or \"duration\", or a vector of\n",
        "  them named by transition label."))
    }
  }
  reestimate_rates(x$pop, x$model, width = width, max_groups = groups,
                   axis = axis)
}


#' The engine of reestimate()
#'
#' @param pop The simulated population.
#' @param model The `micsim_model` that produced it.
#' @param width Band width, or `NULL`.
#' @param max_groups Maximum number of groups of a second time scale.
#' @param axis As in `reestimate()`, or `NULL`.
#' @return A `micsim_reestimation`.
#' @keywords internal
#' @noRd
reestimate_rates <- function(pop, model, width = NULL, max_groups = 3L,
                             axis = NULL) {
  spells <- build_spells(pop, model)
  registry <- model$registry
  decls <- unclass(model$transitions)
  decl_of <- vapply(registry$edges, function(e) {
    if (is.null(e$decl)) NA_integer_ else as.integer(e$decl)
  }, integer(1L))

  out <- lapply(seq_along(decls), function(d) {
    edges <- registry$edges[which(decl_of == d)]
    wanted <- if (is.null(axis)) NULL else if (is.null(names(axis))) axis else {
      unname(axis[describe_declaration(decls[[d]], edges)$label])
    }
    if (length(wanted) == 0L || all(is.na(wanted))) wanted <- NULL
    reestimate_one(decls[[d]], edges, spells, model, width, max_groups,
                   wanted)
  })
  names(out) <- vapply(out, `[[`, character(1L), "label")
  structure(out, class = "micsim_reestimation",
            method = "occurrence/exposure by band; Poisson natural-spline smooth")
}


# ---------------------------------------------------------------------
# 1. Spells
# ---------------------------------------------------------------------

#' Cut the simulated life courses into spells
#'
#' @description
#' One entry per spell, which is an uninterrupted stretch of time an individual
#' spent in one transient state. Each entry holds when the spell began and
#' ended (day counts), the state, how it ended (the destination, or `NA`
#' if it was censored), the individual's birth date, and, for each state
#' variable, the day its duration clock was last reset.
#'
#' @details
#' Entry is the start of the simulation for the starting population, the
#' arrival date for immigrants and the birth date for children born during
#' the run. Observation ends at whichever comes first: the end of the horizon,
#' or the day the individual reaches `max_age`, same as in the engine. 
#' Event times are taken from the dates recorded in `pop`, which are whole
#' days, so a spell can be off by up to a day. Over the many person-years in 
#' a run, this makes no practical difference.
#'
#' Duration clocks start at birth and are reset by a transition that
#' changes the variable, as in the person store.
#'
#' @param pop,model As in `reestimate_rates()`.
#' @return A list of equal-length vectors, plus `clock`, a matrix with one
#'   column per state variable.
#' @keywords internal
#' @noRd
build_spells <- function(pop, model) {
  codec <- model$codec
  start_days <- getInDays(model$horizon[1L])
  end_days   <- getInDays(model$horizon[2L])

  first <- !duplicated(pop$ID)
  people_id    <- pop$ID[first]
  people_birth <- getInDays(pop$birthDate[first])
  people_state <- as.character(pop$initState[first])
  entry <- rep(start_days, length(people_id))
  if (!is.null(model$immigrants)) {
    hit <- match(people_id, model$immigrants$ID)
    entry[!is.na(hit)] <- getInDays(model$immigrants$immigrDate[hit[!is.na(hit)]])
  }
  supplied <- c(model$population$ID,
                if (!is.null(model$immigrants)) model$immigrants$ID)
  newborn <- !(people_id %in% supplied)
  entry[newborn] <- people_birth[newborn]
  end_obs <- pmin(end_days, people_birth + (model$max_age - MAX_AGE_MARGIN) * DAYS_PER_YEAR)

  moved <- !is.na(pop$From)
  ev_person <- match(pop$ID[moved], people_id)
  ev_time   <- getInDays(pop$transitionTime[moved])
  ev_from   <- as.character(pop$From[moved])
  ev_to     <- as.character(pop$To[moved])

  # Marks: one entry mark per person, then one per event, in order. pop is
  # ordered by ID and, within an individual, chronologically.
  n_p <- length(people_id)
  person <- c(seq_len(n_p), ev_person)
  time   <- c(entry, ev_time)
  state  <- c(people_state, ev_to)
  is_ev  <- c(rep(FALSE, n_p), rep(TRUE, length(ev_person)))
  from   <- c(rep(NA_character_, n_p), ev_from)
  ord <- order(person, is_ev, seq_along(person))   # entry mark first
  person <- person[ord]; time <- time[ord]; state <- state[ord]
  is_ev <- is_ev[ord]; from <- from[ord]

  n <- length(person)
  nxt_same <- c(person[-1L] == person[-n], FALSE)
  stop  <- ifelse(nxt_same, c(time[-1L], NA), end_obs[person])
  exit  <- ifelse(nxt_same, c(state[-1L], NA), NA_character_)

  # Duration clocks: for each variable, the latest mark at or before this
  # one that changed it, else the birth day.
  vars <- names(codec$grid)
  clock <- matrix(people_birth[person], nrow = n, ncol = length(vars),
                  dimnames = list(NULL, vars))
  if (any(is_ev) && length(vars) > 0L) {
    from_parts <- split_labels(from, length(vars))
    to_parts   <- split_labels(state, length(vars))
    shift <- min(time) - 1
    big <- (max(time) - shift + 1) * 2
    for (v in seq_along(vars)) {
      changed <- is_ev & !is.na(from_parts[, v]) &
        (is.na(to_parts[, v]) | from_parts[, v] != to_parts[, v])
      key <- person * big + ifelse(changed, time - shift, 0)
      last <- cummax(key) - person * big
      clock[, v] <- ifelse(last > 0, last + shift, people_birth[person])
    }
  }

  keep <- !(state %in% codec$absorbing) & !is.na(stop) & stop >= time
  list(
    person = person[keep],
    start  = time[keep],
    stop   = pmin(stop[keep], end_obs[person[keep]]),
    state  = state[keep],
    exit   = exit[keep],
    birth  = people_birth[person[keep]],
    clock  = clock[keep, , drop = FALSE]
  )
}


#' Split state labels into their components
#' @param labels Character vector. A label without exactly `k` parts, for
#'  example `NA` or an absorbing state like `"dead"`, gives a row of `NA`.
#' @param k Number of state variables.
#' @return A character matrix with `k` columns.
#' @keywords internal
#' @noRd
split_labels <- function(labels, k) {
  parts <- strsplit(ifelse(is.na(labels), "", labels), STATE_SEP, fixed = TRUE)
  t(vapply(parts, function(p) {
    if (length(p) != k) rep(NA_character_, k) else p
  }, character(k)))
}


# ---------------------------------------------------------------------
# 2. One transition
# ---------------------------------------------------------------------

#' Reestimate the rate of one declared transition
#' @param decl The transition as declared.
#' @param edges Its resolved edges.
#' @param spells From `build_spells()`.
#' @param model The model.
#' @param width,max_groups As in `reestimate_rates()`.
#' @param wanted_axis The axis requested for this transition, or `NULL`.
#' @return One element of a `micsim_reestimation`.
#' @keywords internal
#' @noRd
reestimate_one <- function(decl, edges, spells, model, width, max_groups,
                           wanted_axis = NULL) {
  describe <- describe_declaration(decl, edges)
  rate <- decl$rate
  argnames <- names(formals(rate))
  if (is.null(argnames)) argnames <- character(0)
  to_absorbing <- length(edges) > 0L &&
    all(vapply(edges, function(e) model$codec$is_absorbing(e$to_idx), logical(1L)))
  # A transition into an absorbing state sees duration = Inf, so duration
  # is not a time scale it can be shown on.
  scales <- intersect(c("age", "duration", "calTime"), argnames)
  plottable <- if (to_absorbing) setdiff(scales, "duration") else scales
  axis <- if (length(plottable) == 0L) "age" else
    if (!is.null(wanted_axis) && wanted_axis %in% plottable) wanted_axis else
      plottable[1L]
  second <- setdiff(plottable, axis)
  group_by <- if (length(second) == 0L) NA_character_ else
    if ("calTime" %in% second) "calTime" else second[1L]

  # --- the risk set ---------------------------------------------------
  from_labels <- vapply(edges, `[[`, character(1L), "from_label")
  to_labels   <- vapply(edges, `[[`, character(1L), "to_label")
  in_set <- which(spells$state %in% from_labels)
  dest <- to_labels[match(spells$state[in_set], from_labels)]
  event <- !is.na(spells$exit[in_set]) & spells$exit[in_set] == dest
  affected <- vapply(edges, function(e) {
    if (is.na(e$affected_var)) NA_integer_ else as.integer(e$affected_var)
  }, integer(1L))[match(spells$state[in_set], from_labels)]

  start <- spells$start[in_set]; stop <- spells$stop[in_set]
  len <- (stop - start) / DAYS_PER_YEAR
  clocks <- rep(NA_real_, length(in_set))
  has_clock <- !is.na(affected) & !to_absorbing
  if (any(has_clock)) {
    clocks[has_clock] <- spells$clock[cbind(in_set[has_clock], affected[has_clock])]
  }
  at_start <- list(
    age      = (start - spells$birth[in_set]) / DAYS_PER_YEAR,
    calTime  = caltime_in_years(start),
    duration = if (to_absorbing) rep(Inf, length(start)) else (start - clocks) / DAYS_PER_YEAR
  )

  x0 <- at_start[[axis]]
  domain <- if (length(x0) > 0L) range(c(x0, x0 + len)) else default_domain(axis, model)
  if (!all(is.finite(domain)) || diff(domain) <= 0) {
    domain <- default_domain(axis, model)
  }
  breaks <- band_breaks(domain, width)

  # --- groups along a second time scale --------------------------------
  probe_at <- function(sec_value, held) {
    grid <- truth_grid(breaks)
    vals <- evaluate_rate(rate, argnames, axis, grid, group_by, sec_value,
                          held, to_absorbing)
    list(x = grid, rate = vals)
  }
  held <- held_values(argnames, axis, group_by, at_start, len, to_absorbing)
  # Time scales other than the axis that the rate actually varies along
  # (duration is fixed at Inf for a move into an absorbing state).
  varying <- setdiff(intersect(argnames, c("age", "calTime", "duration")),
                     c(axis, if (to_absorbing) "duration"))
  n_events <- sum(event)
  groups <- make_groups(group_by, at_start, len, n_events, max_groups,
                        probe_at, held)

  # --- bands, smooth and truth, per group ------------------------------
  bands <- list(); smooth <- list(); truth <- list(); info <- list()
  for (g in seq_len(nrow(groups))) {
    piece <- clip_to_group(groups[g, ], group_by, at_start, len, event)
    event_x <- x0[piece$event_rows] + len[piece$event_rows]
    b <- band_table(x0[piece$rows] + piece$from, x0[piece$rows] + piece$to,
                    event_x, breaks)
    tr <- if (length(varying) == 0L || length(piece$rows) == 0L) {
      probe_at(groups$at[g], held)
    } else {
      conditional_truth(rate, argnames, axis, truth_grid(breaks), piece,
                        at_start, to_absorbing)
    }
    if (all(is.na(tr$rate)) && length(varying) > 0L) {
      tr <- probe_at(groups$at[g], held)
    }
    b$expected <- expected_by_band(rate, argnames, piece$rows,
                                   piece$from, piece$to, at_start, breaks,
                                   to_absorbing, axis)
    b$group <- g
    bands[[g]] <- b
    truth[[g]] <- data.frame(
      group = g, x = tr$x, rate = tr$rate,
      lower = if (is.null(tr$lower)) NA_real_ else tr$lower,
      upper = if (is.null(tr$upper)) NA_real_ else tr$upper)
    fit <- smooth_bands(b, forced = any(is.infinite(tr$rate)))
    info[[g]] <- fit$info
    if (!is.null(fit$curve)) smooth[[g]] <- cbind(group = g, fit$curve)
  }
  bands  <- do.call(rbind, bands)
  bands  <- bands[, c("group", "lower", "upper", "mid", "events", "exposure",
                      "rate", "ci_lower", "ci_upper", "expected")]
  truth  <- do.call(rbind, truth)
  smooth <- if (length(smooth) > 0L) do.call(rbind, smooth) else
    data.frame(group = integer(0), x = numeric(0), rate = numeric(0),
               lower = numeric(0), upper = numeric(0))
  forced <- any(is.infinite(truth$rate))
  expected <- if (forced || anyNA(bands$expected)) NA_real_ else sum(bands$expected)
  exposure <- sum(bands$exposure)

  edge_ids <- vapply(edges, function(e) as.integer(e$id), integer(1L))
  timed_by <- names(Filter(function(s) {
    any(unlist(s$claims[c("timing")], use.names = FALSE) %in% edge_ids)
  }, model$subsystems))

  notes <- character(0)
  if (length(in_set) == 0L || exposure == 0) {
    notes <- c(notes, "Nobody was at risk of this transition.")
  } else if (n_events == 0L) {
    notes <- c(notes, "No events were observed.")
  }
  if (forced) {
    notes <- c(notes, paste0(
      "Forced transition: the rate is infinite where shaded, so there is no rate to smooth",
      if (length(timed_by) > 0L) paste0("; the timing is set by subsystem ",
                                        toString(paste0("`", timed_by, "`"))) else "",
      ". The points show the rate observed in the simulation."))
  } else if (length(timed_by) > 0L) {
    notes <- c(notes, paste0(
      "The timing is adjusted by subsystem ", toString(paste0("`", timed_by, "`")),
      ", so the estimated rate need not follow the supplied one."))
  }
  if (!forced) {
    why <- unique(unlist(lapply(info, `[[`, "reason")))
    why <- why[!is.na(why)]
    if (length(why) > 0L && n_events > 0L) notes <- c(notes, why)
  }
  if (length(varying) == 0L && anyNA(truth$rate[is.finite(truth$x)]) && !forced) {
    notes <- c(notes, "The supplied rate could not be evaluated everywhere.")
  }

  structure(
    list(
      label     = describe$label,
      title     = describe$title,
      where     = describe$where,
      rate_name = describe$rate_name,
      arguments = argnames,
      axis      = axis,
      group_by  = group_by,
      groups    = groups,
      held      = held,
      varying   = varying,
      timed_by  = timed_by,
      n_edges   = length(edges),
      bands     = bands,
      smooth    = smooth,
      truth     = truth,
      fits      = info,
      events    = as.integer(n_events),
      exposure  = exposure,
      expected  = expected,
      ratio     = if (is.na(expected) || expected <= 0) NA_real_ else n_events / expected,
      ratio_ci  = ratio_interval(n_events, expected),
      forced    = forced,
      notes     = notes
    ),
    class = "micsim_rate_check"
  )
}


#' Exact 95% interval for observed / expected
#'
#' The observed count is Poisson with mean `expected` times the ratio, so
#' the interval is the exact (Garwood) interval for a Poisson mean, 
#' divided by `expected`.
#' @param observed Events.
#' @param expected Expected events.
#' @return `c(lower, upper)`, or `NA`s when there is nothing expected.
#' @keywords internal
#' @noRd
ratio_interval <- function(observed, expected) {
  if (is.na(expected) || expected <= 0) return(c(NA_real_, NA_real_))
  c(if (observed > 0) stats::qchisq(0.025, 2 * observed) / 2 else 0,
    stats::qchisq(0.975, 2 * observed + 2) / 2) / expected
}


#' Expected events per band under the supplied rate
#'
#' @description
#' Integrates the supplied rate along every piece of every spell (the part
#' of a spell that falls in the current group), band by band, using the
#' midpoint of each piece's segment inside the band. All time scales
#' advance together along a spell, so the rate is evaluated at the actual
#' age, calendar time and duration at that midpoint. With bands one year
#' wide this is exact for rates that step at whole years of the axis
#' scale, and very close otherwise: the error shrinks with the square of
#' the band width.
#'
#' @details
#' Segments are processed in blocks, so memory stays bounded however
#' large the population.
#'
#' @return Numeric vector, one value per band; `NA` everywhere if the
#'   rate is infinite or cannot be evaluated on some segment.
#' @keywords internal
#' @noRd
expected_by_band <- function(rate, argnames, rows, from, to, at_start,
                             breaks, to_absorbing, axis) {
  nb <- length(breaks) - 1L
  out <- numeric(nb)
  if (length(rows) == 0L) return(out)
  x0 <- at_start[[axis]][rows]
  a <- x0 + from; b <- x0 + to
  ka <- findInterval(a, breaks, all.inside = TRUE)
  kb <- findInterval(b - 1e-12, breaks, all.inside = TRUE)
  kb <- pmax(kb, ka)
  n_seg <- kb - ka + 1L
  block <- 1e6
  ends <- cumsum(n_seg)
  first <- 1L
  while (first <= length(rows)) {
    base <- if (first > 1L) ends[first - 1L] else 0
    last <- max(first, findInterval(base + block, ends))
    idx <- first:last
    piece <- rep(idx, n_seg[idx])
    k <- sequence(n_seg[idx], from = ka[idx])
    lo <- pmax(a[piece], breaks[k]); hi <- pmin(b[piece], breaks[k + 1L])
    keep <- hi > lo
    piece <- piece[keep]; k <- k[keep]; lo <- lo[keep]; hi <- hi[keep]
    mid <- (lo + hi) / 2
    if (length(mid) > 0L) {
      shift <- mid - x0[piece]
      args <- list()
      for (v in argnames) {
        args[[v]] <- if (v == "duration" && to_absorbing) rep(Inf, length(mid)) else
          at_start[[v]][rows[piece]] + shift
      }
      r <- if (length(argnames) == 0L) rep(as.numeric(rate())[1L], length(mid)) else
        tryCatch(as.numeric(do.call(rate, args)), error = function(e) NULL)
      if (is.null(r) || length(r) != length(mid) || anyNA(r) || any(is.infinite(r))) {
        return(rep(NA_real_, nb))
      }
      out <- out + sum_by_band(k, r * (hi - lo), nb)
    }
    first <- last + 1L
  }
  out
}


#' Sum values by integer band index
#' @keywords internal
#' @noRd
sum_by_band <- function(k, v, nb) {
  s <- rowsum(v, k, reorder = TRUE)
  out <- numeric(nb)
  out[as.integer(rownames(s))] <- s[, 1L]
  out
}


#' Text and plot title for a declared transition
#' @param decl The declaration.
#' @param edges Its edges.
#' @return `list(label, title, where, rate_name)`; `title` is a plotmath
#'   expression.
#' @keywords internal
#' @noRd
describe_declaration <- function(decl, edges) {
  arrow_text <- " -> "
  rate_name <- if (is.null(decl$rate_name)) NA_character_ else decl$rate_name
  where <- NA_character_
  if (!is.null(decl$move)) {
    e <- decl$move
    if (!is.null(e$where)) where <- deparse1(e$where$expr)
    label <- paste0(e$var, ": ", e$from, arrow_text, e$to)
    if (!is.na(where)) label <- paste0(label, " (where ", where, ")")
    title <- bquote(bold(.(e$var) * ":" ~ .(e$from) %->% .(e$to)))
  } else {
    from <- if (length(decl$from) == 1L && is.na(decl$from)) {
      "any transient state"
    } else {
      paste(decl$from, collapse = ", ")
    }
    label <- paste0(from, arrow_text, decl$to)
    title <- bquote(bold(.(from) %->% .(as.character(decl$to))))
  }
  list(label = label, title = title, where = where, rate_name = rate_name)
}


#' The default range of a time scale when nobody was at risk
#' @keywords internal
#' @noRd
default_domain <- function(axis, model) {
  switch(axis,
    age      = c(0, model$max_age),
    duration = c(0, model$max_age),
    calTime  = caltime_in_years(getInDays(model$horizon)))
}


#' Band boundaries over a range
#' @param domain `c(lo, hi)`.
#' @param width Band width, or `NULL` for automatic.
#' @return Increasing numeric vector of boundaries.
#' @keywords internal
#' @noRd
band_breaks <- function(domain, width = NULL) {
  # A range a day short of whole years (a horizon from 1 January to 31
  # December) counts as whole years.
  span <- diff(domain) + 0.003
  if (is.null(width)) {
    nice <- c(1, 0.5, 0.25, 0.2, 0.1, 0.05, 0.025, 0.01)
    width <- nice[c(span >= 10, span / nice[-1L] >= 10)][1L]
    if (is.na(width)) width <- span / 10
  }
  # A boundary within about a day of a whole band is that band's edge:
  # calendar time on 365.25-day years puts 1 January up to 18 hours away
  # from the whole year, which must not open a sliver band of its own.
  tol <- min(0.25, 0.003 / width)
  lo <- floor(domain[1L] / width + tol) * width
  hi <- ceiling(domain[2L] / width - tol) * width
  if (hi <= lo) hi <- lo + width
  n <- round((hi - lo) / width)
  breaks <- lo + (0:n) * width
  # ... and the outer bands stretch to hold that sliver, so no exposure is
  # lost.
  breaks[1L] <- min(breaks[1L], domain[1L])
  breaks[n + 1L] <- max(breaks[n + 1L], domain[2L])
  breaks
}


#' The grid the supplied rate is drawn on
#'
#' A regular grid plus every band boundary and a point just before it,
#' so that a rate that steps at whole years is drawn with sharp steps and
#' a rate that is infinite at a single age (forced transition) is found.
#' @keywords internal
#' @noRd
truth_grid <- function(breaks) {
  lo <- breaks[1L]; hi <- breaks[length(breaks)]
  eps <- (hi - lo) * 1e-6
  inner <- breaks[breaks > lo & breaks < hi]
  sort(unique(c(seq(lo, hi, length.out = 401L), breaks, inner - eps)))
}


#' Evaluate a supplied rate along the x axis
#' @param rate,argnames The rate and its argument names.
#' @param axis The time scale on the x axis.
#' @param x Values of that scale.
#' @param group_by,sec_value The grouping scale and the value it is held at.
#' @param held Named list of values for any other scale.
#' @param to_absorbing Whether duration is `Inf` by definition.
#' @return Numeric vector; `NA` where the rate could not be evaluated.
#' @keywords internal
#' @noRd
evaluate_rate <- function(rate, argnames, axis, x, group_by, sec_value, held,
                          to_absorbing) {
  if (length(x) == 0L) return(numeric(0))
  args <- list()
  for (a in argnames) {
    args[[a]] <- if (a == axis) x else
      if (!is.na(group_by) && a == group_by) rep(sec_value, length(x)) else
        if (a == "duration" && to_absorbing) rep(Inf, length(x)) else
          rep(held[[a]], length(x))
  }
  if (length(argnames) == 0L) {
    v <- tryCatch(as.numeric(rate()), error = function(e) NA_real_)
    return(rep(v[1L], length(x)))
  }
  v <- tryCatch(as.numeric(do.call(rate, args)), error = function(e) NULL)
  if (is.null(v) || length(v) != length(x)) {
    # not vectorised after all, or failing somewhere: point by point
    v <- vapply(seq_along(x), function(i) {
      tryCatch(as.numeric(do.call(rate, lapply(args, `[`, i)))[1L],
               error = function(e) NA_real_)
    }, numeric(1L))
  }
  v[!is.na(v) & v < 0] <- NA_real_
  v
}


#' The supplied rate averaged over the people at risk along the main time scale
#'
#' @description
#' When a rate depends on more than the time scale on the axis, the
#' band estimate at axis value x estimates the supplied rate averaged over
#' the person-years at risk at x, over the calendar times, durations or
#' ages those people have there. This function computes that from the spells,
#' on a grid of x. Every spell at risk at x contributes the
#' supplied rate at its own values of the other time scales. It is the
#' curve the estimated rate should be compared with, and it equals the
#' supplied rate when the rate depends on the axis alone.
#'
#' @details
#' The 10th and 90th percentiles over the people at risk are returned too:
#' how much the supplied rate varies along the other time scales at each x.
#'
#' To limit memory and time, at most `max_pairs` (spell, grid point) pairs
#' are evaluated; above that, every k-th spell is used. The selection is
#' systematic, so no random numbers are drawn.
#'
#' @param rate,argnames The supplied rate and its arguments.
#' @param axis The time scale on the x axis.
#' @param grid Values of the axis scale.
#' @param piece `list(rows, from, to)` from `clip_to_group()`.
#' @param at_start Age, calendar time and duration at the start of every spell
#'  in the risk set.
#' @param to_absorbing Whether duration is `Inf` by definition.
#' @param max_pairs Cap on the number of evaluations.
#' @return `list(x, rate, lower, upper)`; `NA` where nobody is at risk.
#' @keywords internal
#' @noRd
conditional_truth <- function(rate, argnames, axis, grid, piece, at_start,
                              to_absorbing, max_pairs = 3e5) {
  rows <- piece$rows
  x0 <- at_start[[axis]][rows]
  a <- x0 + piece$from; b <- x0 + piece$to
  ja <- findInterval(a - 1e-12, grid) + 1L     # first grid point >= a
  jb <- findInterval(b - 1e-12, grid)          # last grid point <  b
  n <- pmax(0L, jb - ja + 1L)
  if (sum(n) > max_pairs) {
    keep <- seq(1L, length(n), by = ceiling(sum(n) / max_pairs))
    rows <- rows[keep]; x0 <- x0[keep]; ja <- ja[keep]; n <- n[keep]
  }
  na <- rep(NA_real_, length(grid))
  if (sum(n) == 0L) return(list(x = grid, rate = na, lower = na, upper = na))
  i <- rep(seq_along(n), n)
  j <- sequence(n, from = ja)
  shift <- grid[j] - x0[i]
  args <- list()
  for (v in argnames) {
    args[[v]] <- if (v == "duration" && to_absorbing) rep(Inf, length(j)) else
      at_start[[v]][rows[i]] + shift
  }
  val <- tryCatch(as.numeric(do.call(rate, args)), error = function(e) NULL)
  if (is.null(val) || length(val) != length(j)) {
    return(list(x = grid, rate = na, lower = na, upper = na))
  }
  ok <- !is.na(val) & is.finite(val) & val >= 0
  forced_at <- unique(j[is.infinite(val)])
  j <- j[ok]; val <- val[ok]
  mean_rate <- na; lower <- na; upper <- na
  if (length(j) > 0L) {
    cnt <- tabulate(j, nbins = length(grid))
    tot <- sum_by_band(j, val, length(grid))
    has <- cnt > 0L
    mean_rate[has] <- tot[has] / cnt[has]
    o <- order(j, val)
    sv <- val[o]
    start <- c(0L, cumsum(cnt))[seq_along(grid)]
    pick <- function(p) {
      out <- na
      pos <- start[has] + floor(p * (cnt[has] - 1L)) + 1L
      out[has] <- sv[pos]
      out
    }
    lower <- pick(0.1); upper <- pick(0.9)
  }
  # A forced point stays visible.
  mean_rate[forced_at] <- Inf
  list(x = grid, rate = mean_rate, lower = lower, upper = upper)
}


#' Values for a time scale that is neither on the axis nor grouped by
#'
#' To draw the supplied rate as one curve along the axis, any other time
#' scale it depends on has to be fixed at a single value. This uses the
#' average of that scale over the risk set, taken at the middle of each
#' spell and weighted by the spell's length.
#' @keywords internal
#' @noRd
held_values <- function(argnames, axis, group_by, at_start, len, to_absorbing) {
  out <- list()
  for (a in setdiff(argnames, c(axis, group_by))) {
    if (a == "duration" && to_absorbing) { out[[a]] <- Inf; next }
    v <- at_start[[a]] + len / 2
    w <- len
    out[[a]] <- if (length(v) == 0L || sum(w) <= 0 || !any(is.finite(v))) {
      if (a == "calTime") 2000 else 0
    } else {
      ok <- is.finite(v)
      stats::weighted.mean(v[ok], w[ok] + 1e-12)
    }
  }
  out
}


#' Periods of the second time scale
#'
#' Splits the second time scale into up to `max_groups` periods with equal 
#' person-years at risk.
#' Fewer periods are used when there are few events 
#' (at least 150 events per period), and only one when the supplied
#' rate changes by less than 25% across the observed range. This is because a 
#' difference that small cannot be told apart from Monte Carlo error 
#' band by band, and extra colours would only clutter the plot. 
#' With one period, the comparison uses the supplied rate
#' averaged over the exposure in each band, which still takes the second
#' scale fully into account.
#' @return Data frame with `lo`, `hi`, `at` (the person-year-weighted average
#'   of the second scale in the period, used for the label of a single period
#'   and when a rate is evaluated directly) and `label`.
#' @keywords internal
#' @noRd
make_groups <- function(group_by, at_start, len, n_events, max_groups,
                        probe_at, held) {
  one <- data.frame(lo = -Inf, hi = Inf, at = NA_real_, label = "",
                    stringsAsFactors = FALSE)
  if (is.na(group_by)) return(one)
  s0 <- at_start[[group_by]]
  ok <- is.finite(s0)
  s0 <- s0[ok]; l <- len[ok]
  rng <- if (length(s0) > 0L) range(c(s0, s0 + l)) else c(0, 1)
  if (diff(rng) <= 0) rng <- rng + c(0, 1)
  # Cumulative exposure along the second scale, exactly, on a fine grid.
  grid <- seq(rng[1L], rng[2L], length.out = 401L)
  cum <- if (length(s0) > 0L) {
    G <- function(x, v) {
      v <- sort(v); cnt <- findInterval(x, v); cs <- c(0, cumsum(v))
      cnt * x - cs[cnt + 1L]
    }
    G(grid, s0) - G(grid, s0 + l)
  } else grid - rng[1L]
  total <- cum[length(cum)]
  mean_between <- function(lo, hi) {
    # exposure-weighted mean of the scale on [lo, hi]
    k <- which(grid >= lo & grid <= hi)
    if (length(k) < 2L || total <= 0) return((max(lo, rng[1L]) + min(hi, rng[2L])) / 2)
    w <- diff(cum[k]); x <- (grid[k][-1L] + grid[k][-length(k)]) / 2
    if (sum(w) <= 0) return(mean(x))
    sum(w * x) / sum(w)
  }
  k <- min(max_groups, max(1L, floor(n_events / 150)))
  at_all <- mean_between(rng[1L], rng[2L])
  # Does the rate vary materially along this scale?
  a <- probe_at(rng[1L], held)$rate; b <- probe_at(rng[2L], held)$rate
  m <- probe_at(at_all, held)$rate
  spread <- function(u, v) {
    ok <- is.finite(u) & is.finite(v) & u > 0 & v > 0
    if (!any(ok)) return(0)
    max(abs(log(u[ok] / v[ok])))
  }
  material <- max(spread(a, b), spread(a, m), spread(m, b)) > log(1.25) ||
    !identical(is.finite(a) & a > 0, is.finite(b) & b > 0)
  if (!isTRUE(material)) k <- 1L
  if (k == 1L) {
    one$at <- at_all
    one$label <- paste0(scale_name(group_by), " = ", format_scale(at_all))
    return(one)
  }
  # Periods with equal person-years, so that each is equally well estimated.
  cuts <- if (total > 0) {
    stats::approx(cum, grid, xout = total * seq_len(k - 1L) / k, ties = "ordered")$y
  } else seq(rng[1L], rng[2L], length.out = k + 1L)[2:k]
  edges <- c(rng[1L], cuts, rng[2L])
  data.frame(
    lo = c(-Inf, cuts),
    hi = c(cuts, Inf),
    at = vapply(seq_len(k), function(g) mean_between(edges[g], edges[g + 1L]), numeric(1L)),
    label = paste0(scale_name(group_by), " ", format_scale(edges[-(k + 1L)]),
                   " to ", format_scale(edges[-1L])),
    stringsAsFactors = FALSE
  )
}


#' Readable names and number formats of the time scales
#' @keywords internal
#' @noRd
scale_name <- function(scale) {
  switch(scale, age = "age", calTime = "calendar time",
         duration = "duration")
}

#' Format a value on a time scale for a group label: one decimal
#' @keywords internal
#' @noRd
format_scale <- function(v) format(round(v, 1), nsmall = 1, trim = TRUE)


#' The part of each spell that falls in one group
#' @return `list(rows, from, to, event_rows)`: spell rows and offsets
#'   (years from spell start) of the part in the group, and the spells
#'   whose event occurred in the group.
#' @keywords internal
#' @noRd
clip_to_group <- function(group, group_by, at_start, len, event) {
  n <- length(len)
  if (is.na(group_by) || n == 0L) {
    rows <- seq_len(n)
    return(list(rows = rows, from = rep(0, n), to = len,
                event_rows = which(event)))
  }
  s0 <- at_start[[group_by]]
  a <- pmax(0, group$lo - s0); b <- pmin(len, group$hi - s0)
  rows <- which(b > a)
  s_end <- s0 + len
  ev <- which(event & s_end > group$lo & s_end <= group$hi)
  list(rows = rows, from = a[rows], to = b[rows], event_rows = ev)
}


#' Exposure and events per band
#'
#' Exposure is summed exactly. Let G(x) be the total time the pieces spent
#' below x: G(x) = sum((x - a)+) - sum((x - b)+), where (y)+ means
#' max(y, 0). The exposure of a band [l, u) is then G(u) - G(l). After
#' sorting this takes O(n log n) time, without splitting pieces into bands.
#'
#' @param a,b Starts and ends of the pieces on the x scale.
#' @param event_x The x values of the events. An event is counted in the
#'   band the person was at risk in just before it, so an event exactly on
#'   a boundary belongs to the band below.
#' @param breaks Band boundaries.
#' @return Data frame of bands.
#' @keywords internal
#' @noRd
band_table <- function(a, b, event_x, breaks) {
  G <- function(x, v) {
    v <- sort(v)
    cnt <- findInterval(x, v)
    cs <- c(0, cumsum(v))
    cnt * x - cs[cnt + 1L]
  }
  lower <- breaks[-length(breaks)]; upper <- breaks[-1L]
  exposure <- if (length(a) == 0L) rep(0, length(lower)) else {
    g <- G(breaks, a) - G(breaks, b)
    pmax(0, diff(g))
  }
  events <- if (length(event_x) == 0L) integer(length(lower)) else {
    k <- findInterval(event_x - 1e-9, breaks, all.inside = TRUE)
    tabulate(k, nbins = length(lower))
  }
  rate <- ifelse(exposure > 0, events / exposure, NA_real_)
  # Exact (Garwood) 95% interval for a Poisson count, divided by exposure.
  ci_lower <- ifelse(exposure > 0 & events > 0,
                     stats::qchisq(0.025, 2 * events) / (2 * exposure),
                     ifelse(exposure > 0, 0, NA_real_))
  ci_upper <- ifelse(exposure > 0,
                     stats::qchisq(0.975, 2 * events + 2) / (2 * exposure),
                     NA_real_)
  data.frame(lower = lower, upper = upper, mid = (lower + upper) / 2,
             events = as.integer(events), exposure = exposure, rate = rate,
             ci_lower = ci_lower, ci_upper = ci_upper)
}


#' Fit the smooth curve to the bands of one group
#'
#' @description
#' Fits a smooth curve through the band rates of one group. The curve
#' runs from the midpoint of the first band with an event to the midpoint of
#' the last one, and is not extended beyond them.
#'
#' @details
#' The curve is a Poisson regression of each band's event count on a
#' natural cubic spline of the band midpoint, with `log(exposure)` as the
#' offset, so the rate rather than the count is described. Only bands between
#' the first and last band with an event and with positive exposure are used.
#' 
#' The flexibility of the curve depends on the number of interior knots.
#' The function tries every number from 0 (a straight line on the log scale)
#' up to the smallest of three limits: 8, one knot per 10 events,
#' and the number of bands used minus 3. It keeps the fit with the lowest AIC and
#' on a tie the one with fewer knots. With fewer than 80 events, the number of
#' events is usually the binding limit, so more events allow more knots. Beyond
#' that the limit is fixed at 8 or the number of bands used
#' minus 3. The knots are placed at quantiles of the band midpoints,
#' weighted by the number of events, so they fall where events are dense.
#' Fits that fail or do not converge are skipped.
#' 
#' No curve is fitted for a forced transition, one with fewer than 5 events,
#' or when the events span fewer than 3 bands. If there were some events
#' but too few, `info$reason` says why.
#'
#' @param b The bands of one group.
#' @param forced Whether the supplied rate is infinite somewhere.
#' @return `list(curve, info)`. `curve` is `NULL` or a data frame with `x`,
#'   `rate` and the pointwise 95% confidence limits `lower` and `upper`.
#'   `info` has `df`, `knots`, `aic` and `reason`.
#' @keywords internal
#' @noRd
smooth_bands <- function(b, forced = FALSE) {
  none <- function(reason) list(curve = NULL,
                                info = list(df = NA_integer_, knots = numeric(0),
                                            reason = reason))
  if (forced) return(none(NA_character_))
  D <- sum(b$events)
  if (D == 0L) return(none(NA_character_))
  if (D < 5L) return(none(paste0("Too few events for a smooth curve (", D, ").")))
  hit <- which(b$events > 0L)
  support <- seq(min(hit), max(hit))
  d <- b[support, , drop = FALSE]
  d <- d[d$exposure > 0, , drop = FALSE]
  if (nrow(d) < 3L) {
    return(none("Events in too few bands for a smooth curve."))
  }
  boundary <- range(d$mid)
  x_events <- rep(d$mid, d$events)
  max_knots <- max(0L, min(8L, floor(D / 10), nrow(d) - 3L))
  grid <- seq(boundary[1L], boundary[2L], length.out = 200L)

  # Fit every candidate number of interior knots and keep the one with the
  # smallest AIC: as flexible as the events support, and no more.
  best <- NULL
  for (k in 0:max_knots) {
    knots <- if (k == 0L) numeric(0) else {
      q <- unique(stats::quantile(x_events, probs = seq_len(k) / (k + 1L),
                                  names = FALSE, type = 7))
      q[q > boundary[1L] & q < boundary[2L]]
    }
    if (k > 0L && length(knots) < k) next      # ties: same as a smaller k
    basis <- splines::ns(d$mid, knots = knots, Boundary.knots = boundary)
    fit <- tryCatch(
      suppressWarnings(stats::glm.fit(
        x = cbind(1, basis), y = d$events, offset = log(d$exposure),
        family = stats::poisson())),
      error = function(e) NULL)
    if (is.null(fit) || !isTRUE(fit$converged) || anyNA(fit$coefficients) ||
          fit$rank < ncol(basis) + 1L) next
    if (is.null(best) || fit$aic < best$fit$aic - 1e-8) {
      best <- list(fit = fit, basis = basis, knots = knots)
    }
  }
  if (is.null(best)) return(none("The smooth curve could not be fitted."))

  fit <- best$fit
  X <- cbind(1, stats::predict(best$basis, grid))
  eta <- drop(X %*% fit$coefficients)
  V <- tryCatch({
    R <- fit$qr$qr[seq_len(fit$rank), seq_len(fit$rank), drop = FALSE]
    R[lower.tri(R)] <- 0
    Vp <- chol2inv(R)
    piv <- fit$qr$pivot[seq_len(fit$rank)]
    out <- matrix(0, ncol(X), ncol(X)); out[piv, piv] <- Vp; out
  }, error = function(e) NULL)
  if (is.null(V) || any(!is.finite(eta)) || max(eta) > 50) {
    return(none("The smooth curve could not be fitted."))
  }
  se <- sqrt(pmax(0, rowSums((X %*% V) * X)))
  curve <- data.frame(x = grid, rate = exp(eta),
                      lower = exp(eta - 1.96 * se), upper = exp(eta + 1.96 * se))
  list(curve = curve,
       info = list(df = ncol(best$basis), knots = best$knots,
                   aic = fit$aic, reason = NA_character_))
}


# ---------------------------------------------------------------------
# 3. Printing and conversion
# ---------------------------------------------------------------------

#' Print a reestimation
#'
#' One line per declared transition including events, person-years at risk,
#' expected events under the supplied rate, and observed / expected.
#'
#' @param x A `micsim_reestimation`.
#' @param ... Not used.
#' @return `x`, invisibly.
#' @export
print.micsim_reestimation <- function(x, ...) {
  cat("<micsim_reestimation> ", length(x), " transition(s)\n", sep = "")
  if (length(x) == 0L) return(invisible(x))
  tab <- summary_table(x)
  big <- function(v) format(v, big.mark = ",", scientific = FALSE, trim = TRUE)
  # A ratio on less than one expected event says nothing; leave it out.
  shown <- !is.na(tab$ratio) & tab$expected >= 1
  oe <- ifelse(shown, sprintf("%.3f [%.2f, %.2f]", tab$ratio, tab$ratio_lower,
                              tab$ratio_upper),
               ifelse(tab$forced, "forced", "-"))
  cols <- list(
    transition         = tab$transition,
    events             = big(tab$events),
    `person-years`     = big(round(tab$exposure)),
    expected           = ifelse(is.na(tab$expected), "-",
                                big(round(tab$expected, 1))),
    `obs/exp [95% CI]` = oe
  )
  left <- c(TRUE, FALSE, FALSE, FALSE, TRUE)
  width <- mapply(function(v, n) max(nchar(c(v, n))), cols, names(cols))
  pad <- function(v, w, l) formatC(v, width = if (l) -w else w)
  header <- paste(mapply(pad, names(cols), width, left), collapse = "  ")
  cat("  ", header, "\n", sep = "")
  body <- do.call(paste, c(Map(pad, cols, width, left), sep = "  "))
  cat(paste0("  ", body, "\n"), sep = "")
  cat("Observed / expected near 1 means the simulation reproduces the supplied\n",
      "rate; plot() draws each supplied rate beside the estimated one.\n", sep = "")
  invisible(x)
}


#' One row per transition
#' @keywords internal
#' @noRd
summary_table <- function(x) {
  data.frame(
    transition = vapply(x, `[[`, character(1L), "label"),
    axis       = vapply(x, `[[`, character(1L), "axis"),
    events     = vapply(x, `[[`, integer(1L), "events"),
    exposure   = vapply(x, `[[`, numeric(1L), "exposure"),
    expected   = vapply(x, `[[`, numeric(1L), "expected"),
    ratio      = vapply(x, `[[`, numeric(1L), "ratio"),
    ratio_lower = vapply(x, function(e) e$ratio_ci[1L], numeric(1L)),
    ratio_upper = vapply(x, function(e) e$ratio_ci[2L], numeric(1L)),
    forced     = vapply(x, `[[`, logical(1L), "forced"),
    stringsAsFactors = FALSE
  )
}


#' The bands, smooth or supplied rate of every transition as one data frame
#'
#' @param x A `micsim_reestimation`.
#' @param row.names,optional Ignored.
#' @param what "bands" (default), "smooth"
#' (the Poisson-spline fit: x, rate, lower, upper) or
#' "truth" (the supplied rate as drawn).
#' @param ... Not used.
#' @return A data frame with the columns `transition`, `axis` and `group` 
#' (the group's label), followed by the columns of each element's `what` table.
#' @export
as.data.frame.micsim_reestimation <- function(x, row.names = NULL,
                                              optional = FALSE,
                                              what = c("bands", "smooth",
                                                       "truth"), ...) {
  # `what` picks the layer of the figure: the band estimates (the
  # default, as before), the Poisson-spline smooth, or the supplied rate
  # as drawn in the left panel.
  what <- match.arg(what)
  parts <- lapply(x, function(e) {
    tbl <- e[[what]]
    if (is.null(tbl) || nrow(tbl) == 0L) return(NULL)
    cbind(transition = e$label, axis = e$axis,
          group = e$groups$label[tbl$group],
          tbl[, setdiff(names(tbl), "group"), drop = FALSE],
          stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, parts)
  if (is.null(out)) return(data.frame())
  row.names(out) <- NULL
  out
}
