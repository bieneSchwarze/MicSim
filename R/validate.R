# =====================================================================
# validate.R
# ---------------------------------------------------------------------
# WHAT THIS FILE DOES
#   The checks behind validate(). Unlike model construction, which stops
#   at the first structural fault because later checks depend on earlier
#   ones, validate() runs every check independently and reports
#   everything it finds at once.
#
# WHAT IS CHECKED
#   Errors (the run would fail or be wrong), per declared transition:
#     evaluation  the rate raises an error somewhere on the grid
#     type        it returns something other than numbers
#     length      it returns a different number of values than inputs
#     NA          it returns NA
#     negative    it returns a negative value
#   and per fertility subsystem:
#     newborn     a table-form newborn rule that inherits attributes has
#                 no candidate for some state a mother is in after a birth
#
#   Warnings (the run works, but the model is probably not what was
#   meant):
#     zero        the rate is zero everywhere on the grid
#     magnitude   the rate exceeds MAX_PLAUSIBLE_RATE per person-year
#                 (a rate per 1,000 given as a rate per person)
#     unreachable no individual can ever be in any origin state of the
#                 transition
#     stuck       individuals can reach a transient state that has no
#                 outgoing transition, which they leave only at max_age
#
#   Checks that depend on a failed one are skipped: a rate that cannot be
#   evaluated is not also reported as zero.
#
# THE REPORT
#   One `micsim_invalid_model` error listing every error (and, with
#   strict = TRUE, every warning), or one `micsim_model_warning` listing
#   the warnings. Both carry `$problems`, a data frame with one row per
#   finding: transition, check, severity, where, n_points, detail.
# =====================================================================


#' A finite rate above this, per person-year, is reported as implausible.
#' Ten per year means the event is expected within about five weeks.
#'
#' @keywords internal
#' @noRd
MAX_PLAUSIBLE_RATE <- 10

#' At most this many findings are listed in a message; the rest are
#' counted. All of them are in `$problems`.
#'
#' @keywords internal
#' @noRd
MAX_LISTED_PROBLEMS <- 20L


#' One finding, as a one-row data frame.
#'
#' @param transition The transition label.
#' @param check The check type.
#' @param severity The severity of the finding.
#' @param detail The detail of the finding.
#' @param where The location of the finding.
#' @param n_points The number of points at which the finding was found.
#'
#' @return A one-row data frame.
#'
#' @keywords internal
#' @noRd
validation_problem <- function(transition, check, severity, detail,
                               where = NA_character_, n_points = NA_integer_) {
  data.frame(transition = as.character(transition), check = check,
             severity = severity, where = as.character(where),
             n_points = as.integer(n_points), detail = detail,
             stringsAsFactors = FALSE)
}


#' Run every check on a model and act on the findings.
#'
#' @param model A `micsim_model`.
#' @param strict If TRUE, warnings are treated as errors.
#'
#' @return The model invisibly when there is no error; raises
#' `micsim_invalid_model` otherwise.
#'
#' @keywords internal
#' @noRd
validate_model <- function(model, strict = FALSE) {
  if (!isTRUE(strict) && !isFALSE(strict)) {
    micsim_abort("micsim_model_error", "`strict` must be TRUE or FALSE.")
  }
  problems <- rbind(
    check_declared_rates(model, reachable_grid(model)),
    check_reachability(model),
    check_newborn_coverage(model)
  )
  if (is.null(problems) || nrow(problems) == 0L) return(invisible(model))
  row.names(problems) <- NULL

  failing <- problems$severity == "error" |
    (strict & problems$severity == "warning")
  if (any(failing)) {
    abort_invalid_model(problems[failing, , drop = FALSE])
  }
  warning(micsim_warning(
    "micsim_model_warning",
    format_problems(problems, "warning(s) about the simulation model"),
    problems = problems
  ))
  invisible(model)
}


#' The findings as a message: one line per finding, at most `MAX_LISTED_PROBLEMS` of them.
#'
#' @param problems A data frame with one row per finding.
#' @param what The type of finding, problem(s) found in the simulation model.
#'
#' @return A string.
#'
#' @keywords internal
#' @noRd
format_problems <- function(problems, what) {
  shown <- utils_head(seq_len(nrow(problems)), MAX_LISTED_PROBLEMS)
  lines <- vapply(shown, function(i) {
    p <- problems[i, ]
    where <- if (is.na(p$transition)) "" else paste0("[", p$transition, "] ")
    paste0("  - ", where, gsub("\n", "\n    ", p$detail, fixed = TRUE))
  }, character(1L))
  more <- nrow(problems) - length(shown)
  paste0(
    nrow(problems), " ", what, ":\n",
    paste(lines, collapse = "\n"),
    if (more > 0L) paste0("\n  ... and ", more, " more (see `$problems` of the condition).") else ""
  )
}


#' Rate checks, one set per declared transition.
#'
#' @description
#'  Each declaration is checked once on the validation grid, however many
#'  origin states it expands into, and reported under the label the user
#'  wrote (as in the reestimation).
#'
#' @param model A `micsim_model`.
#' @param probe A `micsim_rate_probe`.
#'
#' @return A data frame with one row per finding.
#'
#' @keywords internal
#' @noRd
check_declared_rates <- function(model, probe) {
  registry <- model$registry
  decls <- unclass(model$transitions)
  decl_of <- vapply(registry$edges, function(e) {
    if (is.null(e$decl)) NA_integer_ else as.integer(e$decl)
  }, integer(1L))

  out <- list()
  for (d in seq_along(decls)) {
    edges <- registry$edges[which(decl_of == d)]
    if (length(edges) == 0L) next
    label <- describe_declaration(decls[[d]], edges)$label
    out[[length(out) + 1L]] <- check_one_rate(edges[[1L]], label, probe)
  }
  do.call(rbind, out)
}


#' The checks of one rate function on the grid; NULL if it passes.
#'
#' @param edge The edge of the declaration.
#' @param label The label of the declaration.
#' @param probe The probe of the declaration.
#'
#' @return A data frame with one row per finding.
#'
#' @keywords internal
#' @noRd
check_one_rate <- function(edge, label, probe) {
  grid_size <- function() {
    n <- 1L
    if ("age" %in% edge$argnames)      n <- n * length(probe$ages)
    if ("calTime" %in% edge$argnames)  n <- n * length(probe$cal_years)
    if ("duration" %in% edge$argnames) n <- n * length(probe$durations)
    n
  }
  span <- paste0("ages ", min(probe$ages), " to ", max(probe$ages),
                 ", calendar years ", min(probe$cal_years), " to ",
                 max(probe$cal_years))

  values <- tryCatch(evaluate_rate_on_grid(edge$rate, edge$argnames, probe),
                     error = function(e) e)
  if (inherits(values, "condition")) {
    return(validation_problem(label, "evaluation", "error", paste0(
      "The rate raises an error on the validation grid (", span, "):\n",
      conditionMessage(values))))
  }
  if (!is.numeric(values)) {
    return(validation_problem(label, "type", "error", paste0(
      "The rate returns ", class(values)[1L], " instead of numbers.")))
  }
  expected <- if (length(edge$argnames) == 0L) 1L else grid_size()
  if (length(values) != expected) {
    return(validation_problem(label, "length", "error", paste0(
      "The rate returns ", length(values), " value(s) for ", expected,
      " grid point(s). Rate functions must return one value per element\n",
      "of their input, e.g. function(age) rep(0.1, length(age)).")))
  }

  found <- list()
  is_na <- is.na(values)
  if (any(is_na)) {
    at <- describe_grid_point(edge$argnames, probe, which(is_na)[1L])
    found$na <- validation_problem(label, "NA", "error", paste0(
      "The rate returns NA at ", sum(is_na), " of ", length(values),
      " grid point(s), first at ", at, "."), where = at, n_points = sum(is_na))
  }
  negative <- !is_na & values < 0
  if (any(negative)) {
    at <- describe_grid_point(edge$argnames, probe, which(negative)[1L])
    found$neg <- validation_problem(label, "negative", "error", paste0(
      "The rate returns a negative value at ", sum(negative), " of ",
      length(values), " grid point(s), first at ", at, "."),
      where = at, n_points = sum(negative))
  }
  if (length(found) > 0L) return(do.call(rbind, found))

  # The rate is well defined; now the warnings.
  if (all(values == 0)) {
    found$zero <- validation_problem(label, "zero", "warning", paste0(
      "The rate is zero at every grid point (", span, "),\n",
      "so this transition can never happen."))
  }
  large <- is.finite(values) & values > MAX_PLAUSIBLE_RATE
  if (any(large)) {
    at <- describe_grid_point(edge$argnames, probe, which(large)[1L])
    found$big <- validation_problem(label, "magnitude", "warning", paste0(
      "The rate exceeds ", MAX_PLAUSIBLE_RATE, " per person-year at ",
      sum(large), " grid point(s) (largest ", signif(max(values[large]), 3),
      ", first at ", at, ").\n",
      "Rates are per person-year; a rate per 1,000 must be divided by 1,000."),
      where = at, n_points = sum(large))
  }
  if (length(found) == 0L) NULL else do.call(rbind, found)
}


#' States individuals can enter, and the states reachable from them.
#'
#' @description
#'  Entry states are the starting population's states, the immigrants'
#'  states on arrival and, if a subsystem adds individuals, the states its
#'  newborn rule can assign.
#'
#' @param model A `micsim_model`.
#'
#' @return A logical vector over the codec's state indices.
#'
#' @keywords internal
#' @noRd
reachable_states <- function(model) {
  codec <- model$codec
  entry <- c(as.character(model$population$initState),
             if (!is.null(model$immigrants))
               as.character(model$immigrants$immigrInitState))
  for (s in model$subsystems) {
    if (isTRUE(s$spawns) && !is.null(s$params$newborn)) {
      rule <- tryCatch(resolve_newborn_rule(s$params$newborn, codec, s$name),
                       error = function(e) NULL)
      if (!is.null(rule)) {
        entry <- c(entry, apply(rule$states, 1L, paste, collapse = STATE_SEP))
      }
    }
  }
  reached <- logical(codec$n)
  frontier <- unique(stats::na.omit(match(entry, codec$states)))
  reached[frontier] <- TRUE
  while (length(frontier) > 0L) {
    nxt <- integer(0)
    for (i in frontier) {
      for (e in model$registry$feasible(i)) {
        if (!reached[e$to_idx]) {
          reached[e$to_idx] <- TRUE
          nxt <- c(nxt, e$to_idx)
        }
      }
    }
    frontier <- nxt
  }
  reached
}


#' Transitions nobody can be at risk of, and states nobody can leave.
#'
#' @param model A `micsim_model`.
#'
#' @return A data frame with one row per finding.
#'
#' @keywords internal
#' @noRd
check_reachability <- function(model) {
  codec <- model$codec
  registry <- model$registry
  reached <- reachable_states(model)
  out <- list()

  decls <- unclass(model$transitions)
  decl_of <- vapply(registry$edges, function(e) {
    if (is.null(e$decl)) NA_integer_ else as.integer(e$decl)
  }, integer(1L))
  for (d in seq_along(decls)) {
    edges <- registry$edges[which(decl_of == d)]
    if (length(edges) == 0L) next
    origins <- vapply(edges, `[[`, integer(1L), "from_idx")
    if (!any(reached[origins])) {
      label <- describe_declaration(decls[[d]], edges)$label
      out[[length(out) + 1L]] <- validation_problem(
        label, "unreachable", "warning", paste0(
          "No individual can reach any of its ", length(unique(origins)),
          " origin state(s) from the starting\n",
          "population, the immigrants or the newborn rule, so this transition\n",
          "can never happen."))
    }
  }

  transient <- seq_len(codec$K)
  stuck <- transient[reached[transient] &
                       vapply(transient, function(i) is.null(registry$feasible(i)),
                              logical(1L))]
  if (length(stuck) > 0L) {
    labels <- codec$decode(stuck)
    out[[length(out) + 1L]] <- validation_problem(
      NA_character_, "stuck", "warning", paste0(
        "Individuals can reach ", length(stuck), " transient state(s) with no ",
        "outgoing transition: ", toString(utils_head(labels, 5L)),
        if (length(labels) > 5L) paste0(" and ", length(labels) - 5L, " more") else "",
        ".\nThey stay there until max_age (", model$max_age, ") ends their ",
        "life course."),
      n_points = length(stuck))
  }
  do.call(rbind, out)
}


#' A newborn rule that inherits attributes must have a candidate for every
#' state a mother can be in after a birth. Otherwise the run stops at the
#' first such birth.
#'
#' @param model A `micsim_model`.
#'
#' @return A data frame with one row per finding.
#'
#' @keywords internal
#' @noRd
check_newborn_coverage <- function(model) {
  codec <- model$codec
  out <- list()
  for (s in model$subsystems) {
    if (!isTRUE(s$spawns) || is.null(s$params$newborn) ||
          is.null(s$claims$births)) next
    rule <- tryCatch(resolve_newborn_rule(s$params$newborn, codec, s$name),
                     error = function(e) NULL)
    if (is.null(rule) || length(rule$inherit) == 0L) next
    mothers <- unique(vapply(model$registry$edges[s$claims$births], `[[`,
                             character(1L), "to_label"))
    keys <- apply(rule$states[, rule$inherit, drop = FALSE], 1L, paste,
                  collapse = STATE_SEP)
    uncovered <- mothers[!vapply(mothers, function(m) {
      parts <- strsplit(m, STATE_SEP, fixed = TRUE)[[1L]]
      paste(parts[rule$inherit], collapse = STATE_SEP) %in% keys
    }, logical(1L))]
    if (length(uncovered) > 0L) {
      out[[length(out) + 1L]] <- validation_problem(
        paste0("subsystem `", s$name, "`"), "newborn", "error", paste0(
          "The newborn rule has no candidate state for a mother in state(s) ",
          toString(utils_head(uncovered, 5L)),
          if (length(uncovered) > 5L) paste0(" and ", length(uncovered) - 5L, " more") else "",
          ".\nEvery combination of inherited attributes a mother can hold after a\n",
          "birth must appear in the newborn-state table."),
        n_points = length(uncovered))
    }
  }
  do.call(rbind, out)
}
