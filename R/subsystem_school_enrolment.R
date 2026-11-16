# =====================================================================
# subsystem_school_enrolment.R
# ---------------------------------------------------------------------
# WHAT THIS FILE DOES
#   Adjusts the timing of selected transitions to the date a calendar
#   convention says it must happen.
#
# PLACE IN THE PIPELINE
#   school_enrolment(move("edu", "none", "primary"), month = 9)
#       --add_subsystem()--> compiled against the declared transitions
#   engine.R --(fires on_candidates with every destination's drawn time)-->
#       subsystem_school_enrolment.R --(enrolment candidates on new dates)-->
#   engine.R --(takes the earliest, pushes it)--> engine_queue.R
#
# WHAT YOU NEED TO KNOW  (concepts in ?MicSim-concepts)
#   forced transition.
#
#   Children do not start school on their birthday. They start on the enrolment
#   month of the year in which they are eligible to start school. Those are 
#   two different things, and the model expresses them in two different places:
#
#     the RATE says   "this happens at age seven"      (an Inf rate)
#     the CALENDAR says "on 1 September"               (this file)
#
#   The rate produces a date; this file replaces it. That is the whole
#   mechanism, and it is the second of the two kinds of hook: an ADJUSTING
#   hook, given a value and returning one of the same shape.
#
#   The engine compares candidate times after the adjustment. A death or
#   emigration scheduled before the adjusted enrolment can occur first.
#
# GUARANTEES
#   * It changes only the date, never the destination.
#   * It leaves every non-enrolment transition exactly as drawn.
# =====================================================================


#' Adjust transitions to the school enrolment calendar
#'
#' @description
#' A subsystem for [add_subsystem()]: every time one of the given
#' transitions is drawn, its date is moved to the first day of `month` --
#' in the same year if the drawn date falls in or before that month,
#' otherwise in the next year.
#'
#' @details
#' Elementary-school enrolment is written in two halves. The rate says
#' WHEN in a person's life the transition happens, e.g.,
#' `function(age) ifelse(age == 6, Inf, 0)`, an `Inf` rate that makes the
#' transition *forced*: it happens on that date rather than at a random
#' one. This subsystem then says on WHICH calendar day it happens.
#'
#' So with `month = 9`, a child whose rate fires in March 2018 starts
#' school on 1 September 2018, and one whose rate fires in October 2018 on
#' 1 September 2019.
#'
#' Which transition is the enrolment is up to the model: name it, e.g.
#' `move("edu", "none", "primary")`. No education label has a meaning of
#' its own, so a model can call its levels anything, and can have more
#' than one education variable.
#'
#' The date is moved before competing risks are resolved: the enrolment
#' competes with every other transition at its calendar date, so a child
#' who would die or emigrate before 1 September does so. A date that the
#' rule would move into the past -- a child scheduled after the 1st of the
#' enrolment month with an enrolment due now -- is placed at the present
#' instead.
#'
#' @section Report:
#' After a run, `result$extensions[[name]]` holds
#' \describe{
#'   \item{`n_enrolments`}{the number of enrolments during the run;}
#'   \item{`enrolments`}{one row per enrolment: `ID`, `date` (`yyyymmdd`,
#'     the enrolment date), `age` (years, two decimals), `transition`
#'     (`"from -> to"`), `drawnDate` (the date the rate alone produced) and
#'     `shiftDays` (enrolment date minus drawn date: positive when the
#'     enrolment was postponed to the next enrolment month, negative when it
#'     was brought forward to the first day of the month it was drawn in);}
#'   \item{`n_at_present`}{enrolments the rule would have moved into the
#'     past, placed on the day the child was scheduled instead;}
#'   \item{`by_transition`}{enrolments per enrolment transition;}
#'   \item{`by_month`}{enrolments by calendar month of the drawn date,
#'     `"01"` to `"12"`: how much the calendar rule moved.}
#' }
#'
#' @param transition Which declared transition(s) the rule applies to: an
#'   [move()] or [moves_of()] selector, a list of them, or transition
#'   labels.
#' @param month The enrolment month, a whole number from 1 to 12.
#' @param name The subsystem's name in the model and in
#'   `result$extensions`.
#'
#' @return A subsystem, for [add_subsystem()].
#'
#' @seealso [add_subsystem()], [move()]
#'
#' @examples
#' school_enrolment(move("edu", "none", "primary"), month = 9)
#'
#' @export
school_enrolment <- function(transition, month, name = "school_enrolment") {

  if (missing(transition)) {
    abort_enrolment(paste0(
      "school_enrolment() needs `transition`: the transition that is the\n",
      "  enrolment, e.g. move(\"edu\", \"none\", \"primary\")."
    ))
  }
  month <- check_enrolment_month(month)

  new_subsystem(
    name   = name,
    type   = "school_enrolment",
    params = list(transition = transition, month = month),
    compile = function(ctx) {
      what <- paste0("school_enrolment(transition = ) of subsystem `", name, "`")
      ids <- resolve_claim(transition, ctx$codec, ctx$registry, what)
      compile_school_enrolment(ids, month, ctx$registry)
    }
  )
}


#' Prepare the school enrolment subsystem for a model
#'
#' @description
#' Records which transitions the subsystem controls. Sets up hooks to adjust
#' candidate times before the engine chooses the next event and to log
#' enrolment events when they occur.
#' 
#' @param ids Positions of the claimed edges in the registry.
#' @param month The enrolment month.
#' @param registry The model's `rate_registry`.
#' @return The compiled form; see `R/subsystem.R`.
#' @keywords internal
#' @noRd
compile_school_enrolment <- function(ids, month, registry) {

  n_edges <- length(registry$edges)
  claimed <- logical(n_edges)
  claimed[ids] <- TRUE
  labels <- vapply(registry$edges, `[[`, character(1L), "label")

  list(
    claims   = list(timing = ids),
    describe = paste0("school enrolment (", length(ids),
                      " transition(s), month ", month, ")"),
    merge    = merge_reports,
    instantiate = function(run) {
      st <- new.env(parent = emptyenv())
      st$n_enrolments <- 0L
      st$n_at_present <- 0L
      # The date each person's enrolment candidate had before the calendar
      # rule moved it, by slot: the candidate that later fires is the one
      # set at the person's most recent scheduling.
      st$drawn <- numeric(0)
      st$at_present <- logical(0)
      st$log <- new_event_log(list(ID = integer(), time_days = numeric(),
                                   drawn_days = numeric(),
                                   birth_days = numeric(),
                                   edge_id = integer()))

      on_candidates <- function(cand, ctx) {
        move <- which(claimed[cand$edge_id] & is.finite(cand$time_days))
        for (k in move) {
          if (ctx$slot > length(st$drawn)) {
            n <- max(2L * length(st$drawn), ctx$slot, 64L)
            length(st$drawn) <- n
            length(st$at_present) <- n
          }
          rule_date <- enrolment_date(cand$time_days[k], month)
          st$drawn[ctx$slot] <- cand$time_days[k]
          st$at_present[ctx$slot] <- rule_date < ctx$now_days
          cand$time_days[k] <- max(rule_date, ctx$now_days)
        }
        cand
      }

      on_event <- function(ctx) {
        if (claimed[ctx$edge$id]) {
          st$n_enrolments <- st$n_enrolments + 1L
          drawn <- st$drawn[ctx$slot]
          if (isTRUE(st$at_present[ctx$slot])) {
            st$n_at_present <- st$n_at_present + 1L
          }
          st$log$add(ID = ctx$persons$id_of(ctx$slot),
                     time_days = ctx$now_days, drawn_days = drawn,
                     birth_days = ctx$persons$birth_days_of(ctx$slot),
                     edge_id = ctx$edge$id)
        }
        invisible(NULL)
      }

      list(hooks  = new_engine_hooks(on_candidates = on_candidates,
                                     on_event = on_event),
           finish = function() enrolment_report(st, labels, ids))
    }
  )
}


#' Summarise the enrolments recorded during a run
#'
#' @description
#' Converts the event log into a table of enrolments and counts by transition
#' and original candidate month. The report fields are described in
#' [school_enrolment()].
#'
#' @param st The run's state environment.
#' @param labels The labels of all edges.
#' @param ids The claimed edge ids.
#' @return A named list.
#' @keywords internal
#' @noRd
enrolment_report <- function(st, labels, ids) {
  cols <- st$log$columns()
  drawn_ymd <- as.character(getInDateFormat(cols$drawn_days))
  enrolments <- data.frame(
    ID         = as.integer(cols$ID),
    date       = as.character(getInDateFormat(cols$time_days)),
    age        = round((cols$time_days - cols$birth_days) / DAYS_PER_YEAR, 2),
    transition = labels[cols$edge_id],
    drawnDate  = drawn_ymd,
    shiftDays  = as.integer(floor(cols$time_days) - floor(cols$drawn_days)),
    stringsAsFactors = FALSE
  )
  list(
    n_enrolments  = st$n_enrolments,
    enrolments    = enrolments,
    n_at_present  = st$n_at_present,
    by_transition = count_by(enrolments$transition, labels[ids]),
    by_month      = count_by(substr(drawn_ymd, 5L, 6L), sprintf("%02d", 1:12))
  )
}


#' Check an enrolment month
#' 
#' @param month A single numeric value representing a whole month from 1 to 12.
#' @return The month as an integer. Invalid inputs raise an enrolment error.
#' @keywords internal
#' @noRd
check_enrolment_month <- function(month) {
  whole <- is.numeric(month) && length(month) == 1L && !is.na(month) &&
    month == round(month)
  month <- if (whole) as.integer(month) else NA_integer_
  if (is.na(month) || month < 1L || month > 12L) {
    abort_enrolment("The school-enrolment month must be a whole number from 1 to 12.")
  }
  month
}


#' Find the target date for an enrolment candidate
#'
#' @description
#' Finds the first day of the enrolment month in the candidate year if the
#' candidate falls in or before that month. Uses the following year if it
#' falls later. 
#'
#' @param when_days The drawn event date, as a day count.
#' @param month Enrolment month, 1 to 12.
#' @return The target date as a numeric day count.
#' @keywords internal
#' @noRd
enrolment_date <- function(when_days, month) {
  year <- getYear(when_days)
  if (getMonth(when_days) > month) {
    year <- year + 1L
  }
  as.numeric(getFirstDayOfMonth(year, month))
}


#' Raise a school-enrolment error
#' @param problem The message.
#' @return Never returns.
#' @keywords internal
#' @noRd
abort_enrolment <- function(problem) {
  micsim_abort("micsim_enrolment_error", problem)
}
