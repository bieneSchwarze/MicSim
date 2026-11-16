# =====================================================================
# class_newborn_rule.R
# ---------------------------------------------------------------------
# WHAT THIS FILE DOES
#   Says what state a child is in at the moment of birth.
#
# PLACE IN THE PIPELINE
#   the model spec --> class_newborn_rule.R --> subsystem_fertility.R
#
#   This file is pure specification. It holds no simulation state, knows
#   nothing about the engine, and does exactly one thing when asked: given
#   the mother's state, return the newborn's.
#
# WHAT YOU NEED TO KNOW  (concepts in ?MicSim-concepts)
#   composite state, state variable.
#
#   A newborn's state is settled in two parts. Some of its state variables
#   are INHERITED from the mother, because they are not free to vary: a
#   child born to an Austrian mother is Austrian. The rest are DRAWN, from
#   a list of admissible states with given probabilities: the sex of a
#   newborn is drawn roughly 0.515 male, 0.485 female.
#
#   The two interact, and that is the whole subtlety of this file. The
#   probabilities are CONDITIONAL on what was inherited. If nationality is
#   inherited and sex is drawn, the sex ratio must sum to one separately
#   within each nationality, not once across the whole table.
#
# GUARANTEES
#   * A rule that has been constructed is valid: the conditional
#     probabilities have been checked to sum to one for every combination
#     of inherited values that appears in the table.
#   * Drawing a newborn's state uses exactly one uniform, which the
#     fertility subsystem takes from the mother's birth block.
# =====================================================================


#' Define how a newborn's state is assigned
#'
#' @description
#' Builds the rule that decides which state a child occupies at birth.
#' The rule defines which state variables the child takes from its mother
#' and how the remaining ones are drawn.
#'
#' @details
#' `states` is the table of states a newborn can take. It has one row per
#' possible combination and one column per state variable of the model.
#' The columns follow the model's state variables in order. 
#' `prob` assigns a probability between 0 and 1 for each row to be selected.
#' `inherit` names the column positions the child takes from its mother.
#'
#' At a birth the rule keeps only the rows that agree with the mother on
#' every inherited position, and draws one of them with the corresponding
#' probability. When nothing is inherited, it draws from the whole table to 
#' decide a newborn's state.
#'
#' In the example from `?micSim`, the state variables are sex,
#' nationality and parity. Nationality is inherited and sex is drawn with a
#' sex ratio that differs by nationality; parity
#' is `"0"`.
#'
#' ```
#' states = rbind(c("m","DE","0"), c("f","DE","0"),
#'                c("m","AT","0"), c("f","AT","0"),
#'                c("m","IT","0"), c("f","IT","0"))
#' prob    = c(0.515, 0.485,  0.515, 0.485,  0.515, 0.485)
#' inherit = 2
#' ```
#'
#' A child of an Austrian mother can only be row 3 or row 4, and those two
#' weights sum to one. This is checked by the validation separately 
#' for each unique inherited attribute.
#'
#' @section Assumptions:
#' The weights within each inherited group must sum to one within
#' `PROB_TOL`. A table whose weights do not sum to one in some group returns
#' an error, because a mother in that group would otherwise have no
#' well-defined newborn distribution.
#'
#' A state variable that is neither inherited nor allowed to vary in
#' `states` takes the single value it has in every row. Parity for example will
#' usually have every newborn starting at `"0"`.
#'
#' @param states Character matrix of admissible newborn states. 
#' A data frame is accepted and converted.
#' @param prob Numeric weights, one per row of `states`.
#' @param inherit Integer column positions the newborn takes from its
#'   mother. When empty nothing is inherited.
#'
#' @return An object of class `newborn_rule`.
#'
#' @seealso [newborn_rule()] for the per-variable form, [fertility()],
#'   `assign_newborn_state()`,
#'   [MicSim-concepts]
#'
#' @examples
#' rule <- new_newborn_rule(
#'   states = rbind(c("m", "DE", "0"), c("f", "DE", "0"),
#'                  c("m", "AT", "0"), c("f", "AT", "0")),
#'   prob    = c(0.515, 0.485, 0.515, 0.485),
#'   inherit = 2L
#' )
#' rule
#'
#' @export
new_newborn_rule <- function(states, prob, inherit = integer(0)) {

  if (is.data.frame(states)) {
    states <- as.matrix(states)
  }
  if (is.character(states) && is.null(dim(states))) {
    # A single candidate written as a plain vector is a common and
    # harmless thing to do; treat it as a one-row table.
    states <- matrix(states, nrow = 1L)
  }
  if (!is.matrix(states) || nrow(states) == 0L) {
    abort_newborn_rule(
      "`states` must be a matrix of admissible newborn states, with at\n  least one row."
    )
  }
  states <- matrix(as.character(states), nrow = nrow(states))

  prob <- as.numeric(prob)
  if (length(prob) != nrow(states)) {
    abort_newborn_rule(paste0(
      "`prob` has ", length(prob),
      if (length(prob) == 1) " probability" else " probabilities",
      " but `states` has ", nrow(states), " row(s). ",
      "There must be exactly one probability per row."
    ))
  }
  if (anyNA(prob) || any(prob < 0)) {
    abort_newborn_rule("`prob` must be non-negative and free of NA.")
  }

  inherit <- as.integer(inherit)
  if (length(inherit) > 0L &&
        (anyNA(inherit) || any(inherit < 1L) || any(inherit > ncol(states)))) {
    abort_newborn_rule(paste0(
      "`inherit` must name column positions of `states`, i.e. whole\n",
      "  numbers between 1 and ", ncol(states), "."
    ))
  }

  check_newborn_probabilities(states, prob, inherit)

  structure(
    list(states = states, prob = prob, inherit = inherit,
         n_vars = ncol(states)),
    class = "newborn_rule"
  )
}


#' Say what state a child is born into
#'
#' @description
#' Gives, for every state variable of the model, how a newborn's value is
#' settled: fixed, drawn, or taken from the mother. For [fertility()].
#'
#' @details
#' Each state variable gets it's own argument of the same name:
#'
#' \describe{
#'   \item{a single value}{something every child gets, e.g. `parity = "0"`}
#'   \item{a named vector of probabilities}{the value is drawn, e.g.
#'     `sex = c(male = 0.515, female = 0.485)`. The probabilities must sum
#'     to one.}
#'   \item{`from_mother()`}{the child takes the mother's value, e.g.
#'     `nat = from_mother()`. The mother's value is the one she holds
#'     right after the birth transition.}
#' }
#'
#' Drawn variables are drawn independently of one another. When a
#' probability must depend on an inherited value -- e. g. a sex ratio that
#' differs by nationality -- use the table form, [new_newborn_rule()].
#'
#' The rule is checked against the model when the subsystem is added. To pass,
#' every state variable needs an entry, and every value must exist.
#'
#' @param ... One named entry per state variable.
#'
#' @return An object of class `micsim_newborn_spec`.
#'
#' @seealso [fertility()], [new_newborn_rule()]
#'
#' @examples
#' newborn_rule(sex = c(male = 0.515, female = 0.485),
#'              nat = from_mother(),
#'              parity = "0")
#'
#' @export
newborn_rule <- function(...) {
  entries <- list(...)
  nms <- names(entries)
  if (length(entries) == 0L || is.null(nms) || any(!nzchar(nms)) ||
        anyDuplicated(nms) > 0L) {
    abort_newborn_rule(paste0(
      "newborn_rule() takes one NAMED entry per state variable, e.g.\n",
      "  newborn_rule(sex = c(male = 0.515, female = 0.485), parity = \"0\")."
    ))
  }
  for (v in nms) {
    e <- entries[[v]]
    if (inherits(e, "micsim_from_mother")) next
    if (is.character(e) && length(e) == 1L && !is.na(e)) next
    if (is.numeric(e) && length(e) >= 1L && !is.null(names(e)) &&
          all(nzchar(names(e))) && !anyNA(e) && all(e >= 0)) {
      if (!isTRUE(all.equal(sum(e), 1, tolerance = PROB_TOL))) {
        abort_newborn_rule(paste0(
          "The probabilities for `", v, "` sum to ", signif(sum(e), 6),
          ", not 1."
        ))
      }
      next
    }
    abort_newborn_rule(paste0(
      "The entry for `", v, "` must be a single value, a named vector of\n",
      "  probabilities, or from_mother()."
    ))
  }
  structure(list(entries = entries), class = "micsim_newborn_spec")
}


#' @rdname newborn_rule
#' @export
from_mother <- function() {
  structure(list(), class = "micsim_from_mother")
}


#' Print a newborn rule specification
#' @param x A `micsim_newborn_spec`.
#' @param ... Not used.
#' @return `x` invisible.
#' @exportS3Method print micsim_newborn_spec
#' @noRd
print.micsim_newborn_spec <- function(x, ...) {
  cat("<micsim_newborn_spec>\n")
  for (v in names(x$entries)) {
    e <- x$entries[[v]]
    how <- if (inherits(e, "micsim_from_mother")) {
      "from the mother"
    } else if (is.numeric(e)) {
      paste0("drawn: ", paste(names(e), signif(e, 4), sep = " ", collapse = ", "))
    } else {
      paste0("always ", e)
    }
    cat("  ", v, ": ", how, "\n", sep = "")
  }
  invisible(x)
}


#' Check a newborn rule against a model and turn it into the table form
#'
#' @description
#' Accepts either form of newborn rule and returns a validated table-form
#' `newborn_rule` in which every candidate state exists. Its columns follow
#' the model's state variables in order (`codec$vars`).
#'
#' @details
#' A [newborn_rule()] specification becomes the table of every
#' combination of its drawn values, with each inherited variable taking
#' every value it has in the model. A row's probability is the product of its
#' drawn probabilities, so within every combination of inherited values
#' the probabilities sum to one.
#'
#' @param rule A `micsim_newborn_spec` or a `newborn_rule`.
#' @param codec A `state_codec`.
#' @param name The subsystem's name, for messages.
#' @return A `newborn_rule`.
#' @keywords internal
#' @noRd
resolve_newborn_rule <- function(rule, codec, name = "fertility") {

  if (inherits(rule, "micsim_newborn_spec")) {
    entries <- rule$entries
    missing_vars <- setdiff(codec$vars, names(entries))
    extra_vars <- setdiff(names(entries), codec$vars)
    if (length(missing_vars) > 0L || length(extra_vars) > 0L) {
      abort_newborn_rule(paste0(
        "Subsystem `", name, "`: the newborn rule must give exactly one entry\n",
        "  per state variable (", toString(codec$vars), ").",
        if (length(missing_vars)) paste0("\n  Missing: ", toString(missing_vars), ".") else "",
        if (length(extra_vars)) paste0("\n  Not state variables: ", toString(extra_vars), ".") else ""
      ))
    }
    values <- list()
    weights <- list()
    inherit <- integer(0)
    for (i in seq_along(codec$vars)) {
      v <- codec$vars[i]
      e <- entries[[v]]
      levels <- codec$levels[[i]]
      if (inherits(e, "micsim_from_mother")) {
        values[[v]]  <- levels
        weights[[v]] <- rep(1, length(levels))
        inherit <- c(inherit, i)
      } else if (is.numeric(e)) {
        values[[v]]  <- names(e)
        weights[[v]] <- as.numeric(e)
      } else {
        values[[v]]  <- e
        weights[[v]] <- 1
      }
      unknown <- setdiff(values[[v]], levels)
      if (length(unknown) > 0L) {
        abort_newborn_rule(paste0(
          "Subsystem `", name, "`: the newborn rule gives `", v, "` the value(s) ",
          toString(unknown), ",\n  which it does not take. Values of `", v,
          "`: ", toString(levels), "."
        ))
      }
    }
    combos <- expand.grid(rev(values), stringsAsFactors = FALSE,
                          KEEP.OUT.ATTRS = FALSE)[, rev(seq_along(values)),
                                                  drop = FALSE]
    wgrid <- expand.grid(rev(weights), KEEP.OUT.ATTRS = FALSE)[
      , rev(seq_along(weights)), drop = FALSE]
    rule <- new_newborn_rule(
      states  = as.matrix(combos),
      prob    = apply(wgrid, 1L, prod),
      inherit = inherit
    )
  }

  if (!inherits(rule, "newborn_rule")) {
    abort_newborn_rule(
      "The newborn rule must come from newborn_rule() or new_newborn_rule()."
    )
  }
  if (rule$n_vars != length(codec$vars)) {
    abort_newborn_rule(paste0(
      "Subsystem `", name, "`: the newborn rule describes ", rule$n_vars,
      " state variable(s) but the model has ", length(codec$vars), ".\n",
      "  The candidate-state table needs one column per state variable,\n",
      "  in the model's variable order (", toString(codec$vars), ")."
    ))
  }

  # Every candidate newborn state must exist in the model, or the first
  # birth would fail deep inside the loop.
  candidate_labels <- apply(rule$states, 1L, paste, collapse = STATE_SEP)
  unknown <- setdiff(candidate_labels, codec$transient)
  if (length(unknown) > 0L) {
    abort_newborn_rule(paste0(
      "Subsystem `", name, "`: newborn state(s) not in the state space: ",
      toString(utils_head(unknown, 5L)), "."
    ))
  }
  rule
}


#' Check that the probabilities sum to one within every inherited group
#'
#' @description
#' At construction this function checks that all probabilities in either the
#' whole table or each distinct combination of inherited values sum to 1.
#'
#' @details
#' With nothing inherited there is one group, the whole table, and its
#' probabilities must sum to one. With something inherited, each distinct
#' combination of inherited values forms its own group. Each group's
#' probabilities must sum to one on their own, because a mother in that group
#' can only produce a child from that group.
#'
#' Every failing group is reported together, so all errors are visible at once.
#'
#' @param states,prob,inherit As in `new_newborn_rule()`.
#' @return `invisible(TRUE)`, or raises `micsim_newborn_rule_error`.
#' @keywords internal
#' @noRd
check_newborn_probabilities <- function(states, prob, inherit) {

  group_of <- if (length(inherit) == 0L) {
    rep("all candidates", nrow(states))
  } else {
    apply(states[, inherit, drop = FALSE], 1L, paste, collapse = "/")
  }

  sums <- tapply(prob, group_of, sum)
  wrong <- names(sums)[
    !vapply(sums, function(s) isTRUE(all.equal(s, 1, tolerance = PROB_TOL)),
            logical(1L))
  ]

  if (length(wrong) > 0L) {
    detail <- paste0("    ", wrong, ": ", signif(sums[wrong], 6),
                     collapse = "\n")
    abort_newborn_rule(paste0(
      "The probabilities of assigning a state to a newborn must sum to 1\n",
      if (length(inherit) > 0L) {
        "  separately for each combination of inherited attributes.\n"
      } else {
        "  across all candidate states.\n"
      },
      "  These group(s) do not:\n", detail
    ))
  }

  invisible(TRUE)
}


#' Draw a newborn's state
#'
#' @description
#' Applies the rule: keep the candidates that match the mother on every
#' inherited position, then draw one of them.
#'
#' @details
#'
#' The draw inverts the candidates' cumulative probabilities with one uniform `u`.
#' Candidates with probability zero are dropped first; the first candidate whose
#' cumulative probability exceeds `u` times the total is chosen.
#'
#' @param rule A `newborn_rule`.
#' @param mother_state The mother's state label at the moment of birth,
#'   e.g. `"f/AT/1"`. This is the state she moves INTO, not the one she
#'   leaves.
#' @param u A uniform in (0, 1). The fertility subsystem passes one from the
#'   mother's birth block; the default draws one from the current stream.
#'
#' @return The newborn's state label, e.g. `"m/AT/0"`.
#' @seealso `new_newborn_rule()`
#'
#' @keywords internal
#' @noRd
assign_newborn_state <- function(rule, mother_state, u = stats::runif(1L)) {

  candidates <- seq_len(nrow(rule$states))

  if (length(rule$inherit) > 0L) {
    mother_parts <- strsplit(mother_state, "/", fixed = TRUE)[[1L]]
    if (length(mother_parts) != rule$n_vars) {
      abort_newborn_rule(paste0(
        "The mother is in state `", mother_state, "`, which has ",
        length(mother_parts), " component(s),\n",
        "  but the newborn rule is written for ", rule$n_vars,
        " state variable(s)."
      ))
    }
    wanted <- mother_parts[rule$inherit]
    matches <- apply(rule$states[, rule$inherit, drop = FALSE], 1L,
                     function(row) all(row == wanted))
    candidates <- candidates[matches]
  }

  if (length(candidates) == 0L) {
    abort_newborn_rule(paste0(
      "No candidate newborn state matches a mother in state `",
      mother_state, "`.\n",
      "  Every combination of inherited attributes that a mother can hold\n",
      "  must appear in the newborn-state table."
    ))
  }

  # Inversion of the candidates' cumulative probabilities with one
  # uniform `u`. Candidates with probability zero are dropped first, so
  # they can never be chosen, even through rounding at an interval edge.
  candidates <- candidates[rule$prob[candidates] > 0]
  cum <- cumsum(rule$prob[candidates])
  pick <- findInterval(u * cum[length(cum)], cum) + 1L
  chosen <- candidates[min(pick, length(candidates))]

  paste(rule$states[chosen, ], collapse = "/")
}


#' Raise a newborn-rule error
#' @param problem Message reporting the error.
#' @return Never returns.
#' @keywords internal
#' @noRd
abort_newborn_rule <- function(problem) {
  micsim_abort("micsim_newborn_rule_error", problem)
}


#' Print a newborn rule
#'
#' Shows which state variables a child inherits and which are drawn, with
#' the candidate table small enough to read. The most likely
#' mistake is inheriting the wrong column, which is visible here at a
#' glance.
#'
#' @param x A `newborn_rule`.
#' @param ... Not used; present for compatibility with the [print()] generic.
#' @return `x`, invisibly. Called for the side effect of printing.
#' @exportS3Method print newborn_rule
#' @noRd
print.newborn_rule <- function(x, ...) {
  cat("<newborn_rule>\n")
  cat("  state variables: ", x$n_vars, "\n", sep = "")
  cat("  inherited from the mother: ",
      if (length(x$inherit) == 0L) {
        "none"
      } else {
        paste0("position(s) ", toString(x$inherit))
      },
      "\n", sep = "")
  cat("  candidate states: ", nrow(x$states), "\n", sep = "")
  shown <- seq_len(min(8L, nrow(x$states)))
  for (i in shown) {
    cat("    ", paste(x$states[i, ], collapse = "/"),
        "   p = ", signif(x$prob[i], 4), "\n", sep = "")
  }
  if (nrow(x$states) > length(shown)) {
    cat("    ... and ", nrow(x$states) - length(shown), " more\n", sep = "")
  }
  invisible(x)
}
