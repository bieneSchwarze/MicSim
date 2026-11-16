# =====================================================================
# selectors.R
# ---------------------------------------------------------------------
# WHAT THIS FILE DOES
#   Selects transitions by state-variable names and values, so users do
#   not have to list every composite state label, such as "female/0".
#
# PLACE IN THE PIPELINE
#   the user's move("edu", "none", "primary") -->  selectors.R
#       --(resolved ONCE, before the run, into integer state indices)-->
#   registry.R (declaring transitions)  and  subsystem_*.R (claiming them)
#
# WHAT YOU NEED TO KNOW  (concepts in ?MicSim-concepts)
#   composite state, state variable.
#
#   move(var, from, to) describes a change in one variable, keeping the
#   others fixed. It can declare transitions or select existing ones.
#
#   moves_of(var) selects existing transitions that change only var,
#   in either direction. It cannot declare new transitions.
#
#   For both [move()] and [moves_of()] where is an additional argument that
#   restricts the origin states using a condition such as sex == "female". 
#   A captured condition is called a role internally.
#
# GUARANTEES
#   * Resolving selectors consumes no random numbers.
#   * A selector that matches no states or declared transitions raises
#     an error. Each explicit transition label must also exist.
# =====================================================================


#' Select the transitions that change one state variable
#'
#' @description
#' `move()` describes a change from one value of a state variable to another,
#' keeping every other variable fixed. Use it in [transition()] to declare
#' transitions, or in a subsystem to select which transitions it acts on,
#' as in `school_enrolment(move("edu", "none", "primary"))`.
#'
#' `moves_of()` selects declared transitions in which `var` is the only
#' state variable that changes. Subsystems use it to identify relevant
#' transitions, for example `fertility(births = moves_of("parity"))`.
#' It does not declare new transitions.
#'
#' @details
#' Suppose the state variables are `sex = c("male", "female")` and
#' `edu = c("none", "primary")`. Then `move("edu", "none", "primary")`
#' describes both `"male/none" -> "male/primary"` and
#' `"female/none" -> "female/primary"`.
#'
#' `where` restricts the origin states. For example,
#' `move("parity", "0", "1", where = sex == "female")` selects that parity
#' change only for female origin states. The condition is evaluated against
#' the model's table of transient states when the selector is resolved.
#'
#' @param var Name of the state variable that changes, as a string.
#' @param from,to Its origin and destination values, as different strings.
#' @param where Optional condition on the origin state, written in terms of
#'   the state variables, e.g. `sex == "female" & nat == "DE"`.
#'
#' @return A selector object. `move()` can be used in [transition()] or by
#'   a subsystem; `moves_of()` can only select declared transitions.
#'
#' @seealso [transition()], [add_subsystem()], [MicSim-concepts]
#'
#' @examples
#' # Entry into primary education, regardless of other state variables.
#' move("edu", from = "none", to = "primary")
#'
#' # The change from parity 0 to 1 applies only to female origin states.
#' move("parity", from = "0", to = "1", where = sex == "female")
#'
#' # Treat all declared changes in parity alone as births.
#' fertility(
#'   births = moves_of("parity"),
#'   newborn = newborn_rule(sex = c(female = 0.5, male = 0.5),
#'                          parity = "0")
#' )
#'
#' @export
move <- function(var, from, to, where = NULL) {
  check_selector_string(var, "var")
  check_selector_string(from, "from")
  check_selector_string(to, "to")
  if (identical(from, to)) {
    abort_selector(paste0(
      "move(\"", var, "\", \"", from, "\", \"", to, "\") changes nothing:\n",
      "  `from` and `to` are the same value."
    ))
  }
  structure(
    list(var = var, from = from, to = to,
         where = new_role(substitute(where), parent.frame())),
    class = c("micsim_move", "micsim_selector")
  )
}


#' @rdname move
#' @export
moves_of <- function(var, where = NULL) {
  check_selector_string(var, "var")
  structure(
    list(var = var, where = new_role(substitute(where), parent.frame())),
    class = c("micsim_moves_of", "micsim_selector")
  )
}


#' Print a selector
#' @param x A selector from [move()] or [moves_of()].
#' @param ... Not used.
#' @return `x`, invisibly.
#' @exportS3Method print micsim_selector
#' @noRd
print.micsim_selector <- function(x, ...) {
  cat("<micsim_selector> ", describe_selector(x), "\n", sep = "")
  invisible(x)
}


#' Describe a selector for printing or error messages
#' @param x A selector from [move()] or [moves_of()].
#' @return One string describing the variable, change and any condition.
#' @keywords internal
#' @noRd
describe_selector <- function(x) {
  where <- if (is.null(x$where)) "" else
    paste0(", where ", deparse1(x$where$expr))
  if (inherits(x, "micsim_move")) {
    paste0("move(", x$var, ": ", x$from, " -> ", x$to, where, ")")
  } else {
    paste0("moves_of(", x$var, where, ")")
  }
}


#' Capture a condition for later evaluation against states
#'
#' Stores an expression such as `sex == "female"` with the environment where
#' it was written. When evaluated, state-variable names come from the model;
#' other names can refer to objects in that environment.
#'
#' @param expr An unevaluated expression, or `NULL` for no restriction.
#' @param env Environment used to look up names absent from the state table.
#' @return A `micsim_role` list with `expr` and `env`, or `NULL`.
#' @keywords internal
#' @noRd
new_role <- function(expr, env) {
  if (is.null(expr)) return(NULL)
  structure(list(expr = expr, env = env), class = "micsim_role")
}


#' Test which transient states satisfy a condition
#'
#' @description
#' Evaluates a condition such as `sex == "female"` against the table of
#' transient states, once, and returns the matching state indices.
#'
#' @param codec A `state_codec`.
#' @param role A `micsim_role` from `new_role()`, or `NULL` for "every
#'   transient state".
#' @param what A short description of where the role came from, for the
#'   error message.
#' @return Logical vector over the transient states, in index order.
#' @keywords internal
#' @noRd
states_where <- function(codec, role, what = "condition") {
  if (is.null(role)) return(rep(TRUE, codec$K))
  hit <- tryCatch(eval(role$expr, codec$grid, role$env),
                  error = function(e) e)
  label <- deparse1(role$expr)
  if (inherits(hit, "condition")) {
    abort_selector(paste0(
      "The ", what, " `", label, "` could not be evaluated: ",
      conditionMessage(hit), "\n",
      "  It is evaluated against the state variables ",
      toString(codec$vars), "."
    ))
  }
  if (is.logical(hit) && length(hit) == 1L) hit <- rep(hit, codec$K)
  if (!is.logical(hit) || length(hit) != codec$K || anyNA(hit)) {
    abort_selector(paste0(
      "The ", what, " `", label, "` is not a yes/no condition on the\n",
      "  state variables ", toString(codec$vars), "."
    ))
  }
  if (!any(hit)) {
    abort_selector(paste0(
      "The ", what, " `", label, "` matches no state.\n",
      describe_values(codec)
    ))
  }
  hit
}


#' Resolve a move() selector into state-label pairs
#' 
#' @description
#' Checks the variable and values, then selects origin states with the
#' specified `from` value that satisfy `where`. Each destination is the
#' same state with that variable replaced by `to`. No matches is an error.
#' 
#' @param sel A `micsim_move` selector.
#' @param codec A `state_codec`.
#' @return A two-column character matrix, `from` and `to`, one row per
#'   origin state, in the codec's state order.
#' @keywords internal
#' @noRd
resolve_move <- function(sel, codec) {
  pos <- match(sel$var, codec$vars)
  if (is.na(pos)) {
    abort_selector(paste0(
      describe_selector(sel), " names the state variable `", sel$var,
      "`, which the model does not have.\n",
      "  State variables: ", toString(codec$vars), "."
    ))
  }
  levels <- codec$levels[[pos]]
  unknown <- setdiff(c(sel$from, sel$to), levels)
  if (length(unknown) > 0L) {
    abort_selector(paste0(
      describe_selector(sel), " uses ", toString(unknown),
      ", which `", sel$var, "` does not take.\n",
      "  Values of `", sel$var, "`: ", toString(levels), "."
    ))
  }

  keep <- codec$grid[[pos]] == sel$from &
    states_where(codec, sel$where, paste0("`where` of ", describe_selector(sel)))
  origins <- which(keep)
  if (length(origins) == 0L) {
    abort_selector(paste0(describe_selector(sel), " matches no state."))
  }

  destination <- codec$grid[origins, , drop = FALSE]
  destination[[pos]] <- sel$to
  to_labels <- do.call(paste, c(unname(as.list(destination)), sep = STATE_SEP))

  cbind(from = codec$transient[origins], to = to_labels)
}


#' Find the declared transitions selected by a subsystem
#'
#' @description
#' Converts selectors or transition labels into positions in the registry's
#' edge list. This selection is called a claim.
#'
#' @details
#' A claim selects among transitions the model declares. A `move()` that
#' also covers pairs the model does not declare is fine (they simply cannot
#' happen); a claim that matches no declared transition at all is an
#' error, because it would install a subsystem that never fires.
#'
#' Plain labels are taken literally, and each one must be a declared
#' transition.
#'
#' @param claim The claim.
#' @param codec A `state_codec`.
#' @param registry A `rate_registry`.
#' @param what A short description, for error messages, e.g.
#'   `"fertility(births = )"`.
#' @return Integer vector of edge positions (`edge$id`), sorted, unique.
#' @keywords internal
#' @noRd
resolve_claim <- function(claim, codec, registry, what) {

  pieces <- if (inherits(claim, "micsim_selector") || is.character(claim)) {
    list(claim)
  } else if (is.list(claim) && length(claim) > 0L) {
    claim
  } else {
    abort_selector(paste0(
      what, " must be a move() or moves_of() selector, a list of them,\n",
      "  or a character vector of transition labels such as \"f/0->f/1\"."
    ))
  }

  from_idx <- vapply(registry$edges, `[[`, integer(1L), "from_idx")
  to_idx   <- vapply(registry$edges, `[[`, integer(1L), "to_idx")
  ids <- integer(0)

  for (piece in pieces) {
    if (inherits(piece, "micsim_move")) {
      pairs <- resolve_move(piece, codec)
      key_want <- paste(codec$encode(pairs[, "from"]),
                        codec$encode(pairs[, "to"]))
      found <- which(paste(from_idx, to_idx) %in% key_want)
      if (length(found) == 0L) {
        abort_selector(paste0(
          what, ": ", describe_selector(piece), " matches no transition the\n",
          "  model declares. Declare the transitions with transition() first."
        ))
      }
      ids <- c(ids, found)

    } else if (inherits(piece, "micsim_moves_of")) {
      pos <- match(piece$var, codec$vars)
      if (is.na(pos)) {
        abort_selector(paste0(
          what, ": ", describe_selector(piece), " names the state variable `",
          piece$var, "`, which the model does not have.\n",
          "  State variables: ", toString(codec$vars), "."
        ))
      }
      origin_ok <- states_where(codec, piece$where,
                                paste0("`where` of ", describe_selector(piece)))
      found <- which(vapply(registry$edges, function(e) {
        identical(e$changed_vars, pos) && !codec$is_absorbing(e$to_idx) &&
          origin_ok[e$from_idx]
      }, logical(1L)))
      if (length(found) == 0L) {
        abort_selector(paste0(
          what, ": ", describe_selector(piece), " matches no declared\n",
          "  transition, i.e. the model has no transition that changes `",
          piece$var, "` alone."
        ))
      }
      ids <- c(ids, found)

    } else if (is.character(piece)) {
      labels <- as_edge_labels(piece)
      known <- vapply(registry$edges, `[[`, character(1L), "label")
      missing_labels <- setdiff(labels, known)
      if (length(missing_labels) > 0L) {
        abort_selector(paste0(
          what, " names transition(s) the model does not declare: ",
          toString(missing_labels), "."
        ))
      }
      ids <- c(ids, match(labels, known))

    } else {
      abort_selector(paste0(
        what, " contains a ", class(piece)[1L], "; expected move(),\n",
        "  moves_of() or a transition label."
      ))
    }
  }

  sort(unique(ids))
}


#' Normalise transition labels to the registry's form
#'
#' @description
#' Accepts `"A->B"`, `"A -> B"` or anything in between, and 
#' returns `"A -> B"`.
#'
#' @param labels Character vector of transition labels.
#' @return Character vector in `"from -> to"` form.
#' @keywords internal
#' @noRd
as_edge_labels <- function(labels) {
  if (length(labels) == 0L) return(character(0))
  parts <- strsplit(as.character(labels), "->", fixed = TRUE)
  vapply(parts, function(p) {
    if (length(p) != 2L) {
      abort_selector(
        paste0("`", paste(p, collapse = "->"),
               "` is not a transition label of the form 'from->to'.")
      )
    }
    paste(trimws(p[1L]), "->", trimws(p[2L]))
  }, character(1L))
}


#' List the values of every state variable, for an error message
#' @param codec A `state_codec`.
#' @return A string.
#' @keywords internal
#' @noRd
describe_values <- function(codec) {
  paste0("  Values present: ",
         paste(vapply(seq_along(codec$vars), function(i) {
           paste0(codec$vars[i], " = ", toString(codec$levels[[i]]))
         }, character(1L)), collapse = "; "), ".")
}


#' Check that a selector argument is one non-empty string
#' @param x The argument to check.
#' @param name Argument name for error messages.
#' @return `invisible(TRUE)`.
#' @keywords internal
#' @noRd
check_selector_string <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    abort_selector(paste0("`", name, "` must be a single, non-empty string."))
  }
  invisible(TRUE)
}


#' Raise a selector error
#' @param problem The message.
#' @return Never returns.
#' @keywords internal
#' @noRd
abort_selector <- function(problem) {
  micsim_abort("micsim_selector_error", problem)
}
