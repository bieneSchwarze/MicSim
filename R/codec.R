# =====================================================================
# codec.R  --  State codec: labels <-> integer indices
# ---------------------------------------------------------------------
# WHAT THIS FILE DOES
#   Gives every state a unique integer index and translates both ways.
#   The engine then moves integers instead of text labels such as
#   "f/AT/1". This keeps repeated state comparisons simple; the mapping
#   relies on the calling interface supplying distinct transient states.
#
# PLACE IN THE PIPELINE
#   User-defined state variables --> state_space() --> micsim_model()
#       --> new_state_codec() [this file] --> registry and engine lookups
#
# CONCEPTS  (all defined in ?MicSim-concepts)
#   composite state, state variable, absorbing state
# =====================================================================

#' Assign numbers to state labels
#'
#' @description
#' Builds the one-to-one map between composite state labels and integer
#' indices that the rest of the engine works with. `encode()` turns a
#' label such as `"f/AT/1"` into its index; `decode()` turns an index back
#' into its label.
#'
#' @details
#' Each transient state combines one value from every state variable.
#' For example, two sex categories, three nationalities and two parity
#' values give twelve combinations.
#'
#' Values are joined with `"/"` in the order the variables were supplied.
#' With sex first, nationality second and parity third, `"f/AT/1"` could
#' represent a woman of Austrian nationality with one child.
#'
#' The order of the variables and their values determines the numbering.
#' Keep this order consistent throughout the model. Absorbing states,
#' such as death or emigration, are added after the transient states.
#'
#' Absorbing states are identified by the `absorbing` argument, not by
#' their names. Their labels must be distinct and must not duplicate a
#' transient-state label. State-variable values cannot contain `"/"`,
#' because it separates the parts of a label.
#'
#' @param state_space Named list of state variables and their values, e.g.
#'   `list(sex = c("m", "f"), nat = c("DE", "AT", "IT"))`.
#' @param absorbing Character vector of absorbing state labels, e.g.
#'   `c("dead", "abroad")`. Empty by default.
#' @return A `state_codec`: a list holding `encode()`, `decode()`,
#'   `is_absorbing()`, the state vectors they operate on, and `grid`, a
#'   data frame with one row per transient state (in index order) and one
#'   character column per state variable, which is what role conditions
#'   such as `sex == "f"` are evaluated against.
#' @keywords internal
#' @noRd
new_state_codec <- function(state_space, absorbing = character()) {
  if (!is.list(state_space) || length(state_space) < 1L) {
    abort_codec("The state space must be a list of at least one state variable.")
  }
  state_space <- lapply(state_space, as.character)
  absorbing <- as.character(absorbing)

  values <- unlist(state_space, use.names = FALSE)
  if (anyNA(values) || any(!nzchar(values))) {
    abort_codec("State values must be non-empty strings, not NA or \"\".")
  }
  bad <- unique(values[grepl(STATE_SEP, values, fixed = TRUE)])
  if (length(bad) > 0L) {
    abort_codec(paste0(
      "State values may not contain \"", STATE_SEP, "\", which separates the\n",
      "  components of a composite state. Offending: ", toString(bad), "."
    ))
  }

  # Enumerate the transient state space = Cartesian product of substates,
  # joined by "/" in the order the variables were declared (the contract).
  grid <- expand.grid(rev(state_space), stringsAsFactors = FALSE,
                      KEEP.OUT.ATTRS = FALSE)
  grid <- grid[, rev(seq_along(state_space)), drop = FALSE]
  if (!is.null(names(state_space))) names(grid) <- names(state_space)
  row.names(grid) <- NULL
  transient <- do.call(paste, c(unname(as.list(grid)), sep = STATE_SEP))

  if (anyNA(absorbing) || any(!nzchar(absorbing)) ||
        anyDuplicated(absorbing) > 0L) {
    abort_codec("Absorbing state labels must be distinct, non-empty strings.")
  }
  clash <- intersect(absorbing, transient)
  if (length(clash) > 0L) {
    abort_codec(paste0(
      "Absorbing state label(s) ", toString(clash), " coincide with a\n",
      "  transient state, so a transition into them would be ambiguous."
    ))
  }

  states <- c(transient, absorbing)          # transient first, then absorbing
  K <- length(transient)
  is_abs <- c(rep(FALSE, K), rep(TRUE, length(absorbing)))

  # match() is an exact, one-to-one lookup for any number of states.
  encode <- function(label) {
    idx <- match(label, states)
    if (anyNA(idx)) {
      abort_codec(paste0(
        "Unknown state label(s): ", toString(unique(label[is.na(idx)])), ".\n",
        "  Known states are the combinations of ", toString(names(state_space)),
        if (length(absorbing)) paste0(", and ", toString(absorbing)) else "",
        "."
      ))
    }
    idx
  }
  decode <- function(index) states[index]
  is_absorbing <- function(index) is_abs[index]

  structure(
    list(encode = encode, decode = decode, is_absorbing = is_absorbing,
         states = states, transient = transient, absorbing = absorbing,
         K = K, n = length(states),
         vars = names(state_space), levels = state_space, grid = grid),
    class = "state_codec"
  )
}


#' Raise a codec error
#' @param problem The message.
#' @return Never returns.
#' @keywords internal
#' @noRd
abort_codec <- function(problem) {
  micsim_abort("micsim_codec_error", problem)
}

#' Print an overview of the state mapping
#'
#' Lists the state variables, the number of transient states formed by
#' combining their values, and the absorbing-state labels. Also shows
#' the total number of states and the range of numbers assigned to them.
#' Useful for checking that a model's state space is the one you meant
#' to write.
#'
#' @param x A `state_codec`, from `new_state_codec()`.
#' @param ... Not used; present for compatibility with the [print()] generic.
#' @return `x`, invisibly. Called for the side effect of printing.
#' @exportS3Method print state_codec
#' @noRd
print.state_codec <- function(x, ...) {
  cat("<state_codec>\n")
  cat("  variables : ", paste(x$vars, collapse = ", "), "\n", sep = "")
  cat("  transient : ", x$K, " states\n", sep = "")
  cat("  absorbing : ", paste(x$absorbing, collapse = ", "), "\n", sep = "")
  cat("  total     : ", x$n, " states (indices 1..", x$n, ")\n", sep = "")
  invisible(x)
}
