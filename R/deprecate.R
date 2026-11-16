# =====================================================================
# deprecate.R
# ---------------------------------------------------------------------
# WHAT THIS FILE DOES
#   Tells users of the classic (MicSim 3.0.0) interface that a function is
#   deprecated and what replaces it.
#
# THE DEPRECATION CYCLE
#   The classic entry points stay exported and keep working, so existing
#   scripts run unchanged. They follow R's usual deprecation cycle:
#
#     4.0.0   deprecated: the function works and shows a message, once per
#             R session, naming its replacement (this file).
#     later   the message becomes a warning, then the function stops with
#             an error naming its replacement (.Defunct()), then it is
#             removed.
#
#   The deprecated functions are thin wrappers. Their bodies live in
#   internal functions (legacy_micsim(), convert_to_long(),
#   convert_to_wide()), and the new API calls those directly, so
#   as_long(), as_wide() and the workers of micSimParallel() never show
#   the message.
#
# GUARANTEES
#   * Each deprecated function shows its message at most once per session.
#   * options(MicSim.show_deprecation = FALSE) turns the messages off.
#   * The message is a classed condition, so it can be muffled precisely
#     with suppressMessages(..., classes = "micsim_deprecation_message").
# =====================================================================


#' Replacements of the deprecated functions
#'
#' Named by the deprecated function; each entry says what to use instead.
#' @keywords internal
#' @noRd
DEPRECATED_REPLACEMENTS <- c(
  micSim = paste0(
    "Build the model with micsim_model() (and add_subsystem() for fertility\n",
    "  or school enrolment) and run it with simulate();\n",
    "  see vignette(\"introduction-to-new-micsim\", package = \"MicSim\")."
  ),
  micSimParallel = paste0(
    "Build the model with micsim_model() and run it with\n",
    "  simulate(model, seed = ..., cores = ...)."
  ),
  buildTransitionMatrix = paste0(
    "Declare the transitions with transition(), move() and transitions(),\n",
    "  and pass them to micsim_model()."
  ),
  convertToLongFormat = "Use as_long() on the result of simulate().",
  convertToWideFormat = "Use as_wide() on the result of simulate().",
  getInDays_my = "Use getFirstDayOfMonth(), which takes the same arguments."
)


#' Which deprecation messages have been shown in this session
#' @keywords internal
#' @noRd
.deprecation_shown <- new.env(parent = emptyenv())


#' Tell the user that a function is deprecated
#'
#' @description
#' Shows a message, at most once per R session for each function, saying
#' that `old()` is deprecated and what replaces it. Nothing is shown when
#' the option `MicSim.show_deprecation` is `FALSE`.
#'
#' @param old Name of the deprecated function, a key of
#'   `DEPRECATED_REPLACEMENTS`.
#' @return `invisible(TRUE)` if the message was shown, otherwise
#'   `invisible(FALSE)`.
#' @keywords internal
#' @noRd
deprecate_message <- function(old) {
  if (isFALSE(getOption("MicSim.show_deprecation", TRUE))) {
    return(invisible(FALSE))
  }
  if (isTRUE(get0(old, envir = .deprecation_shown, inherits = FALSE))) {
    return(invisible(FALSE))
  }
  assign(old, TRUE, envir = .deprecation_shown)
  text <- paste0(
    "`", old, "()` is deprecated as of MicSim 4.0.0 and will be removed in a\n",
    "  future release. ", DEPRECATED_REPLACEMENTS[[old]], "\n",
    "  This message is shown once per session; set\n",
    "  options(MicSim.show_deprecation = FALSE) to turn it off.\n"
  )
  message(structure(
    class = c("micsim_deprecation_message", "deprecatedMessage", "message",
              "condition"),
    list(message = text, call = NULL, old = old)
  ))
  invisible(TRUE)
}


#' Forget which deprecation messages have been shown
#'
#' Used by the tests, so that each can see the message again.
#' @return `invisible(NULL)`.
#' @keywords internal
#' @noRd
reset_deprecation_messages <- function() {
  rm(list = ls(envir = .deprecation_shown), envir = .deprecation_shown)
  invisible(NULL)
}
