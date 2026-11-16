# =====================================================================
# engine_hooks.R
# ---------------------------------------------------------------------
# WHAT THIS FILE DOES
#   Defines the fixed points in the simulation loop at which an extension
#   may observe or adjust what the engine is doing, without the engine
#   knowing anything about that extension.
#
# PLACE IN THE PIPELINE
#   engine.R --(at six named points)--> engine_hooks.R --> a subsystem
#
#   Subsystem instantiation creates hook sets. engine.R combines them and
#   calls them during scheduling, after transitions and on absorption.
#   The fertility subsystem also calls on_birth after adding a child.
#
# WHAT YOU NEED TO KNOW
#   Observing hooks receive a context list; their return values are
#   ignored. They can still change mutable stores supplied in that list.
#   Adjusting hooks receive a value and context, and return the value
#   the engine should use. An absent hook is NULL and adds no behaviour.
#
# GUARANTEES
#   * A hook set is inert by default: no hook, no behaviour change.
#   * OBSERVING hooks cannot alter the simulation. Their return value is
#     discarded, so a mistake in one cannot silently change results.
#   * ADJUSTING hooks are exactly three, and each has a stated contract for
#     what it receives and what it must return.
# =====================================================================


# ---------------------------------------------------------------------
# The six hook points, and the two kinds of hook
# ---------------------------------------------------------------------
# OBSERVING hooks are called for their effect on the subsystem's own
# state. Whatever they return is thrown away.
#
#   on_event   fires after a transition is recorded and the person's state
#              is updated. For example, fertility uses it to add a newborn.
#
#   on_birth   is fired by fertility after registering the child and
#              before drawing the child's first event. Other callbacks
#              can read the child and mother slots from the context.
#
#   on_absorb  fires after on_event when the destination is absorbing.
#
# ADJUSTING hooks are given a value and must return one, of the same
# shape. They are the only three places an extension can change what the
# engine does.
#
#   modify_hazard  given the hazard that bind_hazard() built for one
#                  person and one transition (a list with lambda,
#                  forced_time and crossings), returns a list of the same
#                  shape. It runs just before the waiting time is drawn,
#                  so a subsystem can change one individual's rate, e.g.
#                  multiply it by a factor that depends on that person's
#                  attributes or on their mother's state
#                  (intergenerational transmission). No shipped subsystem
#                  uses it yet.
#
#   on_candidates  given every destination's drawn time, as a list
#                  list(time_days, to_idx, edge_id) with one element per
#                  destination (Inf for "no event"), returns the same
#                  shape. The earliest is chosen afterwards, so a calendar
#                  rule that moves one candidate lets the others compete
#                  over the move. No time may be moved before now. The
#                  school-enrolment subsystem uses it.
#
#   on_schedule    given the chosen event as list(time_days, to_idx,
#                  edge_id), returns an event of the same shape, or NULL
#                  to cancel the scheduling entirely. Reserved.
#
# Every hook receives one argument, `ctx` (context), a named list. Its fields
# are documented on new_hook_context() below. The contract is deliberately
# narrow: a hook gets the engine's collaborators and the facts of the
# current event, and nothing else.
# ---------------------------------------------------------------------

#' The names of the observing hooks
#' @keywords internal
#' @noRd
OBSERVING_HOOKS <- c("on_event", "on_birth", "on_absorb")

#' The names of the adjusting hooks
#' @keywords internal
#' @noRd
ADJUSTING_HOOKS <- c("modify_hazard", "on_candidates", "on_schedule")


#' Create a set of engine hooks
#'
#' @description
#' Collects the callbacks an extension wants the simulation driver to make.
#' Called with no arguments it produces an empty set. Subsystems install
#' callbacks only at the points they need.
#'
#' @details
#' A hook left as `NULL` is absent. `fire_hook()` and `adjust_by_hook()`
#' check for it before using `ctx`, so R's lazy evaluation avoids building
#' a context for an absent callback. `hook_active()` also lets the driver
#' skip collecting candidates when no subsystem needs them.
#'
#' @param on_event,on_birth,on_absorb Observing hooks, or `NULL`. Each is
#'   a function of one argument, `ctx`. Their return values are discarded.
#' @param modify_hazard,on_candidates,on_schedule Adjusting hooks, or
#'   `NULL`. Each is a function of `(value, ctx)` returning a value of 
#'   the same shape as the one it was given; see the contracts at the 
#'   top of this file.
#'
#' @return An object of class `engine_hooks`.
#'
#' @seealso `new_hook_context()`, `fire_hook()`, `adjust_by_hook()`
#'
#' @examples
#' hooks <- new_engine_hooks()
#' hook_active(hooks, "on_birth")     # FALSE: core MicSim sets nothing
#'
#' log <- new.env(parent = emptyenv())
#' log$births <- 0L
#' hooks <- new_engine_hooks(
#'   on_birth = function(ctx) log$births <- log$births + 1L
#' )
#' hook_active(hooks, "on_birth")     # TRUE
#'
#' @keywords internal
#' @noRd
new_engine_hooks <- function(on_event      = NULL,
                             on_birth      = NULL,
                             on_absorb     = NULL,
                             modify_hazard = NULL,
                             on_candidates = NULL,
                             on_schedule   = NULL) {

  hooks <- list(
    on_event      = on_event,
    on_birth      = on_birth,
    on_absorb     = on_absorb,
    modify_hazard = modify_hazard,
    on_candidates = on_candidates,
    on_schedule   = on_schedule
  )

  for (name in names(hooks)) {
    if (!is.null(hooks[[name]]) && !is.function(hooks[[name]])) {
      micsim_abort(
        "micsim_hook_error",
        paste0("Hook `", name, "` must be a function or NULL, not ",
               class(hooks[[name]])[1L], ".")
      )
    }
  }

  structure(
    c(hooks, list(active = names(hooks)[!vapply(hooks, is.null, logical(1L))])),
    class = "engine_hooks"
  )
}


#' Is a hook set?
#'
#' The driver's guard. Kept as a named function rather than an inline
#' `!is.null(...)` so that every hook point in `engine.R` reads as a hook
#' point rather than as a null check.
#'
#' @param hooks An `engine_hooks` object.
#' @param which Name of the hook.
#' @return `TRUE` or `FALSE`.
#' @keywords internal
#' @noRd
hook_active <- function(hooks, which) {
  !is.null(hooks[[which]])
}


#' Build the context passed to a hook
#'
#' @description
#' Packages what a hook is allowed to see. Every field is optional, so one
#' constructor serves all six hook points; each point fills in what is
#' meaningful there.
#'
#' @details
#' The fields, and who uses them:
#'
#' \describe{
#'   \item{`persons`, `queue`, `store`}{the engine's mutable stores.
#'     Fertility uses `persons$next_id()` and `persons$add()` to register
#'     a child. Callbacks with a queue reference can cancel or replace
#'     pending events. Not every hook point supplies all three stores.}
#'   \item{`schedule`}{`schedule(slot, now_days)`: draw this individual's
#'     next event and put it in the queue, exactly as the driver does.
#'     This is the primitive that lets a subsystem bring a NEW individual
#'     into a running simulation. Without it, fertility could not be
#'     written as a subsystem at all, and the seam would be adequate for
#'     observing the simulation but not for extending it.}
#'   \item{`hooks`}{the full hook set, so that a subsystem can fire a hook
#'     of its own. The fertility subsystem uses it to fire `on_birth`,
#'     so other installed callbacks can respond after the child is added.}
#'   \item{`codec`, `registry`}{read-only lookups, for a subsystem that
#'     needs to reason about states or rates.}
#'   \item{`now_days`}{the global clock at the moment the hook fires.}
#'   \item{`slot`}{the individual the event concerns.}
#'   \item{`from_idx`, `to_idx`}{the states left and entered, where the
#'     hook point has both.}
#'   \item{`edge`}{the registry edge being solved, for `modify_hazard`.}
#'   \item{`mother_slot`}{for `on_birth`.}
#' }
#'
#' A hook must treat everything here as borrowed, but not as read-only: a
#' subsystem that adds an individual or cancels an event is doing exactly
#' what the seam exists for. What it must not do is assume it is alone.
#' When several subsystems are installed they fire in the order they were
#' combined, and a subsystem that depends on another's writes is depending
#' on that order; `combine_hooks()` documents what the order is.
#'
#' @param ... Named fields as listed above.
#' @return A named list.
#' @keywords internal
#' @noRd
new_hook_context <- function(...) {
  list(...)
}


#' Install several subsystems side by side
#'
#' @description
#' Merges the hook sets of installed subsystems and any raw engine hooks
#' into the single set used for a run.
#'
#' @details
#' The two kinds of hook merge differently, and the difference follows
#' from what each kind is for.
#'
#' **Observing hooks are chained.** Every set's hook is called, in the
#' order the sets were given, and every return value is discarded as
#' usual. None of them can influence another, so the order affects only
#' the sequence of side effects on the subsystems' own state.
#'
#' **Adjusting hooks are composed.** The first set's hook is given the
#' engine's value, the second is given the first's output, and so on, so
#' the last set installed has the final word. 
#'
#' **Core goes first.** The engine installs its own subsystems before the
#' caller's, so a user hook observes a world in which the built-in
#' behaviour has already happened and gets the last word on any value it
#' adjusts. A user who needs to run *before* a built-in subsystem is
#' asking to override core behaviour, which this seam deliberately does
#' not offer.
#'
#' @param ... Two or more `engine_hooks` objects, in the order they should
#'   fire. `NULL` entries are ignored, so a caller can pass an optional
#'   subsystem without writing a branch around it.
#'
#' @return A single `engine_hooks` object.
#'
#' @examples
#' counter <- new.env(parent = emptyenv())
#' counter$n <- 0L
#' a <- new_engine_hooks(on_event = function(ctx) counter$n <- counter$n + 1L)
#' b <- new_engine_hooks(on_event = function(ctx) counter$n <- counter$n + 10L)
#' fire_hook(combine_hooks(a, b), "on_event", new_hook_context())
#' counter$n   # 11: both fired, a first
#'
#' @keywords internal
#' @noRd
combine_hooks <- function(...) {

  sets <- Filter(Negate(is.null), list(...))
  for (set in sets) {
    if (!inherits(set, "engine_hooks")) {
      micsim_abort(
        "micsim_hook_error",
        "combine_hooks() takes engine_hooks objects, or NULL."
      )
    }
  }
  # A set with nothing in it contributes nothing, and dropping it keeps
  # the common case -- one subsystem, no extra hooks -- free of a wrapper
  # closure on every call.
  if (length(sets) > 1L) {
    sets <- Filter(function(set) length(set$active) > 0L, sets)
  }
  if (length(sets) == 0L) return(new_engine_hooks())
  if (length(sets) == 1L) return(sets[[1L]])

  merged <- list()
  for (which in OBSERVING_HOOKS) {
    merged[[which]] <- chain_observers(hooks_named(sets, which))
  }
  for (which in ADJUSTING_HOOKS) {
    merged[[which]] <- compose_adjusters(hooks_named(sets, which))
  }

  do.call(new_engine_hooks, merged)
}


#' Collect the hook of one name from several hook sets
#'
#' @param sets List of `engine_hooks` objects.
#' @param which Name of the hook.
#' @return List of functions, possibly empty, in the order the sets were
#'   given.
#' @keywords internal
#' @noRd
hooks_named <- function(sets, which) {
  Filter(Negate(is.null), lapply(sets, function(s) s[[which]]))
}


#' Chain observing hooks into one
#'
#' @description
#' Returns a single hook that calls each of `present` in turn and
#' discards every result.
#'
#' @details
#' Each returned callback needs its own saved hook list. A closure
#' created directly inside a loop would refer to the loop's current
#' binding and could later call the wrong list. Passing `present` into
#' this helper and forcing it gives each callback a separate binding.
#'
#' @param present List of observing hooks, in firing order.
#' @return A single hook function, or `NULL` when `present` is empty.
#' @keywords internal
#' @noRd
chain_observers <- function(present) {
  force(present)
  if (length(present) == 0L) return(NULL)
  function(ctx) {
    for (hook in present) hook(ctx)
    invisible(NULL)
  }
}


#' Compose adjusting hooks into one
#'
#' Each hook is handed the previous one's output, so the last one
#' installed has the final word.
#' 
#' @param present List of adjusting hooks, in application order.
#' @return A single hook function, or `NULL` when `present` is empty.
#' @keywords internal
#' @noRd
compose_adjusters <- function(present) {
  force(present)
  if (length(present) == 0L) return(NULL)
  function(value, ctx) {
    for (hook in present) value <- hook(value, ctx)
    value
  }
}


#' Call an observing hook
#'
#' @description
#' Runs a hook for its effect and discards whatever it returns, so that an
#' observer cannot change the simulation by accident.
#'
#' @details
#' The `ctx` argument is promised, not evaluated, until something forces
#' it: R's lazy evaluation means a call such as
#'
#' ```
#' fire_hook(hooks, "on_event", new_hook_context(persons = p, ...))
#' ```
#'
#' builds no context at all when `on_event` is `NULL`. A run of core
#' MicSim therefore pays one `is.null()` test per hook point per event and
#' nothing more.
#'
#' The laziness goes one step further than that, which is worth knowing
#' when reading a test: passing the promise on to `hook(ctx)` does not
#' force it either. A hook that never mentions its argument never causes
#' the context to be built. Any hook that actually uses `ctx` forces it
#' immediately, at the moment of the call, so it always sees the values
#' the driver held when it fired and never a later iteration's.
#'
#' @param hooks An `engine_hooks` object.
#' @param which One of `OBSERVING_HOOKS`.
#' @param ctx The context list, typically built inline.
#' @return `invisible(NULL)`.
#' @keywords internal
#' @noRd
fire_hook <- function(hooks, which, ctx) {
  hook <- hooks[[which]]
  if (is.null(hook)) return(invisible(NULL))
  hook(ctx)
  invisible(NULL)
}


#' Pass a value through an adjusting hook
#'
#' @description
#' Gives a hook the chance to replace a value the engine is about to use.
#' When the hook is not set the value is returned untouched.
#'
#' @details
#' The returned value is not validated here, because what counts as the
#' right shape differs between the adjusting hooks and checking it on
#' every event would cost more than it catches. The contracts are stated
#' at the top of this file, and the driver checks whatever it must at the
#' point of use.
#'
#' @param hooks An `engine_hooks` object.
#' @param which One of `ADJUSTING_HOOKS`.
#' @param value The value the engine would otherwise use.
#' @param ctx The context list.
#' @return The hook's return value, or `value` when no hook is set.
#' @keywords internal
#' @noRd
adjust_by_hook <- function(hooks, which, value, ctx) {
  hook <- hooks[[which]]
  if (is.null(hook)) return(value)
  hook(value, ctx)
}


#' Print a hook set
#'
#' @param x An `engine_hooks` object.
#' @param ... Not used; present for compatibility with the [print()] generic.
#' @return `x`, invisibly. Called for the side effect of printing.
#' @exportS3Method print engine_hooks
#' @noRd
print.engine_hooks <- function(x, ...) {
  cat("<engine_hooks>\n")
  if (length(x$active) == 0L) {
    cat("  none set (core MicSim: the driver runs unhooked)\n")
    return(invisible(x))
  }
  observing <- intersect(x$active, OBSERVING_HOOKS)
  adjusting <- intersect(x$active, ADJUSTING_HOOKS)
  if (length(observing) > 0L) {
    cat("  observing: ", paste(observing, collapse = ", "), "\n", sep = "")
  }
  if (length(adjusting) > 0L) {
    cat("  adjusting: ", paste(adjusting, collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}
