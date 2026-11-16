# =====================================================================
# test-hooks.R
# ---------------------------------------------------------------------
# The hook set is an architectural seam, so the tests are about its
# contract rather than about any computation:
#
#   * with no hooks set, nothing happens and nothing is even built;
#   * an observing hook cannot change the simulation, however hard it
#     tries, because its return value is discarded;
#   * an adjusting hook can, and exactly two of them exist.
#
# Nothing in core MicSim sets a hook. These tests exist so that when the
# linked-lives subsystem starts setting them, the seam is already known to
# behave.
#
# Tests are grouped by what they establish:
#   1. the default set is inert
#   2. observing hooks
#   3. adjusting hooks
#   4. laziness of the context
#   5. validation
# =====================================================================


# A counter that needs no superassignment: assigning into an environment
# is an ordinary assignment, so the test suite stays inside the package's
# own lint rules (.lintr bans `<<-`, defect D8).
new_recorder <- function() {
  state <- new.env(parent = emptyenv())
  state$calls <- list()
  list(
    hook  = function(ctx) {
      state$calls[[length(state$calls) + 1L]] <- ctx
      "this return value must be discarded"
    },
    calls = function() state$calls,
    n     = function() length(state$calls)
  )
}


# ---------------------------------------------------------------------
# 1. The default set is inert
# ---------------------------------------------------------------------

test_that("a hook set with nothing in it reports nothing active", {
  hooks <- new_engine_hooks()
  expect_s3_class(hooks, "engine_hooks")
  expect_identical(hooks$active, character(0))

  for (which in c(OBSERVING_HOOKS, ADJUSTING_HOOKS)) {
    expect_false(hook_active(hooks, which))
  }
})


test_that("firing an unset observing hook does nothing", {
  hooks <- new_engine_hooks()
  expect_null(fire_hook(hooks, "on_event", new_hook_context(slot = 1L)))
})


test_that("an unset adjusting hook returns the value untouched", {
  hooks <- new_engine_hooks()
  value <- list(time_days = 100, to_idx = 4L)
  expect_identical(
    adjust_by_hook(hooks, "on_schedule", value, new_hook_context()),
    value
  )
})


# ---------------------------------------------------------------------
# 2. Observing hooks
# ---------------------------------------------------------------------

test_that("an observing hook receives the context it was given", {
  rec <- new_recorder()
  hooks <- new_engine_hooks(on_event = rec$hook)

  expect_true(hook_active(hooks, "on_event"))
  fire_hook(hooks, "on_event",
            new_hook_context(slot = 7L, from_idx = 3L, to_idx = 4L,
                             now_days = 17351))

  expect_identical(rec$n(), 1L)
  ctx <- rec$calls()[[1L]]
  expect_identical(ctx$slot, 7L)
  expect_identical(ctx$from_idx, 3L)
  expect_identical(ctx$to_idx, 4L)
  expect_identical(ctx$now_days, 17351)
})


test_that("an observing hook's return value is discarded", {
  # The guarantee that makes these hooks safe: a subsystem cannot change
  # the run by returning something from an observer, so a mistake in one
  # can only affect that subsystem's own state.
  rec <- new_recorder()
  hooks <- new_engine_hooks(on_birth = rec$hook)
  expect_null(fire_hook(hooks, "on_birth", new_hook_context(slot = 1L)))
})


test_that("the three observing hooks are independent", {
  event  <- new_recorder()
  birth  <- new_recorder()
  absorb <- new_recorder()
  hooks <- new_engine_hooks(on_event  = event$hook,
                            on_birth  = birth$hook,
                            on_absorb = absorb$hook)

  fire_hook(hooks, "on_event",  new_hook_context(slot = 1L))
  fire_hook(hooks, "on_absorb", new_hook_context(slot = 2L))
  fire_hook(hooks, "on_absorb", new_hook_context(slot = 3L))

  expect_identical(event$n(),  1L)
  expect_identical(birth$n(),  0L)
  expect_identical(absorb$n(), 2L)
  expect_setequal(hooks$active, c("on_event", "on_birth", "on_absorb"))
})


test_that("a hook may cancel a pending event through the queue", {
  # The one write a hook is sanctioned to make. Linked lives needs it for
  # separation and widowhood: the partner's scheduled event is withdrawn
  # and an induced transition put in its place.
  q <- new_event_queue()
  q$push(time_days = 100, slot = 1L, to_idx = 2L)
  q$push(time_days = 200, slot = 2L, to_idx = 2L)

  widow <- function(ctx) ctx$queue$cancel(ctx$partner_slot)
  hooks <- new_engine_hooks(on_absorb = widow)

  fire_hook(hooks, "on_absorb",
            new_hook_context(queue = q, slot = 1L, partner_slot = 2L))

  expect_false(q$has_pending(2L))
  expect_identical(q$size(), 1L)
})


# ---------------------------------------------------------------------
# 3. Adjusting hooks
# ---------------------------------------------------------------------

test_that("an adjusting hook replaces the value", {
  # on_schedule is given the chosen event and returns one of the same
  # shape. Snapping a date to a fixed day of the month is the kind of
  # thing it is for -- though core MicSim's school enrolment does that in
  # the driver, not through a hook, because it is not an extension.
  delay <- function(value, ctx) {
    value$time_days <- value$time_days + 30
    value
  }
  hooks <- new_engine_hooks(on_schedule = delay)

  adjusted <- adjust_by_hook(hooks, "on_schedule",
                             list(time_days = 100, to_idx = 4L),
                             new_hook_context(slot = 1L))

  expect_identical(adjusted$time_days, 130)
  expect_identical(adjusted$to_idx, 4L)
})


test_that("an adjusting hook sees the context as well as the value", {
  scale_by_slot <- function(value, ctx) value * ctx$slot
  hooks <- new_engine_hooks(modify_hazard = scale_by_slot)

  expect_identical(
    adjust_by_hook(hooks, "modify_hazard", 10, new_hook_context(slot = 3L)),
    30
  )
})


test_that("modify_hazard can wrap the pair of closures bind_hazard returns", {
  # The realistic shape of the hook: the engine hands over
  # list(lambda, forced_time) and gets a list of the same shape back. This
  # is how an extension would transform a hazard, e.g. the linked-lives
  # transmission that a later MicSim adds as a subsystem.
  double_it <- function(hazard, ctx) {
    inner <- hazard$lambda
    hazard$lambda <- function(x) 2 * inner(x)
    hazard
  }
  hooks <- new_engine_hooks(modify_hazard = double_it)

  original <- list(lambda = function(x) rep(0.1, length(x)),
                   forced_time = NULL)
  adjusted <- adjust_by_hook(hooks, "modify_hazard", original,
                             new_hook_context())

  expect_identical(adjusted$lambda(c(0, 1)), c(0.2, 0.2))
  expect_null(adjusted$forced_time)
})


# ---------------------------------------------------------------------
# 4. Laziness of the context
# ---------------------------------------------------------------------

test_that("no context is built when the hook is not set", {
  # This is what makes the seam free. R evaluates arguments lazily, so a
  # call such as fire_hook(hooks, "on_event", new_hook_context(...)) never
  # runs new_hook_context() at all when on_event is NULL. A run of core
  # MicSim therefore pays one is.null() test per hook point per event.
  built <- new.env(parent = emptyenv())
  built$n <- 0L
  expensive_context <- function() {
    built$n <- built$n + 1L
    new_hook_context(slot = 1L)
  }

  fire_hook(new_engine_hooks(), "on_event", expensive_context())
  expect_identical(built$n, 0L)

  # A hook that uses its argument forces the promise, so the context is
  # built exactly once, at the moment of the call.
  fire_hook(new_engine_hooks(on_event = function(ctx) ctx$slot), "on_event",
            expensive_context())
  expect_identical(built$n, 1L)
})


test_that("a hook that ignores its argument never builds the context either", {
  # Laziness goes one step further than the guard: passing the promise on
  # to hook(ctx) does not force it. This is a bonus rather than a hazard,
  # because a hook that never mentions ctx has nothing to be stale about,
  # and any hook that does use it forces it immediately and therefore sees
  # the values the driver held when it fired.
  built <- new.env(parent = emptyenv())
  built$n <- 0L
  expensive_context <- function() {
    built$n <- built$n + 1L
    new_hook_context(slot = 1L)
  }

  fire_hook(new_engine_hooks(on_event = function(ctx) "ignores ctx"),
            "on_event", expensive_context())
  expect_identical(built$n, 0L)
})


# ---------------------------------------------------------------------
# 5. Validation
# ---------------------------------------------------------------------

test_that("a hook that is not a function is rejected", {
  expect_error(new_engine_hooks(on_event = "not a function"),
               class = "micsim_hook_error")
  expect_error(new_engine_hooks(modify_hazard = 42),
               class = "micsim_hook_error")
})


test_that("the two hook-name constants do not overlap", {
  # Every hook belongs to exactly one kind. If that ever stops being true,
  # the contract at the top of engine_hooks.R has drifted from the code.
  expect_length(intersect(OBSERVING_HOOKS, ADJUSTING_HOOKS), 0L)
  expect_setequal(
    c(OBSERVING_HOOKS, ADJUSTING_HOOKS),
    setdiff(names(new_engine_hooks()), "active")
  )
})
