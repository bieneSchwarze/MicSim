# =====================================================================
# test-solver.R
# ---------------------------------------------------------------------
# The waiting-time solver is the one component of MicSim that can be
# tested against ANALYTIC truth rather than against itself. A constant
# hazard must produce exponentially distributed waiting times, and we know
# that distribution in closed form. MicSim 3.0.0 had no such test; the
# whole package had one test file, and it asserted the *wrong* values of
# the broken date arithmetic.
#
# Tests are grouped by what they establish:
#   1. closed-form correctness      (the mathematics is right)
#   2. the three early exits        (Inf, all-zero, budget not exhausted)
#   3. forced transitions           (Inf rates)
#   4. contract violations          (classed errors on bad rates)
#   5. reliability reporting        (defect D5)
#   6. RNG discipline               (the parity contract)
# =====================================================================


# A counter that does not need `<<-`. Assigning into an environment is an
# ordinary assignment, so this keeps the test suite inside the package's
# own lint rules (see .lintr: superassignment is banned, defect D8).
new_counting_draw <- function(value = 1) {
  state <- new.env(parent = emptyenv())
  state$n <- 0L
  list(
    draw  = function() { state$n <- state$n + 1L; value },
    count = function() state$n
  )
}


# ---------------------------------------------------------------------
# 1. Closed-form correctness
# ---------------------------------------------------------------------

test_that("a constant hazard inverts to u / lambda exactly", {
  # With lambda(x) = c the integrated hazard is Lambda(tau) = c * tau, so
  # solving Lambda(tau) = u gives tau = u / c. This is an identity the
  # solver must reproduce to numerical precision; no simulation needed.
  constant <- function(x) rep(0.1, length(x))

  res <- draw_waiting_time(constant, horizon_years = 100,
                           draw_u = function() 1)
  expect_equal(res$time, 10, tolerance = 1e-6)
  expect_false(res$forced)
  expect_true(res$reliable)

  res2 <- draw_waiting_time(constant, horizon_years = 100,
                            draw_u = function() 0.5)
  expect_equal(res2$time, 5, tolerance = 1e-6)
})


test_that("a linear hazard inverts to the closed-form root", {
  # lambda(x) = a * x  =>  Lambda(tau) = a * tau^2 / 2  =>
  # tau = sqrt(2 * u / a).
  a <- 0.02
  linear <- function(x) a * x
  u <- 0.9

  res <- draw_waiting_time(linear, horizon_years = 100,
                           draw_u = function() u)
  expect_equal(res$time, sqrt(2 * u / a), tolerance = 1e-6)
})


test_that("waiting times from a constant hazard are Exponential(lambda)", {
  # The distributional check: draws must pass a Kolmogorov-Smirnov test
  # against the distribution they are supposed to follow. Seeded, so the
  # outcome is deterministic.
  skip_on_cran()
  set.seed(20260914L)
  lambda_value <- 0.25
  constant <- function(x) rep(lambda_value, length(x))

  draws <- vapply(
    seq_len(1000L),
    function(i) draw_waiting_time(constant, horizon_years = 500)$time,
    numeric(1L)
  )

  expect_true(all(is.finite(draws)))
  ks <- suppressWarnings(stats::ks.test(draws, "pexp", rate = lambda_value))
  expect_gt(ks$p.value, 0.01)
  expect_equal(mean(draws), 1 / lambda_value, tolerance = 0.15)
})


# ---------------------------------------------------------------------
# 2. The three early exits
# ---------------------------------------------------------------------

test_that("a non-positive horizon yields Inf and consumes exactly one draw", {
  # The engine legitimately asks about individuals who have reached maxAge
  # or whose observation window closes today. That is not an error. The
  # budget is used all the same, as on every path.
  constant <- function(x) rep(0.1, length(x))

  for (h in c(0, -0.004)) {
    counter <- new_counting_draw()
    res <- draw_waiting_time(constant, horizon_years = h,
                             draw_u = counter$draw)
    expect_identical(res$time, Inf)
    expect_identical(counter$count(), 1L)
    expect_identical(res$u, 1)
  }
})


test_that("an identically zero hazard yields Inf", {
  zero <- function(x) rep(0, length(x))
  res <- draw_waiting_time(zero, horizon_years = 50, draw_u = function() 1)
  expect_identical(res$time, Inf)
  expect_false(res$forced)
})


test_that("a budget larger than the total accumulated risk yields Inf", {
  # Total integrated hazard over 10 years at rate 0.01 is 0.1. A budget of
  # 5 is never exhausted, so no event occurs inside the window.
  small <- function(x) rep(0.01, length(x))
  res <- draw_waiting_time(small, horizon_years = 10, draw_u = function() 5)
  expect_identical(res$time, Inf)
  expect_equal(res$u, 5)
})


# ---------------------------------------------------------------------
# 3. Forced transitions
# ---------------------------------------------------------------------

test_that("an Inf rate produces a deterministic time and still uses its budget", {
  # School enrolment is modelled as rate = Inf at a known point. The
  # solver must return that point and flag it as forced. The budget is
  # drawn all the same and discarded, so that a forced transition never
  # moves the draws of anything else.
  constant <- function(x) rep(0.1, length(x))
  counter <- new_counting_draw()

  res <- draw_waiting_time(constant, horizon_years = 10,
                           draw_u = counter$draw,
                           forced_time = function(H) 2)

  expect_true(res$forced)
  expect_identical(counter$count(), 1L)
  expect_identical(res$time, 2)
  expect_identical(res$u, NA_real_)
})


test_that("forced_time = NULL, or no forced point, skips to the inversion", {
  constant <- function(x) rep(0.1, length(x))
  res <- draw_waiting_time(constant, horizon_years = 10,
                           draw_u = function() 1, forced_time = NULL)
  expect_false(res$forced)
  expect_equal(res$time, 10, tolerance = 1e-6)
  res <- draw_waiting_time(constant, horizon_years = 10,
                           draw_u = function() 1,
                           forced_time = function(H) NA_real_)
  expect_false(res$forced)
})


test_that("the all-zero shortcut sees a rate that starts in the last partial year", {
  # MicSim 3.0.0 tested whole years from now only. A rate positive only
  # after 0.5 years, in a window of 0.8 years, was declared identically
  # zero and its event could never happen.
  late <- function(x) ifelse(x > 0.5, 50, 0)
  res <- draw_waiting_time(late, horizon_years = 0.8, draw_u = function() 1)
  expect_true(is.finite(res$time))
  expect_gt(res$time, 0.5)
})


test_that("an Inf the forced finder cannot place is an error, not a miss", {
  sloped <- function(x) ifelse(x > 0.3, Inf, 0.1)
  expect_error(draw_waiting_time(sloped, horizon_years = 1,
                                 draw_u = function() 1),
               class = "micsim_rate_error")
})


test_that("degraded quadrature inside the root search is reported", {
  # The total over the window integrates cleanly; the integrals the root
  # search asks for fall back. MicSim 3.0.0, and 4.0.0 until November
  # 2026, reported only the first.
  local_mocked_bindings(integrate_hazard = function(lambda, lower, upper, ...) {
    value <- 0.1 * (upper - lower)
    list(value = value, reliable = upper >= 10)
  })
  constant <- function(x) rep(0.1, length(x))
  res <- draw_waiting_time(constant, horizon_years = 10,
                           draw_u = function() 0.5)
  expect_equal(res$time, 5, tolerance = 1e-6)
  expect_false(res$reliable)
})


# ---------------------------------------------------------------------
# 4. Contract violations
# ---------------------------------------------------------------------

test_that("a non-vectorised rate is rejected with a classed error", {
  # The single most common user error: function(age) 0.1 instead of
  # function(age) rep(0.1, length(age)). MicSim 3.0.0 surfaced this as an
  # opaque message from integrate(), minutes into a run.
  not_vectorised <- function(x) 0.1
  expect_error(
    draw_waiting_time(not_vectorised, horizon_years = 10,
                      draw_u = function() 1),
    class = "micsim_rate_error"
  )
})


test_that("a negative rate is rejected with a classed error", {
  negative <- function(x) rep(-0.1, length(x))
  expect_error(
    draw_waiting_time(negative, horizon_years = 10, draw_u = function() 1),
    class = "micsim_rate_error"
  )
})


test_that("a rate returning NA is rejected with a classed error", {
  na_rate <- function(x) rep(NA_real_, length(x))
  expect_error(
    draw_waiting_time(na_rate, horizon_years = 10, draw_u = function() 1),
    class = "micsim_rate_error"
  )
})


# ---------------------------------------------------------------------
# 5. Reliability reporting (defect D5)
# ---------------------------------------------------------------------

test_that("integrate_hazard reports success as reliable", {
  constant <- function(x) rep(0.1, length(x))
  res <- integrate_hazard(constant, 0, 10)
  expect_equal(res$value, 1, tolerance = 1e-8)
  expect_true(res$reliable)
})


test_that("integrate_hazard always returns a value and a reliability flag", {
  # A hazard that is zero everywhere except a very narrow spike can defeat
  # adaptive quadrature at the default tolerance over a long range. The
  # point is not the numeric value but that whatever happens is REPORTED
  # rather than swallowed, which is the whole of defect D5.
  spike <- function(x) ifelse(abs(x - 500) < 1e-4, 1e5, 0)
  res <- integrate_hazard(spike, 0, 1000)
  expect_type(res$value, "double")
  expect_type(res$reliable, "logical")
  expect_length(res$reliable, 1L)
})


# ---------------------------------------------------------------------
# 6. RNG discipline
# ---------------------------------------------------------------------

test_that("the default draw is a standard exponential by inversion", {
  # The engine's budgets use the same expression (budgets_from_uniforms()).
  set.seed(1L)
  expected <- -log(1 - stats::runif(1L))
  set.seed(1L)
  expect_identical(default_draw_u(), expected)
})


test_that("exactly one random number is consumed on the ordinary path", {
  constant <- function(x) rep(0.1, length(x))

  set.seed(42L)
  draw_waiting_time(constant, horizon_years = 100)
  after_solver <- .Random.seed

  set.seed(42L)
  invisible(stats::runif(1L))
  after_one_runif <- .Random.seed

  expect_identical(after_solver, after_one_runif)
})


test_that("every path uses exactly one budget", {
  # An all-zero rate and an exhausted window use their budget too, like
  # the ordinary and the forced path.
  zero <- function(x) rep(0, length(x))
  counter <- new_counting_draw()
  res <- draw_waiting_time(zero, horizon_years = 10, draw_u = counter$draw)
  expect_identical(res$time, Inf)
  expect_identical(counter$count(), 1L)

  counter <- new_counting_draw()
  res <- draw_waiting_time(zero, horizon_years = 0, draw_u = counter$draw)
  expect_identical(res$time, Inf)
  expect_identical(counter$count(), 1L)
})


test_that("a supplied budget is used and nothing is drawn", {
  constant <- function(x) rep(0.1, length(x))
  counter <- new_counting_draw()
  res <- draw_waiting_time(constant, horizon_years = 100,
                           draw_u = counter$draw, budget = 1)
  expect_identical(counter$count(), 0L)
  expect_equal(res$time, 10, tolerance = 1e-6)
  expect_identical(res$u, 1)
})
