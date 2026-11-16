# =====================================================================
# test-reestimate.R
# ---------------------------------------------------------------------
# The simulation result as a record of the run, and the reestimation of
# the transition rates from it.
#
#   1. what a result records (meta, diagnostics, subsystem reports)
#   2. reestimation: exposure, events, expected events, invariants
#   3. edge cases (forced, nobody at risk, no events)
#   4. printing, conversion and plotting
# =====================================================================


re_mort <- function(age) 0.0004 * exp(0.08 * age)

re_model <- function(n = 400L, seed = 7L) {
  set.seed(seed)
  micsim_model(
    state_space(sex = c("m", "f"), edu = c("no", "low"), par = c("0", "1")),
    transitions(
      transition(move("edu", "no", "low"), function(age) ifelse(age == 7, Inf, 0)),
      transition(move("par", "0", "1", where = sex == "f"),
                 function(age, calTime) ifelse(age > 15 & age < 45, 0.1, 0)),
      transition(NA, "dead", re_mort)
    ),
    data.frame(ID = seq_len(n),
               birthDate = sprintf("%d0315", sample(1950:1999, n, TRUE)),
               initState = rep(c("f/low/0", "m/low/0"), length.out = n)),
    c(20000101, 20201231)
  ) |>
    add_subsystem(fertility(moves_of("par"),
                            newborn_rule(sex = c(m = 0.5, f = 0.5), edu = "no", par = "0"))) |>
    add_subsystem(school_enrolment(move("edu", "no", "low"), month = 9))
}


# ---------------------------------------------------------------------
# 1. What a result records
# ---------------------------------------------------------------------

test_that("a result records the package version, platform and run", {
  r <- simulate(re_model(), seed = 1L)
  expect_named(r, c("pop", "model", "seed", "diagnostics", "extensions",
                    "reestimation", "meta"))
  expect_identical(r$meta$micsim_version,
                   as.character(utils::packageVersion("MicSim")))
  expect_identical(r$meta$r_version, R.version.string)
  expect_s3_class(r$meta$created, "POSIXct")
  expect_type(r$meta$rng_kind, "character")
  expect_gte(r$meta$reestimation_seconds, 0)
  expect_identical(r$seed, 1L)
})


test_that("diagnostics account for every individual", {
  r <- simulate(re_model(), seed = 2L)
  d <- r$diagnostics
  expect_identical(d$n_starting + d$n_immigrants + d$n_newborns, d$n_individuals)
  expect_identical(d$n_individuals, length(unique(r$pop$ID)))
  expect_identical(d$n_events, sum(!is.na(r$pop$From)))
  expect_identical(d$n_newborns, r$extensions$fertility$n_births)
  expect_identical(d$cores, 1L)
  expect_gte(d$seconds, 0)
})


test_that("the fertility report lists every birth, consistent with pop", {
  r <- simulate(re_model(), seed = 3L)
  f <- r$extensions$fertility
  expect_identical(nrow(f$births), f$n_births)
  expect_identical(sum(f$by_transition), f$n_births)
  expect_identical(sum(f$by_child_state), f$n_births)
  expect_named(f$by_transition, c("f/no/0 -> f/no/1", "f/low/0 -> f/low/1"),
               ignore.order = TRUE)
  people <- r$pop[!duplicated(r$pop$ID), ]
  child <- people[match(f$births$childID, people$ID), ]
  expect_identical(child$motherID, f$births$motherID)
  expect_identical(child$birthDate, f$births$date)
  expect_identical(child$initState, f$births$childState)
  mother <- people[match(f$births$motherID, people$ID), ]
  age <- (getInDays(f$births$date) - getInDays(mother$birthDate)) / DAYS_PER_YEAR
  expect_true(all(abs(age - f$births$motherAge) < 0.01))
})


test_that("the enrolment report lists every enrolment and its shift", {
  r <- simulate(re_model(), seed = 4L)
  e <- r$extensions$school_enrolment
  expect_identical(nrow(e$enrolments), e$n_enrolments)
  expect_gt(e$n_enrolments, 0L)
  moved <- r$pop[!is.na(r$pop$From) & grepl("/no/", r$pop$From) &
                   grepl("/low/", r$pop$To), ]
  expect_setequal(e$enrolments$ID, moved$ID)
  expect_true(all(substr(e$enrolments$date, 5L, 6L) == "09" |
                    e$n_at_present > 0L))
  expect_identical(e$enrolments$shiftDays,
                   as.integer(getInDays(e$enrolments$date) -
                                getInDays(e$enrolments$drawnDate)))
  expect_identical(sum(e$by_month), e$n_enrolments)
})


test_that("parallel results renumber report identifiers with the population", {
  skip_on_cran()
  r <- simulate(re_model(600L), seed = 5L, cores = 2L)
  f <- r$extensions$fertility
  expect_identical(nrow(f$births), f$n_births)
  people <- r$pop[!duplicated(r$pop$ID), ]
  child <- people[match(f$births$childID, people$ID), ]
  expect_false(anyNA(child$ID))
  expect_identical(child$motherID, f$births$motherID)
  # the same fields as a single-core result, in the same order
  one <- simulate(re_model(600L), seed = 5L)
  expect_identical(names(r$diagnostics), names(one$diagnostics))
  expect_identical(names(r$meta), names(one$meta))
  expect_identical(r$diagnostics$cores, 2L)
  expect_s3_class(r$reestimation, "micsim_reestimation")
})


# ---------------------------------------------------------------------
# 2. Reestimation
# ---------------------------------------------------------------------

test_that("one entry per declared transition, named by its label", {
  r <- simulate(re_model(), seed = 6L)
  expect_s3_class(r$reestimation, "micsim_reestimation")
  expect_named(r$reestimation, c("edu: no -> low",
                                 "par: 0 -> 1 (where sex == \"f\")",
                                 "any transient state -> dead"))
  expect_identical(r$reestimation[[3L]]$rate_name, "re_mort")     # by name
  expect_identical(r$reestimation[[1L]]$rate_name, NA_character_) # anonymous
})


test_that("exposure and events equal a direct count from pop", {
  # Mortality only, everyone present from the start: the person-years are
  # the time from the start to death or the end of observation.
  mort <- function(age) rep(0.03, length(age))
  pop <- data.frame(ID = 1:300, birthDate = "19600101", initState = "a")
  r <- simulate(micsim_model(state_space(s = c("a", "b")),
                             transitions(transition(NA, "dead", mort)),
                             pop, c(20000101, 20191231)), seed = 1L)
  e <- r$reestimation[[1L]]
  death <- r$pop[!is.na(r$pop$To) & r$pop$To == "dead", ]
  stop <- rep(getInDays(20191231), nrow(pop))
  stop[match(death$ID, pop$ID)] <- getInDays(death$transitionTime)
  py <- sum(stop - getInDays(20000101)) / DAYS_PER_YEAR
  expect_equal(e$exposure, py, tolerance = 1e-9)
  expect_identical(e$events, nrow(death))
  expect_equal(e$expected, 0.03 * py, tolerance = 1e-9)
  expect_equal(sum(e$bands$exposure), py, tolerance = 1e-9)
})


test_that("the supplied rate is recovered within its confidence interval", {
  mort <- function(age) 0.001 * exp(0.07 * age)
  pop <- data.frame(ID = 1:3000, birthDate = sprintf("%d0101", rep(1930:1979, 60)),
                    initState = "a")
  r <- simulate(micsim_model(state_space(s = c("a", "b")),
                             transitions(transition(NA, "dead", mort)),
                             pop, c(20000101, 20291231)), seed = 3L)
  e <- r$reestimation[[1L]]
  expect_true(e$ratio_ci[1L] <= 1 && e$ratio_ci[2L] >= 1)
  expect_gt(nrow(e$smooth), 0L)
  # the smooth follows the supplied rate where the data are dense
  mid <- e$smooth[e$smooth$x > 50 & e$smooth$x < 90, ]
  expect_lt(max(abs(log(mid$rate / mort(mid$x)))), log(1.25))
})


test_that("reestimation draws no random numbers and can be repeated", {
  r <- simulate(re_model(), seed = 8L)
  set.seed(99L)
  before <- .Random.seed
  again <- reestimate(r)
  expect_identical(.Random.seed, before)
  expect_identical(again, r$reestimation)
})


test_that("reestimate = FALSE skips it; plot() then computes it", {
  r <- simulate(re_model(), seed = 9L, reestimate = FALSE)
  expect_null(r$reestimation)
  expect_true(is.na(r$meta$reestimation_seconds))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_no_error(plot(r, which = 3L))
})


test_that("a rate on two time scales is compared with its average over those at risk", {
  r <- simulate(re_model(), seed = 10L)
  e <- r$reestimation[[2L]]
  expect_identical(e$axis, "age")
  expect_identical(e$varying, "calTime")
  expect_true(all(c("lower", "upper") %in% names(e$truth)))
  ok <- is.finite(e$truth$rate) & is.finite(e$truth$lower)
  expect_true(all(e$truth$lower[ok] <= e$truth$rate[ok] + 1e-12))
  expect_true(all(e$truth$upper[ok] >= e$truth$rate[ok] - 1e-12))
})


test_that("calendar bands start on the whole year", {
  calRate <- function(calTime) 0.02 + 0.001 * (calTime - 2000)
  pop <- data.frame(ID = 1:200, birthDate = "19700101", initState = "a")
  r <- simulate(micsim_model(state_space(s = c("a", "b")),
                             transitions(transition("a", "b", calRate)),
                             pop, c(20000101, 20091231)), seed = 1L)
  b <- r$reestimation[[1L]]$bands
  expect_identical(r$reestimation[[1L]]$axis, "calTime")
  expect_equal(b$lower[2L], 2001)
  expect_identical(nrow(b), 10L)
})


test_that("arguments of reestimate() are checked and take effect", {
  r <- simulate(re_model(), seed = 11L)
  expect_error(reestimate(r, width = -1), class = "micsim_reestimation_error")
  expect_error(reestimate(r, groups = 5), class = "micsim_reestimation_error")
  expect_error(reestimate(r, axis = "height"), class = "micsim_reestimation_error")
  w <- reestimate(r, width = 5)
  expect_true(all(diff(w[[3L]]$bands$lower) == 5))
  cal <- reestimate(r, axis = "calTime")
  expect_identical(cal[[2L]]$axis, "calTime")
  expect_identical(cal[[3L]]$axis, "age")            # mortality has no calTime
})


# ---------------------------------------------------------------------
# 3. Edge cases
# ---------------------------------------------------------------------

test_that("forced transitions, empty risk sets and rare events are reported", {
  pop <- data.frame(ID = 1:50, birthDate = "19950101", initState = "x/0")
  m <- micsim_model(
    state_space(s = c("x", "y", "w"), j = c("0", "1")),
    transitions(transition(move("j", "0", "1"), function(age) ifelse(age == 18, Inf, 0)),
                transition(move("s", "w", "x"), function() 0.1),
                transition(move("s", "x", "y"), function(age) rep(1e-9, length(age)))),
    pop, c(20000101, 20201231), max_age = 60)
  r <- simulate(m, seed = 1L)
  forced <- r$reestimation[[1L]]
  expect_true(forced$forced)
  expect_true(is.na(forced$expected) && is.na(forced$ratio))
  expect_identical(forced$events, 50L)
  expect_identical(nrow(forced$smooth), 0L)
  empty <- r$reestimation[[2L]]
  expect_identical(empty$exposure, 0)
  expect_match(empty$notes, "Nobody was at risk", all = FALSE)
  rare <- r$reestimation[[3L]]
  expect_identical(rare$events, 0L)
  expect_match(rare$notes, "No events", all = FALSE)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_no_warning(plot(r, ask = FALSE))
})


test_that("a timing rule on a finite rate is noted", {
  pop <- data.frame(ID = 1:40, birthDate = "19950101", initState = "no")
  m <- micsim_model(state_space(edu = c("no", "low")),
                    transitions(transition("no", "low", function(age)
                      ifelse(age >= 10 & age < 12, 0.3, 0))),
                    pop, c(20000101, 20201231)) |>
    add_subsystem(school_enrolment(move("edu", "no", "low"), 9))
  r <- simulate(m, seed = 1L)
  expect_identical(r$reestimation[[1L]]$timed_by, "school_enrolment")
  expect_match(r$reestimation[[1L]]$notes, "adjusted by subsystem", all = FALSE)
})


# ---------------------------------------------------------------------
# 4. Printing, conversion and plotting
# ---------------------------------------------------------------------

test_that("print, as.data.frame and summary of a reestimation", {
  r <- simulate(re_model(), seed = 12L)
  expect_output(print(r$reestimation), "obs/exp \\[95% CI\\]")
  expect_output(print(r$reestimation), "forced")
  expect_output(print(r$reestimation[[3L]]), "observed/expected")
  expect_output(print(r), "reestimation")
  df <- as.data.frame(r$reestimation)
  expect_true(all(c("transition", "axis", "group", "lower", "upper", "events",
                    "exposure", "rate", "ci_lower", "ci_upper", "expected") %in%
                    names(df)))
  expect_equal(sum(df$events), sum(vapply(r$reestimation, `[[`, integer(1L), "events")))
})


test_that("plot draws every transition and checks its arguments", {
  r <- simulate(re_model(), seed = 13L)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_no_warning(plot(r, ask = FALSE))
  expect_no_warning(plot(r$reestimation, which = "any transient state -> dead",
                         log = TRUE))
  expect_no_warning(plot(r$reestimation[[2L]], log = FALSE))
  expect_error(plot(r, which = "no such transition"),
               class = "micsim_reestimation_error")
  expect_error(plot(r, which = 9L), class = "micsim_reestimation_error")
  expect_invisible(plot(r, which = 1L))
})
