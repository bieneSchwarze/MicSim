# =====================================================================
# test-immigration.R
# ---------------------------------------------------------------------
# Immigration is the third of the four mechanisms beyond an ordinary
# transition, and the only one that is NOT a subsystem. It is not a
# transition at all: there is no origin state and no rate, only a table of
# who arrives and when. So it is handled where individuals are registered,
# and these tests are about registration and timing rather than about
# hooks.
#
# The thing worth checking hardest: an immigrant's clock starts at their
# ARRIVAL. Their age, the calendar time their rates see, and the window
# they have left are all measured from the day they enter the population,
# not from the start of the simulation.
#
# Tests are grouped by what they establish:
#   1. validation
#   2. registration
#   3. timing
#   4. the output contract
# =====================================================================


mig_codec <- function() {
  new_state_codec(list(sex = c("m", "f")), absorbing = "dead")
}

mig_registry <- function(codec = mig_codec(), rate_value = 0.05) {
  constant <- function(age, calTime) rep(rate_value, length(age))
  new_rate_registry(
    list(list(from = NA, to = "dead", rate = constant)),
    codec,
    probe = micsim_rate_probe(0:120, 1990:2060, 0:1)
  )
}

mig_init_pop <- function(n = 10L) {
  data.frame(
    ID        = seq_len(n),
    birthDate = rep("19700101", n),
    initState = rep(c("m", "f"), length.out = n),
    stringsAsFactors = FALSE
  )
}

#' Immigrants arriving one per year from 2002, all aged 30 on arrival
mig_table <- function(n = 5L, first_year = 2002L) {
  years <- first_year + seq_len(n) - 1L
  arrival <- years * 10000L + 601L            # 1 June of each year
  birth   <- (years - 30L) * 10000L + 601L    # aged exactly 30
  data.frame(
    ID              = 1000L + seq_len(n),
    immigrDate      = arrival,
    birthDate       = birth,
    immigrInitState = rep(c("f", "m"), length.out = n),
    stringsAsFactors = FALSE
  )
}

mig_engine <- function(n = 10L, immigrants = mig_table(), ...) {
  codec <- mig_codec()
  new_engine(codec, mig_registry(codec), mig_init_pop(n),
             horizon = c(20000101, 20201231), max_age = 120,
             immigrants = immigrants, ...)
}


# ---------------------------------------------------------------------
# 1. Validation
# ---------------------------------------------------------------------

test_that("the immigrant table must have the contract columns", {
  codec <- mig_codec()
  build <- function(tbl) {
    new_engine(codec, mig_registry(codec), mig_init_pop(5L),
               horizon = c(20000101, 20201231), max_age = 120,
               immigrants = tbl)
  }
  wrong_names <- mig_table(2L)
  names(wrong_names)[2L] <- "arrival"

  expect_error(build(wrong_names), class = "micsim_engine_error")
  expect_error(build(mig_table(0L)), class = "micsim_engine_error")
})


test_that("arriving before being born is refused", {
  tbl <- mig_table(2L)
  tbl$birthDate[1L] <- 20500101L
  expect_error(mig_engine(immigrants = tbl), class = "micsim_engine_error")
})


test_that("arriving outside the horizon is refused", {
  # Before the start: anyone present when the simulation begins belongs in
  # the starting population. After the end: the arrival can have no effect.
  too_early <- mig_table(2L)
  too_early$immigrDate[1L] <- 19990601L
  expect_error(mig_engine(immigrants = too_early),
               class = "micsim_engine_error")

  too_late <- mig_table(2L)
  too_late$immigrDate[2L] <- 20500601L
  expect_error(mig_engine(immigrants = too_late),
               class = "micsim_engine_error")
})


test_that("arriving older than maxAge is refused", {
  codec <- mig_codec()
  tbl <- mig_table(2L)
  expect_error(
    new_engine(codec, mig_registry(codec), mig_init_pop(5L),
               horizon = c(20000101, 20201231), max_age = 20,
               immigrants = tbl),
    class = "micsim_engine_error"
  )
})


test_that("an identifier that collides with the starting population is refused", {
  # Refused when the engine is built (the person store would also refuse
  # at registration). MicSim 3.0.0 only warned about collisions, and only
  # in its parallel path.
  tbl <- mig_table(2L)
  tbl$ID <- c(1L, 2L)                       # already used by the starters
  expect_error(mig_engine(immigrants = tbl), "both the starting",
               class = "micsim_engine_error")
})


# ---------------------------------------------------------------------
# 2. Registration
# ---------------------------------------------------------------------

test_that("immigrants join the population and appear in the output", {
  set.seed(2L)
  result <- mig_engine(n = 10L, immigrants = mig_table(5L))$run()

  expect_identical(result$diagnostics$n_immigrants, 5L)
  expect_identical(result$diagnostics$n_individuals, 15L)
  expect_setequal(unique(result$pop$ID), c(1:10, 1001:1005))
})


test_that("an immigrant's initial state is the state they arrive in", {
  set.seed(3L)
  pop <- mig_engine(n = 4L, immigrants = mig_table(4L))$run()$pop
  for (i in 1:4) {
    expect_identical(unique(pop$initState[pop$ID == 1000L + i]),
                     c("f", "m")[(i - 1L) %% 2L + 1L])
  }
})


test_that("an immigrant's birth date is returned as supplied", {
  set.seed(4L)
  tbl <- mig_table(3L)
  pop <- mig_engine(n = 4L, immigrants = tbl)$run()$pop
  for (i in seq_len(nrow(tbl))) {
    expect_identical(unique(pop$birthDate[pop$ID == tbl$ID[i]]),
                     as.character(tbl$birthDate[i]))
  }
})


test_that("a run with no immigration is unchanged by the feature existing", {
  # immigrants = NULL must be exactly the stage-1 engine.
  codec <- mig_codec()
  without <- new_engine(codec, mig_registry(codec), mig_init_pop(10L),
                        horizon = c(20000101, 20201231), max_age = 120)
  set.seed(5L)
  a <- without$run()$pop
  set.seed(5L)
  b <- without$run()$pop
  expect_identical(a, b)
  expect_identical(nrow(a), 10L)
})


# ---------------------------------------------------------------------
# 3. Timing
# ---------------------------------------------------------------------

test_that("an immigrant experiences nothing before they arrive", {
  # The single most important property. MicSim 3.0.0 achieves it with a
  # lag term added to the waiting time; we achieve it by scheduling from
  # the arrival date, which needs no correction because the queue stores
  # absolute dates.
  set.seed(6L)
  tbl <- mig_table(5L)
  pop <- mig_engine(n = 20L, immigrants = tbl)$run()$pop

  events <- pop[!is.na(pop$transitionTime) & pop$ID %in% tbl$ID, ,
                drop = FALSE]
  expect_gt(nrow(events), 0L)
  for (i in seq_len(nrow(events))) {
    arrival <- tbl$immigrDate[tbl$ID == events$ID[i]]
    expect_gte(getInDays(events$transitionTime[i]), getInDays(arrival))
  }
})


test_that("an immigrant's age at an event counts from their birth", {
  # Not from arrival. They were alive before they got here.
  set.seed(7L)
  tbl <- mig_table(5L)
  pop <- mig_engine(n = 20L, immigrants = tbl)$run()$pop
  events <- pop[!is.na(pop$transitionAge) & pop$ID %in% tbl$ID, ,
                drop = FALSE]
  expect_gt(nrow(events), 0L)

  for (i in seq_len(nrow(events))) {
    birth <- getInDays(events$birthDate[i])
    when  <- getInDays(events$transitionTime[i])
    # Within a few days, not exactly. Two separate roundings stand between
    # the two numbers, and neither is avoidable:
    #
    #   * `transitionTime` is a calendar DATE, so it is the whole day the
    #     event fell inside -- up to one day, 0.0027 years;
    #   * `transitionAge` is rounded to two decimals of a year -- up to
    #     half a unit, 0.005 years.
    #
    # So the age cannot be recovered from the reported date more finely
    # than about 0.008 years, and the same is true of MicSim 3.0.0's
    # output, for the same two reasons.
    expect_lt(abs(events$transitionAge[i] - (when - birth) / 365.25), 0.01)
    # everyone in this fixture arrives aged thirty
    expect_gte(events$transitionAge[i], 29.99)
  }
})


test_that("immigrants are scheduled after the whole starting population", {
  # MicSim 3.0.0 schedules every member of the starting population before
  # it looks at the first immigrant. The random numbers no longer depend on
  # this order (every individual has their own stream), but the order of
  # events on the same day, and so the numbering of newborns, does.
  codec <- mig_codec()
  order_seen <- new.env(parent = emptyenv())
  order_seen$slots <- integer(0)

  # A no-op modify_hazard hook fires once per destination per scheduling
  # call, which is the only way to observe scheduling order from outside.
  watcher <- new_engine_hooks(
    modify_hazard = function(hazard, ctx) {
      order_seen$slots <- c(order_seen$slots, ctx$slot)
      hazard
    }
  )
  engine <- new_engine(codec, mig_registry(codec), mig_init_pop(6L),
                       horizon = c(20000101, 20201231), max_age = 120,
                       immigrants = mig_table(3L), hooks = watcher)
  set.seed(8L)
  engine$run()

  first_seven <- order_seen$slots[1:9]
  expect_identical(first_seven[1:6], 1:6)     # the starting population
  expect_identical(first_seven[7:9], 7:9)     # then the immigrants
})


# ---------------------------------------------------------------------
# 4. The output contract
# ---------------------------------------------------------------------

test_that("an immigrant who never transitions gets one NA row", {
  set.seed(9L)
  tbl <- mig_table(5L)
  low_risk <- mig_registry(mig_codec(), rate_value = 1e-6)
  codec <- mig_codec()
  engine <- new_engine(codec, low_risk, mig_init_pop(5L),
                       horizon = c(20000101, 20201231), max_age = 120,
                       immigrants = tbl)
  pop <- engine$run()$pop

  expect_identical(nrow(pop), 10L)
  expect_true(all(is.na(pop$From)))
  expect_setequal(pop$ID, c(1:5, tbl$ID))
})


test_that("the run is reproducible with immigration on", {
  engine <- mig_engine(n = 10L, immigrants = mig_table(5L))
  expect_identical(engine$run(seed = 31L)$pop, engine$run(seed = 31L)$pop)
})
