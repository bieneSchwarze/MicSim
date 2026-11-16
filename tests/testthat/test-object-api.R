# =====================================================================
# test-object-api.R
# ---------------------------------------------------------------------
# The object API adds to the legacy one; it does not replace it. So the
# single most important test in this file is the last one: the same model,
# expressed both ways, under the same seed, must give the same population
# down to the last row. If it ever does not, one of the two interfaces is
# lying about what it simulates.
#
# Tests are grouped by what they establish:
#   1. building a model
#   2. what building catches
#   3. validate()
#   4. simulating
#   5. the result object
#   6. conversion needs no arguments
#   7. the two interfaces agree
# =====================================================================


obj_mortality <- function(age, calTime) 0.00003 * exp(0.1 * age)
obj_fertility <- function(age, calTime) ifelse(age > 15 & age < 45, 0.2, 0)

obj_population <- function(n = 30L) {
  data.frame(
    ID = seq_len(n), birthDate = rep("19850101", n),
    initState = rep("f/0", n), stringsAsFactors = FALSE
  )
}

obj_model <- function(n = 30L, ...) {
  micsim_model(
    state_space = state_space(sex = c("m", "f"), fert = c("0", "1")),
    transitions = transitions(
      transition("f/0", "f/1",  obj_fertility),
      transition(NA,    "dead", obj_mortality)
    ),
    population = obj_population(n),
    horizon    = c(20000101, 20201231),
    max_age    = 100,
    ...
  ) |>
    add_subsystem(fertility(
      births  = "f/0->f/1",
      newborn = new_newborn_rule(rbind(c("m", "0"), c("f", "0")),
                                 prob = c(0.515, 0.485))
    ))
}


# ---------------------------------------------------------------------
# 1. Building a model
# ---------------------------------------------------------------------

test_that("a model is a value that can be held, printed and reused", {
  model <- obj_model()
  expect_s3_class(model, "micsim_model")
  expect_identical(model$codec$K, 4L)
  expect_length(model$registry$edges, 5L)   # 4 mortality + 1 fertility
  expect_output(print(model), "micsim_model")
})


test_that("building a model consumes no random numbers", {
  # A model is a description, not a run. Building one must be invisible to
  # the random stream, or a pipeline could not be reasoned about.
  set.seed(1L)
  before <- .Random.seed
  invisible(obj_model())
  expect_identical(.Random.seed, before)
})


test_that("a model can be saved and restored", {
  # Provenance: a result six months old should be traceable to the model
  # that made it. That requires the model to survive serialisation.
  model <- obj_model(n = 5L)
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)

  saveRDS(model, path)
  restored <- readRDS(path)

  expect_s3_class(restored, "micsim_model")
  set.seed(7L); a <- as.data.frame(simulate(model))
  set.seed(7L); b <- as.data.frame(simulate(restored))
  expect_identical(a, b)
})


# ---------------------------------------------------------------------
# 2. What building catches
# ---------------------------------------------------------------------

test_that("state values must be unique across variables", {
  # Transition patterns match on values rather than positions, so a value
  # used by two variables could not be told apart.
  expect_error(state_space(a = c("x", "y"), b = c("y", "z")),
               class = "micsim_model_error")
})


test_that("state variables must be named", {
  expect_error(state_space(c("m", "f")), class = "micsim_model_error")
})


test_that("a model with no transitions is refused", {
  expect_error(transitions(), class = "micsim_model_error")
})


test_that("an unknown state in a transition is caught at construction", {
  expect_error(
    micsim_model(
      state_space = state_space(sex = c("m", "f")),
      transitions = transitions(transition("x", "dead", obj_mortality)),
      population  = data.frame(ID = 1L, birthDate = "19850101",
                               initState = "f", stringsAsFactors = FALSE),
      horizon = c(20000101, 20101231)
    ),
    class = "micsim_registry_error"
  )
})


test_that("a non-vectorised rate is caught at construction", {
  broken <- function(age, calTime) 0.1
  expect_error(
    micsim_model(
      state_space = state_space(sex = c("m", "f")),
      transitions = transitions(transition(NA, "dead", broken)),
      population  = data.frame(ID = 1L, birthDate = "19850101",
                               initState = "f", stringsAsFactors = FALSE),
      horizon = c(20000101, 20101231)
    ),
    class = "micsim_rate_error"
  )
})


test_that("the components must be built with their own constructors", {
  expect_error(
    micsim_model(state_space = list(sex = c("m", "f")),
                 transitions = transitions(transition(NA, "dead",
                                                      obj_mortality)),
                 population = obj_population(2L),
                 horizon = c(20000101, 20101231)),
    class = "micsim_model_error"
  )
})


# ---------------------------------------------------------------------
# 3. validate()
# ---------------------------------------------------------------------

test_that("validate returns the model invisibly, so it sits in a pipe", {
  model <- obj_model(n = 5L)
  expect_identical(validate(model), model)
  expect_invisible(validate(model))
})


test_that("validate catches a rate that misbehaves somewhere on the grid", {
  # Construction checks that a rate is a vectorised function of admissible
  # arguments. validate() checks what it RETURNS, over the whole age and
  # calendar grid the simulation can reach, which is the expensive check.
  goes_negative <- function(age, calTime) age - 50
  model <- micsim_model(
    state_space = state_space(sex = c("m", "f")),
    transitions = transitions(transition(NA, "dead", goes_negative)),
    population  = data.frame(ID = 1L, birthDate = "19850101",
                             initState = "f", stringsAsFactors = FALSE),
    horizon = c(20000101, 20101231), max_age = 100
  )
  expect_error(validate(model), class = "micsim_rate_error")
})


test_that("a rate defined only over the model's own age range is accepted", {
  # REGRESSION. A rate table that starts at age 17 and is indexed as
  # vec[age - 16] fails, correctly, if asked about age 0. Neither the
  # construction-time shape probe nor validate() may ask it: the youngest
  # individual in this model is 18 and nobody is ever born into it, so
  # age 0 is unreachable.
  #
  # Both used to ask anyway -- the probe at a hard-coded c(0, 30, 60), the
  # validator over 0:max_age -- and refused a perfectly good model with
  # "only 0's may be mixed with negative subscripts".
  table_from_17 <- c(rep(0.05, 48))          # ages 17..64
  adult_rate <- function(age, calTime) table_from_17[as.integer(age) - 16]

  model <- micsim_model(
    state_space = state_space(sex = c("m", "f")),
    transitions = transitions(transition(NA, "dead", adult_rate)),
    population  = data.frame(
      ID = 1:10, birthDate = rep("19820101", 10L),
      initState = rep("f", 10L), stringsAsFactors = FALSE
    ),
    horizon = c(20000101, 20101231), max_age = 65
  )

  expect_s3_class(model, "micsim_model")
  expect_identical(validate(model), model)

  set.seed(17L)
  expect_s3_class(simulate(model), "micsim_result")
})


test_that("the validation grid stops below maxAge and at the youngest entrant", {
  # The engine gives an individual a window of (max_age - 0.01) - age, so
  # it never evaluates a rate at or beyond max_age. A table with one row
  # per year of age from 0 to max_age - 1 is therefore complete, and
  # asking it for the max_age-th row would be a demand the simulation was
  # never going to make.
  model <- micsim_model(
    state_space = state_space(sex = c("m", "f")),
    transitions = transitions(
      transition(NA, "dead", function(age, calTime) rep(0.01, length(age)))
    ),
    population = data.frame(
      ID = 1:5, birthDate = rep("19700101", 5L),
      initState = rep("f", 5L), stringsAsFactors = FALSE
    ),
    horizon = c(20000101, 20101231), max_age = 100
  )
  grid <- reachable_grid(model)

  expect_identical(max(grid$ages), 99L)      # not 100
  # Born 1970-01-01, observed from 2000-01-01: 10957 days, which is 29.998
  # years on the 365.25-day scale rate functions are defined on, hence 29.
  expect_identical(min(grid$ages), 29L)
  expect_identical(range(grid$cal_years), c(2000L, 2010L))
})


test_that("a model with fertility is validated from age zero", {
  # Because somebody is born into it during the run.
  model <- obj_model(n = 5L)
  expect_identical(min(reachable_grid(model)$ages), 0L)
})


# ---------------------------------------------------------------------
# 4. Simulating
# ---------------------------------------------------------------------

test_that("simulate returns a result carrying its model and seed", {
  result <- simulate(obj_model(), seed = 234L)
  expect_s3_class(result, "micsim_result")
  expect_s3_class(result$model, "micsim_model")
  expect_identical(result$seed, 234L)
})


test_that("the same seed gives the same result", {
  model <- obj_model()
  expect_identical(as.data.frame(simulate(model, seed = 11L)),
                   as.data.frame(simulate(model, seed = 11L)))
})


test_that("without a seed, the seed is drawn from the caller's stream and recorded", {
  model <- obj_model(n = 10L)
  set.seed(5L); a <- simulate(model, reestimate = FALSE)
  set.seed(5L); b <- simulate(model, reestimate = FALSE)
  expect_identical(as.data.frame(a), as.data.frame(b))
  expect_false(is.na(a$seed))
  expect_true(a$meta$seed_drawn)
  again <- simulate(model, seed = a$seed, reestimate = FALSE)
  expect_identical(as.data.frame(again), as.data.frame(a))
  expect_false(again$meta$seed_drawn)
})


test_that("nsim other than one is refused rather than ignored", {
  expect_error(simulate(obj_model(n = 5L), nsim = 3), 
               class = "micsim_model_error")
})


# ---------------------------------------------------------------------
# 5. The result object
# ---------------------------------------------------------------------

test_that("as.data.frame gives the contract data frame", {
  result <- simulate(obj_model(), seed = 3L)
  pop <- as.data.frame(result)
  expect_identical(
    names(pop),
    c("ID", "birthDate", "initState", "From", "To", "transitionTime",
      "transitionAge", "motherID")
  )
})


test_that("print and summary say what happened", {
  result <- simulate(obj_model(), seed = 4L)
  expect_output(print(result), "micsim_result")
  s <- summary(result)
  expect_s3_class(s, "summary.micsim_result")
  expect_identical(s$n_events, result$diagnostics$n_events)
  expect_output(print(s), "Transitions")
})


test_that("diagnostics record what a bare data frame cannot", {
  result <- simulate(obj_model(), seed = 6L)
  d <- result$diagnostics
  expect_type(d$n_unreliable, "integer")
  expect_type(result$meta$rng_kind, "character")
  expect_gte(d$n_individuals, 30L)
})


# ---------------------------------------------------------------------
# 6. Conversion needs no arguments
# ---------------------------------------------------------------------

test_that("as_long needs nothing, because the result knows its model", {
  # Compare convertToLongFormat(), which has to be told the state space,
  # the horizon, the absorbing states and the immigrants all over again.
  result <- simulate(obj_model(), seed = 8L)
  long <- as_long(result)

  expect_identical(
    names(long),
    c("ID", "birthDate", "Tstart", "Tstop", "statusEntry", "statusExit",
      "OD", "ns", "Episode", "sex", "fert")
  )
  expect_setequal(unique(long$ID), unique(result$pop$ID))
})


test_that("the state-variable columns carry the model's own names", {
  result <- simulate(obj_model(), seed = 9L)
  long <- as_long(result)
  expect_true(all(long$sex %in% c("m", "f")))
  expect_true(all(long$fert %in% c("0", "1")))
})


test_that("as_wide needs nothing either", {
  result <- simulate(obj_model(), seed = 10L)
  wide <- as_wide(result)
  expect_identical(nrow(wide), length(unique(result$pop$ID)))
})


test_that("a whole analysis is one pipeline", {
  # The point of the object API, in one expression.
  set.seed(12L)
  long <- obj_model() |>
    validate() |>
    simulate(seed = 12L) |>
    as_long()
  expect_gt(nrow(long), 0L)
})


# ---------------------------------------------------------------------
# 7. The two interfaces agree
# ---------------------------------------------------------------------

test_that("the object API and micSim() simulate the same model", {
  # THE test of this file. The same model expressed both ways, the same
  # seed, the same population down to the last row. The object API is a
  # different way of SAYING a model, not a different model.
  #
  # Transitions are declared in the order buildTransitionMatrix() lays out
  # the columns -- transient states before absorbing ones -- because that
  # order fixes which random numbers go to which destination.
  space <- expand.grid(sex = c("m", "f"), fert = c("0", "1"),
                       stringsAsFactors = FALSE)
  fertRates <- obj_fertility
  mortRates <- obj_mortality

  tm <- buildTransitionMatrix(
    allTransitions = cbind("f/0->f/1", "fertRates"),
    absTransitions = cbind("dead", "mortRates"),
    stateSpace     = space
  )

  set.seed(99L)
  legacy <- micSim(
    initPop          = obj_population(30L),
    transitionMatrix = tm,
    varInitStates    = rbind(c("m", "0"), c("f", "0")),
    initStatesProb   = c(0.515, 0.485),
    maxAge           = 100,
    simHorizon       = c(20000101, 20201231),
    fertTr           = "f/0->f/1",
    verbose          = FALSE
  )
  attr(legacy, "micsim_spec") <- NULL

  set.seed(99L)
  modern <- as.data.frame(simulate(obj_model(30L)))

  expect_identical(legacy, modern)
})


test_that("update() replaces the populations and checks them", {
  mort <- function(age, calTime) 0.0001 * exp(0.085 * age)
  pop <- data.frame(ID = 1:30, birthDate = "19700101",
                    initState = rep(c("m", "f"), 15))
  small <- micsim_model(state_space(sex = c("m", "f")),
                        transitions(transition(NA, "dead", mort)),
                        pop[1:5, ], c(20000101, 20101231))
  big <- update(small, population = pop)
  expect_s3_class(big, "micsim_model")
  expect_identical(nrow(big$population), 30L)
  expect_identical(big$registry, small$registry)
  bad <- pop; bad$ID[2] <- 1L
  expect_error(update(small, population = bad), "repeated",
               class = "micsim_engine_error")
  expect_error(update(small, max_age = 90), class = "micsim_model_error")
})


# ---------------------------------------------------------------------
# Factor and Date columns are read by their values
# ---------------------------------------------------------------------

test_that("factor identifiers keep their values", {
  # A factor converted with as.numeric() gives its codes: IDs 10, 20, 30
  # would become 1, 2, 3.
  population <- obj_population(3L)
  population$ID <- factor(c("30", "10", "20"))
  model <- obj_model(n = 3L)
  model <- update(model, population = population)
  expect_identical(model$population$ID, c("30", "10", "20"))
  pop <- as.data.frame(simulate(model, seed = 1L, reestimate = FALSE))
  expect_setequal(unique(pop$ID[is.na(pop$motherID)]), c(10L, 20L, 30L))
  # Same individuals, same streams, same lives as with numeric IDs.
  numeric_ids <- update(model, population = transform(population,
                                                      ID = c(30L, 10L, 20L)))
  expect_identical(pop, as.data.frame(simulate(numeric_ids, seed = 1L,
                                               reestimate = FALSE)))
})


test_that("Date columns are read as yyyymmdd dates", {
  population <- obj_population(3L)
  as_text <- obj_model(n = 3L)
  with_dates <- update(as_text, population = transform(
    population, birthDate = as.Date("1985-01-01")))
  expect_identical(with_dates$population$birthDate, rep(19850101L, 3L))
  a <- as.data.frame(simulate(as_text, seed = 2L, reestimate = FALSE))
  b <- as.data.frame(simulate(with_dates, seed = 2L, reestimate = FALSE))
  expect_identical(b$birthDate, a$birthDate)
  expect_true(all(grepl("^[0-9]{8}$", b$birthDate)))
})


test_that("micSim() reads factor identifiers by their values", {
  population <- obj_population(3L)
  population$ID <- factor(c("30", "10", "20"))
  population$birthDate <- as.Date("1985-01-01")
  tm <- buildTransitionMatrix(
    allTransitions = cbind(c("f/0->f/1"), c("obj_fertility")),
    absTransitions = cbind("dead", "obj_mortality"),
    stateSpace = expand.grid(sex = c("m", "f"), fert = c("0", "1"))
  )
  set.seed(4L)
  pop <- micSim(initPop = population, transitionMatrix = tm,
                absStates = "dead", maxAge = 100,
                simHorizon = c(20000101, 20101231), verbose = FALSE)
  expect_setequal(unique(pop$ID), c(10L, 20L, 30L))
  expect_true(all(pop$birthDate == "19850101"))
})
