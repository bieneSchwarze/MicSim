# =====================================================================
# test-subsystems.R
# ---------------------------------------------------------------------
# How subsystems are added to a model: specified, compiled by
# add_subsystem(), instantiated afresh in every run, and reported in
# result$extensions. These are the properties the future subsystem
# library will lean on, so they are tested at model level.
#
#   1. adding a subsystem
#   2. what adding refuses
#   3. per-run state and extensions
#   4. parallel runs
# =====================================================================


sub_mort  <- function(age, calTime) 0.00005 * exp(0.09 * age)
sub_fert  <- function(age, calTime) ifelse(age > 15 & age < 45, 0.15, 0)
sub_enrol <- function(age) ifelse(age == 7, Inf, 0)

sub_model <- function(n = 40L, ids = seq_len(n)) {
  micsim_model(
    state_space = state_space(sex = c("male", "female"),
                              edu = c("none", "primary"),
                              parity = c("0", "1")),
    transitions = transitions(
      transition(move("edu", "none", "primary"), sub_enrol),
      transition(move("parity", "0", "1", where = sex == "female"), sub_fert),
      transition(NA, "deceased", sub_mort)
    ),
    population = data.frame(
      ID = ids, birthDate = rep(c("19850101", "19960401"), length.out = n),
      initState = rep(c("female/primary/0", "female/none/0"), length.out = n),
      stringsAsFactors = FALSE
    ),
    horizon = c(20000101, 20201231)
  )
}

sub_fertility <- function(...) {
  fertility(births = moves_of("parity"),
            newborn = newborn_rule(sex = c(male = 0.515, female = 0.485),
                                   edu = "none", parity = "0"), ...)
}


# ---------------------------------------------------------------------
# 1. Adding a subsystem
# ---------------------------------------------------------------------

test_that("subsystems are added in a pipe and listed by the model", {
  model <- sub_model() |>
    add_subsystem(sub_fertility()) |>
    add_subsystem(school_enrolment(move("edu", "none", "primary"), month = 9))

  expect_s3_class(model, "micsim_model")
  expect_named(model$subsystems, c("fertility", "school_enrolment"))
  expect_output(print(model), "fertility \\(2 birth transition")
  expect_output(print(model), "school enrolment \\(4 transition")
})


test_that("a specification is pure data until it is added", {
  spec <- school_enrolment(move("edu", "none", "primary"), month = 9)
  expect_s3_class(spec, "micsim_subsystem")
  expect_output(print(spec), "school_enrolment")
  # the same specification can be added to two different models
  a <- add_subsystem(sub_model(), spec)
  b <- add_subsystem(sub_model(10L), spec)
  expect_identical(a$subsystems$school_enrolment$claims,
                   b$subsystems$school_enrolment$claims)
})


test_that("adding a subsystem consumes no random numbers", {
  set.seed(1L)
  before <- .Random.seed
  invisible(add_subsystem(sub_model(), sub_fertility()))
  expect_identical(.Random.seed, before)
})


test_that("the motherID column and age-zero validation follow the subsystem", {
  plain <- sub_model()
  with_births <- add_subsystem(plain, sub_fertility())
  expect_false("motherID" %in% names(as.data.frame(simulate(plain, seed = 1L))))
  expect_true("motherID" %in% names(as.data.frame(simulate(with_births, seed = 1L))))
  expect_identical(min(reachable_grid(with_births)$ages), 0L)
})


# ---------------------------------------------------------------------
# 2. What adding refuses
# ---------------------------------------------------------------------

test_that("two subsystems with one name are refused", {
  model <- add_subsystem(sub_model(), sub_fertility())
  expect_error(add_subsystem(model, sub_fertility()),
               class = "micsim_subsystem_error")
  # under its own name a second one is a different matter: here it would
  # claim the same births, which is refused for its own reason
  expect_error(add_subsystem(model, sub_fertility(name = "twins")),
               "same purpose", class = "micsim_subsystem_error")
})


test_that("two calendar rules on one transition are refused", {
  model <- add_subsystem(sub_model(),
                         school_enrolment(move("edu", "none", "primary"), 9))
  expect_error(
    add_subsystem(model, school_enrolment(move("edu", "none", "primary"), 1,
                                          name = "other")),
    class = "micsim_subsystem_error"
  )
})


test_that("only subsystems can be added, and only to a model", {
  expect_error(add_subsystem(sub_model(), list()),
               class = "micsim_subsystem_error")
  expect_error(add_subsystem(list(), sub_fertility()),
               class = "micsim_subsystem_error")
})


test_that("a subsystem's mistakes are reported when it is added", {
  expect_error(add_subsystem(sub_model(),
                             school_enrolment(move("edu", "primary", "none"), 9)),
               class = "micsim_selector_error")
  expect_error(add_subsystem(sub_model(), fertility(
    births = moves_of("parity"),
    newborn = newborn_rule(sex = c(male = 0.5, female = 0.5), parity = "0"))),
    class = "micsim_newborn_rule_error")                # edu missing
})


# ---------------------------------------------------------------------
# 3. Per-run state and extensions
# ---------------------------------------------------------------------

test_that("each subsystem reports into result$extensions, run by run", {
  model <- sub_model() |>
    add_subsystem(sub_fertility()) |>
    add_subsystem(school_enrolment(move("edu", "none", "primary"), month = 9))

  a <- simulate(model, seed = 3L)
  b <- simulate(model, seed = 3L)
  expect_named(a$extensions, c("fertility", "school_enrolment"))
  expect_identical(a$extensions, b$extensions)       # nothing carried over

  pop <- a$pop
  n_born <- length(unique(pop$ID[!is.na(pop$motherID)]))
  expect_identical(a$extensions$fertility$n_births, n_born)
  enrolled <- !is.na(pop$To) & grepl("/primary/", pop$To) &
    grepl("/none/", pop$From)
  expect_identical(a$extensions$school_enrolment$n_enrolments, sum(enrolled))
  expect_output(print(a), "extensions")
})


test_that("the object API and micSim() still simulate the same model", {
  # With both subsystems on, through the two front ends.
  ss <- expand.grid(sex = c("m", "f"), edu = c("no", "low"),
                    parity = c("0", "1"), stringsAsFactors = FALSE)
  enrolRates <- sub_enrol
  fertRates  <- sub_fert
  mortRates  <- sub_mort
  tm <- buildTransitionMatrix(
    rbind(c("no->low", "enrolRates"), c("f/0->f/1", "fertRates")),
    cbind("dead", "mortRates"), ss
  )
  pop <- data.frame(ID = 1:40,
                    birthDate = rep(c("19850101", "19960401"), 20),
                    initState = rep(c("f/low/0", "f/no/0"), 20),
                    stringsAsFactors = FALSE)
  set.seed(21L)
  legacy <- micSim(pop, transitionMatrix = tm,
                   varInitStates = rbind(c("m", "no", "0"), c("f", "no", "0")),
                   initStatesProb = c(0.515, 0.485), maxAge = 100,
                   simHorizon = c(20000101, 20201231), fertTr = "f/0->f/1",
                   monthSchoolEnrol = 9, verbose = FALSE)
  attr(legacy, "micsim_spec") <- NULL

  # transitions declared in the matrix's row-by-row, column-by-column order
  trs <- list()
  for (i in seq_len(nrow(tm))) for (j in which(tm[i, ] != "0")) {
    trs[[length(trs) + 1L]] <- transition(rownames(tm)[i], colnames(tm)[j],
                                          get(tm[i, j]))
  }
  model <- micsim_model(
    state_space(sex = c("m", "f"), edu = c("no", "low"), parity = c("0", "1")),
    do.call(transitions, trs), pop, c(20000101, 20201231), max_age = 100
  ) |>
    add_subsystem(fertility(moves_of("parity"),
                            newborn_rule(sex = c(m = 0.515, f = 0.485),
                                         edu = "no", parity = "0"))) |>
    add_subsystem(school_enrolment(move("edu", "no", "low"), 9))
  set.seed(21L)
  modern <- as.data.frame(simulate(model))
  expect_identical(legacy, modern)
})


# ---------------------------------------------------------------------
# 4. Parallel runs
# ---------------------------------------------------------------------

test_that("parallel results merge the subsystems' reports and print", {
  skip_on_cran()
  model <- add_subsystem(sub_model(60L), sub_fertility())
  r <- simulate(model, seed = 5L, cores = 2L)
  n_born <- length(unique(r$pop$ID[!is.na(r$pop$motherID)]))
  expect_identical(as.integer(r$extensions$fertility$n_births), n_born)
  expect_type(r$diagnostics$n_unreliable, "integer")
  expect_output(print(r), "micsim_result")
})


test_that("parallel newborns never take a supplied identifier", {
  skip_on_cran()
  # IDs 21..80: renumbering from N + 1 = 61 would collide with 61..80.
  model <- add_subsystem(sub_model(60L, ids = 21:80), sub_fertility())
  r <- simulate(model, seed = 6L, cores = 2L)
  people <- unique(r$pop[, c("ID", "birthDate", "initState")])
  expect_false(anyDuplicated(people$ID) > 0L)
  newborn_ids <- unique(r$pop$ID[!is.na(r$pop$motherID)])
  expect_true(all(newborn_ids > 80L))
})


test_that("recombination renumbers newborns above every supplied identifier", {
  # The same property without starting workers.
  piece <- function(ids, born) {
    data.frame(ID = c(ids, born), birthDate = "19900101", initState = "f",
               From = NA, To = NA, transitionTime = NA, transitionAge = NA,
               motherID = c(rep(NA, length(ids)), rep(ids[1], length(born))))
  }
  work <- list(list(initPop = data.frame(ID = 21:30), immigrPop = NULL),
               list(initPop = data.frame(ID = 31:40), immigrPop = NULL))
  pop <- recombine_populations(list(piece(21:30, 41L), piece(31:40, 41L)),
                               work, n_supplied = 20L)
  expect_false(anyDuplicated(pop$ID) > 0L)
  expect_true(all(pop$ID[!is.na(pop$motherID)] > 40L))
})
