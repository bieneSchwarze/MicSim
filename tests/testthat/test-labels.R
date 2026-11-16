# =====================================================================
# test-labels.R
# ---------------------------------------------------------------------
# No label has a built-in meaning in the object API. These tests pin that
# down from the outside: a model written with the modeller's own labels
# behaves exactly like the same model written with MicSim 3.0.0's, the
# structure of the model -- not its vocabulary -- decides what is
# absorbing, and the classic interface keeps 3.0.0's conventions, and only
# it does.
#
#   1. absorbing states come from the structure
#   2. move() declares transitions
#   3. renaming every label changes nothing
#   4. the classic interface's conventions
#   5. fixes that came with the change
# =====================================================================


lab_mort <- function(age, calTime) 0.00005 * exp(0.09 * age)
lab_emig <- function(age) rep(0.01, length(age))
lab_fert <- function(age, calTime) ifelse(age > 15 & age < 45, 0.15, 0)

lab_pop <- function(n = 40L, state = "f/0") {
  data.frame(ID = seq_len(n), birthDate = rep("19850101", n),
             initState = rep(state, n), stringsAsFactors = FALSE)
}


# ---------------------------------------------------------------------
# 1. Absorbing states come from the structure
# ---------------------------------------------------------------------

test_that("any destination that is not a state is absorbing, whatever its name", {
  model <- micsim_model(
    state_space = state_space(sex = c("m", "f"), fert = c("0", "1")),
    transitions = transitions(
      transition(NA, "tot",         lab_mort),
      transition(NA, "ausgewandert", lab_emig),
      transition(NA, "institution", lab_emig)
    ),
    population = lab_pop(), horizon = c(20000101, 20201231)
  )
  expect_identical(model$codec$absorbing,
                   c("tot", "ausgewandert", "institution"))
  expect_true(all(is.na(model$exits)))
  long <- as_long(simulate(model, seed = 1L))
  expect_true(all(long$OD %in% c("cens", "tot", "ausgewandert", "institution")))
})


test_that("a model may have no absorbing state at all", {
  model <- micsim_model(
    state_space = state_space(sex = c("m", "f"), fert = c("0", "1")),
    transitions = transitions(transition("f/0", "f/1", lab_fert)),
    population = lab_pop(), horizon = c(20000101, 20201231), max_age = 60
  )
  expect_length(model$codec$absorbing, 0L)
  expect_output(print(model), "none")
  result <- simulate(model, seed = 2L)
  expect_gt(result$diagnostics$n_events, 0L)
  expect_true(all(as_long(result)$OD %in% c("cens", "0->1")))
})


test_that("exits say what kind of exit an absorbing state is", {
  model <- micsim_model(
    state_space = state_space(sex = c("m", "f")),
    transitions = transitions(transition(NA, "deceased", lab_mort),
                              transition(NA, "abroad", lab_emig)),
    population = lab_pop(state = "f"), horizon = c(20000101, 20101231),
    exits = c(deceased = "death", abroad = "emigration")
  )
  expect_identical(model$exits, c(deceased = "death", abroad = "emigration"))
  expect_output(print(model), "deceased \\(death\\)")

  build <- function(exits) {
    micsim_model(state_space(sex = c("m", "f")),
                 transitions(transition(NA, "deceased", lab_mort)),
                 lab_pop(state = "f"), c(20000101, 20101231), exits = exits)
  }
  expect_error(build(c(dead = "death")), class = "micsim_model_error")
  expect_error(build(c(deceased = "passed away")), class = "micsim_model_error")
  expect_error(build("death"), class = "micsim_model_error")
})


test_that("a destination that looks like a mistyped state is refused", {
  build <- function(to) {
    micsim_model(state_space(sex = c("m", "f"), fert = c("0", "1")),
                 transitions(transition("f/0", to, lab_fert)),
                 lab_pop(), c(20000101, 20101231))
  }
  expect_error(build("f/2"), class = "micsim_model_error")   # not a state
  expect_error(build("1"), class = "micsim_model_error")     # a value, not a state
})


test_that("state values may not contain the separator", {
  expect_error(state_space(nat = c("DE/AT", "IT")), class = "micsim_model_error")
  expect_error(state_space(nat = c("DE", "")), class = "micsim_model_error")
})


test_that("the arguments of the development API are refused with advice", {
  build <- function(...) {
    micsim_model(state_space(sex = c("m", "f")),
                 transitions(transition(NA, "dead", lab_mort)),
                 lab_pop(state = "f"), c(20000101, 20101231), ...)
  }
  expect_error(build(absorbing = "dead"), "exits", class = "micsim_model_error")
  expect_error(build(fertility = "f/0->f/1"), "add_subsystem",
               class = "micsim_model_error")
  expect_error(build(school_enrolment = 9), "school_enrolment\\(",
               class = "micsim_model_error")
  expect_error(build(maxage = 9), "Unused", class = "micsim_model_error")
})


# ---------------------------------------------------------------------
# 2. move() declares transitions
# ---------------------------------------------------------------------

test_that("move() expands to one transition per matching origin", {
  model <- micsim_model(
    state_space = state_space(sex = c("m", "f"), nat = c("DE", "AT"),
                              fert = c("0", "1")),
    transitions = transitions(
      transition(move("fert", "0", "1", where = sex == "f"), lab_fert),
      transition(move("nat", "DE", "AT"), lab_emig),
      transition(NA, "dead", lab_mort)
    ),
    population = lab_pop(state = "f/DE/0"), horizon = c(20000101, 20101231)
  )
  labels <- vapply(model$registry$edges, `[[`, character(1L), "label")
  expect_setequal(labels[1:2], c("f/DE/0 -> f/DE/1", "f/AT/0 -> f/AT/1"))
  expect_identical(sum(grepl("/DE/", labels) & grepl("/AT/", labels)), 4L)
})


test_that("a misused selector is a clear error", {
  ss <- state_space(sex = c("m", "f"), fert = c("0", "1"))
  build <- function(tr) {
    micsim_model(ss, transitions(tr), lab_pop(), c(20000101, 20101231))
  }
  expect_error(build(transition(move("parity", "0", "1"), lab_fert)),
               class = "micsim_selector_error")             # no such variable
  expect_error(build(transition(move("fert", "0", "2"), lab_fert)),
               class = "micsim_selector_error")             # no such value
  expect_error(build(transition(move("fert", "0", "1", where = sex == "x"),
                                lab_fert)),
               class = "micsim_selector_error")             # matches no state
  expect_error(transition(moves_of("fert"), lab_fert),
               class = "micsim_model_error")                # cannot declare
  expect_error(move("fert", "0", "0"), class = "micsim_selector_error")
})


test_that("a transition declared twice is refused", {
  expect_error(
    micsim_model(
      state_space(sex = c("m", "f"), fert = c("0", "1")),
      transitions(transition("f/0", "f/1", lab_fert),
                  transition(move("fert", "0", "1"), lab_fert)),
      lab_pop(), c(20000101, 20101231)
    ),
    class = "micsim_registry_error"
  )
})


# ---------------------------------------------------------------------
# 3. Renaming every label changes nothing
# ---------------------------------------------------------------------

test_that("a model in the modeller's own labels simulates identically", {
  classic <- micsim_model(
    state_space = state_space(sex = c("m", "f"), fert = c("0", "1")),
    transitions = transitions(transition("f/0", "f/1", lab_fert),
                              transition(NA, "dead", lab_mort),
                              transition(NA, "rest", lab_emig)),
    population = lab_pop(60L), horizon = c(20000101, 20251231)
  ) |>
    add_subsystem(fertility("f/0->f/1",
                            newborn_rule(sex = c(m = 0.515, f = 0.485),
                                         fert = "0")))

  own <- micsim_model(
    state_space = state_space(geschlecht = c("mann", "frau"),
                              kinder = c("keine", "eins")),
    transitions = transitions(
      transition(move("kinder", "keine", "eins", where = geschlecht == "frau"),
                 lab_fert),
      transition(NA, "tot", lab_mort),
      transition(NA, "weg", lab_emig)
    ),
    population = lab_pop(60L, state = "frau/keine"),
    horizon = c(20000101, 20251231),
    exits = c(tot = "death", weg = "emigration")
  ) |>
    add_subsystem(fertility(moves_of("kinder"),
                            newborn_rule(geschlecht = c(mann = 0.515, frau = 0.485),
                                         kinder = "keine"),
                            mothers = geschlecht == "frau"))

  a <- as.data.frame(simulate(classic, seed = 9L))
  b <- as.data.frame(simulate(own, seed = 9L))
  rename <- c(m = "mann", f = "frau", "0" = "keine", "1" = "eins",
              dead = "tot", rest = "weg")
  relabel <- function(x) {
    out <- vapply(strsplit(x, "/", fixed = TRUE), function(p)
      paste(rename[p], collapse = "/"), character(1L))
    out[is.na(x)] <- NA_character_
    out
  }
  expect_identical(b$ID, a$ID)
  expect_identical(b$motherID, a$motherID)
  expect_identical(b$transitionTime, a$transitionTime)
  expect_identical(b$initState, relabel(a$initState))
  expect_identical(b$To, relabel(a$To))
})


# ---------------------------------------------------------------------
# 4. The classic interface's conventions
# ---------------------------------------------------------------------

test_that("buildTransitionMatrix accepts any absorbing label, keeping 3.0.0's order", {
  ss <- expand.grid(sex = c("m", "f"), stringsAsFactors = FALSE)
  tm <- buildTransitionMatrix(
    NULL,
    rbind(c("institution", "r1"), c("rest", "r2"), c("f/dead", "r3"),
          c("m/dead", "r4")),
    ss
  )
  # dead and rest first, as in 3.0.0; everything else after, in order named
  expect_identical(colnames(tm), c("m", "f", "dead", "rest", "institution"))
  expect_identical(unname(tm["f", "dead"]), "r3")
  expect_identical(unname(tm["m", "dead"]), "r4")
  expect_identical(unname(tm["m", "institution"]), "r1")
})


test_that("an absorbing label that cannot be read one way is refused", {
  ss <- expand.grid(sex = c("m", "f"), fert = c("0", "1"),
                    stringsAsFactors = FALSE)
  expect_error(buildTransitionMatrix(NULL, cbind("f/0", "r"), ss),
               class = "micsim_transition_matrix_error")    # no absorbing part
  expect_error(buildTransitionMatrix(NULL, cbind("f/dead/gone", "r"), ss),
               class = "micsim_transition_matrix_error")    # two of them
})


test_that("micSim() runs a matrix with custom absorbing states", {
  ss <- expand.grid(sex = c("m", "f"), stringsAsFactors = FALSE)
  mortRates <- lab_mort
  emigRates <- lab_emig
  tm <- buildTransitionMatrix(NULL, rbind(c("dead", "mortRates"),
                                          c("abroad", "emigRates")), ss)
  set.seed(4L)
  pop <- micSim(lab_pop(50L, "f"), transitionMatrix = tm, maxAge = 100,
                simHorizon = c(20000101, 20301231), verbose = FALSE)
  expect_true(all(pop$To[!is.na(pop$To)] %in% c("dead", "abroad")))
  expect_identical(attr(pop, "micsim_spec")$absStates, c("dead", "abroad"))
})


test_that("classic fertility patterns that match nothing are refused", {
  ss <- expand.grid(sex = c("m", "f"), fert = c("0", "1"),
                    stringsAsFactors = FALSE)
  fertRates <- lab_fert
  mortRates <- lab_mort
  tm <- buildTransitionMatrix(cbind("f/0->f/1", "fertRates"),
                              cbind("dead", "mortRates"), ss)
  expect_error(
    micSim(lab_pop(5L), transitionMatrix = tm, maxAge = 100,
           varInitStates = rbind(c("m", "0"), c("f", "0")),
           initStatesProb = c(0.5, 0.5),
           simHorizon = c(20000101, 20101231), fertTr = "m/0->m/1",
           verbose = FALSE),
    class = "micsim_compat_error"
  )
})


# ---------------------------------------------------------------------
# 5. Fixes that came with the change
# ---------------------------------------------------------------------

test_that("an absorbing transition sees duration Inf in a one-variable model", {
  # The codec, not the shape of the label, says a state is absorbing. A
  # one-variable model used to compare "m" with "dead" component-wise and
  # hand mortality the time since the sex variable last changed.
  codec <- new_state_codec(list(sex = c("m", "f")), absorbing = "dead")
  registry <- new_rate_registry(
    list(list(from = NA, to = "dead",
              rate = function(age, duration) rep(0.01, length(age)))),
    codec
  )
  expect_true(all(is.na(vapply(registry$edges, `[[`, integer(1L),
                               "affected_var"))))
  expect_true(all(lengths(lapply(registry$edges, `[[`, "changed_vars")) == 0L))
})


test_that("the final censored episode carries the state actually entered", {
  model <- micsim_model(
    state_space(a = c("x", "y"), b = c("p", "q")),
    transitions(transition("x/p", "y/q", function(age) rep(1, length(age)))),
    data.frame(ID = 1L, birthDate = "19900101", initState = "x/p"),
    c(20100101, 20201231)
  )
  long <- as_long(simulate(model, seed = 1L))
  expect_identical(nrow(long), 2L)
  expect_identical(long$a, c("x", "y"))
  expect_identical(long$b, c("p", "q"))
})
