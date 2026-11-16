# =====================================================================
# test-compat.R
# ---------------------------------------------------------------------
# The compatibility layer is what makes a script written in 2012 still
# run. Its job is translation, not simulation, so most of these tests are
# about whether the translation is faithful: does the reconstructed state
# space match the matrix, are destinations declared in the order that
# fixes the random stream, is a partial fertility pattern expanded the way
# MicSim 3.0.0 expands it.
#
# The end-to-end tests at the bottom run whole models through the legacy
# interface, which is the only way to check that the pieces fit.
#
# Tests are grouped by what they establish:
#   1. buildTransitionMatrix: absorbing states
#   2. buildTransitionMatrix: transient patterns
#   3. reconstructing the state space
#   4. resolving rate functions
#   5. destination order (the parity contract)
#   6. expanding fertility patterns
#   7. micSim end to end
# =====================================================================


# Rate functions live at the top level of this file, which is also where
# micSim() will look for them: the compatibility layer resolves names in
# the environment it was called from.
mortRates <- function(age, calTime) 0.00003 * exp(0.1 * age)
fertRates <- function(age, calTime) ifelse(age > 15 & age < 45, 0.15, 0)


# ---------------------------------------------------------------------
# 1. buildTransitionMatrix: absorbing states
# ---------------------------------------------------------------------

test_that("a bare absorbing label applies to every transient state", {
  space <- expand.grid(sex = c("m", "f"), fert = c("0", "1"))
  tm <- buildTransitionMatrix(NULL, cbind("dead", "mortRates"), space)

  expect_identical(dim(tm), c(4L, 5L))
  expect_identical(rownames(tm), c("m/0", "f/0", "m/1", "f/1"))
  expect_identical(colnames(tm)[5L], "dead")
  expect_true(all(tm[, "dead"] == "mortRates"))
})


test_that("a compound absorbing label applies only to matching states", {
  # "f/dead" means only women can reach it.
  space <- expand.grid(sex = c("m", "f"), fert = c("0", "1"))
  tm <- buildTransitionMatrix(NULL, cbind("f/dead", "mortRates"), space)

  expect_true(all(tm[c("f/0", "f/1"), "dead"] == "mortRates"))
  expect_true(all(tm[c("m/0", "m/1"), "dead"] == "0"))
})


test_that("the emigration state appears only when a model uses it", {
  space <- expand.grid(sex = c("m", "f"))
  without <- buildTransitionMatrix(NULL, cbind("dead", "mortRates"), space)
  with_it <- buildTransitionMatrix(
    NULL, rbind(c("dead", "mortRates"), c("rest", "emigRates")), space
  )

  expect_false("rest" %in% colnames(without))
  expect_identical(colnames(with_it), c("m", "f", "dead", "rest"))
})


test_that("absorbing states come after the transient ones", {
  # Not cosmetic. The column order is the order in which MicSim 3.0.0
  # considers destinations, and therefore the order in which it consumes
  # random numbers.
  space <- expand.grid(sex = c("m", "f"), fert = c("0", "1"))
  tm <- buildTransitionMatrix(cbind("f/0->f/1", "fertRates"),
                              cbind("dead", "mortRates"), space)
  expect_identical(colnames(tm)[seq_len(4L)], rownames(tm))
  expect_identical(colnames(tm)[5L], "dead")
})


# ---------------------------------------------------------------------
# 2. buildTransitionMatrix: transient patterns
# ---------------------------------------------------------------------

test_that("a partial pattern expands over the variables it does not name", {
  # "0->1" is a first birth, for every sex and nationality.
  space <- expand.grid(sex = c("m", "f"), nat = c("DE", "AT"),
                       fert = c("0", "1"))
  tm <- buildTransitionMatrix(cbind("0->1", "fertRates"),
                              cbind("dead", "mortRates"), space)

  expect_identical(tm["f/DE/0", "f/DE/1"], "fertRates")
  expect_identical(tm["m/AT/0", "m/AT/1"], "fertRates")
  # nationality may not change on the way
  expect_identical(tm["f/DE/0", "f/AT/1"], "0")
})


test_that("a fully specified pattern expands to exactly one cell", {
  space <- expand.grid(sex = c("m", "f"), nat = c("DE", "AT"),
                       fert = c("0", "1"))
  tm <- buildTransitionMatrix(cbind("f/DE/0->f/DE/1", "fertRates"),
                              cbind("dead", "mortRates"), space)

  expect_identical(sum(tm == "fertRates"), 1L)
  expect_identical(tm["f/DE/0", "f/DE/1"], "fertRates")
})


test_that("a single transition may be given as a bare vector", {
  space <- c("no", "low")
  tm <- buildTransitionMatrix(c("no->low", "enrolRates"),
                              c("dead", "mortRates"), space)
  expect_identical(tm["no", "low"], "enrolRates")
})


test_that("a model may have no transient transitions at all", {
  space <- c("m", "f")
  tm <- buildTransitionMatrix(NULL, cbind("dead", "mortRates"), space)
  expect_identical(dim(tm), c(2L, 3L))
  expect_true(all(tm[, c("m", "f")] == "0"))
})


test_that("a matrix without absorbing transitions is refused", {
  space <- c("m", "f")
  expect_error(buildTransitionMatrix(NULL, NULL, space),
               class = "micsim_transition_matrix_error")
})


# ---------------------------------------------------------------------
# 3. Reconstructing the state space
# ---------------------------------------------------------------------

test_that("the codec is rebuilt from the matrix's dimension names", {
  space <- expand.grid(sex = c("m", "f"), nat = c("DE", "AT"))
  tm <- buildTransitionMatrix(NULL, cbind("dead", "mortRates"), space)
  codec <- codec_from_matrix(tm)

  expect_setequal(codec$transient, rownames(tm))
  expect_identical(codec$absorbing, "dead")
  expect_identical(length(codec$vars), 2L)
})


test_that("absorbing states are read from the columns, not from absStates", {
  # MicSim 3.0.0 took them from the `absStates` argument and assumed its
  # order matched the matrix's. A user who wrote c("rest", "dead") while
  # the matrix put "dead" first had every death recorded as an emigration,
  # silently (defect D17). Reading one source removes the possibility.
  space <- c("m", "f")
  tm <- buildTransitionMatrix(
    NULL, rbind(c("dead", "mortRates"), c("rest", "mortRates")), space
  )
  codec <- codec_from_matrix(tm)
  expect_identical(codec$absorbing, c("dead", "rest"))
})


test_that("a matrix with no absorbing state is refused", {
  tm <- matrix("0", 2L, 2L, dimnames = list(c("a", "b"), c("a", "b")))
  expect_error(codec_from_matrix(tm), class = "micsim_compat_error")
})


test_that("a matrix without dimension names is refused", {
  expect_error(codec_from_matrix(matrix("0", 2L, 3L)),
               class = "micsim_compat_error")
})


test_that("states with different numbers of variables are refused", {
  tm <- matrix("0", 2L, 3L,
               dimnames = list(c("m/0", "f"), c("m/0", "f", "dead")))
  expect_error(codec_from_matrix(tm), class = "micsim_compat_error")
})


# ---------------------------------------------------------------------
# 4. Resolving rate functions
# ---------------------------------------------------------------------

test_that("a rate name is resolved to the function itself", {
  space <- c("m", "f")
  tm <- buildTransitionMatrix(NULL, cbind("dead", "mortRates"), space)
  codec <- codec_from_matrix(tm)
  registry <- registry_from_matrix(tm, codec, environment())

  expect_length(registry$edges, 2L)          # one per transient state
  expect_identical(registry$edges[[1L]]$rate, mortRates)
})


test_that("one name gives one shared closure, not one per cell", {
  # A mortality function written once is one object, whatever the size of
  # the state space. MicSim 3.0.0 looked the name up again on every event.
  space <- expand.grid(sex = c("m", "f"), nat = c("DE", "AT", "IT"))
  tm <- buildTransitionMatrix(NULL, cbind("dead", "mortRates"), space)
  codec <- codec_from_matrix(tm)
  registry <- registry_from_matrix(tm, codec, environment())

  expect_length(registry$edges, 6L)
  for (edge in registry$edges) {
    expect_true(identical(edge$rate, registry$edges[[1L]]$rate))
  }
})


test_that("an unresolvable rate name is refused by name", {
  space <- c("m", "f")
  tm <- buildTransitionMatrix(NULL, cbind("dead", "noSuchFunction"), space)
  codec <- codec_from_matrix(tm)
  expect_error(registry_from_matrix(tm, codec, environment()),
               class = "micsim_compat_error")
})


test_that("a matrix in which nothing can happen is refused", {
  tm <- matrix("0", 2L, 3L,
               dimnames = list(c("m", "f"), c("m", "f", "dead")))
  codec <- codec_from_matrix(tm)
  expect_error(registry_from_matrix(tm, codec, environment()),
               class = "micsim_compat_error")
})


# ---------------------------------------------------------------------
# 5. Destination order
# ---------------------------------------------------------------------

test_that("destinations are declared in the matrix's column order", {
  # The registry keeps a state's destinations in the order of the
  # matrix's columns, as MicSim 3.0.0 walks them. The random numbers do
  # not depend on this order (each destination's uniform follows the
  # state space), but competing-risk ties and hook calls do.
  space <- expand.grid(sex = c("m", "f"), fert = c("0", "1"))
  tm <- buildTransitionMatrix(cbind("f/0->f/1", "fertRates"),
                              cbind("dead", "mortRates"), space)
  codec <- codec_from_matrix(tm)
  registry <- registry_from_matrix(tm, codec, environment())

  from_f0 <- registry$feasible(codec$encode("f/0"))
  destinations <- vapply(from_f0, function(e) e$to_label, character(1L))

  # "f/1" is a transient column and therefore comes before "dead"
  expect_identical(destinations, c("f/1", "dead"))
  expect_lt(which(colnames(tm) == "f/1"), which(colnames(tm) == "dead"))
})


# ---------------------------------------------------------------------
# 6. Expanding fertility patterns
# ---------------------------------------------------------------------

test_that("a partial fertility pattern expands to women only", {
  space <- expand.grid(sex = c("m", "f"), nat = c("DE", "AT"),
                       fert = c("0", "1"))
  tm <- buildTransitionMatrix(cbind("0->1", "fertRates"),
                              cbind("dead", "mortRates"), space)
  codec <- codec_from_matrix(tm)
  edges <- expand_fertility_transitions("0->1", codec)

  expect_setequal(edges, c("f/DE/0 -> f/DE/1", "f/AT/0 -> f/AT/1"))
  expect_false(any(grepl("^m/", edges)))
})


test_that("a birth may not change anything but the fertility variable", {
  space <- expand.grid(sex = c("m", "f"), nat = c("DE", "AT"),
                       fert = c("0", "1"))
  codec <- codec_from_matrix(
    buildTransitionMatrix(NULL, cbind("dead", "mortRates"), space)
  )
  edges <- expand_fertility_transitions("0->1", codec)
  expect_false("f/DE/0 -> f/AT/1" %in% edges)
})


test_that("fully specified fertility patterns pass through", {
  space <- expand.grid(sex = c("m", "f"), nat = c("DE", "AT"),
                       fert = c("0", "1"))
  codec <- codec_from_matrix(
    buildTransitionMatrix(NULL, cbind("dead", "mortRates"), space)
  )
  edges <- expand_fertility_transitions(
    c("f/DE/0->f/DE/1", "f/AT/0->f/AT/1"), codec
  )
  expect_setequal(edges, c("f/DE/0 -> f/DE/1", "f/AT/0 -> f/AT/1"))
})


test_that("a malformed fertility pattern is refused", {
  space <- c("m", "f")
  codec <- codec_from_matrix(
    buildTransitionMatrix(NULL, cbind("dead", "mortRates"), space)
  )
  expect_error(expand_fertility_transitions("not a pattern", codec),
               class = "micsim_compat_error")
})


# ---------------------------------------------------------------------
# 7. micSim end to end
# ---------------------------------------------------------------------

compat_init_pop <- function(n = 20L, birth = "19700101", state = "f/0") {
  data.frame(
    ID        = seq_len(n),
    birthDate = rep(birth, n),
    initState = rep(state, n),
    stringsAsFactors = FALSE
  )
}


test_that("the legacy interface runs a mortality-only model", {
  space <- c("m", "f")
  attr(space, "name") <- "sex"
  tm <- buildTransitionMatrix(NULL, cbind("dead", "mortRates"), space)

  init_pop <- data.frame(
    ID = 1:5,
    birthDate = c("19301231", "19990403", "19561015", "19911111", "19650101"),
    initState = c("f", "m", "f", "m", "m"),
    stringsAsFactors = FALSE
  )

  set.seed(234L)
  pop <- micSim(initPop = init_pop, transitionMatrix = tm,
                absStates = "dead", maxAge = 120,
                simHorizon = c(20000101, 21001231), verbose = FALSE)

  expect_identical(
    names(pop),
    c("ID", "birthDate", "initState", "From", "To", "transitionTime",
      "transitionAge")
  )
  expect_identical(nrow(pop), 5L)
  expect_identical(pop$birthDate, init_pop$birthDate[order(init_pop$ID)])
})


test_that("the legacy interface runs a model with fertility", {
  space <- expand.grid(sex = c("m", "f"), fert = c("0", "1"))
  tm <- buildTransitionMatrix(cbind("f/0->f/1", "fertRates"),
                              cbind("dead", "mortRates"), space)

  set.seed(1L)
  pop <- micSim(
    initPop          = compat_init_pop(30L),
    transitionMatrix = tm,
    absStates        = "dead",
    varInitStates    = rbind(c("m", "0"), c("f", "0")),
    initStatesProb   = c(0.515, 0.485),
    maxAge           = 100,
    simHorizon       = c(20000101, 20201231),
    fertTr           = "f/0->f/1",
    verbose          = FALSE
  )

  expect_identical(names(pop)[8L], "motherID")
  expect_gt(nrow(pop), 30L)
  expect_gt(sum(!is.na(pop$motherID)), 0L)
})


test_that("fertility without a newborn specification is refused", {
  space <- expand.grid(sex = c("m", "f"), fert = c("0", "1"))
  tm <- buildTransitionMatrix(cbind("f/0->f/1", "fertRates"),
                              cbind("dead", "mortRates"), space)
  expect_error(
    micSim(initPop = compat_init_pop(5L), transitionMatrix = tm,
           maxAge = 100, simHorizon = c(20000101, 20201231),
           fertTr = "f/0->f/1", verbose = FALSE),
    class = "micsim_compat_error"
  )
})


test_that("the legacy interface runs a model with immigration", {
  space <- c("m", "f")
  tm <- buildTransitionMatrix(NULL, cbind("dead", "mortRates"), space)

  immigrants <- data.frame(
    ID              = 101:103,
    immigrDate      = c(20050601L, 20080601L, 20120601L),
    birthDate       = c(19750601L, 19780601L, 19820601L),
    immigrInitState = c("f", "m", "f"),
    stringsAsFactors = FALSE
  )

  set.seed(2L)
  pop <- micSim(initPop = compat_init_pop(10L, state = "f"),
                immigrPop = immigrants, transitionMatrix = tm,
                maxAge = 100, simHorizon = c(20000101, 20201231),
                verbose = FALSE)

  expect_setequal(unique(pop$ID), c(1:10, 101:103))
})


test_that("micSim consumes the caller's random stream", {
  # It takes no seed, exactly as it never has, so set.seed() immediately
  # before the call is what makes a run reproducible.
  space <- c("m", "f")
  tm <- buildTransitionMatrix(NULL, cbind("dead", "mortRates"), space)
  run <- function() {
    micSim(initPop = compat_init_pop(10L, state = "f"),
           transitionMatrix = tm, maxAge = 100,
           simHorizon = c(20000101, 20401231), verbose = FALSE)
  }

  set.seed(99L); a <- run()
  set.seed(99L); b <- run()
  expect_identical(a, b)

  set.seed(100L); cc <- run()
  expect_false(identical(a, cc))
})


test_that("verbose = FALSE is silent and verbose = TRUE is not", {
  space <- c("m", "f")
  tm <- buildTransitionMatrix(NULL, cbind("dead", "mortRates"), space)
  call_it <- function(verbose) {
    micSim(initPop = compat_init_pop(5L, state = "f"),
           transitionMatrix = tm, maxAge = 100,
           simHorizon = c(20000101, 20101231), verbose = verbose)
  }

  set.seed(3L)
  expect_silent(call_it(FALSE))

  # capture_messages() rather than expect_message(), which consumes only
  # the first condition and lets the rest print to the console.
  set.seed(3L)
  said <- capture_messages(call_it(TRUE))
  expect_gt(length(said), 0L)
  expect_true(any(grepl("finished", said, fixed = TRUE)))
})
