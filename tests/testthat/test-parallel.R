# =====================================================================
# test-parallel.R
# ---------------------------------------------------------------------
# micSimParallel() ships the model to fresh worker processes. Objects
# that live in the master's GLOBAL environment are not serialised with
# closures (only non-global environments are), so before the fix:
#   1. arguments given as expressions over global objects were lazy
#      promises evaluated on a worker that had never seen those objects;
#   2. rate functions that micSim() finds by inheritance from the caller
#      were exported only if they sat in the caller's own frame.
# Both failed whenever micSimParallel() was called from inside a function.
# Workers call the INSTALLED MicSim, so these tests need it installed.
# =====================================================================

skip_parallel <- function() {
  skip_on_cran()
  skip_if_not(requireNamespace("MicSim", quietly = TRUE) &&
                utils::packageVersion("MicSim") >= "4.0.0",
              "installed MicSim >= 4.0.0 is required by the workers")
  skip_if(parallel::detectCores() < 2L, "needs two cores")
}

parallel_fixture <- function() {
  space <- c("m", "f")
  attr(space, "name") <- "sex"
  list(
    tm = buildTransitionMatrix(NULL, cbind("dead", "parMortRate"), space),
    init_pop = data.frame(
      ID = 1:40,
      birthDate = rep(c("19501231", "19700403", "19851015", "19991111"), 10),
      initState = rep(c("f", "m"), 20),
      stringsAsFactors = FALSE
    )
  )
}

# Put objects into the global environment for the duration of one test.
with_globals <- function(objects, code) {
  g <- globalenv()
  had <- vapply(names(objects), exists, logical(1L), envir = g, inherits = FALSE)
  old <- mget(names(objects)[had], envir = g)
  on.exit({
    rm(list = names(objects)[!had], envir = g)
    for (nm in names(old)) assign(nm, old[[nm]], envir = g)
  }, add = TRUE)
  for (nm in names(objects)) assign(nm, objects[[nm]], envir = g)
  force(code)
}


test_that("global rates and argument expressions work from inside a function", {
  skip_parallel()
  fx <- parallel_fixture()
  runner <- function(tm, init_pop) {
    micSimParallel(initPop = init_pop, transitionMatrix = tm,
                   absStates = "dead", maxAge = pick_max_age(),
                   simHorizon = c(20000101, 20201231),
                   cores = 2L, seeds = 7L, verbose = FALSE)
  }
  environment(runner) <- globalenv()

  pop <- with_globals(
    list(parMortRate  = function(age, calTime) 0.00003 * exp(0.1 * age),
         pick_max_age = function() 100),
    runner(fx$tm, fx$init_pop)
  )
  expect_setequal(unique(pop$ID), fx$init_pop$ID)
})


test_that("rates found by inheritance from the caller reach the workers", {
  skip_parallel()
  fx <- parallel_fixture()
  outer <- new.env(parent = globalenv())
  assign("parMortRate", function(age, calTime) 0.00003 * exp(0.1 * age),
         envir = outer)
  runner <- function(tm, init_pop) {
    micSimParallel(initPop = init_pop, transitionMatrix = tm,
                   absStates = "dead", maxAge = 100,
                   simHorizon = c(20000101, 20201231),
                   cores = 2L, seeds = 7L, verbose = FALSE)
  }
  environment(runner) <- outer

  pop <- runner(fx$tm, fx$init_pop)
  expect_setequal(unique(pop$ID), fx$init_pop$ID)
})


test_that("micSimParallel() gives the same life courses on any number of cores", {
  # Every worker runs with `seeds`, and every individual draws from their
  # own stream, so the split changes nothing.
  skip_parallel()
  fx <- parallel_fixture()
  run <- function(cores) {
    with_globals(
      list(parMortRate = function(age, calTime) 0.002 * exp(0.05 * age)),
      micSimParallel(initPop = fx$init_pop, transitionMatrix = fx$tm,
                     absStates = "dead", maxAge = 100,
                     simHorizon = c(20000101, 20201231),
                     cores = cores, seeds = 7L, verbose = FALSE)
    )
  }
  one <- run(1L)
  two <- run(2L)
  expect_true(any(!is.na(one$To)))
  attr(one, "micsim_spec") <- attr(two, "micsim_spec") <- NULL
  expect_identical(two, one)
})


test_that("micSimParallel() keeps immigrants of a worker without a starting population", {
  skip_parallel()
  fx <- parallel_fixture()
  immigrants <- data.frame(ID = 101:104, immigrDate = "20050101",
                           birthDate = "19800101", immigrInitState = "m",
                           stringsAsFactors = FALSE)
  pop <- with_globals(
    list(parMortRate = function(age, calTime) 0.002 * exp(0.05 * age)),
    micSimParallel(initPopList = list(fx$init_pop, fx$init_pop[0L, ]),
                   immigrPopList = list(NULL, immigrants),
                   transitionMatrix = fx$tm, absStates = "dead",
                   maxAge = 100, simHorizon = c(20000101, 20201231),
                   cores = 2L, seeds = 7L, verbose = FALSE)
  )
  expect_setequal(unique(pop$ID), c(fx$init_pop$ID, immigrants$ID))
})


test_that("micSimParallel() refuses an unusable number of cores", {
  # Checked by the same rule as simulate(cores =), before any worker
  # starts. A fractional value used to be truncated silently.
  fx <- parallel_fixture()
  for (bad in list(1.5, 0, -1, "a", NA, c(1, 2))) {
    expect_error(
      micSimParallel(initPop = fx$init_pop, transitionMatrix = fx$tm,
                     absStates = "dead", maxAge = 100,
                     simHorizon = c(20000101, 20201231),
                     cores = bad, seeds = 7L, verbose = FALSE),
      class = "micsim_parallel_error"
    )
  }
})
