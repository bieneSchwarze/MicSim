# =====================================================================
# test-rng-streams.R
# ---------------------------------------------------------------------
# Every individual draws from their own L'Ecuyer-CMRG stream
# (R/rng_streams.R). These tests pin what that buys:
#   1. the stream arithmetic agrees with R's `parallel` package
#   2. the caller's random-number state is left alone
#   3. common random numbers: what does not change between two scenarios
#      is simulated identically in both
#   4. the result does not depend on processing order, on who else is in
#      the population, or on the order transitions were declared in
#   5. the result does not depend on the number of cores
# =====================================================================


rs_mortality <- function(age, calTime) 0.0001 * exp(0.08 * age)
rs_upgrade <- function(age, calTime) ifelse(age >= 18 & age < 40, 0.08, 0)
rs_fertility <- function(age, calTime) ifelse(age >= 16 & age < 45, 0.12, 0)

# Each individual's attributes depend on their ID only, so that a subset
# of the population describes the same people.
rs_population <- function(n = 60L, ids = seq_len(n)) {
  k <- (ids - 1L) %% 4L + 1L
  data.frame(
    ID = ids,
    birthDate = c("19700101", "19800515", "19900303", "19951224")[k],
    initState = c("f/low", "m/low", "f/high", "m/low")[k],
    stringsAsFactors = FALSE
  )
}

rs_immigrants <- function(n = 8L, first_id = 1001L,
                          states = c("f/low", "m/high")) {
  data.frame(
    ID = first_id + seq_len(n) - 1L,
    immigrDate = rep(c("20030101", "20071001"), length.out = n),
    birthDate = rep(c("19850101", "19750601"), length.out = n),
    immigrInitState = rep(states, length.out = n),
    stringsAsFactors = FALSE
  )
}

# Education upgrade low -> high, mortality by state; `high_mortality`
# multiplies the mortality of the highly educated.
rs_model <- function(population = rs_population(), immigrants = NULL,
                     high_mortality = 1, reverse = FALSE, births = FALSE) {
  mort_low <- rs_mortality
  mort_high <- function(age, calTime) high_mortality * rs_mortality(age, calTime)
  declared <- list(
    transition(move("edu", "low", "high"), rs_upgrade),
    transition(c("f/low", "m/low"), "dead", mort_low),
    transition(c("f/high", "m/high"), "dead", mort_high)
  )
  if (reverse) declared <- rev(declared)
  model <- micsim_model(
    state_space = state_space(sex = c("m", "f"), edu = c("low", "high")),
    transitions = do.call(transitions, declared),
    population = population, immigrants = immigrants,
    horizon = c(20000101, 20301231), max_age = 100
  )
  model
}

# Fertility: women at parity 0 have one child.
rs_fertility_model <- function(population, immigrants = NULL,
                               male_mortality = 1) {
  mort_f <- rs_mortality
  mort_m <- function(age, calTime) male_mortality * rs_mortality(age, calTime)
  micsim_model(
    state_space = state_space(sex = c("m", "f"), parity = c("0", "1")),
    transitions = transitions(
      transition("f/0", "f/1", rs_fertility),
      transition(c("f/0", "f/1"), "dead", mort_f),
      transition(c("m/0", "m/1"), "dead", mort_m)
    ),
    population = population, immigrants = immigrants,
    horizon = c(20000101, 20401231), max_age = 100
  ) |>
    add_subsystem(fertility(
      births = "f/0->f/1",
      newborn = newborn_rule(sex = c(m = 0.515, f = 0.485), parity = "0")
    ))
}

rs_fertility_population <- function(n = 40L) {
  data.frame(
    ID = seq_len(n),
    birthDate = rep(c("19750101", "19820515", "19880303", "19920101"),
                    length.out = n),
    initState = rep(c("f/0", "m/0"), length.out = n),
    stringsAsFactors = FALSE
  )
}

pop_of <- function(model, seed = 11L, cores = 1L) {
  as.data.frame(simulate(model, seed = seed, cores = cores,
                         reestimate = FALSE))
}

rows_of <- function(pop, ids) {
  out <- pop[pop$ID %in% ids, , drop = FALSE]
  row.names(out) <- NULL
  out
}


# ---------------------------------------------------------------------
# 1. Stream arithmetic
# ---------------------------------------------------------------------

test_that("person streams are the nextRNGStream() streams numbered by ID", {
  set.seed(42L, kind = "L'Ecuyer-CMRG")
  master <- .Random.seed[-1L]
  set.seed(NULL, kind = "default")

  nth_stream <- function(n) {
    s <- c(10407L, master)
    for (i in seq_len(n)) s <- parallel::nextRNGStream(s)
    s[-1L]
  }
  # Stream number ID + 2^31, so these IDs own streams 1, 2, 5 and 300.
  ids <- -2^31 + c(1, 2, 5, 300)
  states <- person_stream_states(master, ids)
  for (i in seq_along(ids)) {
    expect_identical(states[i, ], nth_stream(ids[i] + 2^31))
  }
})


test_that("dense and sparse identifiers give the same streams", {
  set.seed(1L, kind = "L'Ecuyer-CMRG")
  master <- .Random.seed[-1L]
  set.seed(NULL, kind = "default")

  dense <- person_stream_states(master, 1:200)
  # Adding one far-away identifier switches to the per-identifier jump.
  sparse <- person_stream_states(master, c(1:200, 2000000000))
  expect_identical(sparse[1:200, ], dense)
  expect_identical(sparse[201L, ], person_stream_states(master, 2000000000)[1L, ])
  # A state is the same whoever else is asked for.
  expect_identical(person_stream_states(master, c(7, 3))[1L, ], dense[7L, ])
})


test_that("seed words of 2^31 and above survive R's signed storage", {
  words <- c(0, 1, 2^31 - 1, 2^31, 2^31 + 1, 2^32 - 1)
  expect_identical(int32_to_unsigned(unsigned_to_int32(words)), words)
  expect_identical(unsigned_to_int32(2^31), NA_integer_)
})


test_that("a newborn's state is a valid generator state", {
  # R silently reseeds from the clock if a state is invalid, so a valid
  # state is what makes newborns reproducible.
  for (u in list(c(1e-12, 1e-12, 1e-12, 0.5, 0.5, 0.5),
                 c(0.999999, 0.5, 0.25, 1e-12, 1e-12, 1e-12),
                 c(2^31 / 4294967087, 0.1, 0.2, 0.3, 0.4, 0.5))) {
    state <- state_from_uniforms(u)
    words <- int32_to_unsigned(state)
    expect_true(all(words[1:3] < 4294967087) && any(words[1:3] > 0))
    expect_true(all(words[4:6] < 4294944443) && any(words[4:6] > 0))
    a <- draw_from_state(state, 3L)
    b <- draw_from_state(state, 3L)
    expect_identical(a, b)
  }
})


test_that("drawing from a state leaves the ambient generator alone", {
  set.seed(9L)
  before <- .Random.seed
  set.seed(5L, kind = "L'Ecuyer-CMRG")
  state <- .Random.seed[-1L]
  set.seed(9L, kind = "default")
  draw_from_state(state, 5L)
  expect_identical(.Random.seed, before)
})


# ---------------------------------------------------------------------
# 2. The caller's random-number state
# ---------------------------------------------------------------------

test_that("a seeded run leaves the caller's state and generator kind alone", {
  old_kind <- RNGkind()
  on.exit({ RNGkind(old_kind[1L], old_kind[2L], old_kind[3L]); set.seed(NULL) },
          add = TRUE)
  RNGkind("Knuth-TAOCP-2002")
  set.seed(3L)
  before <- .Random.seed
  simulate(rs_model(), seed = 1L, reestimate = FALSE)
  expect_identical(.Random.seed, before)
  expect_identical(RNGkind()[1L], "Knuth-TAOCP-2002")
})


test_that("a seeded run does not create a random-number state that was absent", {
  old_kind <- RNGkind()
  on.exit({ RNGkind(old_kind[1L], old_kind[2L], old_kind[3L]); set.seed(NULL) },
          add = TRUE)
  RNGkind("Wichmann-Hill")
  if (exists(".Random.seed", envir = globalenv())) {
    rm(".Random.seed", envir = globalenv())
  }
  simulate(rs_model(), seed = 1L, reestimate = FALSE)
  expect_false(exists(".Random.seed", envir = globalenv()))
  expect_identical(RNGkind()[1L], "Wichmann-Hill")
})


test_that("the state is restored when a run fails", {
  switch <- new.env()
  switch$fail <- FALSE
  failing_rate <- function(age, calTime) {
    if (switch$fail) stop("boom")
    rep(0.01, length(age))
  }
  failing <- micsim_model(
    state_space = state_space(sex = c("m", "f")),
    transitions = transitions(transition(NA, "dead", failing_rate)),
    population = rs_population(5L)[, 1:2] |>
      transform(initState = "f"),
    horizon = c(20000101, 20101231), max_age = 100
  )
  switch$fail <- TRUE
  set.seed(4L)
  before <- .Random.seed
  expect_error(simulate(failing, seed = 1L, reestimate = FALSE))
  expect_identical(.Random.seed, before)
})


test_that("an unusable seed is refused", {
  model <- rs_model()
  for (bad in list("1", 1.5, NA_integer_, c(1L, 2L), Inf, 2^40)) {
    expect_error(simulate(model, seed = bad, reestimate = FALSE),
                 class = "micsim_model_error")
  }
})


# ---------------------------------------------------------------------
# 3. Common random numbers
# ---------------------------------------------------------------------

test_that("individuals a scenario does not touch live identical lives", {
  # The scenario raises the mortality of the highly educated only. Anyone
  # who never becomes highly educated has the same rates in both runs and
  # must have the same life course; anyone else must at least have the
  # same first event, which is drawn before the rates differ for them.
  base <- pop_of(rs_model())
  scenario <- pop_of(rs_model(high_mortality = 3))

  ever_high <- unique(base$ID[grepl("high", base$initState) |
                                grepl("high", base$To)])
  untouched <- setdiff(unique(base$ID), ever_high)
  expect_gt(length(untouched), 10L)
  expect_identical(rows_of(scenario, untouched), rows_of(base, untouched))

  start_low_moved <- unique(base$ID[grepl("low", base$initState) &
                                      base$To %in% c("m/high", "f/high")])
  first_of <- function(pop, ids) {
    p <- rows_of(pop, ids)
    p <- p[!duplicated(p$ID), c("ID", "To", "transitionTime")]
    row.names(p) <- NULL
    p
  }
  expect_gt(length(start_low_moved), 0L)
  expect_identical(first_of(scenario, start_low_moved),
                   first_of(base, start_low_moved))
  # And the scenario did change something.
  expect_false(identical(scenario, base))
})


test_that("births of mothers a scenario does not touch are unchanged", {
  # Raising men's mortality changes nothing for women: their births, and
  # the life courses of their daughters, are the same in both runs.
  pop <- rs_fertility_population()
  run <- function(m) simulate(rs_fertility_model(pop, male_mortality = m),
                              seed = 5L, reestimate = FALSE)
  base <- run(1)
  scenario <- run(4)
  expect_gt(base$extensions$fertility$n_births, 0L)
  expect_identical(scenario$extensions$fertility$births,
                   base$extensions$fertility$births)

  women <- unique(as.data.frame(base)$ID[as.data.frame(base)$initState %in%
                                           c("f/0", "f/1")])
  expect_identical(rows_of(as.data.frame(scenario), women),
                   rows_of(as.data.frame(base), women))
})


test_that("a forced transition does not move anybody else's draws", {
  # In the scenario a forced exit at age 70 is added to an otherwise
  # impossible transition. Whoever dies before 70 in the base run dies on
  # the same day in the scenario: the forced destination uses its own
  # uniform, and the others keep theirs.
  forced_model <- function(forced) {
    exit_rate <- if (forced) function(age, calTime) ifelse(age >= 70, Inf, 0)
                 else function(age, calTime) rep(0, length(age))
    micsim_model(
      state_space = state_space(sex = c("m", "f")),
      transitions = transitions(
        transition(NA, "dead", function(age, calTime) rep(0.08, length(age))),
        transition(NA, "gone", exit_rate)
      ),
      population = data.frame(ID = 1:80, birthDate = "19400101",
                              initState = rep(c("m", "f"), 40),
                              stringsAsFactors = FALSE),
      horizon = c(20000101, 20201231), max_age = 100
    )
  }
  base <- pop_of(forced_model(FALSE))
  scenario <- pop_of(forced_model(TRUE))
  died_young <- base$ID[base$To %in% "dead" & base$transitionAge < 70]
  expect_gt(length(died_young), 5L)
  expect_identical(rows_of(scenario, died_young), rows_of(base, died_young))
  expect_true(any(scenario$To %in% "gone"))
})


# ---------------------------------------------------------------------
# 4. Order and composition do not matter
# ---------------------------------------------------------------------

test_that("the order of the population's rows does not matter", {
  population <- rs_population()
  shuffled <- population[c(seq(2L, 60L, 2L), seq(1L, 59L, 2L)), ]
  row.names(shuffled) <- NULL
  expect_identical(pop_of(rs_model(shuffled)), pop_of(rs_model(population)))
})


test_that("an individual's life does not depend on who else is simulated", {
  full <- pop_of(rs_model(rs_population(60L),
                          immigrants = rs_immigrants()))
  half_ids <- seq(1L, 60L, 3L)
  part <- pop_of(rs_model(rs_population(ids = half_ids),
                          immigrants = rs_immigrants(4L)))
  keep <- c(half_ids, 1001:1004)
  expect_identical(rows_of(part, keep), rows_of(full, keep))
})


test_that("the order transitions are declared in does not matter", {
  expect_identical(pop_of(rs_model(reverse = TRUE)), pop_of(rs_model()))
})


test_that("micSim() and simulate() agree under the same caller seed", {
  # micSim() draws its seed from the caller's stream exactly as
  # simulate(seed = NULL) does.
  model <- rs_model()
  set.seed(77L)
  a <- as.data.frame(simulate(model, reestimate = FALSE))
  set.seed(77L)
  seed <- sample.int(.Machine$integer.max, 1L)
  expect_identical(a, pop_of(model, seed = seed))
})


# ---------------------------------------------------------------------
# 5. Cores
# ---------------------------------------------------------------------

skip_rng_parallel <- function() {
  skip_on_cran()
  skip_if_not(requireNamespace("MicSim", quietly = TRUE) &&
                utils::packageVersion("MicSim") >= "4.0.0",
              "installed MicSim >= 4.0.0 is required by the workers")
  skip_if(parallel::detectCores() < 2L, "needs two cores")
}


test_that("a run gives the same population on any number of cores", {
  skip_rng_parallel()
  model <- rs_fertility_model(rs_fertility_population(30L),
                              immigrants = rs_immigrants(6L, states = c("f/0", "m/0")))
  one <- simulate(model, seed = 3L, cores = 1L, reestimate = FALSE)
  two <- simulate(model, seed = 3L, cores = 2L, reestimate = FALSE)
  # R CMD check --as-cran allows two worker processes at most.
  limited <- nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_")) ||
    parallel::detectCores() < 3L
  three <- if (limited) two else
    simulate(model, seed = 3L, cores = 3L, reestimate = FALSE)

  expect_gt(one$extensions$fertility$n_births, 0L)
  expect_identical(as.data.frame(two), as.data.frame(one))
  expect_identical(as.data.frame(three), as.data.frame(one))

  births <- function(r) {
    b <- r$extensions$fertility$births
    b <- b[order(b$childID), , drop = FALSE]
    row.names(b) <- NULL
    b
  }
  expect_identical(births(three), births(one))
  expect_identical(three$diagnostics[c("n_events", "n_individuals",
                                       "n_newborns", "n_forced")],
                   one$diagnostics[c("n_events", "n_individuals",
                                     "n_newborns", "n_forced")])
})


test_that("immigrants are kept when there are more cores than starting individuals", {
  skip_rng_parallel()
  model <- rs_model(rs_population(2L), immigrants = rs_immigrants(5L))
  one <- pop_of(model, seed = 8L, cores = 1L)
  four <- pop_of(model, seed = 8L, cores = 4L)
  expect_setequal(unique(four$ID), c(1:2, 1001:1005))
  expect_identical(four, one)
})


test_that("a multi-core run without a seed draws one like a single-core run", {
  skip_rng_parallel()
  model <- rs_model()
  set.seed(21L); one <- simulate(model, cores = 1L, reestimate = FALSE)
  set.seed(21L); two <- simulate(model, cores = 2L, reestimate = FALSE)
  expect_identical(two$seed, one$seed)
  expect_true(two$meta$seed_drawn)
  expect_identical(as.data.frame(two), as.data.frame(one))
})
