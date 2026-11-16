# =====================================================================
# test-deprecation.R
# The classic (3.0.0) entry points are exported, still work, and show a
# message once per session naming their replacement.
# =====================================================================

with_messages <- function(code) {
  old <- options(MicSim.show_deprecation = TRUE)
  on.exit(options(old), add = TRUE)
  reset_deprecation_messages()
  force(code)
}

dep_fixture <- function() {
  mort <- function(age, calTime) 0.001 * exp(0.05 * age)
  assign("dep_mort_rates", mort, envir = globalenv())
  ss <- expand.grid(sex = c("m", "f"), stringsAsFactors = FALSE)
  tm <- buildTransitionMatrix(allTransitions = NULL,
                              absTransitions = c("dead", "dep_mort_rates"),
                              stateSpace = ss)
  pop <- data.frame(ID = 1:10, birthDate = "19700101",
                    initState = rep(c("m", "f"), 5))
  list(tm = tm, pop = pop, ss = ss)
}

test_that("every deprecated function is exported", {
  exports <- getNamespaceExports("MicSim")
  expect_true(all(names(DEPRECATED_REPLACEMENTS) %in% exports))
})

test_that("each deprecated function shows its message once per session", {
  with_messages({
    fx <- suppressMessages(dep_fixture())  # buildTransitionMatrix() shows its message here
    expect_message(
      buildTransitionMatrix(NULL, c("dead", "dep_mort_rates"), fx$ss),
      NA
    )
    expect_message(
      pop <- micSim(fx$pop, transitionMatrix = fx$tm,
                    simHorizon = c(20000101, 20101231), verbose = FALSE),
      class = "micsim_deprecation_message"
    )
    expect_message(
      micSim(fx$pop, transitionMatrix = fx$tm,
             simHorizon = c(20000101, 20101231), verbose = FALSE),
      NA
    )
    expect_message(convertToLongFormat(pop, stateSpace = fx$ss),
                   "as_long", class = "micsim_deprecation_message")
    expect_message(convertToWideFormat(pop), "as_wide",
                   class = "micsim_deprecation_message")
    expect_message(getInDays_my(2020, 9), "getFirstDayOfMonth",
                   class = "micsim_deprecation_message")
  })
  rm("dep_mort_rates", envir = globalenv())
})

test_that("buildTransitionMatrix() names its replacement", {
  with_messages(
    expect_message(dep_fixture(), "transition()",
                   class = "micsim_deprecation_message", fixed = TRUE)
  )
  rm("dep_mort_rates", envir = globalenv())
})

test_that("the option turns the messages off", {
  reset_deprecation_messages()
  old <- options(MicSim.show_deprecation = FALSE)
  on.exit(options(old), add = TRUE)
  expect_message(getInDays_my(2020, 9), NA)
})

test_that("deprecated functions behave exactly as before", {
  expect_identical(getInDays_my(c(2020, 1969), c(9, 2)),
                   getFirstDayOfMonth(c(2020, 1969), c(9, 2)))
})

test_that("the new API never shows a deprecation message", {
  with_messages({
    mort <- function(age) 0.001 * exp(0.05 * age)
    model <- micsim_model(
      state_space(sex = c("m", "f")),
      transitions(transition(NA, "dead", mort)),
      data.frame(ID = 1:10, birthDate = "19700101", initState = "f"),
      c(20000101, 20101231)
    )
    expect_message({
      r <- simulate(model, seed = 1, reestimate = FALSE)
      as_long(r)
      as_wide(r)
    }, NA)
  })
})

test_that("micSimParallel() shows one message, not one per worker", {
  skip_on_cran()
  with_messages({
    fx <- suppressMessages(dep_fixture())
    msgs <- character(0)
    withCallingHandlers(
      micSimParallel(initPop = fx$pop, transitionMatrix = fx$tm,
                     simHorizon = c(20000101, 20101231), cores = 2L,
                     verbose = FALSE),
      micsim_deprecation_message = function(m) {
        msgs <<- c(msgs, m$old)
        invokeRestart("muffleMessage")
      }
    )
    expect_identical(msgs, "micSimParallel")
  })
  rm("dep_mort_rates", envir = globalenv())
})
