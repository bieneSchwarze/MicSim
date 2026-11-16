# =====================================================================
# test-converters.R
# ---------------------------------------------------------------------
# The two converters are the legacy OUTPUT surface: what every analysis
# script written against MicSim reads. Long format is by far the more
# intricate of the two, and most of this file is about its three hard
# parts -- left truncation, right censoring, and the episodes that have to
# be manufactured because nothing in `pop` records them.
#
# Tests are grouped by what they establish:
#   1. long format: shape
#   2. long format: episodes are contiguous and complete
#   3. long format: left truncation
#   4. long format: right censoring
#   5. long format: the OD column
#   6. long format: immigrants
#   7. long format: context is never taken from the global environment
#   8. wide format
# =====================================================================


conv_state_space <- function() {
  expand.grid(sex = c("m", "f"), fert = c("0", "1"),
              stringsAsFactors = FALSE)
}

fertRates <- function(age, calTime) ifelse(age > 15 & age < 45, 0.2, 0)
mortRates <- function(age, calTime) 0.00003 * exp(0.1 * age)
neverRates <- function(age, calTime) rep(0, length(age))

conv_matrix <- function() {
  buildTransitionMatrix(
    allTransitions = cbind("f/0->f/1", "fertRates"),
    absTransitions = cbind("dead", "mortRates"),
    stateSpace     = conv_state_space()
  )
}

#' A run with births, deaths and survivors, so that every episode kind is
#' present in one data frame.
conv_pop <- function(n = 40L, seed = 1L, horizon = c(20000101, 20201231),
                     birth = "19800101") {
  set.seed(seed)
  micSim(
    initPop = data.frame(
      ID = seq_len(n), birthDate = rep(birth, n),
      initState = rep("f/0", n), stringsAsFactors = FALSE
    ),
    transitionMatrix = conv_matrix(),
    varInitStates    = rbind(c("m", "0"), c("f", "0")),
    initStatesProb   = c(0.515, 0.485),
    maxAge           = 100,
    simHorizon       = horizon,
    fertTr           = "f/0->f/1",
    verbose          = FALSE
  )
}


# ---------------------------------------------------------------------
# 1. Long format: shape
# ---------------------------------------------------------------------

test_that("long format has the contract columns in the contract order", {
  pop <- conv_pop()
  long <- convertToLongFormat(pop, stateSpace = conv_state_space())

  expect_identical(
    names(long),
    c("ID", "birthDate", "Tstart", "Tstop", "statusEntry", "statusExit",
      "OD", "ns", "Episode", "sex", "fert")
  )
})


test_that("state variables are split into one column each", {
  # `pop` carries "f/0" as one string; an analysis wants two variables.
  pop <- conv_pop()
  long <- convertToLongFormat(pop, stateSpace = conv_state_space())

  expect_true(all(long$sex %in% c("m", "f")))
  expect_true(all(long$fert %in% c("0", "1")))
})


test_that("the state-variable columns take their names from stateSpace", {
  # A transition matrix carries the state space but not the NAMES of the
  # state variables, so without stateSpace they are numbered.
  pop <- conv_pop(n = 10L)
  named <- convertToLongFormat(pop, stateSpace = conv_state_space())
  unnamed <- convertToLongFormat(pop)

  expect_true(all(c("sex", "fert") %in% names(named)))
  expect_true(all(c("var1", "var2") %in% names(unnamed)))
})


test_that("motherID is dropped", {
  # It describes an individual, not an episode, and would be repeated
  # across every one of their rows.
  pop <- conv_pop()
  expect_true("motherID" %in% names(pop))
  long <- convertToLongFormat(pop, stateSpace = conv_state_space())
  expect_false("motherID" %in% names(long))
})


test_that("every individual appears at least once", {
  pop <- conv_pop()
  long <- convertToLongFormat(pop, stateSpace = conv_state_space())
  expect_setequal(unique(long$ID), unique(pop$ID))
})


# ---------------------------------------------------------------------
# 2. Episodes are contiguous and complete
# ---------------------------------------------------------------------

test_that("each individual's episodes join end to end", {
  # The defining property of an episode table: no gaps and no overlaps.
  # One spell ends exactly where the next begins.
  pop <- conv_pop()
  long <- convertToLongFormat(pop, stateSpace = conv_state_space())

  for (id in unique(long$ID)) {
    rows <- long[long$ID == id, , drop = FALSE]
    rows <- rows[order(rows$Episode), , drop = FALSE]
    if (nrow(rows) > 1L) {
      expect_identical(rows$Tstart[-1L], rows$Tstop[-nrow(rows)])
    }
  }
})


test_that("no episode ends before it begins", {
  pop <- conv_pop()
  long <- convertToLongFormat(pop, stateSpace = conv_state_space())
  expect_true(all(getInDays(long$Tstop) >= getInDays(long$Tstart)))
})


test_that("episodes are numbered from one and ns counts them", {
  pop <- conv_pop()
  long <- convertToLongFormat(pop, stateSpace = conv_state_space())

  for (id in unique(long$ID)) {
    rows <- long[long$ID == id, , drop = FALSE]
    expect_identical(sort(rows$Episode), seq_len(nrow(rows)))
    expect_true(all(rows$ns == nrow(rows)))
  }
})


test_that("no episode reaches outside the simulation horizon", {
  horizon <- c(20000101, 20201231)
  pop <- conv_pop(horizon = horizon)
  long <- convertToLongFormat(pop, stateSpace = conv_state_space())

  expect_true(all(getInDays(long$Tstart) >= getInDays(horizon[1L])))
  expect_true(all(getInDays(long$Tstop)  <= getInDays(horizon[2L])))
})


# ---------------------------------------------------------------------
# 3. Left truncation
# ---------------------------------------------------------------------

test_that("someone alive before the simulation began is left truncated", {
  # We know she was childless when observation started, but not for how
  # long, so her first spell's beginning was not observed.
  pop <- conv_pop(birth = "19800101", horizon = c(20000101, 20201231))
  long <- convertToLongFormat(pop, stateSpace = conv_state_space())

  first <- long[long$Episode == 1L & long$ID <= 40L, , drop = FALSE]
  expect_gt(nrow(first), 0L)
  expect_true(all(first$statusEntry == 0L))
  expect_true(all(first$Tstart == "20000101"))
})


test_that("someone born during the simulation is not left truncated", {
  # A newborn's first spell begins at birth, and we saw it begin.
  pop <- conv_pop()
  long <- convertToLongFormat(pop, stateSpace = conv_state_space())

  newborns <- setdiff(unique(long$ID), 1:40)
  expect_gt(length(newborns), 0L)

  first <- long[long$ID %in% newborns & long$Episode == 1L, , drop = FALSE]
  expect_true(all(first$statusEntry == 1L))
  expect_identical(first$Tstart, first$birthDate)
})


test_that("later episodes are always observed at entry", {
  # Only a FIRST episode can be left truncated: every later one begins
  # with a transition we watched happen.
  pop <- conv_pop()
  long <- convertToLongFormat(pop, stateSpace = conv_state_space())
  later <- long[long$Episode > 1L, , drop = FALSE]
  expect_gt(nrow(later), 0L)
  expect_true(all(later$statusEntry == 1L))
})


# ---------------------------------------------------------------------
# 4. Right censoring
# ---------------------------------------------------------------------

test_that("a spell still running at the end is censored there", {
  horizon <- c(20000101, 20201231)
  pop <- conv_pop(horizon = horizon)
  long <- convertToLongFormat(pop, stateSpace = conv_state_space())

  censored <- long[long$statusExit == 0L, , drop = FALSE]
  expect_gt(nrow(censored), 0L)
  expect_true(all(censored$Tstop == as.character(horizon[2L])))
  expect_true(all(censored$OD == "cens"))
})


test_that("an individual has at most one censored episode, and it is last", {
  pop <- conv_pop()
  long <- convertToLongFormat(pop, stateSpace = conv_state_space())

  for (id in unique(long$ID)) {
    rows <- long[long$ID == id, , drop = FALSE]
    rows <- rows[order(rows$Episode), , drop = FALSE]
    n_censored <- sum(rows$statusExit == 0L)
    expect_lte(n_censored, 1L)
    if (n_censored == 1L) {
      expect_identical(rows$statusExit[nrow(rows)], 0L)
    }
  }
})


test_that("a manufactured censoring episode records the state entered", {
  # The row it is copied from describes the spell the individual LEFT, so
  # its state has to be advanced to the one they entered. A woman whose
  # last recorded event was a first birth spends her censored spell at
  # parity one, not parity zero.
  pop <- conv_pop()
  long <- convertToLongFormat(pop, stateSpace = conv_state_space())

  gave_birth <- long$ID[long$OD == "f/0->f/1" | long$OD == "0->1"]
  survivors <- intersect(gave_birth, long$ID[long$statusExit == 0L])
  expect_gt(length(survivors), 0L)

  for (id in survivors) {
    rows <- long[long$ID == id, , drop = FALSE]
    last <- rows[which.max(rows$Episode), ]
    expect_identical(last$fert, "1")
  }
})


test_that("someone who died is not censored", {
  # Their last spell ended, and we saw it end.
  pop <- conv_pop(n = 40L, horizon = c(20000101, 20801231))
  long <- convertToLongFormat(pop, stateSpace = conv_state_space())

  died <- unique(long$ID[long$OD == "dead"])
  expect_gt(length(died), 0L)
  for (id in died) {
    rows <- long[long$ID == id, , drop = FALSE]
    last <- rows[which.max(rows$Episode), ]
    expect_identical(last$statusExit, 1L)
    expect_identical(last$OD, "dead")
  }
})


test_that("someone who never transitioned gets exactly one censored episode", {
  # The simplest case, and the one where `pop` has a single NA row.
  set.seed(5L)
  pop <- micSim(
    initPop = data.frame(ID = 1:5, birthDate = rep("19800101", 5L),
                         initState = rep("f/0", 5L),
                         stringsAsFactors = FALSE),
    transitionMatrix = buildTransitionMatrix(
      NULL, cbind("dead", "neverRates"), conv_state_space()
    ),
    maxAge = 100, simHorizon = c(20000101, 20101231), verbose = FALSE
  )
  long <- convertToLongFormat(pop, stateSpace = conv_state_space())

  expect_identical(nrow(long), 5L)
  expect_true(all(long$Episode == 1L))
  expect_true(all(long$ns == 1L))
  expect_true(all(long$statusEntry == 0L))
  expect_true(all(long$statusExit == 0L))
  expect_true(all(long$OD == "cens"))
})


# ---------------------------------------------------------------------
# 5. The OD column
# ---------------------------------------------------------------------

test_that("an absorbing destination is reported by name", {
  pop <- conv_pop(horizon = c(20000101, 20801231))
  long <- convertToLongFormat(pop, stateSpace = conv_state_space())
  expect_true("dead" %in% long$OD)
})


test_that("a single-variable pattern is reported under the user's name", {
  # The match is made component by component: "f/DE/0" -> "f/DE/1" gives
  # the candidates "f->f", "DE->DE", "0->1". Only a pattern naming ONE
  # state variable can equal one of those, so only such a pattern can ever
  # be reported. MicSim 3.0.0 behaves identically.
  pop <- conv_pop()
  long <- convertToLongFormat(
    pop, stateSpace = conv_state_space(),
    allTransitions = cbind("0->1", "fertRates")
  )
  expect_true("0->1" %in% long$OD)
})


test_that("a multi-variable pattern falls back to the changed variable", {
  # Declaring the transition as "f/0->f/1" does not make OD say so: the
  # full label is not among the candidates, so the changed-variable form
  # is used. This is legacy behaviour and every analysis script written
  # against MicSim expects it.
  pop <- conv_pop()
  long <- convertToLongFormat(
    pop, stateSpace = conv_state_space(),
    allTransitions = cbind("f/0->f/1", "fertRates")
  )
  expect_true("0->1" %in% long$OD)
  expect_false("f/0->f/1" %in% long$OD)
})


test_that("without allTransitions the changed variable is named", {
  pop <- conv_pop()
  long <- convertToLongFormat(pop, stateSpace = conv_state_space())
  expect_true("0->1" %in% long$OD)
  expect_false("f/0->f/1" %in% long$OD)
})


# ---------------------------------------------------------------------
# 6. Immigrants
# ---------------------------------------------------------------------

test_that("an immigrant's first episode begins at their arrival", {
  # Their earlier life happened elsewhere and is outside the study, so the
  # episode that belongs in the data starts at the border -- and its
  # beginning counts as observed.
  immigrants <- data.frame(
    ID              = 101:103,
    immigrDate      = c(20050601L, 20080601L, 20120601L),
    birthDate       = c(19750601L, 19780601L, 19820601L),
    immigrInitState = rep("f/0", 3L),
    stringsAsFactors = FALSE
  )

  set.seed(7L)
  pop <- micSim(
    initPop = data.frame(ID = 1:20, birthDate = rep("19800101", 20L),
                         initState = rep("f/0", 20L),
                         stringsAsFactors = FALSE),
    immigrPop        = immigrants,
    transitionMatrix = conv_matrix(),
    varInitStates    = rbind(c("m", "0"), c("f", "0")),
    initStatesProb   = c(0.515, 0.485),
    maxAge           = 100,
    simHorizon       = c(20000101, 20201231),
    fertTr           = "f/0->f/1",
    verbose          = FALSE
  )

  long <- convertToLongFormat(pop, migr = TRUE,
                              stateSpace = conv_state_space())
  first <- long[long$ID %in% immigrants$ID & long$Episode == 1L, ,
                drop = FALSE]

  expect_identical(nrow(first), 3L)
  expect_identical(first$Tstart,
                   as.character(immigrants$immigrDate[order(immigrants$ID)]))
  expect_true(all(first$statusEntry == 1L))
})


test_that("migr = FALSE leaves immigrants left truncated at the start", {
  # Which is what MicSim 3.0.0 does too: the adjustment is opt-in.
  immigrants <- data.frame(
    ID = 101L, immigrDate = 20050601L, birthDate = 19750601L,
    immigrInitState = "f/0", stringsAsFactors = FALSE
  )
  set.seed(8L)
  pop <- micSim(
    initPop = data.frame(ID = 1:10, birthDate = rep("19800101", 10L),
                         initState = rep("f/0", 10L),
                         stringsAsFactors = FALSE),
    immigrPop = immigrants, transitionMatrix = conv_matrix(),
    varInitStates = rbind(c("m", "0"), c("f", "0")),
    initStatesProb = c(0.515, 0.485), maxAge = 100,
    simHorizon = c(20000101, 20201231), fertTr = "f/0->f/1",
    verbose = FALSE
  )

  long <- convertToLongFormat(pop, migr = FALSE,
                              stateSpace = conv_state_space())
  first <- long[long$ID == 101L & long$Episode == 1L, , drop = FALSE]
  expect_identical(first$Tstart, "20000101")
  expect_identical(first$statusEntry, 0L)
})


# ---------------------------------------------------------------------
# 7. Context is never taken from the global environment
# ---------------------------------------------------------------------

test_that("the converter works inside another function", {
  # MicSim 3.0.0's version read five objects from globalenv() by name, so
  # it could not be called from inside a function and could not be tested
  # (defect D3). This test is the defect, inverted.
  wrapped <- function() {
    pop <- conv_pop(n = 10L)
    convertToLongFormat(pop, stateSpace = conv_state_space())
  }
  long <- wrapped()
  expect_gt(nrow(long), 0L)
})


test_that("missing context is an error, not a silent wrong answer", {
  pop <- conv_pop(n = 5L)
  attr(pop, "micsim_spec") <- NULL
  expect_error(convertToLongFormat(pop, stateSpace = conv_state_space()),
               class = "micsim_convert_error")
})


test_that("explicit arguments override what micSim recorded", {
  pop <- conv_pop(n = 10L, horizon = c(20000101, 20201231))
  long <- convertToLongFormat(pop, stateSpace = conv_state_space(),
                              simHorizon = c(20000101, 20151231))
  censored <- long[long$statusExit == 0L, , drop = FALSE]
  expect_true(all(censored$Tstop == "20151231"))
})


# ---------------------------------------------------------------------
# 8. Wide format
# ---------------------------------------------------------------------

test_that("wide format gives exactly one row per individual", {
  pop <- conv_pop()
  wide <- convertToWideFormat(pop)

  expect_identical(nrow(wide), length(unique(pop$ID)))
  expect_false(any(duplicated(wide$ID)))
})


test_that("transitions are spread across numbered columns", {
  pop <- conv_pop()
  wide <- convertToWideFormat(pop)

  expect_true(all(c("ID", "birthDate", "initState", "ns") %in% names(wide)))
  expect_true(all(c("From.1", "To.1", "transitionTime.1",
                    "transitionAge.1") %in% names(wide)))
})


test_that("ns counts transitions, and is zero for someone with none", {
  # Someone who never transitioned still occupies one row of `pop`, with
  # NA in all four transition columns. Their life course contains no
  # transitions, so ns is 0 and not 1.
  set.seed(9L)
  pop <- micSim(
    initPop = data.frame(ID = 1:5, birthDate = rep("19800101", 5L),
                         initState = rep("f/0", 5L),
                         stringsAsFactors = FALSE),
    transitionMatrix = buildTransitionMatrix(
      NULL, cbind("dead", "neverRates"), conv_state_space()
    ),
    maxAge = 100, simHorizon = c(20000101, 20101231), verbose = FALSE
  )
  wide <- convertToWideFormat(pop)

  expect_identical(nrow(wide), 5L)
  expect_true(all(wide$ns == 0L))
  expect_true(all(is.na(wide$transitionTime.1)))
})


test_that("motherID is dropped from wide format too", {
  pop <- conv_pop()
  expect_false("motherID" %in% names(convertToWideFormat(pop)))
})


test_that("an individual's transitions stay in chronological order", {
  pop <- conv_pop(horizon = c(20000101, 20801231))
  wide <- convertToWideFormat(pop)
  two_or_more <- wide[wide$ns >= 2L, , drop = FALSE]
  expect_gt(nrow(two_or_more), 0L)

  expect_true(all(getInDays(two_or_more$transitionTime.2) >=
                    getInDays(two_or_more$transitionTime.1)))
  expect_identical(two_or_more$To.1, two_or_more$From.2)
})
