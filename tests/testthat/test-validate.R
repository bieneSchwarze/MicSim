# =====================================================================
# test-validate.R
# ---------------------------------------------------------------------
# validate() runs every check and reports every finding at once;
# construction collects all faults of a population table.
# =====================================================================

v_pop <- data.frame(ID = 1:30, birthDate = "19700101", initState = "f/0")
v_ss  <- state_space(sex = c("m", "f"), par = c("0", "1"))
v_ok  <- function(age) rep(0.01, length(age))

v_model <- function(...) {
  micsim_model(v_ss, transitions(...), v_pop, c(20000101, 20101231))
}

test_that("a clean model passes silently and is returned for the pipe", {
  m <- v_model(transition(move("par", "0", "1", where = sex == "f"), v_ok),
               transition(NA, "dead", v_ok))
  expect_silent(out <- validate(m))
  expect_identical(out, m)
})

test_that("all faulty rates are reported in one error", {
  na_rate  <- function(age) ifelse(age > 60, NA, 0.01)
  neg_rate <- function(age) ifelse(age > 70, -1, 0.02)
  m <- v_model(transition(move("par", "0", "1", where = sex == "f"), na_rate),
               transition(NA, "dead", neg_rate))
  err <- tryCatch(validate(m), error = function(e) e)
  expect_s3_class(err, "micsim_invalid_model")
  expect_s3_class(err, "micsim_rate_error")      # still caught as before
  expect_identical(sort(err$problems$check), c("NA", "negative"))
  expect_match(conditionMessage(err), "2 problem")
  expect_match(conditionMessage(err), "age 61", fixed = TRUE)
  expect_match(conditionMessage(err), "age 71", fixed = TRUE)
})

test_that("a rate returning the wrong length or type on the grid is an error", {
  wrong_length <- function(age) {
    if (length(age) > 10) rep(0.01, 5) else rep(0.01, length(age))
  }
  m <- v_model(transition(NA, "dead", wrong_length))
  err <- tryCatch(validate(m), error = function(e) e)
  expect_identical(err$problems$check, "length")

  as_text <- function(age) if (length(age) > 10) rep("a", length(age)) else rep(0.01, length(age))
  m <- v_model(transition(NA, "dead", as_text))
  expect_identical(tryCatch(validate(m), error = function(e) e)$problems$check,
                   "type")
})

test_that("a rate that errors on the grid does not stop the other checks", {
  boom <- function(age) if (any(age > 80)) stop("no data above 80") else rep(0.01, length(age))
  neg  <- function(age) rep(-0.1, length(age))
  m <- v_model(transition(move("par", "0", "1", where = sex == "f"), boom),
               transition(NA, "dead", neg))
  err <- tryCatch(validate(m), error = function(e) e)
  expect_setequal(err$problems$check, c("evaluation", "negative"))
  expect_match(conditionMessage(err), "no data above 80", fixed = TRUE)
})

test_that("suspicious models give one warning listing every finding", {
  zero <- function(age) rep(0, length(age))
  big  <- function(age) rep(25, length(age))
  dur  <- function(age, duration) rep(0.01, length(age))
  m <- v_model(transition(move("par", "0", "1", where = sex == "m"), zero),
               transition(move("par", "0", "1", where = sex == "f"), dur),
               transition(NA, "dead", big))
  w <- tryCatch(validate(m), warning = function(w) w)
  expect_s3_class(w, "micsim_model_warning")
  # A rate that depends on duration is not, by itself, suspicious.
  expect_setequal(w$problems$check, c("zero", "magnitude", "unreachable"))
  # The model is still returned
  expect_identical(suppressWarnings(validate(m)), m)
})

test_that("strict = TRUE turns warnings into errors", {
  zero <- function(age) rep(0, length(age))
  m <- v_model(transition(move("par", "0", "1", where = sex == "f"), zero),
               transition(NA, "dead", v_ok))
  expect_error(validate(m, strict = TRUE), class = "micsim_invalid_model")
  expect_error(validate(m, strict = "yes"), class = "micsim_model_error")
})

test_that("states that can be reached but not left are reported", {
  m <- v_model(transition(move("par", "0", "1", where = sex == "f"), v_ok))
  w <- tryCatch(validate(m), warning = function(w) w)
  expect_identical(w$problems$check, "stuck")
  expect_match(conditionMessage(w), "f/1", fixed = TRUE)
})

test_that("newborn states count as entry states for reachability", {
  m <- v_model(transition(move("par", "0", "1", where = sex == "f"), v_ok),
               transition(NA, "dead", v_ok)) |>
    add_subsystem(fertility(moves_of("par"),
                            newborn_rule(sex = c(m = 0.5, f = 0.5), par = "0")))
  # m/0 is reached only through births; mortality from it is reachable.
  expect_silent(validate(m))
})

test_that("a newborn table without a candidate for a mother's state is an error", {
  rule <- new_newborn_rule(states = matrix(c("m", "0", "f", "0"), ncol = 2,
                                           byrow = TRUE),
                           prob = c(1, 1), inherit = 1L)
  # Mothers are female after the birth, but the table only inherits sex,
  # so a child of an "f" mother has the single candidate f/0: covered.
  m <- v_model(transition(move("par", "0", "1", where = sex == "f"), v_ok),
               transition(NA, "dead", v_ok)) |>
    add_subsystem(fertility(moves_of("par"), rule))
  expect_silent(validate(m))

  # Inheriting parity: a mother has parity 1 after the birth, and the
  # table has no candidate with parity 1.
  rule2 <- new_newborn_rule(matrix(c("m", "0", "f", "0"), ncol = 2, byrow = TRUE),
                            prob = c(0.5, 0.5), inherit = 2L)
  m2 <- v_model(transition(move("par", "0", "1", where = sex == "f"), v_ok),
                transition(NA, "dead", v_ok)) |>
    add_subsystem(fertility(moves_of("par"), rule2))
  err <- tryCatch(validate(m2), error = function(e) e)
  expect_s3_class(err, "micsim_newborn_rule_error")
  expect_identical(err$problems$check, "newborn")
  expect_match(conditionMessage(err), "f/1", fixed = TRUE)
})

test_that("abort_invalid_model() still accepts a character vector", {
  err <- tryCatch(abort_invalid_model(c("first", "second")),
                  error = function(e) e)
  expect_s3_class(err, "micsim_invalid_model")
  expect_identical(err$problems$detail, c("first", "second"))
})

test_that("long reports are capped, the rest counted", {
  rates <- lapply(1:25, function(k) function(age) rep(-k, length(age)))
  pop <- data.frame(ID = 1:5, birthDate = "19700101", initState = "a1")
  vals <- paste0("a", 1:26)
  trs <- lapply(1:25, function(k) transition(vals[k], vals[k + 1], rates[[k]]))
  m <- micsim_model(state_space(x = vals), do.call(transitions, trs), pop,
                    c(20000101, 20101231))
  err <- tryCatch(validate(m), error = function(e) e)
  expect_identical(nrow(err$problems), 25L)
  expect_match(conditionMessage(err), "and 5 more", fixed = TRUE)
})

test_that("all faults of a starting population are reported together", {
  bad <- data.frame(ID = c(1, 1, 3),
                    birthDate = c("19700101", "20300101", "19700230"),
                    initState = c("f/0", "f/0", "x/y"))
  err <- tryCatch(
    micsim_model(v_ss, transitions(transition(NA, "dead", v_ok)), bad,
                 c(20000101, 20101231)),
    error = function(e) e)
  expect_s3_class(err, "micsim_engine_error")
  expect_s3_class(err, "micsim_codec_error")
  expect_length(err$problems, 4L)
  expect_match(conditionMessage(err), "4 problems")
})

test_that("a single population fault keeps its own class and message", {
  bad <- data.frame(ID = 1, birthDate = "19700101", initState = "x/y")
  expect_error(
    micsim_model(v_ss, transitions(transition(NA, "dead", v_ok)), bad,
                 c(20000101, 20101231)),
    "Unknown state label", class = "micsim_codec_error")
})

test_that("all faults of an immigrant table are reported together", {
  imm <- data.frame(ID = c(5, 40), immigrDate = c("19990101", "20050101"),
                    birthDate = c("19800101", "20060101"),
                    immigrInitState = "f/0")
  err <- tryCatch(
    micsim_model(v_ss, transitions(transition(NA, "dead", v_ok)), v_pop,
                 c(20000101, 20101231), immigrants = imm),
    error = function(e) e)
  expect_s3_class(err, "micsim_engine_error")
  expect_match(conditionMessage(err), "before they are born")
  expect_match(conditionMessage(err), "before the simulation starts")
  expect_match(conditionMessage(err), "both the starting and the immigrant")
})
