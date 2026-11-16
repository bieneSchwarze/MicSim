# =====================================================================
# test-time.R
# ---------------------------------------------------------------------
# The time module is the foundation every other number in the package
# stands on: ages, event dates, durations and the arguments handed to
# every rate function are all derived from it. MicSim 3.0.0's version was
# wrong (defect D1) and its one test file asserted the wrong values, so
# there is no legacy test to inherit here.
#
# The tests below are mostly PROPERTY tests rather than table lookups:
# "this round-trips", "this is exact", "this agrees with base R". A table
# of hand-computed dates would only be as trustworthy as the hand that
# computed it.
#
# Tests are grouped by what they establish:
#   1. agreement with base R (the independent reference)
#   2. round-tripping
#   3. fractional day counts          <- regression, see below
#   4. invalid input
#   5. calendar edge cases
#   6. the legacy aliases
# =====================================================================


# ---------------------------------------------------------------------
# 1. Agreement with base R
# ---------------------------------------------------------------------

test_that("getInDays agrees with as.Date over a long span", {
  # base R's Date machinery is an independent implementation of the same
  # calendar, so agreeing with it over tens of thousands of days is strong
  # evidence that the Hinnant kernel is transcribed correctly.
  days <- seq(-40000L, 40000L, by = 37L)          # roughly 1860 to 2079
  dates <- as.Date(days, origin = "1970-01-01")
  ymd <- as.numeric(format(dates, "%Y%m%d"))

  expect_identical(getInDays(ymd), as.integer(days))
})


test_that("the civil decoders agree with as.Date", {
  days <- seq(-40000L, 40000L, by = 37L)
  dates <- as.Date(days, origin = "1970-01-01")

  expect_identical(getYear(days),  as.integer(format(dates, "%Y")))
  expect_identical(getMonth(days), as.integer(format(dates, "%m")))
  expect_identical(getDay(days),   as.integer(format(dates, "%d")))
})


test_that("getCivil matches the three separate decoders", {
  days <- c(-7305, 0, 1, 12345, 20000)
  civil <- getCivil(days)
  expect_identical(civil$year,  getYear(days))
  expect_identical(civil$month, getMonth(days))
  expect_identical(civil$day,   getDay(days))
})


# ---------------------------------------------------------------------
# 2. Round-tripping
# ---------------------------------------------------------------------

test_that("yyyymmdd -> days -> yyyymmdd is the identity", {
  days <- seq(-40000L, 40000L, by = 91L)
  ymd <- getInDateFormat(days)
  expect_identical(getInDays(ymd), as.integer(days))
  expect_identical(getInDateFormat(getInDays(ymd)), ymd)
})


test_that("the conversion is exact rather than approximate", {
  # This is the whole of defect D1. MicSim 3.0.0 returned 18500.44 here,
  # a non-integer, because it treated a year as 365.25 days and a month as
  # 30.42 days. Its unit tests asserted that wrong value.
  expect_identical(getInDays(20200826), as.integer(as.Date("2020-08-26")))
  expect_type(getInDays(20200826), "integer")
})


# ---------------------------------------------------------------------
# 3. Fractional day counts  (regression)
# ---------------------------------------------------------------------

test_that("a fractional day count reports the day it falls inside", {
  # The simulation clock is continuous: an event at day 20315.37 happened
  # part-way through day 20315 and belongs to that calendar date.
  expect_identical(getInDateFormat(20315.37), getInDateFormat(20315))
  expect_identical(getDay(20315.99),          getDay(20315))
})


test_that("fractional day counts never produce an invalid date", {
  # REGRESSION. Before the fix, civil_from_days() was handed fractional
  # input directly, and its integer divisions could tip over a month
  # boundary early and return a day-of-month of 0. getInDateFormat(8278.333)
  # returned 19920900 -- not a date -- and getInDays() then correctly
  # rejected it as NA, so a day count silently failed to round-trip.
  #
  # This is not a hypothetical: the legacy vignette builds its starting
  # population with runif() over a range of day counts and then calls
  # getInDateFormat() on the result, so every user following that idiom
  # feeds fractional day counts in.
  expect_identical(getInDateFormat(8278.333), 19920831L)
  expect_identical(getInDateFormat(15125.556), 20110531L)

  set.seed(4L)
  fractional <- runif(2000L, min = -40000, max = 40000)
  ymd <- getInDateFormat(fractional)
  expect_false(anyNA(getInDays(ymd)))
  expect_identical(getInDays(ymd), as.integer(floor(fractional)))
})


test_that("flooring, not truncation, is used before the epoch", {
  # floor() and trunc() differ for negative day counts. An instant at day
  # -0.5 falls inside day -1, which is 1969-12-31, not 1970-01-01.
  expect_identical(getInDateFormat(-0.5), 19691231L)
  expect_identical(getInDateFormat(0),    19700101L)
})


# ---------------------------------------------------------------------
# 4. Invalid input
# ---------------------------------------------------------------------

test_that("impossible dates are rejected as NA", {
  # Validity is established by inversion: the encode-decode round trip
  # reproduces the input if and only if the input was a real date. No
  # leap-year table is needed.
  expect_true(is.na(getInDays(20210229)))   # 2021 is not a leap year
  expect_true(is.na(getInDays(20211131)))   # November has 30 days
  expect_true(is.na(getInDays(20211301)))   # there is no month 13
  expect_true(is.na(getInDays(20210900)))   # there is no day 0
  expect_true(is.na(getInDays(20210832)))   # there is no 32nd
})


test_that("NA and non-numeric input propagate as NA", {
  expect_true(is.na(getInDays(NA)))
  expect_true(is.na(getInDays("not a date")))
  expect_identical(getInDays(c(20200826, NA)), c(18500L, NA_integer_))
})


test_that("a non-integral yyyymmdd is rejected", {
  # 20200826.5 is not a date, and silently rounding it would hide a
  # user's mistake.
  expect_true(is.na(getInDays(20200826.5)))
})


test_that("character, factor and Date inputs are all accepted", {
  expect_identical(getInDays("20200826"),         getInDays(20200826))
  expect_identical(getInDays(factor("20200826")), getInDays(20200826))
  expect_identical(getInDays(as.Date("2020-08-26")), getInDays(20200826))
})


test_that("empty input gives empty output", {
  expect_length(getInDays(numeric(0)), 0L)
  expect_length(getInDateFormat(numeric(0)), 0L)
})


# ---------------------------------------------------------------------
# 5. Calendar edge cases
# ---------------------------------------------------------------------

test_that("leap years are handled, including the century rules", {
  # 2000 is a leap year (divisible by 400), 1900 is not (divisible by 100
  # but not 400). Getting this wrong is the classic calendar bug.
  expect_identical(getInDays(20000229) - getInDays(20000228), 1L)
  expect_true(is.na(getInDays(19000229)))
  expect_identical(getInDays(20200229) - getInDays(20200228), 1L)
  expect_identical(getInDays(20210301) - getInDays(20210228), 1L)
})


test_that("dates before the epoch give negative day counts", {
  expect_lt(getInDays(19301231), 0L)
  expect_identical(getInDateFormat(getInDays(19301231)), 19301231L)
})


test_that("years before 1000 are zero-padded, not space-padded", {
  # MicSim 3.0.0 built its output with sprintf("%4d%02d%02d", ...), which
  # pads a three-digit year with a leading SPACE, producing " 9990101".
  # Integer arithmetic has no such failure mode.
  expect_identical(getInDateFormat(getInDays(9990101)), 9990101L)
})


# ---------------------------------------------------------------------
# 6. Derived helpers and legacy aliases
# ---------------------------------------------------------------------

test_that("getAgeInDays is the difference of two day counts", {
  expect_identical(getAgeInDays(20200826, 19800605),
                   getInDays(20200826) - getInDays(19800605))
  expect_gt(getAgeInDays(20200826, 19800605), 0L)
})


test_that("getFirstDayOfMonth lands on the first of the month", {
  expect_identical(getInDateFormat(getFirstDayOfMonth(2020, 9)), 20200901L)
  expect_identical(getInDateFormat(getFirstDayOfMonth(1969, 2)), 19690201L)
  expect_identical(getFirstDayOfMonth(2020, 9), getInDays(20200901))
})


test_that("getFirstDayOfMonth gives NA for a month that does not exist", {
  # Like getInDays(), and instead of rolling month 13 into January.
  expect_identical(getFirstDayOfMonth(2001, c(0, 13, 1.5, NA, 12)),
                   c(NA, NA, NA, NA, getInDays(20011201)))
  expect_identical(getFirstDayOfMonth(2001.5, 1), NA_integer_)
})


test_that("getInDays_my is an exact alias of getFirstDayOfMonth", {
  # Kept because MicSim 3.0.0 exported it. Its VALUES have changed: the
  # legacy body was (year - 1970) * 365.25 + month * 30.42 - 30.41, which
  # is neither an integer nor the first day of any month (defect D16).
  years  <- rep(1970:2030, each = 12L)
  months <- rep(1:12, times = 61L)
  expect_identical(getInDays_my(years, months),
                   getFirstDayOfMonth(years, months))
})


test_that("the time helpers are vectorised", {
  ymd <- c(19700101, 20000229, 20200826, 21001231)
  expect_length(getInDays(ymd), 4L)
  expect_length(getInDateFormat(getInDays(ymd)), 4L)
  expect_length(getYear(getInDays(ymd)), 4L)
  expect_identical(getInDateFormat(getInDays(ymd)), as.integer(ymd))
})
