# =====================================================================
# time.R  --  Time model: dates <-> day counts, and the year scale
# ---------------------------------------------------------------------
# WHAT THIS FILE DOES
#   Converts between calendar dates (yyyymmdd) and whole day counts, and
#   provides the fractional-year scale that rate functions are defined on.
#   Day counts are the engine's internal clock. They are exact integers that are
#   easy to compare and to do arithmetic with.
#
# PLACE IN THE PIPELINE
#   user dates (yyyymmdd) --> time.R --(day counts)------------> engine
#   day counts            --> time.R --(yyyymmdd, or years)----> engine
#
# CONCEPTS  (all defined in ?MicSim-concepts)
#   the three time scales (age, calendar time, duration); units in names
#
# HOW DATES ARE CONVERTED
#   Dates map to days with the civil-calendar algorithm of Howard Hinnant
#   (https://howardhinnant.github.io/date_algorithms.html), which counts
#   days from 1970-01-01 using integer arithmetic only. It is exact for
#   every proleptic-Gregorian date and needs no leap-year table.
#
# UNITS
#   Two time scales exist and must never be mixed:
#     * day counts -- exact integers, which are every date, age and duration the
#                     engine stores or reports.
#     * years      -- fractional, these are used ONLY for the age, calendar time
#                     and duration handed to a user's rate function.
#   The day-count-to-years conversion happens only in the rate-argument
#   helpers at the foot of this file.
# =====================================================================

#' Days since 1970-01-01 from a civil (proleptic Gregorian) date
#'
#' @param year,month,day numeric vectors (recycled); month in 1..12,
#'   day in 1..days_in_month. No validity checking is performed.
#' @return numeric vector of days since the epoch; negative before it.
#' @keywords internal
#' @noRd
days_from_civil <- function(year, month, day) {
  y   <- year - (month <= 2)          # shift the year to start in March
  era <- y %/% 400                    # floor division: correct for y < 0
  yoe <- y - era * 400                # year of era,  [0, 399]
  doy <- (153 * (month + 9 - 12 * (month > 2)) + 2) %/% 5 + day - 1  # [0, 365]
  doe <- yoe * 365 + yoe %/% 4 - yoe %/% 100 + doy                   # [0, 146096]
  era * 146097 + doe - 719468         # 719468 = days from 0000-03-01 to 1970-01-01
}

#' Civil (proleptic Gregorian) date from days since 1970-01-01
#'
#' @param days numeric vector of days since the epoch.
#' @return list(year, month, day), each a numeric vector.
#' @keywords internal
#' @noRd
civil_from_days <- function(days) {
  z   <- days + 719468
  era <- z %/% 146097                 # floor division: correct for z < 0
  doe <- z - era * 146097                                            # [0, 146096]
  yoe <- (doe - doe %/% 1460 + doe %/% 36524 - doe %/% 146096) %/% 365  # [0, 399]
  doy <- doe - (365 * yoe + yoe %/% 4 - yoe %/% 100)                 # [0, 365]
  mp  <- (5 * doy + 2) %/% 153                                       # [0, 11]
  d   <- doy - (153 * mp + 2) %/% 5 + 1                              # [1, 31]
  m   <- mp + 3 - 12 * (mp >= 10)                                    # [1, 12]
  list(year = yoe + era * 400 + (m <= 2), month = m, day = d)
}

#' Convert a yyyymmdd date to a day count
#'
#' @description
#' Turns calendar dates written as `yyyymmdd` into whole days since
#' 1970-01-01, which is the scale the engine works on. This is how you usually
#' bring dates into a MicSim model.
#'
#' @details
#' Accepts numeric, character, factor or `Date` input. An input that is
#' not a real calendar date (for example `20210229`, which is not a leap
#' year) becomes `NA`.
#'
#' Validity is checked by inversion. The date is turned into a day count
#' and back again. When that round trip returns the original date, 
#' the date is valid.
#'
#' @param yyyymmdd A vector of dates as `yyyymmdd` (numeric, character or
#'   factor). Alternatively, a `Date` vector, which is passed through.
#' @return Integer vector of days since 1970-01-01, negative before it, and
#'   `NA` for any input that is not a real calendar date.
#' @seealso [getInDateFormat()] for the inverse.
#' @examples
#' getInDays(20200826)
#' getInDays(20210229)   # not a real date: NA
#'
#' @export
getInDays <- function(yyyymmdd) {
  if (inherits(yyyymmdd, "Date")) return(as.integer(floor(unclass(yyyymmdd))))
  if (is.factor(yyyymmdd))        yyyymmdd <- as.character(yyyymmdd)

  z <- suppressWarnings(as.numeric(yyyymmdd))
  y <- z %/% 10000
  m <- (z %/% 100) %% 100
  d <- z %% 100

  days <- days_from_civil(y, m, d)

  # Validity by inversion (see Details), also rejects input that is not a
  # whole number.

  back <- civil_from_days(days)
  ok <- (z == trunc(z)) & back$year == y & back$month == m & back$day == d
  ok[is.na(ok)] <- FALSE
  days[!ok] <- NA_real_

  as.integer(days)
}

#' Age in days between two dates
#'
#' @description
#' The whole number of days between a birth date and a reference date,
#' both written as `yyyymmdd`.
#'
#' @details
#' The result is in **days**. Convert it to years (`/ 365.25`) only when
#' you need the age on the scale a rate function expects. The age reported
#' in the simulation output is also years on that scale (`days / 365.25`,
#' rounded to two decimals). Event *dates*, in contrast, are never put on
#' the year scale: they go through the exact calendar in [getInDateFormat()].
#'
#' @param currDate Reference date(s), `yyyymmdd`.
#' @param birthDate Birth date(s), `yyyymmdd`, recycled against `currDate`.
#' @return Integer vector of days. Negative when `currDate` is before
#'   `birthDate`; this is not treated as an error here.
#' @seealso [getInDays()]
#' @examples
#' getAgeInDays(20200826, 19800605)
#' @export
getAgeInDays <- function(currDate, birthDate) {
  getInDays(currDate) - getInDays(birthDate)
}

#' The whole day that contains a (possibly fractional) day count
#'
#' @description
#' Rounds a day count down to the whole day it falls in. The date decoders
#' need this because a calendar date names a whole day, not an instant
#' part-way through one.
#'
#' @details
#' Since the simulation clock is continuous, an event can happen at day
#' 20315.37, part-way through day 20315. A fractional day count must be
#' floored before decoding, because `civil_from_days()` assumes a whole
#' number. If it is given a fraction, its integer divisions can cross a month
#' boundary early and return a day-of-month of 0, which is not a real date.
#'
#' `floor()` is used instead of `trunc()`, as they differ for dates before 1970
#' (negative counts). An instant at day -0.5 falls inside day -1, which is
#' 1969-12-31, which `floor()` returns.
#'
#' @param days numeric vector of days since 1970-01-01, possibly
#'   fractional.
#' @return numeric vector of whole day counts.
#' @noRd
as_whole_days <- function(days) {
  floor(as.numeric(days))
}

#' Convert a day count back to a yyyymmdd date
#'
#' @description
#' The inverse of [getInDays()]: turns whole days since 1970-01-01 back
#' into a calendar date written as `yyyymmdd`.
#'
#' @details
#' A fractional day count is first floored to the whole day that contains
#' it (see `as_whole_days()`), so the result is always a real date.
#'
#' @param days numeric vector of days since 1970-01-01.
#' @return integer vector of yyyymmdd dates.
#' @seealso [getInDays()] for the inverse.
#' @examples
#' getInDateFormat(18500)
#' 18500 |> getInDateFormat()
#' @export
getInDateFormat <- function(days) {
  ymd <- civil_from_days(as_whole_days(days))
  as.integer(ymd$year * 10000 + ymd$month * 100 + ymd$day)
}

#' Calendar year, month, and day of a day count
#'
#' A fractional day count is floored to the day that contains it; see
#' `as_whole_days()`.
#'
#' @param days numeric vector of days since 1970-01-01.
#' @return integer vector.
#' @examples
#' getYear(18500); getMonth(18500); getDay(18500)
#' @export
#' @rdname civil-components
getYear  <- function(days) as.integer(civil_from_days(as_whole_days(days))$year)

#' @export
#' @rdname civil-components
getMonth <- function(days) as.integer(civil_from_days(as_whole_days(days))$month)

#' @export
#' @rdname civil-components
getDay   <- function(days) as.integer(civil_from_days(as_whole_days(days))$day)

#' All three civil components in one pass
#'
#' @description
#' Returns the calendar year, month and day of each day count in a single
#' decode. Prefer it to calling [getYear()], [getMonth()] and [getDay()]
#' separately, which decode three times.
#'
#' @param days numeric vector of days since 1970-01-01.
#' @return list of three integer vectors: `year`, `month`, `day`.
#' @keywords internal
#' @noRd
getCivil <- function(days) {
  ymd <- civil_from_days(as_whole_days(days))
  list(year  = as.integer(ymd$year),
       month = as.integer(ymd$month),
       day   = as.integer(ymd$day))
}

#' First day of a calendar month, as a day count
#'
#' @description
#' The day count of the first day of `month` in `year`.
#'
#' @details
#' Used to pin a forced school-enrolment event to a chosen month. The
#' engine works out when the event would otherwise fall, then moves it to
#' the first day of the enrolment month.
#'
#' @param year,month integer vectors, recycled against each other; `month`
#'   in 1..12.
#' @return integer vector of days since 1970-01-01 or `NA` where `year` or
#'   `month` is not a whole number, or `month` is not between 1 and 12.
#' @seealso [getInDays()], [getInDateFormat()]
#' @examples
#' getFirstDayOfMonth(2020, 9)
#' getInDateFormat(getFirstDayOfMonth(2020, 9))   # 20200901
#' @export
getFirstDayOfMonth <- function(year, month) {
  year  <- suppressWarnings(as.numeric(year))
  month <- suppressWarnings(as.numeric(month))
  days  <- days_from_civil(year, month, 1)
  # Like getInDays(), give NA for something that is not a date, instead of
  # rolling month 13 over into January of the next year.
  ok <- year == trunc(year) & month == trunc(month) & month >= 1 & month <= 12
  ok[is.na(ok)] <- FALSE
  days[!ok] <- NA_real_
  as.integer(days)
}

#' @rdname getFirstDayOfMonth
#' @details
#' `getInDays_my()` is an old alias of `getFirstDayOfMonth()`, should not be
#' used anymore.
#' @section Deprecated:
#' `getInDays_my()` is deprecated as of MicSim 4.0.0 and will be removed in
#' a future release. It still works, and shows a message once per R session
#' naming its replacement, `getFirstDayOfMonth()`, which takes the same
#' arguments. Set `options(MicSim.show_deprecation = FALSE)` to turn the
#' message off.
#' @export
getInDays_my <- function(year, month) {
  deprecate_message("getInDays_my")
  getFirstDayOfMonth(year, month)
}

# ---- The rate-argument convention ----------------------------------
# Rate functions are defined on fractional YEARS. These helpers do the
# day-count-to-years conversion for the age, calendar time and duration
# that a rate receives. Dates the engine records stay as day counts
# and are turned into yyyymmdd by getInDateFormat(); they are never put on
# the year scale.

#' @keywords internal
#' @noRd
caltime_in_years <- function(now_days)             1970 + now_days / DAYS_PER_YEAR
