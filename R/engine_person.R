# =====================================================================
# engine_person.R
# ---------------------------------------------------------------------
# WHAT THIS FILE DOES
#   Holds everything the simulation knows about each individual.
#   Provides the information the engine asks for at every event of an individual:
#   how old is this person, what state are they in, and how long has it been
#   since the state variable this transition would change last changed.
#
# PLACE IN THE PIPELINE
#   engine.R --(add / move_to)--> engine_person.R --(age, duration)--> registry.R
#                                        |
#                                        +--(id, birth date, states)--> engine_store.R
#
#   The event queue holds SLOTS, not IDs. A slot is simply the row
#   number of an individual in this store. Everything else about that
#   individual is looked up here in constant time.
#
# WHAT YOU NEED TO KNOW  (concepts in ?MicSim-concepts)
#   the three time scales, especially duration; composite state, state
#   variable, absorbing state.
#
#
# GUARANTEES
#   * Slots are handed out consecutively from 1 and never reused, so a slot
#     identifies an individual for the whole run.
#   * Every individual has one clock per state variable. A clock either holds the
#     the date the variable was last seen to change, or, if there hasn't been a
#     change yet, the individual's date of birth. It is never NA.
#   * Every lookup is O(1), the whole population isn't scanned.
#
# UNITS
#   All times are in days on the getInDays() scale: whole days
#   since 1970-01-01 for birth dates, possibly fractional for event dates.
#   The engine converts the fractional years that rate functions expect.
#   This file does not handle conversion.
# =====================================================================


#' Create an empty person store
#'
#' @description
#' Builds the table of individuals in a simulation: who they
#' are, when they were born, what state they are in, and when each of
#' their state variables last changed.
#'
#' @details
#' **Columns, not rows.** Each attribute is its own long vector, and an
#' individual is a position shared across all of them. Adding a person
#' writes one element to each vector. When a vector fills up it doubles in size.
#' Adding N people therefore costs O(N) in total, compared to O(N^2) of
#' growing a data frame one `rbind()` at a time.
#'
#' **IDs and slots.** The user-facing `ID` is whatever the model
#' supplies. The engine works with the *slot*, the individual's position in
#' this store, because that is an index and doesn't need a lookup. Going from
#' an ID to a slot uses a hashed environment, because R
#' hashes environments but not the names of an atomic vector. This means that
#' `env[["4711"]]` is a constant-time lookup, while `vec[["4711"]]` is a
#' scan of every name. The three operations that need the reverse
#' direction are: seeding the starting population, injecting immigrants and
#' recording a newborn's mother.
#'
#' **The duration clocks.** This is the main reason for this file. The store
#' records one date per state variable for each individual, which records when
#' that variable last changed. Asking "how long has this person been a former
#' smoking" is then `now - clock[smoking_status]`. This means that there
#' only needs to be a subtraction not a full search or kept history.
#'
#' The clocks start at the individual's date of birth and are reset to the
#' current date whenever a transition changes the corresponding variable.
#' That single rule handles all three kinds of individual the same way:
#'
#' \describe{
#'   \item{born during the simulation}{duration is measured from the birth
#'     date, which is where the clock starts.}
#'   \item{member of the starting population}{their entry into the state
#'     they start in was not observed, so duration is measured from birth;
#'     the clock, still sitting at the birth date, gives exactly that.}
#'   \item{immigrant}{the same: measured from birth, not from the
#'     immigration date.}
#' }
#'
#' Once an individual is *observed* to change a variable, the clock holds a
#' real date and this starting assumption no longer applies.
#'
#' @section Assumptions:
#' For anyone not born during the simulation, we did not observe when they
#' entered the state they start in, so we assume they have occupied it
#' since birth. This means that a duration-dependent rate applied to the
#' starting population is evaluated at the largest duration the data
#' permit, so rates that rise with duration are overstated at the start of
#' a run and rates that fall with it are understated. The effect fades
#' as simulated transitions replace the assumed ones.
#'
#' @param n_vars Number of state variables in the model, i.e. the number
#'   of components in a composite state label such as `"f/DE/0"`. Each
#'   individual gets this many duration clocks.
#' @param capacity Initial number of rows. The store doubles when full, so
#'   this only affects performance.
#'
#' @return An object of class `person_store`: a list of functions.
#'   \describe{
#'     \item{`add(...)`}{register an individual, returning their slot}
#'     \item{`slot_of(id)`}{slot for a user-facing ID, or `NA_integer_`}
#'     \item{`move_to(slot, to_idx, now_days, changed_vars)`}{record a
#'       transition: set the new state and reset affected clocks}
#'     \item{`age_days(slot, now_days)`}{age in days}
#'     \item{`duration_days(slot, var_idx, now_days)`}{time since that
#'       state variable last changed, in days; `Inf` when `var_idx` is
#'       `NA`}
#'     \item{`state_of`, `id_of`, `birth_days_of`, `birth_ymd_of`,
#'       `init_state_of`, `is_migrant`,
#'       `mother_id_of`}{constant-time attribute reads}
#'     \item{`rng_event_of(slot)`, `rng_birth_of(slot)`,
#'       `set_rng_event(slot, state)`, `set_rng_birth(slot, state)`}{the
#'       individual's random-number stream positions, six integers each;
#'       see `R/rng_streams.R`.}
#'     \item{`next_id()`}{a fresh ID for a newborn}
#'     \item{`size()`}{number of individuals}
#'     \item{`as_columns()`}{the trimmed columns, for building the output}
#'   }
#'
#' @seealso [MicSim-concepts]
#'
#' @keywords internal
#' @noRd
new_person_store <- function(n_vars, capacity = 1024L) {

  n_vars <- as.integer(n_vars)
  if (is.na(n_vars) || n_vars < 1L) {
    micsim_abort(
      "micsim_person_store_error",
      "`n_vars` must be the number of state variables in the model (>= 1)."
    )
  }

  # ------------------------------------------------------------------
  # Private state. `parent = emptyenv()` makes it private: a
  # typo cannot quietly resolve to something else in the package.
  # ------------------------------------------------------------------
  e <- new.env(parent = emptyenv())

  e$n     <- 0L                          # individuals registered so far
  e$cap   <- as.integer(capacity)
  e$nvars <- n_vars

  e$id          <- integer(e$cap)        # user-facing identifier
  e$birth_days  <- numeric(e$cap)        # exact day count
  e$birth_ymd   <- character(e$cap)      # as supplied, for the output contract
  e$state       <- integer(e$cap)        # current state index (codec)
  e$init_state  <- integer(e$cap)        # state on entering the population
  e$is_migrant  <- logical(e$cap)
  e$mother_id   <- integer(e$cap)        # NA unless born during the run

  # One row per individual and one column per state variable: the date on
  # which that variable last changed.
  e$cov_change <- matrix(NA_real_, nrow = e$cap, ncol = n_vars)

  # Each individual's random-number stream (see rng_streams.R): the start
  # of the substream their latest scheduling used, and the current
  # position in their birth block. Six integers each, one row per person.
  e$rng_event <- matrix(NA_integer_, nrow = e$cap, ncol = 6L)
  e$rng_birth <- matrix(NA_integer_, nrow = e$cap, ncol = 6L)

  # ID -> slot. An ENVIRONMENT, because environments are hashed and the
  # names of an atomic vector are not; see the details section.
  e$slot_index <- new.env(parent = emptyenv(), hash = TRUE)

  e$max_id <- 0L                         # largest ID seen, for next_id()

  # ------------------------------------------------------------------
  # Growth by doubling: appending n individuals costs O(n) copying in
  # total rather than O(n^2). Same idea as the event queue's arrays.
  # ------------------------------------------------------------------
  grow <- function() {
    extra <- e$cap
    e$id         <- c(e$id,         integer(extra))
    e$birth_days <- c(e$birth_days, numeric(extra))
    e$birth_ymd  <- c(e$birth_ymd,  character(extra))
    e$state      <- c(e$state,      integer(extra))
    e$init_state <- c(e$init_state, integer(extra))
    e$is_migrant <- c(e$is_migrant, logical(extra))
    e$mother_id  <- c(e$mother_id,  integer(extra))
    e$cov_change <- rbind(e$cov_change,
                          matrix(NA_real_, nrow = extra, ncol = e$nvars))
    e$rng_event  <- rbind(e$rng_event,
                          matrix(NA_integer_, nrow = extra, ncol = 6L))
    e$rng_birth  <- rbind(e$rng_birth,
                          matrix(NA_integer_, nrow = extra, ncol = 6L))
    e$cap <- e$cap + extra
  }

  # Identifiers are used as environment keys, so they have to be turned
  # into strings. Going through as.integer() first avoids the scientific
  # notation that as.character() produces for large doubles, which would
  # make "1e+05" and "100000" different keys for the same person.
  id_key <- function(id) as.character(as.integer(id))

  # ------------------------------------------------------------------
  # Registering an individual
  # ------------------------------------------------------------------

  add <- function(id, birth_days, state_idx, init_state_idx,
                  birth_ymd = NA_character_, is_migrant = FALSE,
                  mother_id = NA_integer_, rng_state = NULL) {

    whole <- suppressWarnings(as.numeric(id))
    if (length(whole) != 1L || is.na(whole) || whole != round(whole) ||
          abs(whole) > .Machine$integer.max) {
      micsim_abort("micsim_person_store_error",
                   paste0("Person identifiers must be whole numbers; got `",
                          paste(id, collapse = ", "), "`."))
    }
    id <- as.integer(whole)

    key <- id_key(id)
    # get0() rather than slot_index[[key]]: `[[` on an environment raises
    # "subscript out of bounds" for a name that is not there, whereas a
    # missing identifier is the normal case here, not an error.
    if (!is.na(slot_of(id))) {
      micsim_abort(
        "micsim_person_store_error",
        paste0("Duplicate person identifier: ", id, ".\n",
               "  Identifiers must be unique across the starting population,\n",
               "  the immigrant population, and individuals born during the run.")
      )
    }

    if (e$n == e$cap) grow()
    slot <- e$n + 1L

    # Coerce every value FIRST, into locals, so the in-place writes that
    # follow contain no call that could fail. A failure between a column's
    # `<- NULL` and its restore would leave that column NULL. Nothing here
    # can error: these are coercions of already-validated scalars.
    id_v  <- id
    bd_v  <- as.numeric(birth_days)
    ymd_v <- if (is.na(birth_ymd)) as.character(getInDateFormat(birth_days))
             else as.character(birth_ymd)
    st_v  <- as.integer(state_idx)
    is_v  <- as.integer(init_state_idx)
    mg_v  <- isTRUE(is_migrant)
    mo_v  <- as.integer(mother_id)
    cc_v  <- as.numeric(birth_days)  # every clock starts at the date of birth
    # The start of the individual's stream: both the event cursor and the
    # birth block begin there (rng_streams.R). Stores built outside a run,
    # as in the unit tests, may leave it out.
    rs_v  <- if (is.null(rng_state)) rep(NA_integer_, 6L) else as.integer(rng_state)

    # WRITE IN PLACE via the "steal-back" idiom. `e$col[slot] <- v` looks
    # like an O(1) write, but the column is reachable from the captured
    # environment `e`, so its reference count is >= 2 and R DUPLICATES the
    # whole column on every write. Loading N individuals is then O(N^2)
    # (measured: ~11s for 40k). Dropping the store's own reference first
    # (`e$col <- NULL`) leaves the local `v` as the sole reference, so
    # `v[slot] <- ...` mutates in place; then we hand it back. This is the
    # documented base-R way to keep an environment-held buffer truly O(1),
    # with no change to the store's structure or API.
    v <- e$id;         e$id         <- NULL; v[slot]  <- id_v;  e$id         <- v
    v <- e$birth_days; e$birth_days <- NULL; v[slot]  <- bd_v;  e$birth_days <- v
    v <- e$birth_ymd;  e$birth_ymd  <- NULL; v[slot]  <- ymd_v; e$birth_ymd  <- v
    v <- e$state;      e$state      <- NULL; v[slot]  <- st_v;  e$state      <- v
    v <- e$init_state; e$init_state <- NULL; v[slot]  <- is_v;  e$init_state <- v
    v <- e$is_migrant; e$is_migrant <- NULL; v[slot]  <- mg_v;  e$is_migrant <- v
    v <- e$mother_id;  e$mother_id  <- NULL; v[slot]  <- mo_v;  e$mother_id  <- v
    m <- e$cov_change; e$cov_change <- NULL; m[slot, ] <- cc_v; e$cov_change <- m
    m <- e$rng_event;  e$rng_event  <- NULL; m[slot, ] <- rs_v; e$rng_event  <- m
    m <- e$rng_birth;  e$rng_birth  <- NULL; m[slot, ] <- rs_v; e$rng_birth  <- m

    assign(key, slot, envir = e$slot_index)
    e$n <- slot
    if (id > e$max_id) e$max_id <- id

    slot
  }

  slot_of <- function(id) {
    get0(id_key(id), envir = e$slot_index, inherits = FALSE,
         ifnotfound = NA_integer_)
  }

  next_id <- function() {
    e$max_id <- e$max_id + 1L
    e$max_id
  }

  # ------------------------------------------------------------------
  # Recording a transition
  # ------------------------------------------------------------------

  move_to <- function(slot, to_idx, now_days, changed_vars) {
    # Steal-back write, for the same reason as add(): a plain
    # `e$state[slot] <- .` copies the whole column each call because the
    # column is reachable from `e`. move_to() runs once per event, so left
    # unfixed it makes the run itself O(E * N), not just the load.
    to_v <- as.integer(to_idx)
    v <- e$state; e$state <- NULL; v[slot] <- to_v; e$state <- v

    # Reset the clock of every state variable this transition disturbed.
    # `changed_vars` is `integer(0)` for a move into an absorbing state,
    # whose label has no components to compare; an absorbed individual
    # never schedules again, so there is nothing to reset.
    if (length(changed_vars) > 0L) {
      now_v <- as.numeric(now_days)
      m <- e$cov_change; e$cov_change <- NULL
      m[slot, changed_vars] <- now_v; e$cov_change <- m
    }
    invisible(NULL)
  }

  # ------------------------------------------------------------------
  # The three questions asked on every event
  # ------------------------------------------------------------------

  age_days <- function(slot, now_days) {
    as.numeric(now_days) - e$birth_days[slot]
  }

  duration_days <- function(slot, var_idx, now_days) {
    # A transition with no single affected state variable has no
    # meaningful duration. That covers moves into absorbing states, where
    # the destination label is atomic, and moves that change two or more
    # variables at once. In both cases the duration is Inf, which is why a
    # mortality rate may declare a `duration` argument and simply ignore it.
    if (is.na(var_idx)) return(Inf)
    as.numeric(now_days) - e$cov_change[slot, var_idx]
  }

  # ------------------------------------------------------------------
  # Attribute reads. Trivial by design: an index into a vector.
  # ------------------------------------------------------------------

  state_of        <- function(slot) e$state[slot]
  id_of           <- function(slot) e$id[slot]
  birth_days_of   <- function(slot) e$birth_days[slot]
  birth_ymd_of    <- function(slot) e$birth_ymd[slot]
  init_state_of   <- function(slot) e$init_state[slot]
  is_migrant      <- function(slot) e$is_migrant[slot]
  mother_id_of    <- function(slot) e$mother_id[slot]
  cov_change_of   <- function(slot) e$cov_change[slot, , drop = TRUE]

  # Random-number stream positions, read and written once per scheduling
  # and once per birth. Writes use the steal-back idiom of add().
  rng_event_of <- function(slot) e$rng_event[slot, ]
  rng_birth_of <- function(slot) e$rng_birth[slot, ]
  set_rng_event <- function(slot, state) {
    m <- e$rng_event; e$rng_event <- NULL; m[slot, ] <- state; e$rng_event <- m
    invisible(NULL)
  }
  set_rng_birth <- function(slot, state) {
    m <- e$rng_birth; e$rng_birth <- NULL; m[slot, ] <- state; e$rng_birth <- m
    invisible(NULL)
  }

  # ------------------------------------------------------------------
  # Output
  # ------------------------------------------------------------------

  as_columns <- function() {
    keep <- seq_len(e$n)
    list(
      ID         = e$id[keep],
      birthDate  = e$birth_ymd[keep],
      birthDays  = e$birth_days[keep],
      initState  = e$init_state[keep],
      state      = e$state[keep],
      isMigrant  = e$is_migrant[keep],
      motherID   = e$mother_id[keep]
    )
  }

  structure(
    list(
      add           = add,
      slot_of       = slot_of,
      next_id       = next_id,
      move_to       = move_to,
      age_days      = age_days,
      duration_days = duration_days,
      state_of      = state_of,
      id_of         = id_of,
      birth_days_of = birth_days_of,
      birth_ymd_of  = birth_ymd_of,
      init_state_of = init_state_of,
      is_migrant    = is_migrant,
      mother_id_of  = mother_id_of,
      cov_change_of = cov_change_of,
      rng_event_of  = rng_event_of,
      rng_birth_of  = rng_birth_of,
      set_rng_event = set_rng_event,
      set_rng_birth = set_rng_birth,
      as_columns    = as_columns,
      size          = function() e$n,
      n_vars        = function() e$nvars,
      max_id        = function() e$max_id,
      state         = e            # exposed for tests and diagnostics only
    ),
    class = "person_store"
  )
}


#' Print a person store
#'
#' Reports how many individuals are registered, how many clocks each one
#' carries, and how many arrived after the simulation started. This is a
#' quick way to check that immigration and fertility did something.
#'
#' @param x A `person_store`.
#' @param ... Not used; present for compatibility with the [print()] generic.
#' @return `x`, invisibly. Called for the side effect of printing.
#' @exportS3Method print person_store
#' @noRd
print.person_store <- function(x, ...) {
  n <- x$size()
  cat("<person_store>\n")
  cat("  individuals    : ", n, "\n", sep = "")
  cat("  duration clocks: ", x$n_vars(), " per individual\n", sep = "")
  if (n > 0L) {
    cols <- x$as_columns()
    cat("  immigrants     : ", sum(cols$isMigrant), "\n", sep = "")
    cat("  born in run    : ", sum(!is.na(cols$motherID)), "\n", sep = "")
    cat("  largest ID     : ", x$max_id(), "\n", sep = "")
  }
  invisible(x)
}
