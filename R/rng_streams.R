# =====================================================================
# rng_streams.R
# ---------------------------------------------------------------------
# WHAT THIS FILE DOES
#   Gives every individual their own random-number stream, so that what
#   happens to one person never depends on what happened to anybody else,
#   in which order the engine processed them, or on which core they ran.
#
# PLACE IN THE PIPELINE
#   engine.R run() --(seed)--> open_run_rng() --(master state)-->
#   person_stream_states() --(one state per person)--> engine_person.R
#   engine.R schedule_next() --> next_event_uniforms() --> solver.R
#   subsystem_fertility.R    --> birth_uniforms() --> a newborn's stream
#
# WHY PER-PERSON STREAMS
#   With one global stream, the k-th random number goes to whoever the
#   event loop happens to schedule k-th. Two runs of the same population
#   under slightly different rates schedule people in a different order,
#   so almost everyone receives different random numbers and individual
#   life courses cannot be compared between scenarios (no common random
#   numbers). With one stream per person, a person whose rates did not
#   change lives exactly the same life in both runs.
#
# THE LAYOUT
#   The generator is L'Ecuyer-CMRG (MRG32k3a), the one R's `parallel`
#   package uses. Its cycle is cut into streams of length 2^127, and every
#   stream into substreams of length 2^76.
#
#   * The run's MASTER state comes from the seed (stream 0). It is the
#     ambient generator during a run; the engine itself draws nothing from
#     it.
#   * A supplied individual (starting population or immigrant) with
#     identifier `ID` owns stream number ID + 2^31, which is between 1
#     and 2^32 - 1 for every admissible identifier. It depends on the seed
#     and the identifier only, not on who else is in the population.
#   * Within a person's stream, substream 0 is the BIRTH block: every child
#     the person has takes the next seven uniforms from it, six to seed the
#     child's own stream and one to draw the child's state. Substream k,
#     k >= 1, holds the draws of the person's k-th scheduling: one uniform
#     per possible destination, assigned to the destinations in the order
#     of the state space (not the order the transitions were declared in).
#   * A newborn's stream starts at a point of the generator's cycle chosen
#     by those six uniforms. Its layout is the same as everybody else's.
#
# GUARANTEES
#   * A given seed always produces the same run, on any number of cores.
#   * An individual's random numbers depend only on the seed, their
#     identifier (or, for a newborn, their mother and birth order), and
#     their own history: never on other people or on processing order.
#   * A forced transition, a zero rate or an exhausted window still uses
#     its uniform, so the draws of later destinations and schedulings do
#     not move.
#   * The caller's random-number state is left exactly as it was when a
#     seed is supplied, including the generator kind and the absence of
#     `.Random.seed`. Without a seed, the run takes exactly one draw from
#     the caller's stream (the seed it then uses) and nothing else.
#
# UNITS
#   States are L'Ecuyer-CMRG seeds without the leading kind code: six
#   integers, three per component, oldest first, in R's signed 32-bit
#   storage. The arithmetic below works on the same values as unsigned
#   numbers held in doubles, where every intermediate result is an exact
#   integer below 2^53.
# =====================================================================


#' Kind code of the generator in `.Random.seed[1]`
#'
#' The code R stores in `.Random.seed[1]` for L'Ecuyer-CMRG with the
#' Inversion normal generator and the Rejection sampler. Assigning a state
#' with this code to `.Random.seed` switches R to that generator on its
#' next use.
#'
#' @keywords internal
#' @noRd
LECUYER_KIND_CODE <- 10407L

#' Moduli of the two MRG32k3a components
#'
#' `LECUYER_M1` = 2^32 - 209 and `LECUYER_M2` = 2^32 - 22853. Every word of
#' a valid state is below the modulus of its component.
#'
#' @keywords internal
#' @noRd
LECUYER_M1 <- 4294967087

#' @rdname LECUYER_M1
#' @noRd
LECUYER_M2 <- 4294944443

#' Jump matrices one stream ahead
#'
#' `A1^(2^127)` and `A2^(2^127)`, the transition matrices of the two
#' components raised to the stream length, from L'Ecuyer, Simard, Chen and
#' Kelton (2002), "An object-oriented random-number package with many long
#' streams and substreams", *Operations Research* 50(6). They are the
#' matrices `parallel::nextRNGStream()` applies.
#'
#' @keywords internal
#' @noRd
LECUYER_A1P127 <- matrix(c(2427906178, 3580155704,  949770784,
                            226153695, 1230515664, 3580155704,
                           1988835001,  986791581, 1230515664),
                         nrow = 3L, byrow = TRUE)

#' @rdname LECUYER_A1P127
#' @noRd
LECUYER_A2P127 <- matrix(c(1464411153,  277697599, 1610723613,
                             32183930, 1464411153, 1022607788,
                           2824425944,   32183930, 2093834863),
                         nrow = 3L, byrow = TRUE)

#' Offset from an identifier to a stream number
#'
#' Individual `ID` owns stream number `ID + STREAM_ID_OFFSET`. Identifiers
#' lie between `-.Machine$integer.max` and `.Machine$integer.max`, so the
#' stream number lies between 1 and 2^32 - 1; stream 0 is the master.
#'
#' @keywords internal
#' @noRd
STREAM_ID_OFFSET <- 2^31


# ---------------------------------------------------------------------
# Opening and closing a run
# ---------------------------------------------------------------------

#' Check a seed
#'
#' @description
#' Accepts `NULL` or a single finite whole number in R's integer range,
#' which is what [set.seed()] can use without changing it.
#'
#' @param seed The `seed` argument of a run.
#' @return `TRUE`, invisibly. Otherwise an error of class
#'   `micsim_model_error`.
#' @keywords internal
#' @noRd
check_seed <- function(seed) {
  if (is.null(seed)) return(invisible(TRUE))
  ok <- is.numeric(seed) && length(seed) == 1L && !is.na(seed) &&
    is.finite(seed) && seed == round(seed) &&
    abs(seed) <= .Machine$integer.max
  if (!ok) {
    micsim_abort("micsim_model_error", paste0(
      "`seed` must be NULL or a single whole number between\n",
      "  -", .Machine$integer.max, " and ", .Machine$integer.max, "."
    ))
  }
  invisible(TRUE)
}


#' Set up the random numbers of one run
#'
#' @description
#' Records the caller's random-number state, settles the seed and makes the
#' run's master state the ambient generator.
#'
#' @details
#' Without a seed, one is drawn from the caller's stream with
#' `sample.int(.Machine$integer.max, 1L)`. This is the run's only draw from
#' that stream; it also creates `.Random.seed` if the caller had none, so
#' the state saved afterwards includes the draw.
#'
#' Restoring assigns the saved `.Random.seed` back. That also restores the
#' generator kinds, because R reads them from `.Random.seed[1]` on its next
#' use. When the caller had no `.Random.seed`, the kinds are reset with
#' [RNGkind()] and the object is removed again, so the caller's next random
#' number is seeded from the clock exactly as it would have been.
#'
#' @param seed `NULL` or a single whole number; see `check_seed()`.
#' @return A list with
#'   \describe{
#'     \item{`seed`}{the integer seed the run uses;}
#'     \item{`drawn`}{`TRUE` when the seed was drawn from the caller's
#'       stream;}
#'     \item{`master`}{the master state, six integers;}
#'     \item{`restore`}{a function of no arguments that puts the caller's
#'       state back. Call it on exit, whatever happens.}
#'   }
#' @keywords internal
#' @noRd
open_run_rng <- function(seed = NULL) {
  check_seed(seed)
  drawn <- is.null(seed)
  if (drawn) {
    # The run's only draw from the caller's stream. It also creates
    # .Random.seed if the caller had none, so it is saved below.
    seed <- sample.int(.Machine$integer.max, 1L)
  }
  seed <- as.integer(seed)

  global <- globalenv()
  had_state <- exists(".Random.seed", envir = global, inherits = FALSE)
  saved_state <- if (had_state) get(".Random.seed", envir = global,
                                    inherits = FALSE)
  saved_kind <- RNGkind()

  restore <- function() {
    if (had_state) {
      assign(".Random.seed", saved_state, envir = globalenv())
    } else {
      # "Rounding" warns that it is non-uniform; the caller chose it.
      suppressWarnings(RNGkind(kind = saved_kind[1L],
                               normal.kind = saved_kind[2L],
                               sample.kind = saved_kind[3L]))
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
        rm(".Random.seed", envir = globalenv())
      }
    }
    invisible(NULL)
  }

  set.seed(seed, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion",
           sample.kind = "Rejection")
  master <- get(".Random.seed", envir = global, inherits = FALSE)[-1L]

  list(seed = seed, drawn = drawn, master = master, restore = restore)
}


# ---------------------------------------------------------------------
# Streams of supplied individuals
# ---------------------------------------------------------------------

#' The starting state of every supplied individual's stream
#'
#' @description
#' Individual `ID` owns stream number `ID + 2^31` counted from the master
#' state, i.e. the state that `parallel::nextRNGStream()` would reach after
#' that many applications. The state depends on the master state and the
#' identifier only.
#'
#' @details
#' Jumping `n` streams ahead multiplies each component by `A^(2^127 n)`.
#' The powers `A^(2^127 2^k)`, `k = 0, ..., 31`, are formed once by repeated
#' squaring. When the identifiers are dense (their range is at most 16
#' times their number, as for the usual `1:N`), the smallest is reached by
#' its bits and the whole range is then filled by doubling, so each stream
#' costs one matrix-vector product. Otherwise every identifier is jumped by
#' its own bits (`jump_by_bits()`). All individuals are handled at once.
#'
#' @param master The master state, six integers (`open_run_rng()$master`).
#' @param ids Whole-number identifiers with absolute value at most
#'   `.Machine$integer.max`.
#' @return An integer matrix with one row per identifier and six columns,
#'   in R's `.Random.seed` storage (without the kind code).
#' @keywords internal
#' @noRd
person_stream_states <- function(master, ids) {
  n <- length(ids)
  out <- matrix(0L, nrow = n, ncol = 6L)
  if (n == 0L) return(out)

  index <- as.numeric(ids) + STREAM_ID_OFFSET
  first <- min(index)
  span <- max(index) - first + 1
  state <- int32_to_unsigned(master)

  for (comp in 1:2) {
    cols <- if (comp == 1L) 1:3 else 4:6
    base <- if (comp == 1L) LECUYER_A1P127 else LECUYER_A2P127
    mod <- if (comp == 1L) LECUYER_M1 else LECUYER_M2
    # powers[[k + 1]] = A^(2^127 * 2^k): k streams ahead, for k = 0..31.
    powers <- vector("list", 32L)
    powers[[1L]] <- base
    for (k in 2:32) powers[[k]] <- matmul_mod(powers[[k - 1L]], powers[[k - 1L]], mod)

    if (span <= 16 * n) {
      # Identifiers are dense (the usual 1..N): reach the smallest one by
      # its bits, then fill the whole range by doubling, so that every
      # stream in the range costs one matrix-vector product.
      s <- jump_by_bits(matrix(state[cols], nrow = 1L), first, powers, mod)
      k <- 1L
      while (nrow(s) < span) {
        s <- rbind(s, matvec_mod(powers[[k]], s, mod))
        k <- k + 1L
      }
      s <- s[index - first + 1, , drop = FALSE]
    } else {
      # Sparse identifiers: jump every one by its own bits.
      s <- jump_by_bits(matrix(state[cols], nrow = n, ncol = 3L, byrow = TRUE),
                        index, powers, mod)
    }
    out[, cols] <- unsigned_to_int32(s)
  }
  out
}


#' Jump states ahead by a number of streams
#'
#' @description
#' Advances each row of `S` by `index` streams, applying the powers whose
#' bit is set in the index.
#'
#' @param S Numeric matrix with three columns, one component's state per
#'   row, as unsigned values.
#' @param index Stream counts, one per row of `S` or one for all rows;
#'   whole numbers below 2^32.
#' @param powers List of 32 matrices, `A^(2^127 2^k)` for `k = 0, ..., 31`.
#' @param mod The component's modulus.
#' @return `S` with every row advanced.
#' @keywords internal
#' @noRd
jump_by_bits <- function(S, index, powers, mod) {
  for (k in 0:31) {
    bit <- floor(index / 2^k) %% 2 == 1
    if (length(bit) == 1L) bit <- rep(bit, nrow(S))
    if (any(bit)) S[bit, ] <- matvec_mod(powers[[k + 1L]], S[bit, , drop = FALSE], mod)
  }
  S
}


# ---------------------------------------------------------------------
# Drawing
# ---------------------------------------------------------------------

#' Draw uniforms from a given state
#'
#' @description
#' Draws `n` uniforms with R's own L'Ecuyer-CMRG generator started at
#' `state`, then puts the ambient `.Random.seed` back (or removes it if
#' there was none).
#'
#' @param state Six integers, a generator state without the kind code.
#' @param n Number of uniforms.
#' @return A list with `u`, the uniforms, and `state`, the state after them.
#' @keywords internal
#' @noRd
draw_from_state <- function(state, n) {
  global <- globalenv()
  ambient <- get0(".Random.seed", envir = global, inherits = FALSE)
  assign(".Random.seed", c(LECUYER_KIND_CODE, state), envir = globalenv())
  u <- stats::runif(n)
  after <- get(".Random.seed", envir = global, inherits = FALSE)[-1L]
  if (is.null(ambient)) {
    rm(".Random.seed", envir = globalenv())
  } else {
    assign(".Random.seed", ambient, envir = globalenv())
  }
  list(u = u, state = after)
}


#' The uniforms of one scheduling
#'
#' @description
#' Moves the individual to their next event substream, records it in the
#' person store, and draws one uniform per destination there. The k-th
#' scheduling of an individual therefore always starts at substream k,
#' however many uniforms earlier schedulings used.
#'
#' @param persons A `person_store`.
#' @param slot The individual's slot.
#' @param n Number of destinations.
#' @return `n` uniforms.
#' @keywords internal
#' @noRd
next_event_uniforms <- function(persons, slot, n) {
  start <- next_substream(persons$rng_event_of(slot))
  persons$set_rng_event(slot, start)
  draw_from_state(start, n)$u
}


#' The next uniforms of an individual's birth block
#'
#' @description
#' Draws from substream 0 of the individual's stream, continuing where the
#' previous birth stopped, and records the new position.
#'
#' @param persons A `person_store`.
#' @param slot The mother's slot.
#' @param n Number of uniforms; seven per birth.
#' @return `n` uniforms.
#' @keywords internal
#' @noRd
birth_uniforms <- function(persons, slot, n) {
  d <- draw_from_state(persons$rng_birth_of(slot), n)
  persons$set_rng_birth(slot, d$state)
  d$u
}


#' The start of the next substream
#'
#' @param state Six integers, the start of a substream.
#' @return Six integers, the start of the substream after it
#'   (`parallel::nextRNGSubStream()` without the kind code).
#' @keywords internal
#' @noRd
next_substream <- function(state) {
  parallel::nextRNGSubStream(c(LECUYER_KIND_CODE, state))[-1L]
}


#' A generator state built from six uniforms
#'
#' @description
#' Gives a newborn their own stream, starting at a point of the
#' generator's cycle chosen by six uniforms from the mother's birth block.
#'
#' @details
#' The state must be valid: three words below `LECUYER_M1`, three below
#' `LECUYER_M2`, neither triple all zero. R silently reseeds an invalid
#' state from the clock, which would make the run irreproducible. A triple
#' that comes out all zero (probability about 2^-96) gets a 1 in its first
#' word.
#'
#' @param u Six uniforms in (0, 1).
#' @return Six integers in `.Random.seed` storage.
#' @keywords internal
#' @noRd
state_from_uniforms <- function(u) {
  x <- floor(u[1:3] * LECUYER_M1)
  y <- floor(u[4:6] * LECUYER_M2)
  if (all(x == 0)) x[1L] <- 1
  if (all(y == 0)) y[1L] <- 1
  unsigned_to_int32(c(x, y))
}


#' Exponential budgets for one scheduling
#'
#' @description
#' Turns the uniforms of a scheduling into one standard exponential budget
#' per feasible edge, by inversion (`-log(1 - u)`, as `default_draw_u()`).
#'
#' @param u The uniforms, in the order of the destinations in the state
#'   space.
#' @param position `position[j]` is the rank, in that order, of the
#'   destination of the j-th feasible edge.
#' @return One budget per feasible edge, in the registry's order.
#' @keywords internal
#' @noRd
budgets_from_uniforms <- function(u, position) {
  -log(1 - u[position])
}


# ---------------------------------------------------------------------
# Exact modular arithmetic on 32-bit values held in doubles
# ---------------------------------------------------------------------

#' Exact modular product of 32-bit values
#'
#' @description
#' `(a * b) mod m` for `0 <= a, b < 2^32`, without leaving the range in
#' which doubles hold integers exactly: `a` is split into 16-bit halves, so
#' no product exceeds 2^48.
#'
#' @param a,b Whole numbers below 2^32; vectors are recycled.
#' @param m The modulus, below 2^32.
#' @return `(a * b) mod m`.
#' @keywords internal
#' @noRd
mulmod <- function(a, b, m) {
  a_hi <- floor(a / 65536)
  a_lo <- a - a_hi * 65536
  ((a_hi * b) %% m * 65536 + a_lo * b) %% m
}

#' Exact modular product of 3 x 3 matrices
#'
#' @param A,B 3 x 3 matrices of whole numbers below `m`.
#' @param m The modulus, below 2^32.
#' @return `A %*% B` mod `m`.
#' @keywords internal
#' @noRd
matmul_mod <- function(A, B, m) {
  C <- matrix(0, 3L, 3L)
  for (i in 1:3) for (j in 1:3) {
    C[i, j] <- (mulmod(A[i, 1L], B[1L, j], m) + mulmod(A[i, 2L], B[2L, j], m) +
                  mulmod(A[i, 3L], B[3L, j], m)) %% m
  }
  C
}

#' Exact modular matrix-vector product, row by row
#'
#' @description
#' Replaces each row `s` of `S` by `(A %*% s) mod m`. Every term of the sum
#' is below 2^49, so the sum of three is exact before the one reduction.
#'
#' @param A A 3 x 3 matrix of whole numbers below `m`.
#' @param S A matrix with three columns of whole numbers below `m`.
#' @param m The modulus, below 2^32.
#' @return `S` with every row transformed.
#' @keywords internal
#' @noRd
matvec_mod <- function(A, S, m) {
  A_hi <- floor(A / 65536)
  A_lo <- A - A_hi * 65536
  out <- S
  for (i in 1:3) {
    acc <- 0
    for (j in 1:3) {
      acc <- acc + (A_hi[i, j] * S[, j]) %% m * 65536 + A_lo[i, j] * S[, j]
    }
    out[, i] <- acc %% m
  }
  out
}

#' Convert seed words between R's storage and unsigned values
#'
#' @description
#' R stores the unsigned 32-bit words of a generator state as signed
#' integers: a word of 2^31 or more appears as a negative number, and the
#' word 2^31 itself as `NA_integer_` (the bit pattern of `INT_MIN`).
#' `int32_to_unsigned()` turns stored words into their unsigned values,
#' held in doubles; `unsigned_to_int32()` does the reverse and keeps the
#' dimensions of a matrix.
#'
#' @param x Integer vector or matrix in `.Random.seed` storage.
#' @param v Numeric vector or matrix of whole numbers in `[0, 2^32)`.
#' @return The converted values.
#' @keywords internal
#' @noRd
int32_to_unsigned <- function(x) {
  v <- as.numeric(x)
  v[is.na(x)] <- 2^31
  v[!is.na(x) & v < 0] <- v[!is.na(x) & v < 0] + 2^32
  v
}

#' @rdname int32_to_unsigned
#' @noRd
unsigned_to_int32 <- function(v) {
  dims <- dim(v)
  is_min <- v == 2^31
  signed <- ifelse(v >= 2^31, v - 2^32, v)
  signed[is_min] <- 0
  out <- as.integer(signed)
  out[is_min] <- NA_integer_
  dim(out) <- dims
  out
}
