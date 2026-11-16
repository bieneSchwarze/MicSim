# =====================================================================
# test-fertility.R
# ---------------------------------------------------------------------
# Fertility is the first subsystem, so this file tests two things at once:
# that births happen correctly, and that the seam they happen through
# works. The second is the reason fertility was written as a subsystem
# rather than as a branch of the driver; see the header of
# R/subsystem_fertility.R.
#
# Tests are grouped by what they establish:
#   1. the newborn rule: construction and validation
#   2. the newborn rule: drawing a state
#   3. the subsystem in isolation
#   4. births end to end
#   5. the output contract with fertility on
#   6. the seam itself: composition, ordering, user subsystems
# =====================================================================


# ---------------------------------------------------------------------
# Fixtures
#
# Two state variables, sex and parity, plus an absorbing "dead". States
# are "m/0", "m/1", "f/0", "f/1" (indices 1..4) and "dead" (5). Women at
# parity 0 give birth; nobody dies unless a test says so.
# ---------------------------------------------------------------------

fert_codec <- function() {
  new_state_codec(list(sex = c("m", "f"), fert = c("0", "1")),
                  absorbing = "dead")
}

#' Every newborn is born childless; sex is drawn with a realistic ratio
simple_newborn_rule <- function() {
  new_newborn_rule(
    states = rbind(c("m", "0"), c("f", "0")),
    prob   = c(0.515, 0.485)
  )
}

fert_registry <- function(codec = fert_codec(), rate_value = 0.2) {
  births <- function(age, calTime) rep(rate_value, length(age))
  new_rate_registry(
    list(list(from = "f/0", to = "f/1", rate = births)),
    codec,
    probe = micsim_rate_probe(0:120, 2000:2040, 0:1)
  )
}

fert_init_pop <- function(n = 50L) {
  data.frame(
    ID        = seq_len(n),
    birthDate = rep("19800101", n),
    initState = rep("f/0", n),
    stringsAsFactors = FALSE
  )
}

fert_subsystem <- function(codec, registry, rule = simple_newborn_rule(),
                           births = "f/0->f/1") {
  compiled(fertility(births = births, newborn = rule), codec, registry)
}

fert_engine <- function(n = 50L, horizon = c(20000101, 20101231), ...) {
  codec <- fert_codec()
  registry <- fert_registry(codec)
  new_engine(
    codec, registry, fert_init_pop(n),
    horizon = horizon, max_age = 120,
    subsystems = list(fert_subsystem(codec, registry)),
    ...
  )
}


# ---------------------------------------------------------------------
# 1. The newborn rule: construction and validation
# ---------------------------------------------------------------------

test_that("a well-formed rule is accepted", {
  rule <- simple_newborn_rule()
  expect_s3_class(rule, "newborn_rule")
  expect_identical(rule$n_vars, 2L)
  expect_identical(rule$inherit, integer(0))
})


test_that("probabilities must sum to one when nothing is inherited", {
  expect_error(
    new_newborn_rule(states = rbind(c("m", "0"), c("f", "0")),
                     prob = c(0.5, 0.4)),
    class = "micsim_newborn_rule_error"
  )
})


test_that("probabilities must sum to one within each inherited group", {
  # This is the check MicSim 3.0.0 got wrong. Its version looked only at
  # the FIRST value of the inherited column, so it validated the "DE"
  # group three times and never looked at "AT" or "IT". The table below,
  # whose AT group sums to 0.8, passed there and is rejected here.
  bad <- new_newborn_rule
  expect_error(
    bad(
      states = rbind(c("m", "DE", "0"), c("f", "DE", "0"),
                     c("m", "AT", "0"), c("f", "AT", "0")),
      prob    = c(0.515, 0.485, 0.400, 0.400),
      inherit = 2L
    ),
    class = "micsim_newborn_rule_error"
  )
})


test_that("a correct conditional table is accepted", {
  rule <- new_newborn_rule(
    states = rbind(c("m", "DE", "0"), c("f", "DE", "0"),
                   c("m", "AT", "0"), c("f", "AT", "0"),
                   c("m", "IT", "0"), c("f", "IT", "0")),
    prob    = c(0.515, 0.485, 0.515, 0.485, 0.515, 0.485),
    inherit = 2L
  )
  expect_identical(nrow(rule$states), 6L)
  expect_identical(rule$inherit, 2L)
})


test_that("malformed inputs are rejected with a classed error", {
  expect_error(new_newborn_rule(rbind(c("m", "0")), prob = c(1, 1)),
               class = "micsim_newborn_rule_error")          # wrong length
  expect_error(new_newborn_rule(rbind(c("m", "0")), prob = -1),
               class = "micsim_newborn_rule_error")          # negative
  expect_error(new_newborn_rule(rbind(c("m", "0")), prob = 1, inherit = 5L),
               class = "micsim_newborn_rule_error")          # no such column
})


# ---------------------------------------------------------------------
# 2. The newborn rule: drawing a state
# ---------------------------------------------------------------------

test_that("a drawn state comes from the candidate table", {
  rule <- simple_newborn_rule()
  set.seed(1L)
  drawn <- replicate(50L, assign_newborn_state(rule, "f/0"))
  expect_true(all(drawn %in% c("m/0", "f/0")))
})


test_that("the sex ratio of newborns matches the weights", {
  skip_on_cran()
  rule <- simple_newborn_rule()
  set.seed(20260917L)
  drawn <- replicate(4000L, assign_newborn_state(rule, "f/0"))
  expect_equal(mean(drawn == "m/0"), 0.515, tolerance = 0.03)
})


test_that("inherited attributes are taken from the mother", {
  # An Austrian mother can only have an Austrian child, whatever the draw.
  rule <- new_newborn_rule(
    states = rbind(c("m", "DE", "0"), c("f", "DE", "0"),
                   c("m", "AT", "0"), c("f", "AT", "0")),
    prob    = c(0.515, 0.485, 0.515, 0.485),
    inherit = 2L
  )
  set.seed(2L)
  drawn <- replicate(50L, assign_newborn_state(rule, "f/AT/1"))
  expect_true(all(drawn %in% c("m/AT/0", "f/AT/0")))
  expect_true(all(grepl("/AT/", drawn, fixed = TRUE)))
})


test_that("a mother whose inherited combination has no candidate is an error", {
  rule <- new_newborn_rule(
    states = rbind(c("m", "DE", "0"), c("f", "DE", "0")),
    prob = c(0.515, 0.485), inherit = 2L
  )
  expect_error(assign_newborn_state(rule, "f/IT/1"),
               class = "micsim_newborn_rule_error")
})


test_that("the newborn's state is the inversion of one uniform", {
  # The fertility subsystem passes the uniform from the mother's birth
  # block; the candidates are taken in table order, with probability mass
  # 0.515 and 0.485 here.
  rule <- simple_newborn_rule()
  states <- vapply(c(0.001, 0.514, 0.516, 0.999), function(u)
    assign_newborn_state(rule, "f/0", u = u), character(1L))
  first <- paste(rule$states[1L, ], collapse = "/")
  second <- paste(rule$states[2L, ], collapse = "/")
  expect_identical(states, c(first, first, second, second))

  # Without `u`, exactly one uniform is drawn from the ambient stream.
  set.seed(3L)
  assign_newborn_state(rule, "f/0")
  after_rule <- .Random.seed
  set.seed(3L)
  invisible(stats::runif(1L))
  expect_identical(after_rule, .Random.seed)
})


test_that("a candidate with probability zero is never drawn", {
  rule <- new_newborn_rule(rbind(c("m", "0"), c("x", "0"), c("f", "0")),
                           prob = c(0.5, 0, 0.5))
  drawn <- vapply(c(1e-12, 0.5, 0.5 + 1e-12, 1 - 1e-12), function(u)
    assign_newborn_state(rule, "f/0", u = u), character(1L))
  expect_false(any(drawn == "x/0"))
})


# ---------------------------------------------------------------------
# 3. The subsystem in isolation
# ---------------------------------------------------------------------

test_that("the compiled subsystem is an ordinary hook set", {
  # Nothing about it is privileged. It is exactly what a user would write.
  codec <- fert_codec()
  sub <- fert_subsystem(codec, fert_registry(codec))
  expect_s3_class(sub, "micsim_compiled_subsystem")
  expect_true(sub$spawns)
  expect_true(sub$parallel_safe)
  inst <- sub$instantiate(list(codec = codec))
  expect_s3_class(inst$hooks, "engine_hooks")
  expect_identical(inst$hooks$active, "on_event")
})


test_that("transition labels are normalised to the registry's spacing", {
  # MicSim 3.0.0 writes "f/0->f/1"; the registry builds "f/0 -> f/1".
  expect_identical(as_edge_labels("f/0->f/1"), "f/0 -> f/1")
  expect_identical(as_edge_labels("f/0 -> f/1"), "f/0 -> f/1")
  expect_identical(as_edge_labels(character(0)), character(0))
  expect_error(as_edge_labels("not a transition"),
               class = "micsim_selector_error")
})


test_that("births can be named with a selector instead of labels", {
  codec <- fert_codec()
  registry <- fert_registry(codec)
  by_label <- fert_subsystem(codec, registry)
  by_edges <- fert_subsystem(codec, registry, births = moves_of("fert"))
  by_edge  <- fert_subsystem(codec, registry,
                             births = move("fert", "0", "1"))
  expect_identical(by_edges$claims, by_label$claims)
  expect_identical(by_edge$claims, by_label$claims)
})


test_that("a newborn state outside the state space is caught at build time", {
  rule <- new_newborn_rule(states = rbind(c("x", "0")), prob = 1)
  codec <- fert_codec()
  expect_error(fert_subsystem(codec, fert_registry(codec), rule),
               class = "micsim_newborn_rule_error")
})


test_that("a rule with the wrong number of state variables is caught", {
  rule <- new_newborn_rule(states = rbind(c("m")), prob = 1)
  codec <- fert_codec()
  expect_error(fert_subsystem(codec, fert_registry(codec), rule),
               class = "micsim_newborn_rule_error")
})


test_that("fertility() needs both halves", {
  expect_error(fertility(births = "f/0->f/1"),
               class = "micsim_fertility_error")
  expect_error(fertility(newborn = simple_newborn_rule()),
               class = "micsim_fertility_error")
})


test_that("births that are not declared transitions are refused", {
  codec <- fert_codec()
  registry <- fert_registry(codec)
  expect_error(fert_subsystem(codec, registry, births = "m/0->m/1"),
               class = "micsim_selector_error")
  expect_error(fert_subsystem(codec, registry, births = moves_of("sex")),
               class = "micsim_selector_error")
})


test_that("an optional `mothers` condition checks every birth transition", {
  codec <- fert_codec()
  registry <- fert_registry(codec)
  rule <- simple_newborn_rule()
  ok <- compiled(fertility(moves_of("fert"), rule, mothers = sex == "f"),
                 codec, registry)
  expect_length(ok$claims$births, 1L)
  expect_error(
    compiled(fertility(moves_of("fert"), rule, mothers = sex == "m"),
             codec, registry),
    class = "micsim_fertility_error"
  )
  expect_error(
    compiled(fertility(moves_of("fert"), rule, mothers = sex == "female"),
             codec, registry),
    class = "micsim_selector_error"         # matches no state at all
  )
})


test_that("the per-variable newborn rule is checked and expanded", {
  codec <- fert_codec()
  rule <- resolve_newborn_rule(
    newborn_rule(sex = c(m = 0.515, f = 0.485), fert = "0"), codec
  )
  expect_s3_class(rule, "newborn_rule")
  expect_identical(apply(rule$states, 1L, paste, collapse = "/"),
                   c("m/0", "f/0"))
  expect_equal(rule$prob, c(0.515, 0.485))

  expect_error(resolve_newborn_rule(newborn_rule(sex = "m"), codec),
               class = "micsim_newborn_rule_error")              # fert missing
  expect_error(resolve_newborn_rule(
    newborn_rule(sex = "x", fert = "0"), codec),
    class = "micsim_newborn_rule_error")                          # unknown value
  expect_error(newborn_rule(sex = c(m = 0.5, f = 0.4), fert = "0"),
               class = "micsim_newborn_rule_error")              # sums to 0.9
  expect_error(newborn_rule(c(m = 1)), class = "micsim_newborn_rule_error")
})


test_that("from_mother() inherits a variable, conditionally on each value", {
  codec <- new_state_codec(list(sex = c("m", "f"), nat = c("DE", "AT"),
                                fert = c("0", "1")))
  rule <- resolve_newborn_rule(
    newborn_rule(sex = c(m = 0.5, f = 0.5), nat = from_mother(), fert = "0"),
    codec
  )
  expect_identical(rule$inherit, 2L)
  set.seed(1L)
  kids <- vapply(1:50, function(i) assign_newborn_state(rule, "f/AT/1"),
                 character(1L))
  expect_true(all(grepl("^[mf]/AT/0$", kids)))
})


test_that("the per-variable rule reproduces the equivalent table", {
  # Same candidates, same order, same weights: same newborns.
  codec <- fert_codec()
  table_rule <- simple_newborn_rule()
  spec_rule <- resolve_newborn_rule(
    newborn_rule(sex = c(m = 0.515, f = 0.485), fert = "0"), codec
  )
  expect_identical(spec_rule$states, table_rule$states)
  expect_equal(spec_rule$prob, table_rule$prob)
})


# ---------------------------------------------------------------------
# 4. Births end to end
# ---------------------------------------------------------------------

test_that("a birth adds an individual to the population", {
  set.seed(4L)
  result <- fert_engine(n = 50L)$run()

  expect_gt(result$diagnostics$n_individuals, 50L)
  expect_setequal(unique(result$pop$ID), seq_len(result$diagnostics$n_individuals))
})


test_that("newborns get identifiers above the starting population's", {
  set.seed(5L)
  result <- fert_engine(n = 50L)$run()
  pop <- result$pop
  newborns <- pop$ID[!is.na(pop$motherID)]
  expect_true(all(newborns > 50L))
})


test_that("every newborn has a mother who gave birth", {
  set.seed(6L)
  pop <- fert_engine(n = 50L)$run()$pop

  newborn_rows <- pop[!is.na(pop$motherID), , drop = FALSE]
  expect_gt(nrow(newborn_rows), 0L)

  mothers <- unique(newborn_rows$motherID)
  for (mother in mothers) {
    her_rows <- pop[pop$ID == mother, , drop = FALSE]
    expect_true(any(her_rows$From == "f/0" & her_rows$To == "f/1",
                    na.rm = TRUE))
  }
})


test_that("a newborn's birth date is the date of the birth event", {
  set.seed(7L)
  pop <- fert_engine(n = 50L)$run()$pop
  newborn_rows <- pop[!is.na(pop$motherID), , drop = FALSE]
  expect_gt(nrow(newborn_rows), 0L)

  for (i in seq_len(nrow(newborn_rows))) {
    child <- newborn_rows[i, ]
    mother_events <- pop[pop$ID == child$motherID &
                           !is.na(pop$To) & pop$To == "f/1", ]
    expect_identical(child$birthDate, mother_events$transitionTime[1L])
  }
})


test_that("a newborn is at risk from birth, not from the start of the run", {
  # The child's own life course begins when they are born. With a parity
  # transition available to female newborns, some of them should go on to
  # give birth themselves if the horizon is long enough -- and none of
  # their events may predate their birth.
  set.seed(8L)
  pop <- fert_engine(n = 40L, horizon = c(20000101, 20601231))$run()$pop

  newborns <- unique(pop$ID[!is.na(pop$motherID)])
  expect_gt(length(newborns), 0L)
  for (id in newborns) {
    rows <- pop[pop$ID == id, , drop = FALSE]
    events <- rows$transitionTime[!is.na(rows$transitionTime)]
    if (length(events) > 0L) {
      expect_true(all(getInDays(events) >= getInDays(rows$birthDate[1L])))
    }
  }
})


test_that("no newborn appears when no fertility transition fires", {
  # The same model with the birth rate set to zero, so that the assertion
  # is exact rather than probable.
  codec <- fert_codec()
  never <- function(age, calTime) rep(0, length(age))
  registry <- new_rate_registry(
    list(list(from = "f/0", to = "f/1", rate = never)),
    codec,
    probe = micsim_rate_probe(0:120, 2000:2040, 0:1)
  )
  engine <- new_engine(codec, registry, fert_init_pop(20L),
                       horizon = c(20000101, 20101231), max_age = 120,
                       subsystems = list(fert_subsystem(codec, registry)))

  set.seed(9L)
  result <- engine$run()
  expect_identical(result$diagnostics$n_individuals, 20L)
  expect_identical(result$diagnostics$n_events, 0L)
  expect_true(all(is.na(result$pop$motherID)))
})


# ---------------------------------------------------------------------
# 5. The output contract with fertility on
# ---------------------------------------------------------------------

test_that("motherID appears only when fertility is modelled", {
  set.seed(10L)
  with_fertility <- fert_engine(n = 20L)$run()$pop
  expect_identical(names(with_fertility)[8], "motherID")

  codec <- fert_codec()
  without <- new_engine(codec, fert_registry(codec), fert_init_pop(20L),
                        horizon = c(20000101, 20101231), max_age = 120)
  set.seed(10L)
  expect_false("motherID" %in% names(without$run()$pop))
})


test_that("a newborn who never transitions still gets an NA row", {
  # They are part of the population, so the contract lists them once with
  # NA in every transition column, exactly like a survivor of the starting
  # population.
  set.seed(11L)
  pop <- fert_engine(n = 30L)$run()$pop

  quiet_newborns <- pop[!is.na(pop$motherID) & is.na(pop$From), ,
                        drop = FALSE]
  expect_gt(nrow(quiet_newborns), 0L)
  expect_true(all(is.na(quiet_newborns$To)))
  expect_true(all(is.na(quiet_newborns$transitionTime)))
  expect_false(anyNA(quiet_newborns$initState))
})


test_that("motherID is repeated across all of a newborn's rows", {
  set.seed(12L)
  pop <- fert_engine(n = 40L, horizon = c(20000101, 20601231))$run()$pop

  newborns <- unique(pop$ID[!is.na(pop$motherID)])
  expect_gt(length(newborns), 0L)
  for (id in newborns) {
    rows <- pop[pop$ID == id, , drop = FALSE]
    expect_identical(length(unique(rows$motherID)), 1L)
  }
})


test_that("the run is reproducible with fertility on", {
  engine <- fert_engine(n = 40L)
  expect_identical(engine$run(seed = 77L)$pop, engine$run(seed = 77L)$pop)
})


test_that("each run starts the subsystem's state afresh", {
  # The births counter lives in per-run state: a second run of the same
  # engine must not inherit the first run's count.
  engine <- fert_engine(n = 40L)
  a <- engine$run(seed = 5L)
  b <- engine$run(seed = 5L)
  expect_identical(a$extensions, b$extensions)
  expect_identical(a$extensions$fertility$n_births,
                   a$diagnostics$n_individuals - 40L)
})


# ---------------------------------------------------------------------
# 6. The seam itself
# ---------------------------------------------------------------------

test_that("a user subsystem runs alongside the built-in one", {
  # The point of writing fertility through the seam: a researcher can add
  # a mechanism of their own without touching the package. Here a user
  # hook simply counts events, but it could equally well add people,
  # cancel events, or maintain a relationship graph.
  log <- new.env(parent = emptyenv())
  log$events <- 0L
  log$births <- 0L

  watcher <- new_engine_hooks(
    on_event = function(ctx) log$events <- log$events + 1L,
    on_birth = function(ctx) log$births <- log$births + 1L
  )

  set.seed(13L)
  result <- fert_engine(n = 40L, hooks = watcher)$run()

  expect_identical(log$events, result$diagnostics$n_events)
  expect_identical(log$births,
                   result$diagnostics$n_individuals - 40L)
})


test_that("on_birth fires once per birth, from inside the subsystem", {
  # The fertility subsystem fires on_birth itself, which is how linked
  # lives will learn about a birth without having to recognise birth
  # transitions for itself.
  seen <- new.env(parent = emptyenv())
  seen$mothers <- integer(0)

  watcher <- new_engine_hooks(
    on_birth = function(ctx) {
      seen$mothers <- c(seen$mothers, ctx$mother_slot)
    }
  )

  set.seed(14L)
  result <- fert_engine(n = 30L, hooks = watcher)$run()

  expect_length(seen$mothers, result$diagnostics$n_individuals - 30L)
  expect_true(all(seen$mothers >= 1L))
})


test_that("combine_hooks chains observing hooks in order", {
  order_log <- new.env(parent = emptyenv())
  order_log$seen <- character(0)

  first  <- new_engine_hooks(
    on_event = function(ctx) order_log$seen <- c(order_log$seen, "first"))
  second <- new_engine_hooks(
    on_event = function(ctx) order_log$seen <- c(order_log$seen, "second"))

  fire_hook(combine_hooks(first, second), "on_event", new_hook_context())
  expect_identical(order_log$seen, c("first", "second"))
})


test_that("combine_hooks composes adjusting hooks, last one winning", {
  add_one    <- new_engine_hooks(on_schedule = function(v, ctx) v + 1)
  double_it  <- new_engine_hooks(on_schedule = function(v, ctx) v * 2)

  # (10 + 1) * 2, not 10 * 2 + 1
  expect_identical(
    adjust_by_hook(combine_hooks(add_one, double_it), "on_schedule", 10,
                   new_hook_context()),
    22
  )
})


test_that("combine_hooks ignores NULL and passes a lone set through", {
  only <- new_engine_hooks(on_event = function(ctx) NULL)
  expect_identical(combine_hooks(NULL, only, NULL), only)
  expect_identical(combine_hooks()$active, character(0))
  expect_error(combine_hooks(only, "not a hook set"),
               class = "micsim_hook_error")
})


test_that("a subsystem can bring a new individual into the simulation", {
  # ctx$schedule is the primitive that makes fertility expressible as a
  # subsystem at all. This test uses it directly, the way a user-written
  # mechanism would: on the first transition, conjure one extra person.
  codec <- fert_codec()
  conjure <- new_engine_hooks(
    on_event = function(ctx) {
      if (ctx$persons$size() > 60L) return(invisible(NULL))
      id <- ctx$persons$next_id()
      slot <- ctx$persons$add(
        id = id, birth_days = ctx$now_days,
        state_idx = ctx$codec$encode("f/0"),
        init_state_idx = ctx$codec$encode("f/0")
      )
      ctx$schedule(slot, ctx$now_days)
    }
  )

  engine <- new_engine(codec, fert_registry(codec), fert_init_pop(50L),
                       horizon = c(20000101, 20101231), max_age = 120,
                       hooks = conjure)
  set.seed(15L)
  result <- engine$run()

  expect_gt(result$diagnostics$n_individuals, 50L)
})
