# =====================================================================
# test-release-prep.R
# ---------------------------------------------------------------------
# Version stamp of models, summary() by dimension, the clock guarantee
# for on_schedule, and the reestimation plot's options and robustness.
# =====================================================================

rp_mort <- function(age) 0.0004 * exp(0.08 * age)

rp_model <- function(n = 300L, seed = 7L) {
  set.seed(seed)
  micsim_model(
    state_space(sex = c("m", "f"), edu = c("no", "low"), par = c("0", "1", "2")),
    transitions(
      transition(move("edu", "no", "low"), function(age) ifelse(age == 7, Inf, 0)),
      transition(move("par", "0", "1", where = sex == "f"),
                 function(age) ifelse(age > 15 & age < 45, 0.1, 0)),
      transition(move("par", "1", "2", where = sex == "f"),
                 function(age) ifelse(age > 15 & age < 45, 0.08, 0)),
      transition(NA, "dead", rp_mort)
    ),
    data.frame(ID = seq_len(n),
               birthDate = sprintf("%d0315", sample(1950:1999, n, TRUE)),
               initState = rep(c("f/low/0", "m/low/0"), length.out = n)),
    c(20000101, 20201231)
  ) |>
    add_subsystem(fertility(moves_of("par"),
                            newborn_rule(sex = c(m = 0.5, f = 0.5), edu = "no", par = "0"))) |>
    add_subsystem(school_enrolment(move("edu", "no", "low"), month = 9))
}

# ---------------------------------------------------------------------
# Version stamp
# ---------------------------------------------------------------------

test_that("a model records the MicSim version that built it", {
  m <- rp_model(n = 20L)
  expect_identical(m$micsim_version,
                   as.character(getNamespaceVersion("MicSim")))
  expect_silent(simulate(m, seed = 1, reestimate = FALSE))
})

test_that("simulating a model built by another version warns", {
  m <- rp_model(n = 20L)
  m$micsim_version <- "3.9.0"
  expect_warning(simulate(m, seed = 1, reestimate = FALSE),
                 "MicSim 3.9.0", class = "micsim_version_warning")
  m$micsim_version <- NULL
  expect_warning(simulate(m, seed = 1, reestimate = FALSE),
                 class = "micsim_version_warning")
})

test_that("bind_hazard() no longer takes genparams", {
  expect_false("genparams" %in% names(formals(bind_hazard)))
  m <- rp_model(n = 5L)
  expect_false(any(vapply(m$registry$edges,
                          function(e) "genparams" %in% names(e), logical(1L))))
})

test_that("subsystems no longer declare moves_events_back", {
  expect_false("moves_events_back" %in% names(formals(new_subsystem)))
  m <- rp_model(n = 5L)
  for (s in m$subsystems) expect_null(s$moves_events_back)
})

test_that("an on_schedule hook may not move an event into the past", {
  codec <- new_state_codec(list(sex = c("m", "f")), absorbing = "dead")
  reg <- new_rate_registry(list(list(from = NA, to = "dead", rate = rp_mort)),
                           codec)
  back <- function(event, ctx) { event$time_days <- ctx$now_days - 10; event }
  eng <- new_engine(codec, reg,
                    data.frame(ID = 1:3, birthDate = "19500101",
                               initState = "f"),
                    horizon = c(20000101, 20101231), max_age = 100,
                    hooks = new_engine_hooks(on_schedule = back))
  expect_error(eng$run(seed = 1), "into the past",
               class = "micsim_engine_error")
})

# ---------------------------------------------------------------------
# summary() by dimension of the state space
# ---------------------------------------------------------------------

test_that("summary() gives one age summary per dimension of the state space", {
  r <- simulate(rp_model(), seed = 3, reestimate = FALSE)
  s <- summary(r)
  ages <- s$ages
  expect_true(is.matrix(ages))
  expect_identical(colnames(ages),
                   c("n", "Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max."))
  # sex never changes: no row. edu and par change; dead is an exit.
  expect_setequal(rownames(ages), c("edu", "par", "exit: dead"))

  pop <- r$pop
  moved <- pop[!is.na(pop$From), ]
  par_of <- function(x) vapply(strsplit(x, "/", fixed = TRUE),
                                function(p) p[3], character(1L))
  alive <- moved[moved$To != "dead", ]
  births <- alive[par_of(alive$From) != par_of(alive$To), ]
  # All births are pooled, first and second alike
  expect_identical(ages["par", "n"], as.numeric(nrow(births)))
  expect_equal(ages["par", "Median"], stats::median(births$transitionAge))
  expect_identical(ages["exit: dead", "n"], as.numeric(sum(moved$To == "dead")))
  expect_output(print(s), "by dimension of the state space")
})

test_that("summary() of a run without events has no age table", {
  m <- micsim_model(state_space(sex = c("m", "f")),
                    transitions(transition(NA, "dead", function(age) rep(0, length(age)))),
                    data.frame(ID = 1:3, birthDate = "19700101", initState = "f"),
                    c(20000101, 20011231))
  s <- summary(simulate(m, seed = 1, reestimate = FALSE))
  expect_null(s$ages)
  expect_output(print(s), "State at the end")
})

test_that("a transition changing two variables counts in both dimensions", {
  m <- micsim_model(state_space(a = c("a0", "a1"), b = c("b0", "b1")),
                    transitions(transition("a0/b0", "a1/b1", function() 0.5)),
                    data.frame(ID = 1:50, birthDate = "19700101", initState = "a0/b0"),
                    c(20000101, 20201231))
  s <- summary(simulate(m, seed = 1, reestimate = FALSE))
  expect_setequal(rownames(s$ages), c("a", "b"))
  expect_identical(s$ages["a", "n"], s$ages["b", "n"])
})

# ---------------------------------------------------------------------
# The reestimation plot
# ---------------------------------------------------------------------

rp_result <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) cache <<- simulate(rp_model(), seed = 5)
    cache
  }
})

on_png <- function(width, height, code) {
  f <- tempfile(fileext = ".png")
  grDevices::png(f, width = width, height = height)
  on.exit({ grDevices::dev.off(); unlink(f) }, add = TRUE)
  grDevices::dev.control("enable")
  force(code)
}

test_that("a device too small for any layout is refused before drawing", {
  r <- rp_result()
  on_png(200, 150, {
    expect_error(plot(r, which = 1), class = "micsim_plot_error")
    expect_false(graphics::par("new"))
    expect_identical(graphics::par("mfrow"), c(1L, 1L))
    expect_silent(plot(1:3))          # the device is still usable
  })
})

test_that("a small device gets the compact layout, with one message", {
  r <- rp_result()
  on_png(320, 220, {
    expect_message(plot(r, which = 1:2, ask = FALSE),
                   "2 of 2 figures", class = "micsim_plot_message")
    expect_false(graphics::par("new"))
  })
})

test_that("a figure replays on a smaller device", {
  r <- rp_result()
  on_png(1000, 450, {
    plot(r, which = 2)
    p <- grDevices::recordPlot()
  })
  on_png(400, 250, expect_error(grDevices::replayPlot(p), NA))
})

test_that("plot() restores the caller's graphics settings", {
  r <- rp_result()
  on_png(1000, 450, {
    graphics::par(mar = c(1, 2, 3, 4), las = 0)
    plot(r, which = 1)
    expect_identical(graphics::par("mar"), c(1, 2, 3, 4))
    expect_identical(graphics::par("las"), 0L)
    expect_identical(graphics::par("mfrow"), c(1L, 1L))
  })
})

test_that("titles, labels and limits are taken from `...`", {
  r <- rp_result()
  on_png(1000, 450, {
    expect_silent(plot(r, which = 1:2, ask = FALSE,
                       main = c("A", "B"), xlab = "Age", ylim = c(0, 0.3),
                       xlim = c(10, 50), cex = 0.9))
  })
  on_png(1000, 450, {
    expect_warning(plot(r, which = 1, col = "red"), "col",
                   class = "micsim_plot_warning")
  })
  on_png(1000, 450, {
    expect_error(plot(r, which = 1, xlim = 5), class = "micsim_plot_error")
    expect_error(plot(r, which = 1, "unnamed"), class = "micsim_plot_error")
  })
})

test_that("`style` overrides colours and is checked", {
  r <- rp_result()
  on_png(1000, 450, {
    expect_silent(plot(r, which = 1, style = list(smooth = "darkred")))
    expect_error(plot(r, which = 1, style = list(foo = "red")),
                 class = "micsim_plot_error")
    expect_error(plot(r, which = 1, style = list(smooth = "not-a-colour")),
                 class = "micsim_plot_error")
  })
  expect_identical(plot_style(), PLOT_STYLE)   # reset afterwards
})

test_that("single panels and layout = FALSE draw into the caller's grid", {
  r <- rp_result()
  on_png(1000, 800, {
    graphics::par(mfrow = c(2, 2))
    expect_silent(plot(r, which = 1:4, panels = "recovered", layout = FALSE))
    expect_identical(graphics::par("mfrow"), c(2L, 2L))
  })
  on_png(1000, 450, {
    expect_silent(plot(r, which = 1, panels = "supplied", legend = FALSE))
  })
  expect_error(plot(r, panels = "left"), class = "micsim_plot_error")
})

test_that("as.data.frame() returns every layer of the figure", {
  x <- rp_result()$reestimation
  bands <- as.data.frame(x)
  expect_identical(bands, as.data.frame(x, what = "bands"))
  sm <- as.data.frame(x, what = "smooth")
  tr <- as.data.frame(x, what = "truth")
  expect_true(all(c("transition", "axis", "group", "x", "rate") %in% names(sm)))
  expect_true(all(c("transition", "axis", "group", "x", "rate") %in% names(tr)))
  expect_error(as.data.frame(x, what = "nothing"))
})
