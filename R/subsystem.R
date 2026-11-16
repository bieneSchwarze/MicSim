# =====================================================================
# subsystem.R
# ---------------------------------------------------------------------
# WHAT THIS FILE DOES
#   Defines what a subsystem IS, and how one is added to a model.
#
# PLACE IN THE PIPELINE
#   fertility(), school_enrolment(), ...      (a SPECIFICATION: pure data)
#        |
#        v  add_subsystem(model, spec)        (COMPILED against the model:
#        |                                     selectors -> integer edge ids,
#        |                                     everything validated)
#        v
#   engine.R, at the start of every run        (INSTANTIATED: fresh hooks
#                                                and fresh state, per run)
#        |
#        v
#   result$extensions[[name]]                  (what the subsystem reports)
#
# WHAT YOU NEED TO KNOW
#   A SUBSYSTEM is a mechanism that one person's own transition rates
#   cannot express, such as a transition that creates a person (fertility), or a
#   calendar rule that decides exactly when something happens (school
#   enrolment). The engine knows none of them. 
#   It knows only the named hook points of engine_hooks.R,
#   and a subsystem is a set of callbacks for some of those points.
#
#   The three stages are the point of this file:
#
#     specification   what the user writes. No model needed, nothing
#                     resolved, safe to build, print and reuse.
#     compilation     happens once, in add_subsystem(). Every selector is
#                     resolved against the model's state space and declared
#                     transitions. Mistakes are reported here, naming
#                     the subsystem where they happened.
#     instantiation   happens at the start of EVERY run. The hooks and any
#                     state they keep are built fresh, so nothing leaks
#                     from one run into the next, and nothing is shared
#                     between parallel workers.
#
#  WHAT A SUBSYSTEM MUST PROVIDE (new_subsystem() checks this)
#   name               unique within the model; its output appears in
#                      result$extensions[[name]]
#   compile(ctx)       ctx = list(codec, registry). Returns a list with
#     instantiate(run)   -> list(hooks = <engine_hooks>, finish = function)
#     claims             named list: kind -> integer edge ids. Two
#                        subsystems may not claim the same edge for the
#                        same kind (two calendar rules on one transition,
#                        two birth mechanisms on one birth).
#     describe           one line for print(model)
#     merge(outputs)     optional: combine the outputs of parallel workers
#                        (default merge_reports(): scalars and named
#                        counts are added, tables are stacked)
#   parallel_safe      may the population be split across cores?
#   spawns             does it add people? (then the output has motherID
#                      and the validation grid starts at age zero)
#
# WHAT A SUBSYSTEM REPORTS  (result$extensions[[name]])
#   A named list of: single numbers (counts), named integer vectors
#   (counts by category) and data frames (one row per event). Person
#   identifiers in a data frame live in columns named ID, motherID or
#   childID. Those are the columns a parallel run
#   renumbers together with the population, so that a table and `pop`
#   always agree on who is who.
#
#   Users never write a subsystem directly. The presets
#   fertility() and school_enrolment() come with the package.
# =====================================================================


#' Define a subsystem
#'
#' @description
#' The internal constructor every subsystem specification goes through,
#' including presets. It checks the declarations and returns a
#' specification. Resolving only happens once [add_subsystem()] compiles it
#' against a model.
#'
#' @param name Unique name; the key of `result$extensions`.
#' @param type The kind of subsystem, e.g. `"fertility"`.
#' @param compile Function of one argument, `ctx` (a list with `codec` and
#'   `registry`), returning the compiled form described in the file header.
#' @param params Named list of the user's parameters, kept for printing.
#' @param parallel_safe `TRUE` if splitting the population across cores
#'   leaves the subsystem's behaviour unchanged.
#' @param spawns `TRUE` if the subsystem adds individuals to the run.
#' @return An object of class `micsim_subsystem`.
#' @keywords internal
#' @noRd
new_subsystem <- function(name, type, compile, params = list(),
                          parallel_safe = TRUE, spawns = FALSE) {
  if (!is.character(name) || length(name) != 1L || is.na(name) ||
        !nzchar(name)) {
    abort_subsystem("A subsystem's `name` must be a single, non-empty string.")
  }
  if (!is.function(compile)) {
    abort_subsystem(paste0("Subsystem `", name, "`: `compile` must be a function."))
  }
  flags <- list(parallel_safe = parallel_safe, spawns = spawns)
  for (flag in names(flags)) {
    if (!isTRUE(flags[[flag]]) && !isFALSE(flags[[flag]])) {
      abort_subsystem(paste0("Subsystem `", name, "`: `", flag,
                             "` must be TRUE or FALSE."))
    }
  }
  structure(
    c(list(name = name, type = type, compile = compile, params = params),
      flags),
    class = "micsim_subsystem"
  )
}


#' Add a subsystem to a model
#'
#' @description
#' Switches on a mechanism that ordinary transitions cannot express, such
#' as transitions that create a newborn ([fertility()]) or a calendar rule
#' for when a transition takes place ([school_enrolment()]).
#'
#' @details
#' Subsystems are added one at a time, and the call returns the model.
#' Thus, a workflow can look like this:
#'
#' ```
#' model <- micsim_model(...) |>
#'   add_subsystem(fertility(births = moves_of("parity"), newborn = nb)) |>
#'   add_subsystem(school_enrolment(move("edu", "none", "primary"), month = 9))
#' ```
#'
#' Adding a subsystem compiles it against the model. That means every selector 
#' is resolved into the model's declared transitions and every parameter is
#' checked. If there is a mistake, it is reported and names the subsystem. Two
#' subsystems with the same name, or two that claim the same transition
#' for the same purpose, will give an error.
#'
#' Subsystems take effect in the order they were added. Their per-run
#' state is created fresh at the start of every [simulate()] call, and
#' whatever they report is returned in `result$extensions`, under their
#' name.
#'
#' @param model A [micsim_model()].
#' @param subsystem A subsystem, from either
#'  [fertility()] or [school_enrolment()].
#'
#' @return The model with the subsystem installed.
#'
#' @seealso [fertility()], [school_enrolment()], [move()], [MicSim-concepts]
#'
#' @examples
#' mortality <- function(age) 0.00003 * exp(0.1 * age)
#' first_birth <- function(age) ifelse(age > 15 & age < 45, 0.1, 0)
#'
#' model <- micsim_model(
#'   state_space = state_space(sex = c("male", "female"), parity = c("0", "1")),
#'   transitions = transitions(
#'     transition(move("parity", "0", "1", where = sex == "female"), first_birth),
#'     transition(NA, "deceased", mortality)
#'   ),
#'   population = data.frame(ID = 1:20, birthDate = "19850101",
#'                           initState = "female/0"),
#'   horizon = c(20140101, 20241231)
#' ) |>
#'   add_subsystem(fertility(
#'     births  = moves_of("parity"),
#'     newborn = newborn_rule(sex = c(male = 0.515, female = 0.485),
#'                            parity = "0")
#'   ))
#' model
#'
#' @export
add_subsystem <- function(model, subsystem) {
  if (!inherits(model, "micsim_model")) {
    abort_subsystem("`model` must be built with micsim_model().")
  }
  if (!inherits(subsystem, "micsim_subsystem")) {
    abort_subsystem(paste0(
      "`subsystem` must be built by a subsystem constructor such as\n",
      "  fertility() or school_enrolment(), not a ", class(subsystem)[1L], "."
    ))
  }
  if (subsystem$name %in% names(model$subsystems)) {
    abort_subsystem(paste0(
      "The model already has a subsystem named `", subsystem$name, "`.\n",
      "  Give the second one its own name, e.g. ", subsystem$type,
      "(..., name = \"", subsystem$name, "_2\")."
    ))
  }

  compiled <- compile_subsystem(subsystem, model$codec, model$registry)
  check_claim_conflicts(compiled, model$subsystems)

  model$subsystems[[subsystem$name]] <- compiled

  # As in micsim_model(): a model that exists is a model that can run.
  invisible(micsim_engine_for(model, verbose = FALSE))
  model
}


#' Compile a subsystem specification against a model
#'
#' @description
#' Runs the specification's `compile()` and checks that what comes back
#' satisfies what a subsystem has to provide (in the header).
#' Used by [add_subsystem()].
#'
#' @param subsystem A `micsim_subsystem`.
#' @param codec A `state_codec`.
#' @param registry A `rate_registry`.
#' @return A `micsim_compiled_subsystem`: the specification's declarations
#'   plus `instantiate`, `claims`, `describe` and `merge`.
#' @keywords internal
#' @noRd
compile_subsystem <- function(subsystem, codec, registry) {
  out <- subsystem$compile(list(codec = codec, registry = registry))
  where <- paste0("Subsystem `", subsystem$name, "`")
  if (!is.list(out) || !is.function(out$instantiate)) {
    abort_subsystem(paste0(where, ": compile() must return a list with an\n",
                           "  `instantiate` function."))
  }
  claims <- if (is.null(out$claims)) list() else out$claims
  if (!is.list(claims) || (length(claims) > 0L && is.null(names(claims)))) {
    abort_subsystem(paste0(where, ": `claims` must be a named list of edge ids."))
  }
  structure(
    list(
      name              = subsystem$name,
      type              = subsystem$type,
      params            = subsystem$params,
      parallel_safe     = subsystem$parallel_safe,
      spawns            = subsystem$spawns,
      instantiate       = out$instantiate,
      claims            = claims,
      describe          = if (is.null(out$describe)) subsystem$type else out$describe,
      merge             = if (is.function(out$merge)) out$merge else merge_reports
    ),
    class = "micsim_compiled_subsystem"
  )
}


#' Refuse two subsystems that claim the same transition for one purpose
#'
#' @param new A compiled subsystem about to be added.
#' @param existing The model's compiled subsystems.
#' @return `invisible(TRUE)`, or raises `micsim_subsystem_error`.
#' @keywords internal
#' @noRd
check_claim_conflicts <- function(new, existing) {
  for (old in existing) {
    for (kind in intersect(names(new$claims), names(old$claims))) {
      shared <- intersect(new$claims[[kind]], old$claims[[kind]])
      if (length(shared) > 0L) {
        abort_subsystem(paste0(
          "Subsystems `", old$name, "` and `", new$name, "` both claim ",
          length(shared), " transition(s)\n",
          "  for the same purpose (", kind, "). Each transition can have one ",
          kind, " rule."
        ))
      }
    }
  }
  invisible(TRUE)
}


#' Build the hooks of every subsystem for one run
#'
#' @description
#' Called by the engine at the start of every run, so that each run gets
#' fresh hooks and a fresh subsystem state.
#'
#' @param subsystems List of compiled subsystems in the order they were
#'   added.
#' @param run A list describing the run: `codec`, `registry`.
#' @return A list of instances, each a list with `hooks` (an
#'   `engine_hooks` object) and `finish` (a function, or `NULL`).
#' @keywords internal
#' @noRd
instantiate_subsystems <- function(subsystems, run) {
  out <- lapply(subsystems, function(s) {
    inst <- s$instantiate(run)
    if (!is.list(inst) || !inherits(inst$hooks, "engine_hooks")) {
      abort_subsystem(paste0(
        "Subsystem `", s$name, "`: instantiate() must return a list whose\n",
        "  `hooks` element is an engine_hooks object."
      ))
    }
    inst
  })
  names(out) <- vapply(subsystems, `[[`, character(1L), "name")
  out
}


#' Collect what each subsystem reports at the end of a run
#' @param instances The instances from `instantiate_subsystems()`.
#' @return A named list; subsystems that report nothing are left out.
#' @keywords internal
#' @noRd
finish_subsystems <- function(instances) {
  out <- lapply(instances, function(inst) {
    if (is.function(inst$finish)) inst$finish() else NULL
  })
  out[!vapply(out, is.null, logical(1L))]
}


#' Combine the outputs of parallel workers, subsystem by subsystem
#' @param subsystems Compiled subsystems.
#' @param per_worker List, one element per worker, of `extensions` lists.
#' @return A named list.
#' @keywords internal
#' @noRd
merge_extensions <- function(subsystems, per_worker) {
  out <- list()
  for (s in subsystems) {
    pieces <- lapply(per_worker, function(ext) ext[[s$name]])
    pieces <- pieces[!vapply(pieces, is.null, logical(1L))]
    if (length(pieces) > 0L) out[[s$name]] <- s$merge(pieces)
  }
  out
}


#' Merge the reports of parallel workers
#'
#' @description
#' The default `merge` of a subsystem, which combines the workers' reports
#' element by element, following the conventions in the file header.
#'
#' @details
#' * a single number is summed;
#' * a named numeric vector (counts by category) is summed by name, so a
#'   category seen by one worker only is kept;
#' * a data frame is stacked and ordered by its first column if that is a
#'   person identifier.
#'
#' Anything else is taken from the first worker that has it.
#'
#' @param pieces List of reports, one per worker.
#' @return One report.
#' @keywords internal
#' @noRd
merge_reports <- function(pieces) {
  keys <- unique(unlist(lapply(pieces, names)))
  out <- lapply(keys, function(k) {
    parts <- lapply(pieces, `[[`, k)
    parts <- parts[!vapply(parts, is.null, logical(1L))]
    first <- parts[[1L]]
    if (is.data.frame(first)) {
      df <- do.call(rbind, parts)
      if (ncol(df) > 0L && names(df)[1L] %in% ID_COLUMNS) {
        df <- df[order(df[[1L]]), , drop = FALSE]
      }
      row.names(df) <- NULL
      df
    } else if (is.numeric(first) && is.null(names(first)) &&
                 length(first) == 1L) {
      total <- sum(vapply(parts, as.numeric, numeric(1L)))
      if (is.integer(first)) as.integer(total) else total
    } else if (is.numeric(first) && !is.null(names(first))) {
      nm <- unique(unlist(lapply(parts, names)))
      total <- vapply(nm, function(n) {
        sum(vapply(parts, function(p) if (n %in% names(p)) as.numeric(p[[n]]) else 0,
                   numeric(1L)))
      }, numeric(1L))
      if (all(vapply(parts, is.integer, logical(1L)))) {
        total <- stats::setNames(as.integer(total), nm)
      }
      total
    } else {
      first
    }
  })
  stats::setNames(out, keys)
}


#' Columns of a report table that hold person identifiers
#' @keywords internal
#' @noRd
ID_COLUMNS <- c("ID", "motherID", "childID")


#' Renumber the identifiers in one worker's reports
#'
#' @description
#' Applies the renumbering that `recombine_populations()` applied to the
#' worker's population to every identifier column of every table in the
#' worker's `extensions`.
#'
#' @param extensions One worker's `extensions` list.
#' @param id_map `list(old =, new =)`, or `NULL` when nothing was
#'   renumbered.
#' @return The `extensions` list with identifiers replaced.
#' @keywords internal
#' @noRd
renumber_extensions <- function(extensions, id_map) {
  if (is.null(id_map) || length(id_map$old) == 0L) return(extensions)
  lapply(extensions, function(report) {
    if (!is.list(report)) return(report)
    for (k in seq_along(report)) {
      tbl <- report[[k]]
      if (!is.data.frame(tbl)) next
      for (col in intersect(names(tbl), ID_COLUMNS)) {
        hit <- match(tbl[[col]], id_map$old)
        tbl[[col]][!is.na(hit)] <- id_map$new[hit[!is.na(hit)]]
      }
      report[[k]] <- tbl
    }
    report
  })
}


#' A growable table for one run's events
#'
#' @description
#' Columns are kept as preallocated vectors in an environment and doubled
#' when full, same as in the transition store. 
#' Used by subsystems to keep one row per event they handle.
#'
#' @param columns Named list of zero-length prototypes, e.g.
#'   `list(ID = integer(), time = numeric())`.
#' @return A list with `add(...)` (one value per column, by name) and
#'   `columns()` (the filled part, as a list of vectors).
#' @keywords internal
#' @noRd
new_event_log <- function(columns) {
  e <- new.env(parent = emptyenv())
  e$n <- 0L
  e$cap <- 64L
  e$cols <- lapply(columns, function(proto) {
    v <- proto; length(v) <- e$cap; v
  })
  add <- function(...) {
    row <- list(...)
    if (e$n == e$cap) {
      e$cap <- 2L * e$cap
      cols <- e$cols; e$cols <- NULL
      for (k in names(cols)) length(cols[[k]]) <- e$cap
      e$cols <- cols
    }
    e$n <- e$n + 1L
    cols <- e$cols; e$cols <- NULL
    for (k in names(row)) cols[[k]][e$n] <- row[[k]]
    e$cols <- cols
    invisible(NULL)
  }
  list(
    add     = add,
    size    = function() e$n,
    columns = function() lapply(e$cols, function(v) v[seq_len(e$n)])
  )
}


#' Print a subsystem specification
#' @param x A `micsim_subsystem`.
#' @param ... Not used.
#' @return `x`, invisibly.
#' @exportS3Method print micsim_subsystem
#' @noRd
print.micsim_subsystem <- function(x, ...) {
  cat("<micsim_subsystem> ", x$type, " (name: ", x$name, ")\n", sep = "")
  for (p in names(x$params)) {
    cat("  ", p, ": ", describe_param(x$params[[p]]), "\n", sep = "")
  }
  cat("  add it to a model with add_subsystem()\n")
  invisible(x)
}


#' Describe one subsystem parameter on one line
#' @param x A parameter value.
#' @return A string.
#' @keywords internal
#' @noRd
describe_param <- function(x) {
  if (inherits(x, "micsim_selector")) return(describe_selector(x))
  if (inherits(x, "micsim_role")) return(deparse1(x$expr))
  if (inherits(x, "micsim_newborn_spec") || inherits(x, "newborn_rule")) {
    return("a newborn rule")
  }
  if (is.list(x)) {
    return(paste(vapply(x, describe_param, character(1L)), collapse = ", "))
  }
  if (is.atomic(x)) return(toString(x))
  class(x)[1L]
}


#' Raise a subsystem error
#' @param problem The message.
#' @return Never returns.
#' @keywords internal
#' @noRd
abort_subsystem <- function(problem) {
  micsim_abort("micsim_subsystem_error", problem)
}
