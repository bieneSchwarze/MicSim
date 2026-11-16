# MicSim release notes

### Changes in version 4.0.0

MicSim 4.0.0 replaces the simulation engine and adds an object-based interface for building, checking and running models. Scripts written for 3.0.0 still run through the classic functions, which are now deprecated. Simulated results differ from 3.0.0, also with the same seed due to numerous bug fixes (see "Changes to simulated results"). New vignette "Introduction to MicSim 4.0.0" guides the user through the 4.0.0 interface.

#### New model interface

- Models are now built as objects. `state_space()` declares the state variables, `transition()` and `transitions()` declare the transitions with their rate functions (supplied as functions, not by name), and `move()` and `moves_of()` select transitions that change one state variable. `micsim_model()` combines them with the population, immigrants, simulation horizon, maximum age and MicSim version. `update()` changes parts of a model, `validate()` checks it and `simulate()` runs it. The workflow supports the native pipe (`|>`).
- Optional mechanisms are added with `add_subsystem()`: `fertility()` creates newborns, whose states are set with `newborn_rule()` and `from_mother()`, and `school_enrolment()` moves school entry to a fixed calendar month. The subsystem approach will be the backbone of any future extensions of MicSim.
- `simulate()` returns a `micsim_result` with the simulated population, the model, the seed and run diagnostics. It has `print()`, `summary()`, `as.data.frame()`, `as_long()` and `as_wide()` methods. `reestimate()` estimates the transition rates back from the result and compares them with the supplied rates and `plot()` shows the comparison. `simulate()` does this by default (`reestimate = FALSE` turns it off).
- A new help page, `?MicSim-concepts`, explains the terms used throughout the documentation.

#### Deprecated

- The classic interface is deprecated: `micSim()`, `micSimParallel()`, `buildTransitionMatrix()`, `convertToLongFormat()`, `convertToWideFormat()` and `getInDays_my()`. They remain exported and run on the new engine, but show a message once per R session naming the replacement (`micsim_model()`, `add_subsystem()` and `simulate()`; `as_long()` and `as_wide()`; `getFirstDayOfMonth()`). Set `options(MicSim.show_deprecation = FALSE)` to turn the messages off. They will be removed in a future release.

#### Changes to simulated results

- Calendar dates are exact. They follow the Gregorian calendar, using Howard Hinnant's date algorithms [1], so transition dates and ages can differ by a few days from 3.0.0. Invalid dates give `NA` and are rejected as simulation inputs. Rate functions still receive age, calendar time and duration in years of 365.25 days.
- Every individual draws from their own random-number stream (L'Ecuyer-CMRG), derived from the seed and their ID. A newborn's stream is derived from the mother's. A seed therefore gives the same result on any number of cores, and scenarios simulated with the same seed share their random numbers.
- Without a seed, `simulate()` and `micSim()` draw one from R's current random-number state, so `set.seed()` before the call still makes a run reproducible; `simulate()` records the seed in its result. Apart from that one draw, the caller's random-number state, including the generator kind, is left unchanged.
- Forced transitions (infinite rates) are detected at every whole-number age, calendar year and duration. 3.0.0 found them only at some ages, so, for example, children born during the simulation were never enrolled in school.
- School enrolment is scheduled before competing risks are resolved. 3.0.0 moved the event to the enrolment date afterwards, which could block other transitions or move the simulation clock backwards.
- In the long format, episodes are censored at the maximum age as well as at the end of the simulation horizon.

#### Bug fixes

- States of a state variable with ten or more values were coded wrongly, so different states could be confused.
- Parallel runs seeded only one random-number stream, whatever the number of cores.
- When the numerical integration of a rate failed, it was silently repeated at lower accuracy. Such cases are now counted in the run diagnostics, and an infinite rate met during integration is an error naming the transition.
- `convertToLongFormat()` no longer reads `stateSpace`, `simHorizon` and other objects from the global environment. `micSim()` stores what the conversion needs with its result. Pass `stateSpace` to give the state-variable columns their names.
- Absorbing states are taken from the transition matrix, so an `absStates` argument in a different order no longer mislabels them. `absStates` is ignored.
- Newborn-state probabilities must sum to one within every group of inherited attributes. 3.0.0 checked only the first group.
- In the long format, the last episode keeps the state entered by the last transition, also when that transition changed two state variables.
- `micSimParallel()` splits the population into exactly `cores` parts, and the identifiers of newborns no longer collide with supplied identifiers.
- Faulty inputs are reported with clear errors before the simulation starts: identifiers that are not whole numbers or are duplicated (also between the starting and the immigrant population), absorbing initial states, rate functions that are not vectorised, `fertTr` patterns that match no transition, and duration-dependent rates on transitions that change several state variables at once.
- Various other bug fixes.

#### Other changes

- Simulation is much faster for large populations: transitions are recorded without copying the history, pending events are kept in a priority queue, and rate functions are resolved once when the model is built.
- Named arguments of the date helpers have changed: `currDate` becomes `yyyymmdd` in `getInDays()`, and `daysSince01011970` becomes `days` in `getYear()`, `getMonth()`, `getDay()` and `getInDateFormat()`. Calls that pass the arguments by position are unaffected.
- `micsim_model()` uses `max_age = 100` by default `maxAge` in the classic functions keeps its default of 99.
- Models record the MicSim version that built them, and `simulate()` warns when it differs from the installed version.
- Parallel execution uses R's `parallel` package. `snowfall` and `rlecuyer` are no longer needed. R >= 4.1 is required.
- The vignettes have been rewritten for the new interface: an introduction, an application to educational inequalities in depression, and the simulation of migration flows.
 
 [1] Howard Hinnant, *chrono-Compatible Low-Level Date Algorithms*, page dated 2021-09-01: <https://howardhinnant.github.io/date_algorithms.html>. MicSim adapts `days_from_civil` and `civil_from_days`.

### Changes in version 2.0.0

- The simulation engine has been speed up by making its source code more efficient.
- The dependency to the chron package of R has been removed to get rid of unhandy format definitions of dates in the input data to the simulation. Instead now a new date format has been defined that is based on numerical values. It is `yyyymmdd`, e.g. `20221017` for October 17th, 2022. New functions have been designed and included into the package facilitating the handling of dates within the simulation.
- There is now an initial check in the package to see whether the input data is in the right format and the rates required to run the simulation are available for the complete simulation horizon and along all ages considered.
- A vignette has been added to show the capabilities of the package for running large-scale simulation applications. The example of the vignette has been designed by researchers of the Eureopean Commission. It is on the modelling and simulation of migration flows within Europe. They are now added as contributors to the MicSim package. They have also tested and checked the code of the package to remove bugs.
- The example section of the `micSim` function has been extended by another more complex example to highlight the versality of the package.
- It is possible with this new version to define which parts of the initial population and/or migrant population go to which CPUs/cores when conducting parallel computing.

### Changes in version 1.1.0

- Fix two bugs related to the creation of newborns during simulation: (1) Before that version all transitions of the type "1+"->"1+" led to newborns; now only transitions where no other attributes of the state change cause newborns, e.g. "f/M/+1" -> "f/D/1+" caused a newborn before this version. This is not any longer the case. Now only transitions of the form "f/M/1+" -> "f/M/+1" cause a newborn. (2) There was a bug in transforming the fertility matrix to detect birth events; instead of rows, columns were taken as transitions. This has been fixed now.
- Fix an inconsistency in the computation of transition ages (in years) related to leap years (now the related function is vectorized).
- A further functionality was added to the package: now it is possible to link individuals born during simulation to their mothers (see argument 'reportMothers' in the `micSim` or `micSimParallel` function).

### Changes in version 1.0.15

- Fix bug in parallel computing function concerning case when there are migrants in the model (this bug resulted from fix to version 1.0.14).

### Changes in version 1.0.14

- Fix bug in parallel computing function concerning case when there are no migrants in the model.

### Changes in version 1.0.13

- Update of chron package caused an error due to new date format. This problem has been fixed.

### Changes in version 1.0.12

- Fix bug concerning the transition age of immigrants.

### Changes in version 1.0.11

- Up to version 1.0.10 deterministic events could only be defined to occur after whole numbered years. Now, they can also be defined to happen after whole numbered days.

### Changes in version 1.0.10

- Fix bug concerning the uniqueness of IDs of newborns when simulating natives and migrants on several cores, i.e., when using the function `micSimParallel`.
- Slight changes in the description of the (life-course) model used in MicSim.

### Changes in version 1.0.9

- Add error message arising if negative rate values appear.
- Besides not defining impossible transitions in the transition matrix, they can now also be specified by giving only zero rate values for the respective rate function/s.

### Changes in version 1.0.8

- Missing package dependency has been added. MicSim requires `rlecuyer`.

### Changes in version 1.0.7

- Add a new reference.
- Modify function `micSim` that it can also deal with a single life course.

### Changes in version 1.0.6

- Fix bug in function `convertToLongFormat` leading to wrong ordering of the event sequences of individuals born during simulation.

### Changes in version 1.0.5

- Fix minor conversion problem in function `convertToWideFormat`, related to request for NA values.
- Now function `buildTransitionMatrix` can handle covariate dependent mortality rates of every kind, e.g., mortality rates depending on health status and gender.

### Changes in version 1.0.4

- Add function `convertToWideFormat` to convert simulation output into wide format.
- Slight changes in the examples to make them more feasible. Concretely, in several rates definitions I adjusted the year in which dynamics are assumed to change.
- Fix a bug in the optional features of the `micSim` function: Children who are part of the base population had been enrolled one year to late. Now, enrollment dates are OK.

### Changes in version 1.0.3

- Refine funtionality of `micSim` function to handle deterministc events (marked by Inf rates).
- Fix bug in function `buildTransitionMatrix`.

### Changes in version 1.0.2

- Add News file.

### Changes in version 1.0.1

- Extend function `buildTransitionMatrix` to facilitate using transition rates which depend on covariates such as gender and marital status.
- Correct typing error in Details section of help page of function `buildTransitionMatrix`.
- Fix a problem in function `micSim` occuring if the simulation is asked not to produce newborns (concerns argument `fertTr=c()`).
