## Submission

This is a major release, from 3.0.0 to 4.0.0.

* The simulation engine has been rewritten, and models are now built as
  objects: `micsim_model()`, `add_subsystem()`, `validate()` and
  `simulate()`.
* The 3.0.0 functions (`micSim()`, `micSimParallel()`,
  `buildTransitionMatrix()`, `convertToLongFormat()`,
  `convertToWideFormat()`, `getInDays_my()`) remain exported and keep
  working, but are deprecated. Each shows a message once per R session
  naming its replacement.
* Simulated results differ from 3.0.0 because of bug fixes, which are
  listed in NEWS.md.
* The package no longer depends on `snowfall` and `rlecuyer`; parallel
  runs use the base package `parallel`. R >= 4.1 is now required.
* `R.rsp` is suggested for one static PDF vignette, a diagram of how the
  package builds, checks and simulates a model.

## Test environments

* local macOS 26 (aarch64), R 4.3.2
* TODO: win-builder, R-devel and R-release
* TODO: macOS builder, R-release
* TODO: R-hub, Linux R-devel

## R CMD check results

0 errors | 0 warnings | 0 notes

The local check showed two NOTEs that come from the local setup and not
from the package:

* "checking for future file timestamps ... unable to verify current time":
  the time server could not be reached.
* "checking HTML version of manual": validation warnings from the
  outdated `tidy` shipped with macOS.

TODO: replace with the results from the remote checks.

## Tests

No test, example or vignette uses more than two cores. The migration
vignette runs one simulation on two cores.

## Reverse dependencies

There are currently no reverse dependencies on CRAN.
