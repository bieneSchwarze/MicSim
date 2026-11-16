# MicSim

MicSim is an R package for **continuous-time microsimulation**. You can use it
to simulate how individual life courses unfold under different demographic,
social or economic conditions.

A life course consists of the states a person occupies and the times at which
they move between them. In continuous time, these transitions can occur at any
point during the simulation period. You specify which transitions are possible
and supply a rate, also called a hazard, for each one. MicSim uses those
rates to draw when events occur and which state a person enters next.

- **State space.** Choose the characteristics that matter for your research
  question, such as health status, employment or country of residence. A state
  combines one value for each characteristic.
- **Transition rates.** Use constant rates or let them vary with age,
  calendar time or both. For transitions that change a single state variable,
  rates can also depend on duration, the time since that variable last changed.
  Rates can differ between states.
- **Population processes.** Model births, deaths, immigration and emigration.
  You can also specify how newborns acquire their starting characteristics and
  schedule enrolment according to a school calendar.
- **Event histories.** Follow individual transitions and use the simulated
  histories to calculate population measures such as age-specific prevalence.
- **Large populations.** Simulate life courses for large synthetic populations,
  with the option to use multiple processor cores.
  
MicSim builds on the methods developed in Sabine Zinn's doctoral dissertation
at the University of Rostock in 2011.
[A Continuous-Time Microsimulation and First Steps Towards a Multi-Level Approach in Demography](https://doi.org/10.18453/rosdok_id00000951)

## Installation

MicSim 4.0.0 requires R 4.1 or later. Install MicSim from CRAN and load it in R.

```r
install.packages("MicSim")
library(MicSim)
```

## Getting started

Use `micsim_model()` to bring together the starting population, simulation
period, state space and transitions. You can check the rates with `validate()`
before running the model with `simulate()`. The result contains the model and
the simulated life courses. Use `as_long()` for one row per spell in a state
or `as_wide()` for one row per person.

- Start here for a general introduction to MicSim 4.0.0's model-object design
  and new functions. The vignette shows how to build a model, add subsystems,
  run the simulation and check its results.
  [Introduction to MicSim 4.0.0](vignettes/introduction-to-new-micsim.Rmd)
- Model moves between the Netherlands, Spain and Sweden alongside births, deaths,
  immigration and emigration.
  [Simulating migration flows](vignettes/simulating-migration-flows.Rmd)
- Use rates estimated from data to simulate major depressive disorder and
  compare prevalence across education groups under two scenarios.
  [Educational inequalities in depression](vignettes/educational-inequalities-in-depression.Rmd)

Changes in version 4.0.0
are described in the [release notes](NEWS.md).

## Citations

Run `citation("MicSim")` to get the citation for your installed version.

Zinn, S. (2014). The MicSim package of R: An entry-level toolkit for
continuous-time microsimulation. *International Journal of Microsimulation*,
7(3), 3–32. [Article](https://doi.org/10.34196/IJM.00105)

## Licence

MicSim is licensed under the GNU General Public License, version 2 or later
(`GPL (>= 2)`). [CRAN package page](https://CRAN.R-project.org/package=MicSim)