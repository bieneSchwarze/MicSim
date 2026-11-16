# =====================================================================
# concepts.R
# ---------------------------------------------------------------------
# WHAT THIS FILE DOES
#   A glossary. It holds no code: it gives `?MicSim-concepts` a place to
#   look up every technical term the rest of the documentation uses, so
#   the other help pages can link here instead of re-explaining.
#
# WHO IT IS FOR
#   Researchers who may be new to microsimulation or the way MicSim
#   processes events.
# =====================================================================

#' Concepts and notation used throughout MicSim
#'
#' @description
#' A glossary of the terms and symbols that appear in MicSim's
#' documentation and source code.
#'
#' @section State space:
#' \describe{
#'   \item{**state variable** (substate)}{A categorical attribute used to
#'     describe an individual in the model, with a set of possible values.
#'     For example, an education variable might take the values `"none"`,
#'     `"primary"` and `"secondary"`. Another variable might record the
#'     number of children ever born, also called parity. A state variable
#'     need not change during a person's simulated life course.}
#'   \item{**composite state**}{The combination of one value from each
#'     state variable, describing an individual's current situation.
#'     With variables for sex, nationality and parity, `"f/DE/1"` could
#'     describe a German woman with one child. MicSim joins the values
#'     with `"/"`, in the order the variables were supplied.}
#'   \item{**transient state**}{A state that can be left again.}
#'   \item{**absorbing state**}{A state from which the model allows no
#'     further transitions. Death is a common example; leaving the study
#'     population permanently is another. In [micsim_model()], a
#'     transition destination outside the transient state space is
#'     treated as absorbing. Its label need not be `"dead"`.}
#'   \item{**role**}{A condition on the state variables, such as
#'     `sex == "female"`, by which a model or a subsystem says what a state
#'     means. No label has a built-in meaning; see [move()].}
#'   \item{**subsystem**}{An additional rule for events that ordinary
#'     transitions alone cannot describe.  A birth changes the
#'     mother's state but also adds a child to the population:
#'     [fertility()] supplies that extra action. [school_enrolment()]
#'     adjusts event dates to follow a specified school-enrolment calendar.
#'     Use [add_subsystem()] to add these rules to a model.}
#' }
#'
#' @section Random numbers:
#' \describe{
#'   \item{**random number stream**}{Each individual has their own
#'     random-number stream, derived from the run's seed and the individual's
#'     ID; a newborn's stream comes from the mother's. Each time an
#'     individual's next event is scheduled, one uniform is drawn for every
#'     possible destination, also for destinations whose rate is zero
#'     or infinite.}
#'   \item{**common random numbers**}{A consequence of the random number
#'     streams is that results do not depend on the number of cores or on
#'     who else is simulated, and two scenarios simulated with the same seed
#'     use the same random numbers for every individual whose rates
#'     do not differ (common random numbers).}
#' }
#'
#' @section Rates and waiting times:
#' \describe{
#'   \item{**transition rate**, **hazard rate**, **transition intensity**}{
#'     A measure of how quickly a particular transition
#'     can occur, written \eqn{\lambda(t)}. In MicSim, rates are expressed
#'     per year. A rate is not a probability and can exceed one.}
#'   \item{**integrated hazard** (cumulative hazard)}{The rate accumulated
#'     over a period of time. If the rate stays at 2 per year for half a year,
#'     the integrated hazard is `2 * 0.5 = 1`. When the rate changes over
#'     time, the calculation accounts for those changes. For example, a rate
#'     of 0.2 per year for one year, followed by 0.4 per year for another,
#'     gives an integrated hazard of 0.6.
#'
#'   This accumulated risk is written as
#'   \eqn{\Lambda(\tau) = \int_0^{\tau} \lambda(s)\,\mathrm{d}s},
#'   where \eqn{\tau} is the interval's length in years. It represents
#'   the individual's total exposure to the risk of this particular
#'   transition over that interval.}
#'   \item{**waiting time** (sojourn time)}{How long the individual stays
#'     in the current state before a transition occurs.}
#'   \item{**inverse-transform sampling**}{A method for turning a random
#'     draw into a waiting time that follows the supplied rate. For an
#'     ordinary finite rate, MicSim draws a positive random threshold,
#'     then calculates when the integrated hazard reaches it. That time
#'     is the candidate waiting time. For example, with a constant rate
#'     of 2 per year and a threshold of 1, the candidate is half a year
#'     away. If the threshold is not reached within the remaining
#'     observation period, there is no candidate event in that period.
#'
#'   The threshold, written \eqn{u}, is drawn from an exponential
#'   distribution with mean one. The calculation finds \eqn{\tau}
#'   such that \eqn{\Lambda(\tau) = u}. Each draw generally gives a
#'   different threshold and hence a different waiting time.}
#'   \item{**forced transition**}{A transition whose rate is `Inf` at some
#'     point, encoding an event that happens on a known date rather than
#'     at random. Elementary-school enrolment is a common example.}
#' }
#'
#' @section Competing risks:
#' Several transitions may be possible from the same state. For example,
#' someone who is a former smoker might resume smoking or die. These
#' possibilities compete because whichever happens first ends the
#' current stay in that state.
#'
#' MicSim draws a candidate waiting time for each possible destination,
#' as if that transition were the only one possible. If the draw places
#' relapse two years away and death five years away, MicSim schedules
#' relapse. After relapse, it draws new candidates from the person's new
#' state, including a new candidate for death if that transition is allowed.
#' A school-calendar rule adjusts the relevant candidates before the
#' earliest is selected.
#'
#' The random draws for the different candidates are independent. This
#' is part of the sampling method; it does not claim that real-world
#' causes are unrelated.
#'
#' @section Time scales:
#' A transition rate can use three ways of measuring time:
#' \describe{
#'   \item{**age**}{Time since the individual's birth.}
#'   \item{**calendar time**}{When something happens in the study period,
#'     such as 2020 rather than 2010. A rate using `calTime` can therefore
#'     represent changes over historical time, such as falling mortality.}
#'   \item{**duration**}{Time since the last change of the state variable
#'     that the transition would change. If quitting moves someone from
#'     `"current_smoker"` to `"former_smoker"` during the simulation, a
#'     relapse rate can use `duration` to measure time since quitting.}
#' }
#' All three advance together as time passes. Someone aged 35 who quit two years earlier
#' in the simulation will be 36 with three years since quitting after
#' another year of abstinence. Calendar time also advances by one year.
#' If they resume smoking and later quit again, the new abstinent spell
#' starts with duration zero. Until the first simulated change in a state
#' variable, MicSim measures its duration from birth.
#'
#' @section Observation scheme:
#' \describe{
#'   \item{**episode**}{One uninterrupted stay in a state during
#'     observation, recorded with a start and an end date. A person who
#'     starts childless, has a first child and then remains in the
#'     one-child state until observation ends has two episodes.}
#'   \item{**left truncation**}{We know the individual was in a state at
#'     the start of observation but not when they entered it. Every member
#'     of the starting population is left-truncated.}
#'   \item{**right censoring**}{Observation ends before the state is left.
#'     Everyone still alive at the end of the simulation horizon is
#'     right-censored.}
#'   \item{**occurrence-exposure rate**}{The number of relevant events
#'     divided by the total time individuals were able to experience
#'     them. That time is called **exposure** or **person-time at risk**.
#'     For a first-birth rate, for example, count first births and sum
#'     the observed time spent in the eligible childless states.
#'
#'   Ten events during 200 person-years give `10 / 200 = 0.05` events
#'   per person-year. One person observed for two years contributes
#'   two person-years, as do two people observed for one year each.
#'   This estimates a rate, not the proportion of individuals with an
#'   event. It can be compared with the supplied rate when sufficient
#'   events and exposure are available. Sampling variation and how rates
#'   change within the estimation intervals affect this comparison.}
#' }
#'
#' @section Simulation machinery:
#' \describe{
#'   \item{**discrete-event simulation**}{A simulation that processes
#'     events in time order instead of updating everyone at fixed time
#'     steps. If the next scheduled event is on 5 March, the engine
#'     advances directly to it rather than processing each intervening
#'     day. Individuals remain in their current states between events, 
#'     although their ages and transition rates can change.}
#'   \item{**event queue**}{The list of events already scheduled but
#'     not yet processed. Each entry identifies a person, a destination
#'     state and an event time. The queue lets the engine find the next
#'     event across the whole population. A person has at most one
#'     active pending event, chosen from their competing candidates.}
#'   \item{**global clock**}{The date of the soonest pending event. It only
#'     ever moves forward.}
#' }
#'
#' @section Notation:
#' \tabular{lll}{
#'   **Symbol** \tab **Meaning** \tab **Name in the code** \cr
#'   \eqn{\lambda(x)} \tab transition rate `x` years from now \tab `lambda` \cr
#'   \eqn{\Lambda(\tau)} \tab rate accumulated over `tau` years \tab `integrated_hazard` \cr
#'   \eqn{u} \tab standard-exponential draw \tab `u` \cr
#'   \eqn{\tau} \tab waiting time, in years \tab `tau_years` \cr
#'   \eqn{E} \tab number of events in a run \tab - \cr
#'   \eqn{N} \tab number of individuals \tab - \cr
#' }
#'
#' @section Units in variable names:
#' Many internal names end with a suffix that helps identify the value:
#' \describe{
#'   \item{`_days`}{Days since 1970-01-01 for dates, or elapsed days for
#'     ages and durations. Values may include fractions of a day.}
#'   \item{`_years`}{Time in fractional years, using 365.25 days per year.}
#'   \item{`_ymd`}{A date written as year, month and day without separators,
#'     for example `20201103` for 3 November 2020.}
#'   \item{`_idx`}{A state's number in the internal mapping between labels
#'     and numbers.}
#'   \item{`_slot`}{A person's row number in the internal person store.}
#'   \item{`_id`}{A person identifier used in the inputs or output, such as
#'     `ID` or `motherID`. It need not equal their row number in the store.}
#' }
#'
#' @name MicSim-concepts
#' @keywords internal
NULL
