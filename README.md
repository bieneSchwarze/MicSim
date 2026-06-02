# MicSim

MicSim is an R package for **continuous-time microsimulation**, designed to model life-course trajectories of individuals or populations under complex demographic, social, or economic processes. It provides tools to define state spaces, transition rates, event histories, and to simulate individual life paths using flexible hazard-based models.

MicSim is suitable for applications in demography, epidemiology, social sciences, and any domain where event-driven stochastic processes are required.

Microsimulation models represent individual life courses by simulating transitions between discrete states (e.g., employment, health status, marital status). MicSim implements these models in continuous time, allowing:

-   Arbitrary state spaces
-   Time-dependent transition rates
-   Individual-level heterogeneity
-   Event history generation
-   Simulation of large synthetic populations

The MicSim package was developed following the methodological concepts introduced in Zinn's dissertation "Microsimulation of Life Histories: A Continuous-Time Approach" (University of Rostock, 2011). The dissertation outlines the statistical, demographic, and computational principles that form the core of MicSim's continuous-time microsimulation engine.

[Zinn, S. (2011). Microsimulation of Life Histories: A Continuous-Time Approach. University of Rostock.](https://rosdok.uni-rostock.de/file/rosdok_derivate_0000004766/Dissertation_Zinn_2011.pdf)

The vignettes included in this package provide a more detailed explanation of how **MicSim** can be applied in research settings, offering step-by-step examples and methodological guidance.

### Installing the micsimlink-v1.0 Version of MicSim for **EPC 2026**

A dedicated tagged release of the MicSimLink extension is available under the tag `micsimlink-v1.0` in the GitHub repository. It is possible to install this specific version directly from GitHub using the [remotes package](https://remotes.r-lib.org/).

#### Step-by-step installation
1.  Make sure the remotes package is installed: 

```
install.packages("remotes")
library(remotes)
```

2.  Install the `micsimlink-v1.0` version of MicSim:

```
remotes::install_github("bieneSchwarze/MicSim@micsimlink-v1.0", build_vignettes=TRUE)
```

3.  Load the package and check the Vignettes: 

```
library(MicSim)
browseVignettes(package = "MicSim")
```

------------------------------------------------------------------------

#### Advanced Example: Fertility, Mortality, and Maternal Attribute Inheritance

This illustrative example shows how to prepare the input needed by MicSim to simulate a population with: - Mortality (Gompertz model) - Fertility (Hadwiger mixture model) - Inheritance of nationality from the mother - Multiple sub-states (sex x nationality x fertility status)

```         
# Clean workspace
rm(list = ls())

# Simulation horizon
startDate <- 20140101
endDate   <- 20241231
simHorizon <- c(startDate = startDate, endDate = endDate)

set.seed(234)

# Maximal age
maxAge <- 100

# State space
sex <- c("m", "f")
nat <- c("DE", "AT", "IT")
fert <- c("0", "1")
stateSpace <- expand.grid(sex = sex, nat = nat, fert = fert)
absStates <- "dead"

# Initial population (random example)
N <- 100
birthDates <- runif(N, min = getInDays(19500101), max = getInDays(20131231))

getRandInitState <- function(birthDate) {
  age <- trunc((getInDays(simHorizon[1]) - birthDate) / 365.25)
  s1 <- sample(sex, 1)
  s2 <- sample(nat, 1)
  s3 <- ifelse(age <= 18, fert[1], sample(fert, 1))
  paste(c(s1, s2, s3), collapse = "/")
}

initPop <- data.frame(
  ID = 1:N,
  birthDate = birthDates,
  initState = sapply(birthDates, getRandInitState)
)
initPop$birthDate <- getInDateFormat(initPop$birthDate)

# Initial states for newborns
fixInitStates <- 2  # inherit nationality from mother

varInitStates <- rbind(
  c("m", "DE", "0"), c("f", "DE", "0"),
  c("m", "AT", "0"), c("f", "AT", "0"),
  c("m", "IT", "0"), c("f", "IT", "0")
)

initStatesProb <- c(
  0.515, 0.485,
  0.515, 0.485,
  0.515, 0.485
)

# Transition rates
fertRates <- function(age, calTime) {
  b <- ifelse(calTime <= 2020, 3.5, 3.0)
  c <- ifelse(calTime <= 2020, 28, 29)
  rate <- (b / c) * (c / age)^(3/2) * exp(-b^2 * (c/age + age/c - 2))
  rate[age <= 15 | age >= 45] <- 0
  rate
}

mortRates <- function(age, calTime) {
  a <- 0.00003
  b <- ifelse(calTime <= 2020, 0.1, 0.097)
  a * exp(b * age)
}

# Transition matrices
fertTrMatrix <- cbind(
  c("f/DE/0->f/DE/1", "f/AT/0->f/AT/1", "f/IT/0->f/IT/1"),
  rep("fertRates", 3)
)

absTransitions <- cbind(
  c("f/DE/dead", "f/AT/dead", "f/IT/dead",
    "m/DE/dead", "m/AT/dead", "m/IT/dead"),
  rep("mortRates", 6)
)

transitionMatrix <- buildTransitionMatrix(
  allTransitions = fertTrMatrix,
  absTransitions = absTransitions,
  stateSpace = stateSpace
)

# Fertility-triggering transitions
fertTr <- fertTrMatrix[, 1]

# Run microsimulation
pop <- micSim(
  initPop = initPop,
  transitionMatrix = transitionMatrix,
  absStates = absStates,
  varInitStates = varInitStates,
  initStatesProb = initStatesProb,
  fixInitStates = fixInitStates,
  maxAge = maxAge,
  simHorizon = simHorizon,
  fertTr = fertTr
)
```
