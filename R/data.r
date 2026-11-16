#' Immigrant population data for the migration-flow examples
#'
#' Entrants arriving after the simulation start, used in the migration-flow
#' vignette.
#'
#' @format ## `immigrPopMigrExp`
#' A data frame with 3758 rows and 4 columns:
#' \describe{
#'   \item{ID}{individual identifier}
#'   \item{immigrDate}{immigration date}
#'   \item{birthDate}{birth date}
#'   \item{immigrInitState}{initial state of the individual, as per the state space}
#' }

"immigrPopMigrExp"

#' Initial population data for the migration-flow examples
#'
#' Starting population used in the migration-flow vignette.
#'
#' @format ## `initPopMigrExp`
#' A data frame with 72965 rows and 3 columns:
#' \describe{
#'  \item{ID}{individual identifier}
#'  \item{birthDate}{birth date}
#'  \item{initState}{initial state of the individual, as per the state space}
#' }

"initPopMigrExp"

#' Transition-rate data for the migration-flow examples
#'
#' Data used to construct transition rates in the migration-flow vignette.
#'
#' @format ## `migrExpRates`
#' A data frame with 100 rows and 30 columns:
#' \describe{
#'   \item{mort_f_ES}{mortality rate for females from Spain}
#'   \item{mort_f_NL}{mortality rate for females from the Netherlands}
#'   \item{mort_f_SE}{mortality rate for females from Sweden}
#'   ...
#' }

"migrExpRates"

#' Transition-rate data for the MDD example
#'
#' Data used to construct transition rates for the MDD example.
#'
#' @format ## `rates_mdd`
#' A data frame with 49 rows and 17 columns:
#' \describe{
#'   \item{age}{age of the individual}
#'   \item{incidence_lo_m}{incidence rate, low education, men}
#'   \item{cf_incidence_lo_f}{counterfactual incidence rate, low education, women}
#'   ...
#' }
#' @source https://doi.org/10.1093/eurpub/ckae066
"rates_mdd"
