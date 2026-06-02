#' Rate Data for Simulation of MDD example
#'
#' Data used to construct transition rates from @lepe2024
#'
#' @format ## `rates_mdd`
#' A data frame with 49 rows and 17 columns:
#' \describe{
#'   \item{age}{Age of individual}
#'   \item{incidence_lo_m}{Incidence Rates low education men}
#'   \item{cf_incidence_lo_f}{Counterfactual Incidence Rate low female}
#'
#' }
#' @source https://doi.org/10.1093/eurpub/ckae066

"rates_mdd"

#' Initial Population Data for Examples on Simulating Migration Flows
#'
#' Data used to construct initial population for Migration Flow Vignette
#'
#' @format ## `immigrPopMigrExp`
#' A data frame with 3758 rows and 4 columns:
#' \describe{
#'   \item{ID}{ID of Individual}
#'   \item{immigrDate}{Immigration Date}
#'   \item{birthDate}{Birth Date}
#'   \item{immigrInitState}{initial State of Individual as per stateSpace}
#' }

"immigrPopMigrExp"

#' Immigrant Population Data for Examples on Simulating Migration Flows
#'
#' Data used to construct initial immigrant population for Migration Flow Vignette
#'
#' @format ## `initPopMigrExp`
#' A data frame with 72965 rows and 3 columns:
#' \describe{
#'  \item{ID}{Id of Individual}
#'  \item{birthDate}{Birth Date}
#'  \item{immigrInitState}{initial State of Individual as per stateSpace}
#' }

"initPopMigrExp"

#' Rate Data for Examples on Simulating Migration Flows
#'
#' Data used to construct transition rates for Migration Flow Vignette
#'
#' @format ## `migrExpRates`
#' A data frame with 100 rows and 30 columns:
#' \describe{
#'   \item{mort_f_ES}{Mortality Rate for females from Spain}
#'   \item{mort_f_NL}{Mortality Rate for females from the Netherlands}
#'   \item{mort_f_SE}{Mortality Rate for females from Swedem}
#'
#' }

"migrExpRates"

#' Eurostat Population Data (ES) – EPC Extract
#'
#' Population counts by age, sex and year for Germany, Spain and Hungary, extracted from Eurostat
#' (dataset: DEMO_PJAN). Used as an example to construct a synthetic population.
#'
#' @format A data frame with 54,624 rows and 21 variables:
#' \describe{
#'   \item{STRUCTURE}{Metadata field identifying the SDMX structure.}
#'   \item{STRUCTURE_ID}{Eurostat SDMX identifier (e.g., "ESTAT:DEMO_PJAN(1.0)").}
#'   \item{STRUCTURE_NAME}{Full dataset name ("Population on 1 January by age and sex").}
#'   \item{freq}{Frequency code ("A" for annual).}
#'   \item{Time.frequency}{Frequency label.}
#'   \item{unit}{Unit code ("NR").}
#'   \item{Unit.of.measure}{Unit label ("Number").}
#'   \item{age}{Age category (e.g., "TOTAL", "0", "1", ...).}
#'   \item{Age.class}{Age label.}
#'   \item{sex}{Sex code ("F", "M", "T").}
#'   \item{Sex}{Sex label.}
#'   \item{geo}{Country code ("DE").}
#'   \item{Geopolitical.entity..reporting.}{Country label ("Germany").}
#'   \item{TIME_PERIOD}{Calendar year.}
#'   \item{Time}{Unused field (always NA).}
#'   \item{OBS_VALUE}{Population count.}
#'   \item{Observation.value}{Unused field (always NA).}
#'   \item{OBS_FLAG}{Observation flag.}
#'   \item{Observation.status..Flag..V2.structure}{Flag label.}
#'   \item{CONF_STATUS}{Confidentiality flag (NA).}
#'   \item{Confidentiality.status..flag.}{Confidentiality label (NA).}
#' }
#"demo_pjan_custom"

#' Eurostat Housing & Tenure Data (ES) – EPC Extract
#'
#' Distribution of population by tenure status, household type and income group
#' for Spain, extracted from Eurostat (dataset: ILC_LVHO02).
#'
#' @format A data frame with 6,783 rows and 23 variables:
#' \describe{
#'   \item{STRUCTURE}{Metadata field identifying the SDMX structure.}
#'   \item{STRUCTURE_ID}{Eurostat SDMX identifier ("ESTAT:ILC_LVHO02(1.0)").}
#'   \item{STRUCTURE_NAME}{Dataset name.}
#'   \item{freq}{Frequency code ("A").}
#'   \item{Time.frequency}{Frequency label.}
#'   \item{incgrp}{Income group code.}
#'   \item{Income.situation.in.relation.to.the.risk.of.poverty.threshold}{Income group label.}
#'   \item{hhcomp}{Household composition code.}
#'   \item{Household.composition}{Household composition label.}
#'   \item{tenure}{Tenure code ("OWN", "RENT", etc.).}
#'   \item{Tenure.status}{Tenure label.}
#'   \item{unit}{Unit code ("PC").}
#'   \item{Unit.of.measure}{Unit label ("Percentage").}
#'   \item{geo}{Country code ("DE").}
#'   \item{Geopolitical.entity..reporting.}{Country label.}
#'   \item{TIME_PERIOD}{Calendar year.}
#'   \item{Time}{Unused field (NA).}
#'   \item{OBS_VALUE}{Percentage value.}
#'   \item{Observation.value}{Unused field (NA).}
#'   \item{OBS_FLAG}{Observation flag.}
#'   \item{Observation.status..Flag..V2.structure}{Flag label.}
#'   \item{CONF_STATUS}{Confidentiality flag (NA).}
#'   \item{Confidentiality.status..flag.}{Confidentiality label (NA).}
#' }
#"ilc_lvho02_custom"

#' Female Exposure to Risk by Year, Age and Parity (Spain)
#'
#' Raw exposure-to-risk table extracted from the Human Fertility Database.
#' Contains exposure counts by calendar year, age and parity.
#'
#' @format A data frame with 1,454 rows and 3 variables:
#' \describe{
#'   \item{Spain}{Raw text column containing metadata and exposure rows.}
#'   \item{Female.exposure.to.risk.by.calendar.year}{Placeholder column (NA).}
#'   \item{age.and.parity}{Placeholder column (NA).}
#' }
#'
#' @source https://www.humanfertility.org/Country/Country?cntr=ESP
#'
# "ex_year_age_parity"

#' Spanish Mortality Rates (Period 1x1)
#'
#' Raw mortality rates by year, age and sex for Spain, extracted from the Human Mortality Database.
#'
#' @format A data frame with 12,877 rows and 3 variables:
#' \describe{
#'   \item{Spain}{Raw text column containing metadata and mortality rows.}
#'   \item{Death.rates..period.1x1.}{Placeholder column (NA).}
#'   \item{Last.modified..20.Feb.2025...Methods.Protocol..v6..2017.}{Placeholder column (NA).}
#' }
#'
#' @source https://mortality.org/Country/Country?cntr=ESP
#'
# "es_m1x1"

#' Spanish Age-Specific Fertility Rates (ASFR)
#'
#' Raw ASFR values by year and age for Spain, extracted from the Human Fertility Database.
#'
#' @format A data frame with 4,490 rows and 3 variables:
#' \describe{
#'   \item{Spain}{Raw text column containing metadata and ASFR rows.}
#'   \item{Period.fertility.rates.by.calendar.year.and.age..Lexis.squares}{Placeholder column (NA).}
#'   \item{age.in.completed.years..ACY..}{Placeholder column (NA).}
#' }
#'
#' @source https://www.humanfertility.org/Country/Country?cntr=ESP#age
#'
# "es_asfr"

#' Fertility Rates by Year, Age and Parity (Spain)
#'
#' Fertility rates extracted from the HFD database, including
#' age-specific fertility rates (ASFR) and parity-specific fertility rates
#' for Spain. Used in microsimulation examples.
#'
#' @format A data frame with 840 rows and 8 variables:
#' \describe{
#'   \item{year}{Calendar year.}
#'   \item{age}{Age in completed years.}
#'   \item{asfr}{Age-specific fertility rate (all parities combined).}
#'   \item{asfr_0}{Fertility rate for parity 0 (first births).}
#'   \item{asfr_1}{Fertility rate for parity 1.}
#'   \item{asfr_2}{Fertility rate for parity 2.}
#'   \item{asfr_3}{Fertility rate for parity 3.}
#'   \item{asfr_4plus}{Fertility rate for parity 4 or higher.}
#' }
#'
#'
"fert_parity"

#' Mortality Rates by Year, Age and Sex (Spain)
#'
#' Mortality rates extracted from the HMD mortality database for Spain.
#' Contains period 1x1 mortality rates by age and sex.
#'
#' @format A data frame with 12,876 rows and 5 variables:
#' \describe{
#'   \item{year}{Calendar year.}
#'   \item{age}{Age (character-coded).}
#'   \item{female}{Female mortality rate (character-coded numeric).}
#'   \item{male}{Male mortality rate (character-coded numeric).}
#'   \item{total}{Total mortality rate (character-coded numeric).}
#' }
#'
#'
"mort_es"

#' Example Initial Population for Microsimulation
#'
#' Synthetic initial population used in microsimulation examples.
#' Contains individual IDs, birthdates and initial states.
#'
#' @format A tibble with 5,000 rows and 3 variables:
#' \describe{
#'   \item{id}{Individual identifier.}
#'   \item{birthdate}{Birthdate in YYYYMMDD numeric format.}
#'   \item{initState}{Initial state at simulation start (e.g., "f/h/0").}
#' }
"initpop_epc_example"

#' #' Cleaned Population Structure (Spain)
#' #'
#' #' Population counts by year, age and sex, cleaned and harmonized for use
#' #' in demographic simulations.
#' #'
#' #' @format A tibble with 15,942 rows and 4 variables:
#' #' \describe{
#' #'   \item{year}{Calendar year.}
#' #'   \item{age}{Age category (Eurostat Y-coding).}
#' #'   \item{sex}{Sex code ("F", "M").}
#' #'   \item{population}{Population count.}
#' #' }
#' "pop_structure"

#' Housing Tenure Data (Spain)
#'
#' Housing tenure distribution for Spain, extracted from Eurostat dataset
#' ILC_LVHO02. Contains the share of population living in different tenure
#' arrangements over time.
#'
#' @format A tibble with 133 rows and 3 variables:
#' \describe{
#'   \item{year}{Calendar year.}
#'   \item{tenure}{Tenure category (e.g., "OWN").}
#'   \item{obs}{Observed percentage (character-coded numeric).}
#' }
"housing_es"

#' Homeownership Rates by Year and Age (Spain)
#'
#' Cross-sectional prevalence rates extracted from EU SILC (EU-SILC 2004-2023 Cross-sectional).
#' Contains homeownership rates by age and sex back filled from 1950 until 2023.
#'
#' @format A numeric matrix with 101 rows (ages 0-100) and 24 columns (years 1950-2023):
#' \describe{
#'   \item{Rows}{Age (0-100)}
#'   \item{Columns}{Years (1900-2023)}
#'   \item{Values}{Homeownership transition rates (0-1). Ages <20 and >70 are 0.
#'   Years 1950-2004 are backfilled from 2005. Years >2020 use 2020 data.}
#' }
#'
#' @source EU-SILC cross-sectional data for Spain (2004-2023).
#' Years 1950-2004 extrapolated from 2005 data.
#'
#' @examples
#' #Get homeownership rate for age 30 in year 2010
#' ho_matrix["30", "2010"]
#'
"ho_matrix"
