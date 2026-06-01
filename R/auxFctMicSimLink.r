#####################################################################################
#####################################################################################
### AUXILIARY FUNCTIONS                                                            ##
### - TO DEFINE MICROSIMULATION INPUT                                              ##
### - FOR COMPUTATION WITH DATES                                                   ##
### SZ, June 2022                                                                  ##
#####################################################################################
#####################################################################################
#' build_genMatrix
#'
#' Build Generation Matrix
#'
#' Construct a matrix specifying inheritance or intergenerational transition rules.
#' Function to save the user-defined info given to process intergenerational transitions into the standardized output needed by MicSim.
#'
#' @usage build_genMatrix(geninfo)
#'
#' @param geninfo A matrix containing user‑defined inheritance
#'   specifications. Each row should provide a transition identifier, the
#'   relationship attribute to update (e.g. "motherID", "fatherID"), and the
#'   parameters alpha and beta used to compute the inheritance effect.
#'
#' @returns A 4‑column matrix with standardized fields
#'   (`tr`, `inherit`, `alpha`, `beta`) used by MicSim to process
#'   intergenerational transitions.
#'
#' @examples
#' \dontrun{
#'   geninfo <- c("ownhome", "motherID", 0, 0.5)
#'   genMatrix <- build_genMatrix(geninfo)
#' }
#'
#' @export
#'
build_genMatrix <- function(geninfo){ #TODO: check robustness!
  genMatrix <- matrix(0, ncol = 4)
  genMatrix <- rbind(geninfo, deparse.level = 0)
  genMatrix <- as.matrix(genMatrix)
  colnames(genMatrix) <- c("tr", "inherit", "alpha", "beta") # second row: parents alpha = 0, beta = 1
  return(genMatrix)
}
#'
#' builtStatesCodes
#'
#' Assign to all states and subStates numerical codes
#' Build Transition Matrix
#'
#' @usage builtStatesCodes(transitionMatrix)
#' @description The buildStatesCodes function converts a state transition matrix into a numerical coding scheme.
#' Each state in the transition matrix is divided into sub-states (if separated by “/”) and each sub-state is mapped to a unique numerical code.
#'
#' @param transitionMatrix The row and column names represent states.
#' States can be compound, separated by “/”.
#'
#' @returns A data frame containing the coding scheme for the states.
#' @keywords internal
#' @noRd
#'
builtStatesCodes <- function(transitionMatrix){

  transitionMatrixNum <- transitionMatrix
  allStates <- rownames(transitionMatrix)
  allStatesMatrix <- do.call(cbind,sapply(allStates, strsplit, "/"))
  absStates <- setdiff(colnames(transitionMatrix), rownames(transitionMatrix))
  codeList <- vector(length=nrow(allStatesMatrix)+1, mode="list")
  for(j in 1:nrow(allStatesMatrix)){
    usj <- unique(allStatesMatrix[j,])
    codeList[[j]] <- cbind(usj, 1:length(usj))
    codeList[[j]][,2] <- ifelse(nchar(codeList[[j]][,2])%in% 1,paste("0", codeList[[j]][,2], sep=""),nchar(codeList[[j]][,2]))
  }
  codeList[[j+1]] <- cbind(absStates, -c(1:length(absStates)))

  allCodes <- c()
  allCodesSep <- NULL
  for(i in 1:length(allStates)){
    st <- unlist(strsplit(allStates[i], "/"))
    stNum <- c()
    for(k in 1:length(st)){
      subCode <- codeList[[k]][codeList[[k]][,1] %in% st[k],2]
      stNum <- c(stNum, subCode)
    }
    allCodes <- c(allCodes, paste(stNum, collapse = ""))
    allCodesSep <- rbind(allCodesSep, as.numeric(stNum))
  }
  codesSepAbs <- matrix(rep(NA, length(absStates)*ncol(allCodesSep)), ncol=ncol(allCodesSep), nrow=length(absStates))
  codesSepAbs[,1] <- -c(1:length(absStates))
  codingScheme <- data.frame(allStates = c(allStates,absStates),
                             allCodes= as.numeric(c(allCodes, -c(1:length(absStates)))),
                             allCodesSep=rbind(allCodesSep,codesSepAbs))
  return(codingScheme=codingScheme)
}
#'
#' eventExposureMatrices
#'
#' Prepare Age- and Year-Specific Event and Exposure Matrices
#'
#' Constructs matrices of exposures and events by single-year ages and calendar
#' years based on simulated population data. Individuals are included only if
#' alive and part of the population at the start of each year, optionally
#' accounting for immigration. Exposure is counted for individuals whose state
#' at mid-year belongs to the specified risk set. Events are counted when
#' transitions match the user-defined event types.
#'
#' @usage
#' eventExposureMatrices(
#'   pop,
#'   riskSet,
#'   events,
#'   ages = 0:(maxAge - 1),
#'   absStates = "dead",
#'   immigrPop = NULL,
#'   confint = FALSE
#' )
#'
#' @param pop A data frame describing the simulated population, containing at
#'   minimum the variables \code{ID}, \code{birthDate}, \code{From}, \code{To},
#'   \code{transitionTime}, \code{transitionAge}, and \code{initState}. Each row
#'   represents a transition or state episode for an individual.
#'
#' @param riskSet A character vector of states considered to be at risk. Only
#'   individuals in these states at mid-year contribute exposure time.
#'
#' @param events A two-column matrix or data frame specifying the transitions to
#'   be counted as events. Column 1 contains origin states (\code{From}), and
#'   column 2 contains destination states (\code{To}). Only transitions matching
#'   a row in this table are counted.
#'
#' @param ages A vector of integer ages for which exposure and event counts
#'   should be computed. Defaults to all ages from 0 to \code{maxAge - 1}.
#'
#' @param absStates A character vector of absorbing states (e.g., \code{"dead"}).
#'   Individuals who have entered an absorbing state before the start of a given
#'   year are excluded from exposure and event calculations for that year.
#'
#' @param immigrPop Optional data frame describing immigration into the
#'   population. Must contain \code{ID} and \code{immigrDate}. Individuals whose
#'   immigration date is after the current calendar year are excluded from that
#'   year's risk set.
#'
#' @param confint Logical; whether to compute confidence intervals. Currently not
#'   implemented but reserved for future extensions.
#'
#' @returns A list with two matrices:
#'   \itemize{
#'     \item \code{eventMatrix}: counts of events by age (rows) and calendar year (columns)
#'     \item \code{exposureSet}: counts of individuals at risk by age (rows) and calendar year (columns)
#'   }
#'
#' @keywords internal
#' @noRd
#'
eventExposureMatrices <- function(pop, riskSet, events, ages=c(0:(maxAge-1)), absStates = "dead", immigrPop=NULL, confint=FALSE){

  years <- c(trunc(startDate/10000):trunc(endDate/10000))
  expMatrix <- matrix(0, ncol=length(years), nrow=length(ages)) # along ages and calendar years (ages measured at the end of a year)
  evMatrix <- matrix(0, ncol=length(years), nrow=length(ages)) # along ages and calendar years (ages measured at the end of a year)

  for(year in years) {
      # kick out who is already dead or has left population at year start
      popOut <- pop[pop$To %in% absStates & trunc(as.numeric(pop$transitionTime)/10000)<=year,]
      popRed <- pop[!(pop$ID %in% popOut$ID),]

      # kick out who is in year not yet part of population (immigrants and newborns)
      if(!is.null(immigrPop)){
        notYetImmigr <- immigrPop[trunc(as.numeric(immigrPop$immigrDate)/10000)>year, "ID"]
        popRed <- popRed[!(popRed$ID %in% notYetImmigr),]
      }
      notYetBorn <- unique(as.numeric(pop[trunc(as.numeric(pop$birthDate)/10000)>year, "ID"]))
      popRed <- popRed[!(popRed$ID %in% notYetBorn),]
      midyear <- getInDays(as.numeric(paste(year, "0701", sep="")))
      popRed$ageInYear <- trunc(c(midyear - getInDays(popRed$birthDate))/365.25)

      for(age in ages){
        popAge <- popRed[popRed$ageInYear %in% age,]
        ids <- unique(popAge$ID)
        for(i in ids){
          setI <- popAge[popAge$ID %in% i,]
          # if i experienced transition during sim
          if(!is.na(setI$From[1])){
            stateAtAge <- setI[setI$transitionAge > age,][1,"From"]
          }
          # if i has not experienced any transition during sim
          if(is.na(setI$From[1])){
            stateAtAge <- setI$initState
          }
          if(stateAtAge %in% riskSet){
            expMatrix[age+1, year-years[1]+1] <- expMatrix[age+1, year-years[1]+1] + 1 # TODO: for the moment I count whole years as exposure
          }
        }
        popEventAge <- popAge[trunc(popAge$transitionAge) %in% age,,drop=FALSE]

        if(nrow(popEventAge)>0) {
          for(k in c(1:nrow(popEventAge))){ # theoretically possible several events during a year
            rowF <- which(popEventAge[k,"From"] %in% events[,1])
            rowT <- which(popEventAge[k,"To"] %in% events[,2])
            if(length(rowF)==1 & length(rowT)==1){
              if(rowF == rowT) {
                evMatrix[age+1, year-years[1]+1] <- evMatrix[age+1, year-years[1]+1] + 1
              }
            }
          }
        }
      } # end age
    } # end year
  return(list(eventMatrix=evMatrix, exposureSet=expMatrix))
}
#'
#' estimateAgeFreq
#'
#' Estimate Age-Specific Relative Frequencies
#'
#' Computes age-specific relative frequencies (rates) by dividing the total
#' number of events by the total exposure across all calendar years. Event and
#' exposure counts are obtained from \code{eventExposureMatrices()}, and
#' frequencies are calculated for each single-year age group.
#'
#' @usage
#' estimateAgeFreq(
#'   pop,
#'   riskSet,
#'   events,
#'   ages = 0:(maxAge - 1),
#'   absStates = "dead",
#'   immigrPop = NULL,
#'   confint = FALSE
#' )
#'
#' @param pop A data frame describing the simulated population, containing at
#'   minimum the variables \code{ID}, \code{birthDate}, \code{From}, \code{To},
#'   \code{transitionTime}, \code{transitionAge}, and \code{initState}. Each row
#'   represents a transition or state episode for an individual.
#'
#' @param riskSet A character vector of states considered to be at risk. Only
#'   individuals in these states contribute exposure time.
#' @param events A two-column matrix or data frame specifying the transitions to
#'   be counted as events. Column 1 contains origin states (\code{From}), and
#'   column 2 contains destination states (\code{To}).
#' @param ages A vector of integer ages for which age-specific frequencies should
#'   be estimated. Defaults to all ages from 0 to \code{maxAge - 1}.
#'
#' @param absStates A character vector of absorbing states (e.g., \code{"dead"}).
#'   Individuals who have entered an absorbing state before the start of a given
#'   year are excluded from exposure and event calculations.
#'
#' @param immigrPop Optional data frame describing immigration into the
#'   population. Must contain \code{ID} and \code{immigrDate}. Individuals whose
#'   immigration date is after the current calendar year are excluded from that
#'   year's risk set.
#'
#' @param confint Logical; whether to compute confidence intervals. Currently not
#'   implemented but reserved for future extensions.
#'
#' @returns A numeric vector of age-specific relative frequencies (events divided
#'   by exposure) for each age in \code{ages}.
#'
#' @keywords internal
#' @noRd
#'
estimateAgeFreq <- function(pop, riskSet, events, ages=c(0:(maxAge-1)), absStates = "dead", immigrPop=NULL, confint=FALSE){

  evExpList <- eventExposureMatrices(pop, riskSet, events, ages=c(0:(maxAge-1)), absStates = "dead", immigrPop=NULL, confint=FALSE)

  evAges <- apply(evExpList$eventMatrix, 1,sum)
  expAges <- apply(evExpList$exposureMatrix, 1,sum) # whole years
  ageRates <- evAges/expAges
  ageRates[ageRates %in% NaN] <- 0
  return(ageRates)
}
#'
#' rate_cS
#'
#' Construct matrix noting dependencies of functions
#' Build dependency matrix
#'
#' @usage rate_cS(allTr)
#' @description Helper function to track arguments used by transition rate functions and if intergenerational transmission is used for the microsimulation.
#' A function for building a matrix containing the arguments of the transition rate functions.
#'
#' @param allTr A matrix containing names of all transition rate functions.
#'
#' @returns A dependency matrix
#' @keywords internal
#' @noRd
#'
#'
rate_cS <- function(allTr){
  #TODO:  Check if genMatrix exists in the global environment -> print error message otherwise
  genMatrix_exists <- exists("genMatrix", envir = .GlobalEnv)

  rates <-  unique(allTr)

  form <- c("age", "calTime", "duration", "genarg")

  depMatrix <- matrix(0, nrow = length(rates), ncol = 4)

  colnames(depMatrix) <- form
  rownames(depMatrix) <- rates

  col_count <- rep(0, ncol(depMatrix))

  for(i in 1:nrow(depMatrix)){
    args <- names(formals(match.fun(rates[i])))
      for(k in 1:length(args)) {
        col <- match(args[k], colnames(depMatrix))
        col_count[col] <- col_count[col] + 1
        depMatrix[i, col] <- col_count[col]
      }

    if (genMatrix_exists && rates[i] %in% genMatrix[, 1]) {
      col <- match("genarg", colnames(depMatrix))
      col_count[col] <- col_count[col] + 1
      depMatrix[i, col] <- col_count[col]
    }
  }

  depMatrix <- as.matrix(depMatrix)
  # avoid warning
  depMatrix[,3] <- depMatrix[,3] > 0
  depMatrix[,4] <- depMatrix[,4] > 0
  return(depMatrix)

}
