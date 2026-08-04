#####################################################################################
#####################################################################################
### AUXILIARY FUNCTIONS                                                            ##
### - TO DEFINE MICROSIMULATION INPUT                                              ##
### - FOR COMPUTATION WITH DATES                                                   ##
### SZ, June 2022                                                                  ##
#####################################################################################
#####################################################################################
#'
#' builtStatesCodes
#'
#' Assign to all states and subStates numerical codes
#' Build Transition Matrix
#'
#' @usage builtStatesCodes(transitionMatrix)
#' @description The buildStatesCodes function converts a state transition matrix into a numerical coding scheme.
#' Each state in the transition matrix is divided into sub-states (if separated by ???/???) and each sub-state is mapped to a unique numerical code.
#'
#' @param transitionMatrix Is the row and column names represent states.
#' States can be compound, separated by ???/???.
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
#' @export
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
    args <- names(formals(rates[i]))
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
  return(depMatrix)
  
}

