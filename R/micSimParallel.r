# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# II. Execute microsimulation distributed (by executing as many single thread microsimulations in parallel as cores
#     are available)
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
#'
#' micSimParallel
#'
#' Parallel Continuous-Time Microsimulation
#'
#' Executes a continuous-time microsimulation in parallel by distributing the
#' initial population and, if applicable, the immigrant population across
#' multiple CPU cores. Each core runs an independent simulation using
#' \code{micSim()}, and the resulting life histories are merged into a single
#' output population. Newborn IDs are harmonized across cores to ensure global
#' uniqueness.
#'
#' @usage
#' micSimParallel(
#'   initPop = NULL,
#'   immigrPop = NULL,
#'   initPopList = c(),
#'   immigrPopList = c(),
#'   transitionMatrix,
#'   absStates = NULL,
#'   varInitStates = c(),
#'   initStatesProb = c(),
#'   fixInitStates = c(),
#'   maxAge = 99,
#'   simHorizon,
#'   fertTr = c(),
#'   monthSchoolEnrol = c(),
#'   cores = 1,
#'   seeds = 1254
#' )
#'
#' @description
#' Performs a continuous-time microsimulation using multiple CPU cores. The
#' function splits the initial and immigrant populations across cores (either
#' automatically or according to user-provided lists), runs independent
#' simulations in parallel, and merges the results. This enables substantial
#' speed improvements for large populations.
#'
#' @param initPop A data frame describing the initial population at the start of
#'   the simulation. Must contain \code{ID}, \code{birthDate}, and
#'   \code{initState}. Used when no manual split is provided.
#'
#' @param immigrPop Optional data frame describing immigrants entering the
#'   population during the simulation. Must contain \code{ID}, \code{immigrDate},
#'   \code{birthDate}, and \code{immigrInitState}. Used when no manual split is
#'   provided.
#'
#' @param initPopList A list of data frames, each containing a subset of the
#'   initial population to be simulated on a specific core. If supplied, the
#'   function uses this split directly and does not perform automatic
#'   partitioning. Each element must have the same structure as \code{initPop}.
#'
#' @param immigrPopList A list of data frames, each containing a subset of the
#'   immigrant population assigned to a specific core. If supplied, the function
#'   uses this split directly. Each element must have the same structure as
#'   \code{immigrPop}.
#'
#' @param transitionMatrix A matrix defining the transition pattern of the
#'   multi-state model. Row names represent origin states, column names represent
#'   destination states, and entries contain the names of functions returning
#'   age- and/or time-specific transition rates.
#'
#' @param absStates A character vector of absorbing states. If not provided,
#'   absorbing states are inferred as states that appear as column names but not
#'   as row names in \code{transitionMatrix}.
#'
#' @param varInitStates A vector or matrix of substates/attributes that are
#'   assigned to newborns probabilistically according to \code{initStatesProb}.
#'
#' @param initStatesProb A vector of probabilities corresponding to
#'   \code{varInitStates}. If inheritance is used via \code{fixInitStates},
#'   probabilities must sum to one conditional on the inherited attributes;
#'   otherwise they must sum to one overall.
#'
#' @param fixInitStates Indices of substates/attributes that newborns inherit
#'   directly from their mother. If empty, no attributes are inherited.
#'
#' @param maxAge Maximum age (in years) an individual can reach during the
#'   simulation. Must be strictly positive.
#'
#' @param simHorizon A length-two vector giving the start and end dates of the
#'   simulation in \code{"yyyymmdd"} format. The start date must precede the end
#'   date.
#'
#' @param fertTr A character vector specifying transitions that trigger a birth
#'   event. Each element must be of the form \code{"stateA -> stateB"}.
#'
#' @param monthSchoolEnrol Numeric month (1–12) indicating the general school
#'   enrollment month for elementary school. If omitted and no school enrollment
#'   transition is defined, school enrollment is not simulated.
#'
#' @param cores Integer specifying the number of CPU cores to use. Must not
#'   exceed the number of available cores on the machine.
#'
#' @param seeds A numeric seed or vector of seeds used to initialize the random
#'   number generators across cores, ensuring reproducibility of the parallel
#'   simulation.
#'
#' @details
#' If \code{initPopList} and/or \code{immigrPopList} are provided, the function
#' uses these lists directly to assign individuals to cores. Otherwise, the
#' populations are automatically partitioned into approximately equal-sized
#' subsets. Each core runs an independent call to \code{micSim()}, after which
#' the resulting life histories are merged. Newborn IDs are reassigned to ensure
#' global uniqueness across cores. Parallelization is implemented using the
#' \code{snowfall} framework.
#'
#' @returns
#' A data frame containing the merged simulated population from all cores,
#' including all individuals, their transition histories, births, deaths,
#' parental links, and any events generated during the simulation. The output
#' structure matches that of \code{micSim()} but includes harmonized IDs and
#' combined results from all parallel threads.
#'
#' @export
#'
micSimParallel <- function(initPop=NULL, immigrPop=NULL, initPopList = c(), immigrPopList = c(),
                           transitionMatrix, absStates=NULL, varInitStates=c(), initStatesProb=c(),
                           fixInitStates = c(), maxAge=99, simHorizon, fertTr=c(), monthSchoolEnrol=c(),
                           cores=1, seeds=1254){

  cat('Starting at '); print(Sys.time())
  if(!is.null(initPop)) {
    N <- dim(initPop)[1]
  } else {
    if(length(initPopList)>0){
      N <- 0
      for(k in 1:cores){
        N <- N +nrow(initPopList[[k]])
      }
    } else {
      stop("No initial population for parallel computing has been given.\n")
    }
  }

  if(!is.null(immigrPop)) {
    M <- dim(immigrPop)[1]
  } else {
    M <- 0
    if(length(immigrPopList)>0){
      for(k in 1:cores){
        M <- M +nrow(immigrPopList[[k]])
      }
    }
  }

  # Split starting population and (if available) immigrant population according to available cores
  if(length(cores) %in% 0)
    stop("At least one core must be given.\n")

  condSplit <- ((length(initPopList) %in% cores) & is.null(immigrPop)) |
    ((length(initPopList) %in% cores) & (!is.null(immigrPop) & (length(immigrPopList) %in% cores)))

  if(condSplit){
    cat('\nAssign cases to distinct cores according to the split provided.\n')
    cat('Beware: It is not checked whether cases appear twice in the splits.\n')
    cat('If duplicates are in the different splits, this will result in duplicate life histories for the same entities.\n')
    cat('Thus, please check for duplicate IDs in advance.\n')
  }

  if(!condSplit) {

    if(length(initPopList)>0 & !(length(initPopList) %in% cores)){
      cat('\nSplit of initial population given for parallel computing does not match the number of cores determined.\n')
      cat('Therefore, MicSim makes an automated assignment of cases of the initial population to the distinct cores.\n')
      cat('At this, cases are distributed to the cores such that at each core approx. the same number of cases is simulated.\n')
    }
    if(!is.null(immigrPop) & (length(immigrPopList)>0 & !(length(immigrPopList) %in% cores))){
      cat('\nSplit of immigrant population given for parallel computing does not match the number of cores determined.\n')
      cat('Therefore, MicSim makes an automated assignment of cases of the immigrant population to the distinct cores.\n')
      cat('At this, cases are distributed to the cores such that at each core approx. the same number of immigrant cases is simulated.\n')
    }

    widthV <- max(trunc(N/cores), 10)
    widthW <- max(trunc(M/cores), 10)
    intV <- matrix(NA,ncol=2,nrow=cores)
    intW <- matrix(NA,ncol=2,nrow=cores)
    nI <- trunc(N/widthV)
    nIM <- trunc(M/widthW)
    ni <- 1
    for(i in 1:(nI-1)){
      intV[i,1] <- ni
      intV[i,2] <- ni+widthV-1
      ni <- ni+widthV
    }
    intV[nI,1] <- ni
    intV[nI,2] <- N
    ni <- 1
    if(nIM>1){
      for(i in 1:(nIM-1)){
        intW[i,1] <- ni
        intW[i,2] <- ni+widthW-1
        ni <- ni+widthW
      }
    }
    intW[nIM,1] <- ni
    intW[nIM,2] <- M
    initPopList <- list()
    immigrPopList <- list()
    for(core in 1:cores){
      if(!is.na(intV[core,1])){
        initPopList[[core]] <- initPop[intV[core,1]:intV[core,2],]
      } else {
        initPopList[[core]] <- NA
      }
      if(!is.na(intW[core,1])){
        immigrPopList[[core]] <- immigrPop[intW[core,1]:intW[core,2],]
      } else {
        immigrPopList[[core]] <- NA
      }
    }
  }

  nL <- cores - sum(unlist((lapply(initPopList, is.na))))
  mL <- cores - sum(unlist((lapply(immigrPopList, is.na))))

  sfInit(parallel=T,cpus=cores,slaveOutfile=NULL)
  sfLibrary("MicSim")
  sfExportAll(debug=FALSE)
  sfClusterSetupRNGstream(seed=(rep(seeds,35)[1:length(cores)]))
  myPar <- function(itt){
    #cat('Starting thread: ',itt,'\n')
    if(itt<=mL){
      immigrPopL <- immigrPopList[[itt]]
    } else {
      immigrPopL <- NULL
    }
    if (itt<=nL) {
      initPopL <- initPopList[[itt]]
    } else {
      initPopL <- NULL
      stop("\nCompared to the number of migrants, the starting population is too small to justify running a distributed simulation on several cores.")
    }
    popIt <- micSim(initPop=initPopL, immigrPop=immigrPopL, transitionMatrix=transitionMatrix,
                    absStates=absStates, varInitStates=varInitStates, initStatesProb=initStatesProb,
                    fixInitStates=fixInitStates, maxAge=maxAge, simHorizon=simHorizon, fertTr=fertTr,
                    monthSchoolEnrol=monthSchoolEnrol)
    #cat('Thread: ',itt,' has stopped.\n')
    return(popIt)
  }
  pop <- sfLapply(1:max(nL,mL), myPar)
  # create unique IDs for newborns
  refID <- 0
  replaceID <- function(rr){
    pop[[i]][which(as.numeric(pop[[i]][,1])==rr[1]),1] <<- rr[2]
    return(NULL)
  }
  for(i in 1:length(pop)){
    if(!is.na(immigrPopList[[i]])[1]){
      allIDs <- c(initPopList[[i]]$ID, immigrPopList[[i]]$ID)
    } else {
      allIDs <- initPopList[[i]]$ID
    }
    exIDs <- unique(as.numeric(pop[[i]][,1]))
    repl <- setdiff(exIDs, allIDs)
    if(length(repl)>0) {
      newIDs <- cbind(repl,-(refID+(1:length(repl))))
      idch <- apply(newIDs,1,replaceID)
      refID <- max(abs(newIDs[,2]))
    }
  }
  pop <- do.call(rbind,pop)
  pop[as.numeric(pop[,1])<0,1]  <- abs(as.numeric(pop[as.numeric(pop[,1])<0,1]))+N+M
  pop <- pop[order(as.numeric(pop[,1])),]
  sfStop()
  cat('Stopped at '); print(Sys.time())
  return(pop)
}



