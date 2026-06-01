####################################################################################
####################################################################################
## FUNCTION EXECUTING MICROSIMULATION WITH LINKED LIVES                           ##
## SZ, December 2023                                                              ##
####################################################################################
####################################################################################
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# I. Execute microsimulation as single thread
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
#'
#' micSimLink
#'
#' Microsimulation With Linked Lives
#'
#' Performs continuous-time microsimulation of individual life courses with
#' linked lives, allowing for fertility, partnership formation, separation,
#' widowhood, immigration, school enrollment, and duration-dependent transition
#' rates. Attributes can also influence following generations through intergenerational transmission. The simulation is driven by a user-specified transition matrix and
#' supports inheritance of attributes from mothers to newborns as well as
#' probabilistic assignment of non-inherited characteristics.
#'
#' @usage
#' micSimLink(
#'   initPop,
#'   immigrPop = NULL,
#'   transitionMatrix,
#'   absStates = NULL,
#'   fixInitStates = c(),
#'   varInitStates = c(),
#'   initStatesProb = c(),
#'   maxAge = 99,
#'   simHorizon,
#'   fertTr = c(),
#'   partTr = c(),
#'   rule = 0,
#'   ageDiffDistr = NULL,
#'   sepTr = c(),
#'   probSepTr = c(),
#'   absPartTr = c(),
#'   monthSchoolEnrol = c(),
#'   genMatrix = c()
#' )
#'
#' @description
#' Simulates life courses in continuous time for a population in which
#' individuals may form partnerships, dissolve them, have children, migrate into
#' the population, and experience transitions defined by a multi-state model.
#' The function supports linked lives, meaning that transitions of one
#' individual can affect the state trajectories of others either vertically (e.g., partnership formation, separation, widowhood) or horizontally (e.g., a mother's state influencing a child's transition rates).
#'
#' @param initPop A data frame describing the initial population at the start of
#'   the simulation. Must contain the variables \code{ID}, \code{birthDate}, and
#'   \code{initState}.
#' @param immigrPop Optional data frame describing immigrants entering the
#'   population during the simulation. Must contain \code{ID}, \code{immigrDate},
#'   \code{birthDate}, and \code{immigrInitState}.
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
#' @param fixInitStates Indices of substates/attributes that newborns inherit
#'   directly from their mother. If empty, no attributes are inherited.
#'
#' @param varInitStates A matrix or vector of substates/attributes that are
#'   assigned to newborns probabilistically according to \code{initStatesProb}.
#'
#' @param initStatesProb A vector of probabilities corresponding to
#'   \code{varInitStates}. If inheritance is used, probabilities must sum to one
#'   *conditional on the inherited attributes*; otherwise they must sum to one
#'   overall.
#'
#' @param maxAge Maximum age (in years) an individual can reach during the
#'   simulation. Must be strictly positive.
#'
#' @param simHorizon A length-two vector giving the start and end dates of the
#'   simulation in \code{"yyyymmdd"} format. The start date must precede the end date.
#'
#' @param fertTr A character vector specifying transitions that trigger a birth
#'   event. Each element must be of the form \code{"stateA -> stateB"}.
#'
#' @param partTr A character vector specifying transitions that represent the
#'   onset of a partnership. Each element must be of the form
#'   \code{"stateA -> stateB"}. Used to link two individuals as partners.
#'
#' @param rule Integer specifying the partner selection rule:
#'   \itemize{
#'     \item \code{0}: random partner selection among eligible individuals aged ≥18
#'     \item \code{1}: partner selection based on an age-difference distribution
#'       provided via \code{ageDiffDistr}
#'     \item \code{2}: reserved for future extensions (currently defaults to random)
#'   }
#'
#' @param ageDiffDistr A function returning probabilities for selecting a partner
#'   based on age differences. Used only when \code{rule = 1}. The function must
#'   accept a vector of age differences (male age minus female age) and return a
#'   vector of probabilities of equal length.
#'
#' @param sepTr A character vector specifying transitions that dissolve a
#'   partnership (separation). Each element must be of the form \code{"stateA -> stateB"}.
#'
#' @param probSepTr A numeric vector of probabilities corresponding to
#'   \code{sepTr}. If multiple separation transitions are possible and no valid
#'   probabilities are supplied, equal probabilities are assumed.
#'
#' @param absPartTr A character vector specifying transitions triggered when a
#'   partner enters an absorbing state (e.g., widowhood transitions). Each
#'   element must be of the form \code{"stateA -> stateB"}.
#'
#' @param monthSchoolEnrol Numeric month (1–12) indicating the general school
#'   enrollment month for elementary school. If omitted and no school enrollment
#'   transition is defined, school enrollment is not simulated.
#'
#' @param genMatrix A matrix defining which family or partnership relations to record.
#'   It must include four columns: transition rate function to be considered, relation type
#'   (e.g., "motherID", "fatherID", "partnerID"), alpha,
#'   and beta.
#'
#' @details
#' The simulation proceeds in continuous time using integrated hazard rates. For
#' each individual, the next event time is drawn based on the transition rates
#' defined in \code{transitionMatrix}. Linked lives are handled explicitly:
#' partnership formation links two individuals, separation breaks the link, and
#' transitions of one partner (e.g., death) may trigger transitions in the
#' other. Fertility events create new individuals whose attributes may be
#' inherited or assigned probabilistically. Immigration introduces individuals
#' after the simulation start. Duration-dependent transition rates are supported.
#'
#' @returns
#' A data frame containing the full simulated population, including all
#' individuals ever present, their transition histories, parental links,
#' partnership links, birth and death information, and all simulated events.
#'
#' @examples
#'
#' \dontrun{
#' # Clean workspace
#' rm(list=ls())
#'
#' # Defining simulation horizon
#' startDate  <- 20140101
#' endDate    <- 20401231
#' simHorizon <- c(startDate=startDate, endDate=endDate)
#'
#' # Seed for random number generator
#' set.seed(123)
#'
#' #Definition of maximal age
#' maxAge <- 100
#'
#' # Definition of non-absorbing and absorbing states
#' sex    <- c("m", "f")
#' #PH = Parental Home, PA = Partnership and A = Alone
#' livArr <- c("PH", "PA", "A")
#' stateSpace <- expand.grid(sex = sex, livArr = livArr)
#' absStates  <- c("dead")
#'
#' # Definition of an initial population
#' N <- 300
#' birthDates <- runif(N, min = getInDays(19500101), max = getInDays(20131231))
#'
#' getInit <- function(bd){
#'   age <- trunc((getInDays(simHorizon[1]) - bd)/365.25)
#'   s1  <- sample(sex, 1)
#'   s2  <- ifelse(age <= 18, "PH", sample(livArr, 1))
#'   paste(s1, s2, sep = "/")
#' }
#' initPop <- data.frame(
#'   ID        = 1:N,
#'   birthDate = getInDateFormat(birthDates),
#'   initState = sapply(birthDates, getInit)
#' )
#' # Definition of mortality rates
#' mortRates <- function(age, calTime){
#'   a <- 0.00003
#'   b <- 0.1
#'   a * exp(b * age)
#' }
#'
#' moveOut <- function(age, calTime){
#'   ifelse(age > 18, 0.02, 0)
#' }
#'
#' movePH <- function(age, calTime){
#'   ifelse(age < 35, 0.01, 0.002)
#' }
#'
#' startPA <- function(age, calTime){
#'   r <- dnorm(age, mean = 30, sd = 8)
#'   r[age <= 18] <- 0
#'   r
#' }
#'
#' stopPA_A  <- function(age, calTime) dnorm(age, mean = 45, sd = 10)
#' stopPA_PH <- function(age, calTime) dnorm(age, mean = 30, sd = 8) / 5
#'
#' ageDiffDistr <- function(ageDiff){
#'   dnorm(ageDiff, sd = 5)
#' }
#'
#' baseTransitions <- cbind(
#'   c(
#'     "m/PH->m/A", "m/A->m/PH",
#'     "f/PH->f/A", "f/A->f/PH"
#'   ),
#'   c(
#'     "moveOut", "movePH",
#'     "moveOut", "movePH"
#'   )
#' )
#'
#' partnerTransitions <- cbind(
#'   c(
#'     "f/PH->f/PA", "f/A->f/PA",
#'     "f/PA->f/A", "f/PA->f/PH",
#'     "m/PA->m/A", "m/PA->m/PH"
#'   ),
#'   c(
#'     "startPA", "startPA",
#'     "stopPA_A", "stopPA_PH",
#'     "stopPA_A", "stopPA_PH"
#'   )
#' )
#'
#' # Transition pattern and assignment of functions specifying transition rates
#' allTransitions <- rbind(baseTransitions, partnerTransitions)
#'
#' absTransitions <- cbind(
#'   c("f/dead", "m/dead"),
#'   c("mortRates", "mortRates")
#' )
#'
#' transitionMatrix <- buildTransitionMatrix(
#'   allTransitions = allTransitions,
#'   absTransitions = absTransitions,
#'   stateSpace     = stateSpace
#' )
#'
#' partTr    <- c("PH->PA", "A->PA")
#' sepTr     <- c("PA->A", "PA->PH")
#' absPartTr <- c("dead -> A")
#' probSepTr <- c(0.7, 0.3)
#'
#' # Execute microsimulation
#' poplink <- micSimLink(initPop=initPop, immigrPop=NULL, transitionMatrix=transitionMatrix,
#'                       absStates=absStates, varInitStates=NULL, initStatesProb=NULL,
#'                       maxAge=maxAge, simHorizon=simHorizon, fertTr=c(),
#'                       partTr=partTr, rule=0, ageDiffDistr=ageDiffDistr,
#'                       sepTr=sepTr, probSepTr=probSepTr, absPartTr=c())
#' }
#' @export
#'
  micSimLink <- function(initPop, immigrPop=NULL, transitionMatrix, absStates=NULL, fixInitStates = c(),
                         varInitStates=c(), initStatesProb=c(), maxAge=99, simHorizon,
                         fertTr=c(), partTr=c(), rule=0, ageDiffDistr=NULL, sepTr=c(),
                         probSepTr = c(), absPartTr = c(), monthSchoolEnrol=c(), genMatrix = c()) {

  # --------------------------------------------------------------------------------------------------------------------
  # --------------------------------------------------------------------------------------------------------------------
  # A. CHECK INPUT FOR CONSISTENCY
  # --------------------------------------------------------------------------------------------------------------------
  # --------------------------------------------------------------------------------------------------------------------
  if(is.null(initPop))
    stop('No starting population has been defined.')
  if(!is.null(initPop)){
    if(paste(colnames(initPop),collapse='/')!='ID/birthDate/initState')
      stop('Matrix specifying the starting population has not been defined properly.')
  }
  if(!is.null(immigrPop)){
    if(paste(colnames(immigrPop),collapse='/')!='ID/immigrDate/birthDate/immigrInitState')
      stop('Matrix specifying immigrants has not been defined properly.')
  }
  if(is.null(transitionMatrix))
    stop('Matrix defining transition pattern und functions has not been defined properly.')
  if(maxAge<=0)
    stop('The maximal age until which individual life courses are simulated should exceed zero.')
  if(length(simHorizon)!=2)
    stop('The simulation horizon has not been defined properly.')
  if(is.null(absStates))
    absStates <- setdiff(colnames(transitionMatrix),rownames(transitionMatrix))
  if(length(fertTr)>0){
    if((is.null(varInitStates) & is.null(initStatesProb)))
      stop('For children potentially born during simulation no inital state(s) and/or corresponding occurrence probabilities have been defined.')
    if(length(fixInitStates)>0){
      for(i in 1:length(fixInitStates)){
        ssum <- sum(initStatesProb[apply(varInitStates,1, function(rr){varInitStates[,fixInitStates[i]][1] %in% rr})])
        if(ssum!=1)
          stop('The sum of the probabilities to assign initial states to newborns must equal 1.')
      }
    } else {
      if(sum(initStatesProb)!=1)
        stop('The sum of the probabilities to assign initial states to newborns must equal 1.')
    }
  }
  # check whether all functions delivering transition rates deliver vectors of rates as output (necessary for integration procedure later on)
  allTr <- unique(as.vector(transitionMatrix)[as.vector(transitionMatrix) !="0"])
  simStartInDays <- getInDays(simHorizon[1])
  simStopInDays  <- getInDays(simHorizon[2])
  ranYear <- c(getYear(simStartInDays), getYear(simStopInDays))
  if(length(fertTr)>0){
    minAge <- 0
  } else {
    minAge <- min(trunc(getAgeInDays(simHorizon[1], initPop$birthDate)/365.25))
  }
  ranAge <- c(minAge,maxAge)
  ran <- min(c(diff(ranYear), diff(ranAge)))
  depMatrix <- rate_cS(allTr)
  for (tr in 1:length(depMatrix)) {
    if(ncol(depMatrix) >= 3 && any(depMatrix[,3] == 1)){
      tr_dur <- rownames(depMatrix)[depMatrix[,3] != 0]
      if(length(tr_dur) != 0){
        for(tr in 1:length(tr_dur)){ # check whether all functions delivering transition rates deliver vectors of rates as output (necessary for integration procedure later on)
          for(cal in c(getYear(simStartInDays): getYear(simStopInDays))){
            for(age in minAge:(maxAge-1)){
              for(dur in 0:ran){
                res <- eval(do.call(tr_dur[tr], args=list(age=age,calTime=cal,duration=dur)))
                if(anyNA(res)){
                  cat("The rates function for ", allTr[i], " does not deliver a vector of rates for an input vector of age, calendar time, and/or duration (all in years).\n")
                  cat("The missing rate occurs at year ",cal, " for age ", age, " and duration ", dur, "\n.")
                  cat("This is a requirement for the simulation procedure to run since it is based on integrated hazard rates.\n")
                  stop('Incorrect definition of input rates function!')
                }
              }
            }
          }
        }
      }
      else{
        tr_nd <- rownames(depMatrix)[depMatrix[,3] == 0]
        for(tr in 1:length(tr_nd)){ # check whether all functions delivering transition rates deliver vectors of rates as output (necessary for integration procedure later on)
          for(cal in c(getYear(simStartInDays): getYear(simStopInDays))){
            for(age in minAge:(maxAge-1)){
              res <- eval(do.call(tr_nd[tr], args=list(age=age,calTime=cal)))
              if(anyNA(res)){
                cat("The rates function for ", allTr[i], " does not deliver a vector of rates for an input vector of age and/ or calendar time (all in years).\n")
                cat("The missing rate occurs at year ",cal, " and for age ", age, "\n.")
                cat("This is a requirement for the simulation procedure to run since it is based on integrated hazard rates.\n")
                stop('Incorrect definition of input rates function!')
              }
            }
          }
        }
      }
    }
    if (ncol(depMatrix) == 2) {
      tr_nd <- rownames(depMatrix)
      for (tr in 1:length(tr_nd)) {
        for (cal in getYear(simStartInDays):getYear(simStopInDays)) {
          for (age in minAge:(maxAge - 1)) {
            res <- eval(do.call(tr_nd[tr], args = list(age = age, calTime = cal)))
            if (anyNA(res)) {
              cat("The rates function for ", tr_nd[tr], " does not deliver a vector of rates for an input vector of age and calendar time.\n")
              cat("The missing rate occurs at year ", cal, " for age ", age, "\n.")
              cat("This is a requirement for the simulation procedure to run since it is based on integrated hazard rates.\n")
              stop('Incorrect definition of input rates function!')
            }
          }
        }
      }
    }
  }
  if(length(monthSchoolEnrol)==0){
    schoolEnrol <- FALSE
  } else {
    schoolEnrol <- TRUE
  }

  if(length(sepTr)>1){
    if(length(probSepTr)!=length(sepTr)) {
      cat("There is more than one option to trigger a separation event for couples (linked lives), but no (correctly specified) related probabilites (`probSepTr') to decide which to choose\n.")
      cat("Therefore during simulation a coin toss (pure random draw) decides which option to choose from.")
      probSepTr <- rep(1/length(sepTr), length(sepTr))
    }
  }

  if(length(absPartTr)>1){
    if(length(partTr)==0){
      cat("You try to model dissoluation of partnerships without having defined a transition causing a partnership, i.e. the object `partTr' is missing. \n")
      stop("Miss partnership indicator.")
    }
  }

  # --------------------------------------------------------------------------------------------------------------------
  # --------------------------------------------------------------------------------------------------------------------
  # B. DEFINITION OF GLOBAL PARAMETERS
  # --------------------------------------------------------------------------------------------------------------------
  # --------------------------------------------------------------------------------------------------------------------

  # Simulation horizon
  simStartInDays <- getInDays(simHorizon[1])
  simStopInDays  <- getInDays(simHorizon[2])
  # Assign to each string state and substate an unique numerical code
  codingScheme <- builtStatesCodes(transitionMatrix)
  nSubStates <- ceiling(max(floor(log10(abs(codingScheme[,2]))) + 1)/2) # number of substates (without absorbing states)
  # Put numerical codes to transition matrix as well
  transitionMatrixNum <- transitionMatrix
  indCodesR <- match(rownames(transitionMatrix), codingScheme[,1])
  rownames(transitionMatrixNum) <- as.numeric(codingScheme[indCodesR,2])
  indCodesC <- match(colnames(transitionMatrix), codingScheme[,1])
  colnames(transitionMatrixNum) <- as.numeric(codingScheme[indCodesC,2])
  absStatesNum <- -c(1:length(absStates))
  if(!is.null(varInitStates)){
    varInitStatesStr <- apply(varInitStates, 1, paste, collapse="/")
    varInitStatesNum <- codingScheme[codingScheme[,1] %in% varInitStatesStr, c(2,2+c(1:nSubStates))]
  }
  # Event queue
  queue <- matrix(NA,ncol=11,nrow=0) # columns: 'ID','currTime','currState','currAge','nextState','timeToNextState', 'birthtime', 'initState', 'motherID', 'fatherID', 'partnerID'
  # Global time
  t.clock <- simStartInDays  # counts in days since 01-01-1970
  # Recording transitions performed
  transitions <- matrix(NA,ncol=8,nrow=0) # columns: ID, From, To, transitionTime, transitionAge, motherID, fatherID, partnerID
  # Record linkage of mothers to newborns via their IDs
  if(length(fertTr)>0){
    mothers <- matrix(NA,ncol=2,nrow=0) # 'motherID', 'childID' #this could become unneccessary with relMatrix?
  }
  #Recording of Relationships
  if(length(genMatrix)>0){
    relMatrix <- matrix(NA, ncol= 8, nrow = 0) #ID, birthtime,  motherID, motherbirth, fatherID, fatherbirth, partnerID, partnerbirth #where do I get births -> initpop (yikes)
    mapgenMatrix_ID <- c("ID" = 1, "motherID" = 3, "fatherID" = 5, "partnerID" = 7) #refers to columns of relMatrix
    # Function to recode tr in genMatrix
    mapgenMatrix_R <- function(genMatrix, depMatrix) {
      # Extract the tr values from genMatrix
      tr <- genMatrix[, "tr"]
      crows <- match(tr, rownames(depMatrix))
      gen <- depMatrix[crows, "genarg"]
      genMatrix[, "tr"] <- gen

      return(genMatrix)
    }

    genMatrix[, 2] <- mapgenMatrix_ID[match(genMatrix[, 2], names(mapgenMatrix_ID))] #TODO: add to function above
    genMatrix <- mapgenMatrix_R(genMatrix, depMatrix) #basically build genMatrix expanded

    genMatrix <- matrix(as.numeric(genMatrix), ncol = 4) #TODO: ideally all of this goes in the function above

  }

  # Maximal Id / counter for individuals in the simulation
  maxId <- 0

  # ----------------------------------------------------------------------------------------------------------------------
  # ----------------------------------------------------------------------------------------------------------------------
  # C. FUNCTIONS REQUIRED FOR SIMULATION
  # ----------------------------------------------------------------------------------------------------------------------
  # ----------------------------------------------------------------------------------------------------------------------

  # Function building matrix indicating the transitions between states causing a newborn
  buildFertTrExpanded <- function(){
    allStates <- matrix(rownames(transitionMatrix), ncol=1)
    allStatesSplit <- apply(allStates,1,strsplit, split="/")
    fert <- do.call(rbind,(strsplit(fertTr,split='->')))
    fertTrExpandedStr <- NULL
    for(i in 1:nrow(allStates)){
      cS <- allStatesSplit[[i]][[1]]
      for(j in 1:nrow(allStates)){
        dS <- allStatesSplit[[j]][[1]]
        if(("f" %in% cS) & ("f" %in% dS)){
          for(k in 1:nrow(fert)){
            ff <- fert[k,]
            oS <- strsplit(ff[1],'/')[[1]]
            bS <- strsplit(ff[2],'/')[[1]]
            cond1 <- !(F %in% (oS %in% cS)) & !(F %in% (bS %in% dS))
            cond2 <- paste((cS[!(cS %in% oS)]),collapse="/") == paste((dS[!(dS %in% bS)]),collapse="/") # if there are a fertility event only one substate can change, namely that one belonging to the fertility attribute
            if(cond1 & cond2){
              fertTrExpandedStr <- rbind(fertTrExpandedStr, c(paste0(cS,collapse="/"), paste0(dS,collapse="/")))
            }
          }
        }
      }
    }
    indCodes1 <- match(fertTrExpandedStr[,1], codingScheme[,1]) # transform numerical state codes to strings according to codingScheme
    indCodes2 <- match(fertTrExpandedStr[,2], codingScheme[,1])
    fertTrExpanded <- cbind(codingScheme[indCodes1,2],codingScheme[indCodes2,2])
    return(fertTrExpanded)
  }

  # Function checks whether a transition causes a newborn.
  # (Demands `fertTrExpanded': matrix indicating the transitions between states causing a newborn (defined by `fertTr').)
  isBirthEvent <- function(currState, destState){
    oS <- which(fertTrExpanded[,1] %in% currState)
    dS <- which(fertTrExpanded[,2] %in% destState)
    if(length(intersect(oS,dS))>0)
      return(TRUE)
    return(FALSE)
  }

  # Function adds to simulation population a newborn (using combined sim. step with duration dep. trans. functions)
  addNewNewborn_dur <- function(birthTime=birthTime, motherID=motherID, motherState=motherState, motherbirth=motherbirth, fatherID=fatherID){
    if(length(fixInitStates)>0){
      motherStatePart <-codingScheme[codingScheme[,2] %in% motherState, 2+c(1:nSubStates)][fixInitStates]
      if(length(fixInitStates)>1){
       notSuitedStatesInd <- which(is.na(apply(varInitStatesNum[,-1,drop=F][,fixInitStates,drop=F],1,match, x=motherStatePart)), arr.ind=TRUE)[,2]
      } else{
       notSuitedStatesInd <- which(is.na(apply(varInitStatesNum[,-1,drop=F][,fixInitStates,drop=F],1,match, x=motherStatePart)), arr.ind=TRUE)
      }
      inSt <- setdiff(1:nrow(varInitStatesNum), notSuitedStatesInd)
      varInitStatesR <- varInitStatesNum[inSt,,drop=F]
      initStatesProbR <- initStatesProb[inSt]
      birthState <- varInitStatesR[sample(1:nrow(varInitStatesR),size=1,replace=T,prob=initStatesProbR),1]
    } else {
      birthState <- varInitStatesNum[sample(1:nrow(varInitStatesNum),size=1,replace=T,prob=initStatesProb),1]
    }
    maxId <<- maxId + 1
    newInd <- c(maxId,getInDateFormat(birthTime),codingScheme[codingScheme[,2] %in% birthState,1], motherID, fatherID, NA)
    #cat('NewBorn: ',newInd,'\n')
    initPop <<- rbind(initPop,newInd)
    mothers <<- rbind(mothers, c(motherID, maxId))

    if(length(genMatrix)>0){
    relMatrix <<- rbind(relMatrix, c(maxId, birthTime, motherID, motherbirth, fatherID, NA, NA, NA)) #TODO: add fatherbirthtime + partnership (likely in partnership creation!)
    #cat('RelationMatrix:', tail(relMatrix, n = 3))
    }

    nE <- getNextStep(c(maxId,birthState,0,birthTime,birthTime,birthState,motherID, fatherID, NA))

    if(length(partTr) > 0){
    matePool <- rbind(matePool, c(maxId, birthState, birthTime, birthState, motherID, fatherID, birthTime+(18*365.25))) # can become partner when turning 18
    }
    #cat('\n------------n')
  }

  # Function building matrix indicating the transitions between states causing a school enrollment
  buildEduTrExpanded <- function(){
    allStates <- matrix(rownames(transitionMatrix), ncol=1)
    eduTrExpandedStr <- NULL
    for(i in 1:nrow(allStates)){
      for(j in 1:nrow(allStates)){
        cond <- all( c(any('no' %in% strsplit(allStates[i],'/')[[1]]), any('low' %in% strsplit(allStates[j],'/')[[1]])))
        if(cond){
          eduTrExpandedStr <- rbind(eduTrExpandedStr,c(allStates[i], allStates[j]))
        }
      }
    }
    eduInd <- gsub(x=eduTrExpandedStr[,1], pattern="no", replacement="")==gsub(x=eduTrExpandedStr[,2], pattern="low", replacement="")
    eduTrExpandedStr <- eduTrExpandedStr[eduInd,]
    indCodes1 <- match(eduTrExpandedStr[,1], codingScheme[,1]) # transform numerical state codes to strings according to codingScheme
    indCodes2 <- match(eduTrExpandedStr[,2], codingScheme[,1])
    eduTrExpanded <- cbind(codingScheme[indCodes1,2],codingScheme[indCodes2,2])
    return(eduTrExpanded)
  }

  # Function checks whether a transition implies a school enrollment (in the year when child turns seven).
  # (If state 1 comprises value `no' and state 2 comprises value `low', the transition is marked as `school enrollment'.)
  isSchoolEnrolment <- function(currState,destState){
    oS <- which(eduTrExpanded[,1] %in% currState)
    dS <- which(eduTrExpanded[,2] %in% destState)
    if(length(intersect(oS,dS))>0)
      return(TRUE)
    return(FALSE)
  }

  # Function building matrix indicating the transitions between states causing a partnership event
  buildPartTrExpanded <- function(){
    allStates <- matrix(rownames(transitionMatrix), ncol=1)
    allStatesSplit <- apply(allStates,1,strsplit, split="/")
    part <- do.call(rbind,(strsplit(partTr,split='->')))
    partTrExpandedStr <- NULL
    for(i in 1:nrow(allStates)){
      cS <- allStatesSplit[[i]][[1]]
      for(j in 1:nrow(allStates)){
        dS <- allStatesSplit[[j]][[1]]
        if((("f" %in% cS) & ("f" %in% dS)) | (("m" %in% cS) & ("m" %in% dS))){
          for(k in 1:nrow(part)){
            ff <- part[k,]
            oS <- strsplit(ff[1],'/')[[1]]
            bS <- strsplit(ff[2],'/')[[1]]
            cond1 <- !(F %in% (oS %in% cS)) & !(F %in% (bS %in% dS))
            cond2 <- paste((cS[!(cS %in% oS)]),collapse="/") == paste((dS[!(dS %in% bS)]),collapse="/") # if there are a partnership event only one substate can change, namely that one belonging to the partnership attribute
            if(cond1 & cond2){
              partTrExpandedStr <- rbind(partTrExpandedStr, c(paste0(cS,collapse="/"), paste0(dS,collapse="/")))
            }
          }
        }
      }
    }
    indCodes1 <- match(partTrExpandedStr[,1], codingScheme[,1]) # transform numerical state codes to strings according to codingScheme
    indCodes2 <- match(partTrExpandedStr[,2], codingScheme[,1])
    partTrExpanded <- cbind(codingScheme[indCodes1,2],codingScheme[indCodes2,2])
    return(partTrExpanded)
  }

  # Function checks whether a transition causes a partnership.
  # (Demands `partTrExpanded': matrix indicating the transitions between states causing a partnership (defined by `partTr').)
  isPartEvent <- function(currState, destState){
    oS <- which(partTrExpanded[,1] %in% currState)
    dS <- which(partTrExpanded[,2] %in% destState)
    if(length(intersect(oS,dS))>0)
      return(TRUE)
    return(FALSE)
  }

  # Decide whether a match between potential partners will happen (based on Bernoulli experiment)
  isPartnMatch <- function(pr){
    if(length(pr)==0) return(0) # nobody available
    if(length(pr)==1) return(1) # only one option, return this one guy
    candidates <- rep(0,length(pr)) # pool of candidates
    for(k in 1:length(pr)){
      candidates[k] <- rbinom(n=1,size=1,prob=pr[k])
    }
    candIndex <- which(candidates %in% 1, arr.ind=TRUE)
    if(length(candIndex)>0) return(sample(x=candIndex,size=1))
    return(which(pr %in% max(pr), arr.ind = TRUE)) # if no match can be found stochastically, take the one with the highest probability
  }

  # If an event means the onset of a partnership: a partner has to be searched for and linked
  createPartLink <- function(partTime, egoID, egoState, egoBirth, egoMothID, egoFathID, rule=0){
    #partTime=t.clock; egoID=indS[1]; egoState=indS[5]; egoBirth=indS[7]; egoMothID=indS[9]; egoFathID=indS[10]; rule=rule

    # if female, take males
    if(egoState %in% femStatesNum){
      matePool_egoM <- matePool[matePool[,7] <= partTime & matePool[,2] %in% maleStatesNum,,drop=F]
      # mate pool: id, currentState, birthday in days, initState, motherID, fatherID, time when available
    }
    if(egoState %in% maleStatesNum){
      matePool_egoM <- matePool[matePool[,7] <= partTime & matePool[,2] %in% femStatesNum,,drop=F]
      # mate pool: id, currentState, birthday in days, initState, motherID, fatherID, time when available
    }
    # no mate available
    if(nrow(matePool_egoM)==0){
      return(NA)
    }
    # incest is not allowed
    if(nrow(matePool_egoM)>=1){
      matePool_egoM <- matePool_egoM[!(matePool_egoM[,1] %in% egoMothID),,drop=F] # sort out mother
      matePool_egoM <- matePool_egoM[!(matePool_egoM[,1] %in% egoFathID),,drop=F] # sort out father
      if(!is.na(egoMothID)) # if mother is known (by ID)
        matePool_egoM <- matePool_egoM[!(matePool_egoM[,5] %in% egoMothID),,drop=F] # sort out siblings (siblings with same mother)
      if(!is.na(egoFathID)) # if father is known (by ID)
        matePool_egoM <- matePool_egoM[!(matePool_egoM[,6] %in% egoFathID),,drop=F] # sort out siblings (siblings with same father)
    }
    # partner choice rule
    if(rule==0) { # random draw from pool of opposite sex individuals older than 18 years
      matePool_egoM <- matePool_egoM[c(partTime-matePool_egoM[,3])/365.25>=18,,drop=F]
    }
    if(rule==1){ # select partner who is currently older than 18 years and then according to age difference (ageMale-ageFem) distribution provided
      #don't implement go directly to GEN rule
      #cat("-----------------------------")
      #matePool_egoM <- matePool_egoM[c(partTime-matePool_egoM[,3])/365.25>=18,,drop=F]
      if(egoState %in% femStatesNum) {
        #cat("Fem search with fem birthtime: ", egoBirth,"\n")
        #cat("Size pool: ", nrow(matePool_egoM), "\n")
        #cat("Male Mate Pool with birthtime: ", matePool_egoM[,3],"\n")
        ageDiff <- (egoBirth-matePool_egoM[,3])/365.25 # ageMale - ageFem TODO Check sequence!
      } else {
        #cat("Male search with male birthtime: ", egoBirth,"\n")
        #cat("Size pool: ", nrow(matePool_egoM), "\n")
        #cat("Fem  Mate Pool with birthtime: ", matePool_egoM[,3])
        ageDiff <- (matePool_egoM[,3]-egoBirth)/365.25
      }
      matePool_egoM <- matePool_egoM[abs(ageDiff)<=15,,drop=F]
      ageDiff <- ageDiff[abs(ageDiff)<=15]
      #cat("AgeDiff: ", ageDiff, "\n")
      prr <- ageDiffDistr(ageDiff)
      #cat("Probs: ",prr,"\n")
      matePool_egoM <- matePool_egoM[isPartnMatch(prr),,drop=F]
      #cat("Mate: ", matePool_egoM, "\n")
    }
    if(rule==2){ # TODO: not implemented yet; default random draw
      ## GENERIC NEW RULE: 2nd col of matePool gives currentState of potential mates & egoState is input argument
      #logit?
      matePool_egoM <- matePool_egoM[c(partTime-matePool_egoM[,3])/365.25>=18,,drop=F]
    }
    # if no partner is available return NA
    if(nrow(matePool_egoM)==0){
      return(NA)
    }
    # if only one mate available, this will be the partner
    if(nrow(matePool_egoM)==1){
      matePool_part <- matePool_egoM
    }
    # if more than potential mate is available, select a partner at random
    if(nrow(matePool_egoM)>1){
      matePool_part <- matePool_egoM[sample(size=1, x=1:nrow(matePool_egoM)),]
    }
    # partner ID
    partID <- matePool_part[1]
    # remove newly formed partners from mate pool
    matePool <<- matePool[!(matePool[,1] %in% c(egoID, partID)),,drop=F]
    return(matePool_part)
  }

  # Function building matrix indicating the transitions between states causing a partnership event
  buildSepTrExpanded <- function(){
    allStates <- matrix(rownames(transitionMatrix), ncol=1)
    allStatesSplit <- apply(allStates,1,strsplit, split="/")
    sep <- do.call(rbind,(strsplit(sepTr,split='->')))
    sepTrExpandedStr <- NULL
    for(i in 1:nrow(allStates)){
      cS <- allStatesSplit[[i]][[1]]
      for(j in 1:nrow(allStates)){
        dS <- allStatesSplit[[j]][[1]]
        if((("f" %in% cS) & ("f" %in% dS)) | (("m" %in% cS) & ("m" %in% dS))){
          for(k in 1:nrow(sep)){
            ff <- sep[k,]
            oS <- strsplit(ff[1],'/')[[1]]
            bS <- strsplit(ff[2],'/')[[1]]
            cond1 <- !(F %in% (oS %in% cS)) & !(F %in% (bS %in% dS))
            cond2 <- paste((cS[!(cS %in% oS)]),collapse="/") == paste((dS[!(dS %in% bS)]),collapse="/") # if there are a partnership event only one substate can change, namely that one belonging to the partnership attribute
            if(cond1 & cond2){
              sepTrExpandedStr <- rbind(sepTrExpandedStr, c(paste0(cS,collapse="/"), paste0(dS,collapse="/")))
            }
          }
        }
      }
    }
    indCodes1 <- match(sepTrExpandedStr[,1], codingScheme[,1]) # transform numerical state codes to strings according to codingScheme
    indCodes2 <- match(sepTrExpandedStr[,2], codingScheme[,1])
    partTrExpanded <- cbind(codingScheme[indCodes1,2],codingScheme[indCodes2,2])
    return(partTrExpanded)
  }  #TODO: septr everything correct??

  # Function checks whether a transition causes a separation
  # (Demands `sepTrExpanded': matrix indicating the transitions between states causing a separation (defined by `sepTr').)
  isSepEvent <- function(currState, destState){
    oS <- which(sepTrExpanded[,1] %in% currState)
    dS <- which(sepTrExpanded[,2] %in% destState)
    if(length(intersect(oS,dS))>0)
      return(TRUE)
    return(FALSE)
  }

  # Function building matrix indicating an event due to the transition of the partner to an absorbent state
  buildAbsPartTrExpanded <- function(){
    allStates <- matrix(rownames(transitionMatrix), ncol=1)
    allStatesSplit <- apply(allStates,1,strsplit, split="/")
    wid <- do.call(rbind,(strsplit(absPartTr,split='->')))
    wid <- gsub(" ", "", unique(wid[,2]))
    part <- do.call(rbind,(strsplit(partTr,split='->')))
    part <- gsub(" ", "", unique(part[,2]))
    absPartTrExpandedStr <- NULL
    for(i in 1:length(part)){
      ll <- part[i]
      for(j in 1:length(wid)){
        ff <- wid[j]
        for(k in 1: length(allStatesSplit)){
          cS <- allStatesSplit[[k]][[1]]
          if(ll %in% cS){
            dS <- cS
            dS[dS %in% ll] <- ff
            absPartTrExpandedStr <- rbind(absPartTrExpandedStr, c(paste0(cS,collapse="/"), paste0(dS,collapse="/")))
          }
        }
      }
    }
    indCodes1 <- match(absPartTrExpandedStr[,1], codingScheme[,1]) # transform numerical state codes to strings according to codingScheme
    indCodes2 <- match(absPartTrExpandedStr[,2], codingScheme[,1])
    absPartTrExpanded <- cbind(codingScheme[indCodes1,2],codingScheme[indCodes2,2])
    absPartTrExpanded <- absPartTrExpanded[which(absPartTrExpanded[,1] %in% partTrExpanded [,2]),] # states of origin can only be states indicating being partnered
    return(absPartTrExpanded)
  }

  # ----------------------------------------------------------------------------------------------------------------------
  # ----------------------------------------------------------------------------------------------------------------------
  # D. SIMULATION STEP
  # ----------------------------------------------------------------------------------------------------------------------
  # ----------------------------------------------------------------------------------------------------------------------

#TODO: check that no ncol == n is inside of getnextstep anymore

  # Function to compute the next transition state and time of an individual who at age `currAge', at time `calTime'
  # entered its current state `currState'. At this, consider possible duration dependencies of transition rates.
  getNextStep <- function(inp, isIMInitEvent=F){

    # Extract input data
    id <- inp[1]
    currState <- inp[2] # current state in numerical code
    currAge <- inp[3] # age in days
    calTime <- inp[4] # calendar time in days since 01-01-1970
    birthTime <- inp[5] # birth time in days since 01-01-1970
    initState <- inp[6] # initial state at sim. start, NA for individuals not yet born at that time
    mothID <- inp[7]
    fathID <- inp[8]
    partID <- inp[9]

    #relMatrix input here?

    # First event of an immigrant: he/she enters the population later than sim. starting time
    if(isIMInitEvent) lagToWaitingTime <- (calTime - simStartInDays)/365.25 # in years
    if(!isIMInitEvent) lagToWaitingTime <- 0 # in years

    #cat('\n-----\nID: ',id,'\n')
    #print(inp)
    ageInYears <- currAge/365.25
    #cat('Age: ',ageInYears,' - CalTime: ',getYear(calTime),'-',getMonth(calTime),'-',getDay(calTime),'\n')
    # Possible destination states
    possTr <- transitionMatrixNum[match(currState, rownames(transitionMatrixNum)),]
    possTr <- possTr[which(possTr !=0)]

    nextEventMatrix <- matrix(0, ncol=2, nrow=length(possTr))

    # How many years (along age scale) remain until `maxAge'?
    ranMaxAge <- (maxAge-0.01)-ageInYears
    # How many years (along cal. time scale) remain until simulation end?
    ranMaxYear <-  (simStopInDays - calTime)/365.25
    # Identify the time range that should be considered.
    ran <- min(ranMaxYear,ranMaxAge)
    #cat('Ran: ',ran,' - ranMaxAge: ',ranMaxAge,' - ranMaxYear: ',ranMaxYear,'\n')
    #ranAge <- c(ageInYears,ageInYears+ranMaxAge) # age range in years
    #ranYear <- c(getYear(calTime), getYear(calTime)+ran) # year range in years
    #cat('RanAge: ',ranAge,' - ranYear: ',ranYear,'\n')

    #duration argument = yes, genarg argument = no
    if(any(depMatrix[,3] == 1) && all(depMatrix[,4] == 0)) {
      tr_dur <- rownames(depMatrix)[depMatrix[,3] == 1]
      # Extract transition history of individual until current cal. time.
      historiesInd <- transitions[transitions[,1] %in% id & transitions[,4] <= calTime,,drop=F]
      # Extract for each state the duration until transition (in days).
      # Here, we have to differ between states of which we do not know when they are entered (i.e., `initial states' of members
      # of the starting population and the states of migrants when they entered the country), and
      # the states we know the `entering date' as well as the `leaving date' (if the state has been left).
      if(birthTime < simStartInDays | isIMInitEvent) {
        dur <- rbind(c(initState,NA),historiesInd[,c(3,4),drop=F]) # extract from historiesInd 'To' and 'transitionTime'
        dur <- cbind(dur,c(diff(dur[,2]),0)) # columns: TransitionTo, AtTime, durUntil
        dur[which(is.na(dur[,2])),3] <- NA
      } else {  # Individual is born during simulation.
        dur <- rbind(c(initState,birthTime),historiesInd[,c(3,4),drop=F])
        dur <- cbind(dur,c(diff(dur[,2]),0)) # columns: TransitionTo, AtTime, durUntil
      }
      # Compute for each possible destination state a waiting time.
      for(i in 1:length(possTr)){
        tr <- possTr[i]
        destState <-  as.numeric(names(tr))
        cS <- codingScheme[codingScheme[,2] %in% currState, 2+c(1:nSubStates)]
        dS <- codingScheme[codingScheme[,2] %in% destState, 2+c(1:nSubStates)]
        # To determine the duration (time elapsed since last transition) that applies for the considered destination state,
        # we have to determine the duration since the last change in the covariate affected.
        # For example, to specify the time being married, we have to determine the duration since (last) marriage.
        covToCh <- which((cS==dS)==F)
        durSinceLastCovCh <- Inf  # For the transition to `dead' so far the time elapsed since the last transition does not play any role.
        if(length(covToCh)==1){
          covHist <- codingScheme[codingScheme[,2] %in% dur[,1], 2+c(1:nSubStates),drop=F][,covToCh]
          idd <- which(covHist==cS[covToCh])
          if(length(idd)>1){
            if(F %in% (diff(idd)==1)){
              y <- rev(idd)[c(-1,diff(rev(idd)))==-1]
              idd <- rev(y)[c(diff(rev(y)),1)==1]
            }
          }
          durSinceLastCovCh <- sum(dur[idd,3]) # If I do not know how long an individual already is in a state: This gives NA.
          if(is.na(durSinceLastCovCh))
            durSinceLastCovCh <- currAge # Then assume the individual is already for his/her whole life in the state.
        }
        if(length(covToCh)>1 & (!destState %in% absStatesNum)){
          cat('Recognized a possible transition implying a change of two or more covariates.',
              'Concerning the derivation of the time being elapsed since the last transition this feature is not yet implemented.',
              'Current State: ',currState,' -> Possible transition to ',destState,'\n')
        }
        tageInYears <- trunc(ageInYears)
        tCalTime <- trunc(1970.001+calTime/365.25)
        tdurSinceLastCovCh <- trunc(durSinceLastCovCh/365.25)
        if(tr %in% tr_dur) {
          indRateFctDET <- function(x){
            res <- eval(do.call(tr,
                                args=list(age=tageInYears+x,calTime=tCalTime+x,duration=tdurSinceLastCovCh+x)))
            return(res)
          }
          ranAccuracyInDays <- (0:(trunc(ran*365.25)+0.99))/365.25
          detE <- indRateFctDET(ranAccuracyInDays)
          daysToTrInYears <- (which(detE == Inf)[1] - 1)/365.25
          if (Inf %in% detE) {
            timeToNext <- daysToTrInYears
          } else {
            u <- -log(1-runif(1))
            #cat('It: ',i,'--u: ',u,'\n')
            # Extract individual transition rate (depending on age, calendar time, and time elapsed)
            indRateFct <- function(x){
              ageIn <- ageInYears+x
              calIn <- 1970.001+calTime/365.25+x
              durIn <- durSinceLastCovCh/365.25+x
              res <- eval(do.call(tr, args=list(age=ageIn,calTime= calIn,duration=durIn)))
              if(TRUE %in% (res<0))
                stop('I have found negative rate value/s for transition: ',tr,'\n
                 This is implausible. Please check this. Simulation has been stopped.\n')
              #cat('x: ',x,' -- res', res,'\n')
              #cat('\n---\n')
              return(res)
            }

            if(sum(indRateFct(0:ran))==0){ # Rate function contains only zeros.
              intHaz <- 0
            } else {
              # Integrated hazard at max. value
              intHaz <- try(integrate(indRateFct, lower=0, upper=ran)$value, silent=TRUE)
              if(inherits(intHaz, 'try-error')){
                intHaz <- integrate(indRateFct, lower=0, upper=ran, stop.on.error = FALSE, rel.tol = 0.01)$value
              }
            }
            # If transformed random variate exceeds max. value of integr. hazard, we will not find a finite random waiting time.
            if(u<=intHaz){
              invHazFct <- function(x){
                #cat('x: ',x,'\n')
                try.res <- try(integrate(indRateFct, lower=0, upper=x)$value-u, silent=TRUE)
                #print(try.res)
                if(inherits(try.res, 'try-error')){
                  #cat('Seemingly, divergent intergral for ID ',id,
                  # ' in state ',currState,' at age ',currAge,' at time ',calTime, ' to state ',destState,
                  #  ' for random number: ',u,'\n')
                  try.res <- integrate(indRateFct, lower=0, upper=x, stop.on.error = FALSE, rel.tol = 0.01)$value-u
                }
                #cat('res: ',try.res,'\n-----\n')
                return(try.res)
              }
              # Find random waiting time.
              timeToNext <- uniroot(invHazFct,interval=c(0,ran))$root
            } else {
              timeToNext <- Inf
            }
          } }else {
            indRateFctDET <- function(x){
              res <- eval(do.call(tr,
                                  args=list(age=tageInYears+x,calTime=tCalTime+x)))
              return(res)
            }
            ranAccuracyInDays <- (0:(trunc(ran*365.25)+0.99))/365.25
            detE <- indRateFctDET(ranAccuracyInDays)
            daysToTrInYears <- (which(detE == Inf)[1] - 1)/365.25
            if (Inf %in% detE) {
              timeToNext <- daysToTrInYears
            } else {
              u <- -log(1-runif(1))
              #cat('It: ',i,'--u: ',u,'\n')
              # Extract individual transition rate (depending on age and calendar time)
              indRateFct <- function(x){
                ageIn <- ageInYears+x
                calIn <- 1970.001+calTime/365.25+x
                res <- eval(do.call(tr, args=list(age=ageIn,calTime= calIn)))
                if(TRUE %in% (res<0))
                  stop('I have found negative rate value/s for transition: ',tr,'\n
                 This is implausible. Please check this. Simulation has been stopped.\n')
                #cat('x: ',x,' -- res', res,'\n')
                #cat('\n---\n')
                return(res)
              }
              if(sum(indRateFct(0:ran))==0){ # Rate function contains only zeros.
                intHaz <- 0
              } else {
                # Integrated hazard at max. value
                intHaz <- try(integrate(indRateFct, lower=0, upper=ran)$value, silent=TRUE)
                if(inherits(intHaz, 'try-error')){
                  intHaz <- integrate(indRateFct, lower=0, upper=ran, stop.on.error = FALSE, rel.tol = 0.01)$value
                }
              }
              # If transformed random variate exceeds max. value of integr. hazard, we will not find a finite random waiting time.
              if(u<=intHaz){
                invHazFct <- function(x){
                  #cat('x: ',x,'\n')
                  try.res <- try(integrate(indRateFct, lower=0, upper=x)$value-u, silent=TRUE)
                  #print(try.res)
                  if(inherits(try.res, 'try-error')){
                    #cat('Seemingly, divergent intergral for ID ',id,
                    # ' in state ',currState,' at age ',currAge,' at time ',calTime, ' to state ',destState,
                    #  ' for random number: ',u,'\n')
                    try.res <- integrate(indRateFct, lower=0, upper=x, stop.on.error = FALSE, rel.tol = 0.01)$value-u
                  }
                  #cat('res: ',try.res,'\n-----\n')
                  return(try.res)
                }
                # Find random waiting time.
                timeToNext <- uniroot(invHazFct,interval=c(0,ran))$root
              } else {
                timeToNext <- Inf
              }
            }
          }
        nextEventMatrix[i,1] <- destState
        nextEventMatrix[i,2] <- (timeToNext+lagToWaitingTime)*365.25    # time to next event in days
      }
      #print(nextEventMatrix)
      nE <- nextEventMatrix[which(nextEventMatrix[,2]==min(as.numeric(nextEventMatrix[,2]))),,drop=F]
      if(dim(nE)[1]>1)
        nE <- nE[1,,drop=F]
      if(nE[1,2]!=Inf){
        # Cal. time of next event of individual. (If there is one.)
        tt <- calTime + as.numeric(nE[1,2])
        #print(tt)
        #cat(nE[1,1],'---',tt,'\n')
        # Check whether next event implies school enrollment. If yes, adjust transition time to ensure that the individual
        # enters school at Sept. 1 in the year he/she turns seven.
        if(schoolEnrol){ # Is school enrollment considered in this simulation model? Yes, then continue; otherwise skip this part
          if(isSchoolEnrolment(currState,nE[1,1])){
            enYear <- getYear(tt)
            if(getMonth(tt) <= monthSchoolEnrol) {
              enDate <- getInDays_my(enYear, monthSchoolEnrol)
            } else {
              enYear <- enYear+1
              enDate <- getInDays_my(enYear, monthSchoolEnrol)
            }
            diffToEn <- as.numeric(enDate-tt)
            nE[1,2] <- as.numeric(nE[1,2]) + diffToEn
          }
        }
        # Enqueue new event (if there is one).
        queue <<- rbind(queue, c(id, t.clock, currState, currAge - lagToWaitingTime*365.25, nE[1,1], nE[1,2], birthTime, initState, mothID, fathID, partID))
      }
      #cat('\n----------\n')
      return(nE)
    }

    #duration argument = no, genarg argument = no
    if(all(depMatrix[,3] == 0) && all(depMatrix[,4] == 0)){
      # Compute for each possible destination state a waiting time.
      for(i in 1:length(possTr)){
        tr <- possTr[i]
        destState <-  as.numeric(names(tr))
        cS <- codingScheme[codingScheme[,2] %in% currState, 2+c(1:nSubStates)]
        dS <- codingScheme[codingScheme[,2] %in% destState, 2+c(1:nSubStates)]

        tageInYears <- trunc(ageInYears)
        tCalTime <- trunc(1970.001+calTime/365.25)

        indRateFctDET <- function(x){
          res <- eval(do.call(tr,
                              args=list(age=tageInYears+x,calTime=tCalTime+x)))
          return(res)
        }
        ranAccuracyInDays <- (0:(trunc(ran*365.25)+0.99))/365.25
        detE <- indRateFctDET(ranAccuracyInDays)
        daysToTrInYears <- (which(detE == Inf)[1] - 1)/365.25
        if (Inf %in% detE) {
          timeToNext <- daysToTrInYears
        } else {
          u <- -log(1-runif(1))
          #cat('It: ',i,'--u: ',u,'\n')
          # Extract individual transition rate (depending on age and calendar time)
          indRateFct <- function(x){
            ageIn <- ageInYears+x
            calIn <- 1970.001+calTime/365.25+x
            res <- eval(do.call(tr, args=list(age=ageIn,calTime= calIn)))
            if(TRUE %in% (res<0))
              stop('I have found negative rate value/s for transition: ',tr,'\n
                   This is implausible. Please check this. Simulation has been stopped.\n')
            #cat('x: ',x,' -- res', res,'\n')
            #cat('\n---\n')
            return(res)
          }
          if(sum(indRateFct(0:ran))==0){ # Rate function contains only zeros.
            intHaz <- 0
          } else {
            # Integrated hazard at max. value
            intHaz <- try(integrate(indRateFct, lower=0, upper=ran)$value, silent=TRUE)
            if(inherits(intHaz, 'try-error')){
              intHaz <- integrate(indRateFct, lower=0, upper=ran, stop.on.error = FALSE, rel.tol = 0.01)$value
            }
          }
          # If transformed random variate exceeds max. value of integr. hazard, we will not find a finite random waiting time.
          if(u<=intHaz){
            invHazFct <- function(x){
              #cat('x: ',x,'\n')
              try.res <- try(integrate(indRateFct, lower=0, upper=x)$value-u, silent=TRUE)
              #print(try.res)
              if(inherits(try.res, 'try-error')){
                #cat('Seemingly, divergent intergral for ID ',id,
                # ' in state ',currState,' at age ',currAge,' at time ',calTime, ' to state ',destState,
                #  ' for random number: ',u,'\n')
                try.res <- integrate(indRateFct, lower=0, upper=x, stop.on.error = FALSE, rel.tol = 0.01)$value-u
              }
              #cat('res: ',try.res,'\n-----\n')
              return(try.res)
            }
            # Find random waiting time.
            timeToNext <- uniroot(invHazFct,interval=c(0,ran))$root
          } else {
            timeToNext <- Inf
          }
        }

        nextEventMatrix[i,1] <- destState
        nextEventMatrix[i,2] <- (timeToNext+lagToWaitingTime)*365.25    # time to next event in days
      }
      #print(nextEventMatrix)
      nE <- nextEventMatrix[which(nextEventMatrix[,2]==min(as.numeric(nextEventMatrix[,2]))),,drop=F]
      if(dim(nE)[1]>1)
        nE <- nE[1,,drop=F]
      if(nE[1,2]!=Inf){
        # Cal. time of next event of individual. (If there is one.)
        tt <- calTime + as.numeric(nE[1,2])
        #print(tt)
        #cat(nE[1,1],'---',tt,'\n')
        # Check whether next event implies school enrollment. If yes, adjust transition time to ensure that the individual
        # enters school at Sept. 1 in the year he/she turns seven.
        if(schoolEnrol){ # Is school enrollment considered in this simulation model? Yes, then continue; otherwise skip this part
          if(isSchoolEnrolment(currState,nE[1,1])){
            enYear <- getYear(tt)
            if(getMonth(tt) <= monthSchoolEnrol) {
              enDate <- getInDays_my(enYear, monthSchoolEnrol)
            } else {
              enYear <- enYear+1
              enDate <- getInDays_my(enYear, monthSchoolEnrol)
            }
            diffToEn <- as.numeric(enDate-tt)
            nE[1,2] <- as.numeric(nE[1,2]) + diffToEn
          }
        }
        # Enqueue new event (if there is one).
        queue <<- rbind(queue, c(id, t.clock, currState, currAge - lagToWaitingTime*365.25, nE[1,1], nE[1,2], birthTime, initState, mothID, fathID, partID))
      }
      #cat('\n----------\n')
      return(nE)
    }

    if(any(depMatrix[,4] == 1 )){ #&& any(depMatrix[,3] == 1)){ ##if dep duration is 0 skip needs to be build in
      #print("test genarg access")
      #TODO: go through added code again and check for stuff leftover
      tr_dur <- rownames(depMatrix)[depMatrix[,3] == 1]

      #extract genarg functions
      tr_genarg <- rownames(depMatrix)[depMatrix[,4] != 0]

      #Exclude transitions for which info isn't available
      if(is.na(match(id, relMatrix[,1]))){
        possTr <- possTr[!possTr %in% tr_genarg]
        nextEventMatrix <- matrix(0, ncol=2, nrow=length(possTr))
      }

      # Extract transition history of individual until current cal. time.
      historiesInd <- transitions[transitions[,1] %in% id & transitions[,4] <= calTime,,drop=F]
      #print(historiesInd)
      # Extract for each state the duration until transition (in days).
      # Here, we have to differ between states of which we do not know when they are entered (i.e., `initial states' of members
      # of the starting population and the states of migrants when they entered the country), and
      # the states we know the `entering date' as well as the `leaving date' (if the state has been left).
      if(birthTime < simStartInDays | isIMInitEvent) {
        dur <- rbind(c(initState,NA),historiesInd[,c(3,4),drop=F]) # extract from historiesInd 'To' and 'transitionTime'
        dur <- cbind(dur,c(diff(dur[,2]),0)) # columns: TransitionTo, AtTime, durUntil
        dur[which(is.na(dur[,2])),3] <- NA
      } else {  # Individual is born during simulation.
        dur <- rbind(c(initState,birthTime),historiesInd[,c(3,4),drop=F])
        dur <- cbind(dur,c(diff(dur[,2]),0)) # columns: TransitionTo, AtTime, durUntil
      }
      # Compute for each possible destination state a waiting time.
      for(i in 1:length(possTr)){
        tr <- possTr[i]
        destState <-  as.numeric(names(tr))
        cS <- codingScheme[codingScheme[,2] %in% currState, 2+c(1:nSubStates)]
        dS <- codingScheme[codingScheme[,2] %in% destState, 2+c(1:nSubStates)]
        # To determine the duration (time elapsed since last transition) that applies for the considered destination state,
        # we have to determine the duration since the last change in the covariate affected.
        # For example, to specify the time being married, we have to determine the duration since (last) marriage.
        covToCh <- which((cS==dS)==F)
        durSinceLastCovCh <- Inf  # For the transition to `dead' so far the time elapsed since the last transition does not play any role.
        if(length(covToCh)==1){
          covHist <- codingScheme[codingScheme[,2] %in% dur[,1], 2+c(1:nSubStates),drop=F][,covToCh]
          idd <- which(covHist==cS[covToCh])
          if(length(idd)>1){
            if(F %in% (diff(idd)==1)){
              y <- rev(idd)[c(-1,diff(rev(idd)))==-1]
              idd <- rev(y)[c(diff(rev(y)),1)==1]
            }
          }
          durSinceLastCovCh <- sum(dur[idd,3]) # If I do not know how long an individual already is in a state: This gives NA.
          if(is.na(durSinceLastCovCh))
            durSinceLastCovCh <- currAge # Then assume the individual is already for his/her whole life in the state.
        }
        if(length(covToCh)>1 & (!destState %in% absStatesNum)){
          cat('Recognized a possible transition implying a change of two or more covariates.',
              'Concerning the derivation of the time being elapsed since the last transition this feature is not yet implemented.',
              'Current State: ',currState,' -> Possible transition to ',destState,'\n')
        }
        tageInYears <- trunc(ageInYears)
        tCalTime <- trunc(1970.001+calTime/365.25)
        tdurSinceLastCovCh <- trunc(durSinceLastCovCh/365.25)
        if(tr %in% tr_dur) {
          indRateFctDET <- function(x){
            res <- eval(do.call(tr,
                                args=list(age=tageInYears+x,calTime=tCalTime+x,duration=tdurSinceLastCovCh+x)))
            return(res)
          }
          ranAccuracyInDays <- (0:(trunc(ran*365.25)+0.99))/365.25
          detE <- indRateFctDET(ranAccuracyInDays)
          daysToTrInYears <- (which(detE == Inf)[1] - 1)/365.25
          if (Inf %in% detE) {
            timeToNext <- daysToTrInYears
          } else {
            u <- -log(1-runif(1))
            #cat('It: ',i,'--u: ',u,'\n')
            # Extract individual transition rate (depending on age, calendar time, and time elapsed)
            indRateFct <- function(x){
              ageIn <- ageInYears+x
              calIn <- 1970.001+calTime/365.25+x
              durIn <- durSinceLastCovCh/365.25+x
              res <- eval(do.call(tr, args=list(age=ageIn,calTime= calIn,duration=durIn)))
              if(TRUE %in% (res<0))
                stop('I have found negative rate value/s for transition: ',tr,'\n
                 This is implausible. Please check this. Simulation has been stopped.\n')
              #cat('x: ',x,' -- res', res,'\n')
              #cat('\n---\n')
              return(res)
            }

            if(sum(indRateFct(0:ran))==0){ # Rate function contains only zeros.
              intHaz <- 0
            } else {
              # Integrated hazard at max. value
              intHaz <- try(integrate(indRateFct, lower=0, upper=ran)$value, silent=TRUE)
              if(inherits(intHaz, 'try-error')){
                intHaz <- integrate(indRateFct, lower=0, upper=ran, stop.on.error = FALSE, rel.tol = 0.01)$value
              }
            }
            # If transformed random variate exceeds max. value of integr. hazard, we will not find a finite random waiting time.
            if(u<=intHaz){
              invHazFct <- function(x){
                #cat('x: ',x,'\n')
                try.res <- try(integrate(indRateFct, lower=0, upper=x)$value-u, silent=TRUE)
                #print(try.res)
                if(inherits(try.res, 'try-error')){
                  #cat('Seemingly, divergent intergral for ID ',id,
                  # ' in state ',currState,' at age ',currAge,' at time ',calTime, ' to state ',destState,
                  #  ' for random number: ',u,'\n')
                  try.res <- integrate(indRateFct, lower=0, upper=x, stop.on.error = FALSE, rel.tol = 0.01)$value-u
                }
                #cat('res: ',try.res,'\n-----\n')
                return(try.res)
              }
              # Find random waiting time.
              timeToNext <- uniroot(invHazFct,interval=c(0,ran))$root
            } else {
              timeToNext <- Inf
            }
          }
          } else {
            if(tr %in% tr_genarg){
              #print("genarg!")
              #TODO: run tests
              genMatrix <- genMatrix[genMatrix[,1] == depMatrix[tr,4],, drop=F]
              #what to do if two access same genMatrix[,1]? -> FORBID to safe ourselves the for loop

              rel_id <- relMatrix[id == relMatrix[,1],, drop=F]

              rel_id_birth <- genMatrix[,2] + 1

              tLinkBirth <- getYear(rel_id[rel_id_birth]) + tageInYears  #trgetInXX(genarg)])/365.25 #
              #TODO:check for robustness

              #TODO:Think about if special consideration for if indiv birthtime is used
              #rel_id[genMatrix[,2] == column relMatrix]
              #if(genMatrix[,2] == 1)
              #tCalTime <- trunc(1970.001+calTime/365.25)

              indRateFctDET <- function(x){
                res <- eval(genMatrix[3] + do.call(tr,
                                    args=list(age=tageInYears+x,calTime=tLinkBirth+x))* genMatrix[4]) #tCalTime = ; #tLinkBirth+tageInYears+x
                #print(res[1]) #check if function works
                return(res)
              }
              ranAccuracyInDays <- (0:(trunc(ran*365.25)+0.99))/365.25
              detE <- indRateFctDET(ranAccuracyInDays)
              daysToTrInYears <- (which(detE == Inf)[1] - 1)/365.25
              if (Inf %in% detE) {
                timeToNext <- daysToTrInYears
              } else {
                u <- -log(1-runif(1))
                #cat('It: ',i,'--u: ',u,'\n')
                # Extract individual transition rate (depending on age and calendar time)
                indRateFct <- function(x){
                  ageIn <- ageInYears+x
                  calIn <- tLinkBirth + ageInYears + x #1970.001+calTime/365.25+x #truncation for linkbirth should be irrelevant? check!
                  res <- eval(genMatrix[3] + do.call(tr, args=list(age=ageIn,calTime= calIn))* genMatrix[4]) #untruncated versions
                  if(TRUE %in% (res<0))
                    stop('I have found negative rate value/s for transition: ',tr,'\n
                 This is implausible. Please check this. Simulation has been stopped.\n')
                  #cat('x: ',x,' -- res', res,'\n')
                  #cat('\n---\n')
                  return(res)
                }
                if(sum(indRateFct(0:ran))==0){ # Rate function contains only zeros.
                  intHaz <- 0
                } else {
                  # Integrated hazard at max. value
                  intHaz <- try(integrate(indRateFct, lower=0, upper=ran)$value, silent=TRUE)
                  if(inherits(intHaz, 'try-error')){
                    intHaz <- integrate(indRateFct, lower=0, upper=ran, stop.on.error = FALSE, rel.tol = 0.01)$value
                  }
                }
                # If transformed random variate exceeds max. value of integr. hazard, we will not find a finite random waiting time.
                if(u<=intHaz){
                  invHazFct <- function(x){
                    #cat('x: ',x,'\n')
                    try.res <- try(integrate(indRateFct, lower=0, upper=x)$value-u, silent=TRUE)
                    #print(try.res)
                    if(inherits(try.res, 'try-error')){
                      #cat('Seemingly, divergent intergral for ID ',id,
                      # ' in state ',currState,' at age ',currAge,' at time ',calTime, ' to state ',destState,
                      #  ' for random number: ',u,'\n')
                      try.res <- integrate(indRateFct, lower=0, upper=x, stop.on.error = FALSE, rel.tol = 0.01)$value-u
                    }
                    #cat('res: ',try.res,'\n-----\n')
                    return(try.res)
                  }
                  # Find random waiting time.
                  timeToNext <- uniroot(invHazFct,interval=c(0,ran))$root
                } else {
                  timeToNext <- Inf
                }
              }
              } else {
            indRateFctDET <- function(x){
              res <- eval(do.call(tr,
                                  args=list(age=tageInYears+x,calTime=tCalTime+x)))
              return(res)
            }
            ranAccuracyInDays <- (0:(trunc(ran*365.25)+0.99))/365.25
            detE <- indRateFctDET(ranAccuracyInDays)
            daysToTrInYears <- (which(detE == Inf)[1] - 1)/365.25
            if (Inf %in% detE) {
              timeToNext <- daysToTrInYears
            } else {
              u <- -log(1-runif(1))
              #cat('It: ',i,'--u: ',u,'\n')
              # Extract individual transition rate (depending on age and calendar time)
              indRateFct <- function(x){
                ageIn <- ageInYears+x
                calIn <- 1970.001+calTime/365.25+x
                res <- eval(do.call(tr, args=list(age=ageIn,calTime= calIn)))
                if(TRUE %in% (res<0))
                  stop('I have found negative rate value/s for transition: ',tr,'\n
                 This is implausible. Please check this. Simulation has been stopped.\n')
                #cat('x: ',x,' -- res', res,'\n')
                #cat('\n---\n')
                return(res)
              }
              if(sum(indRateFct(0:ran))==0){ # Rate function contains only zeros.
                intHaz <- 0
              } else {
                # Integrated hazard at max. value
                intHaz <- try(integrate(indRateFct, lower=0, upper=ran)$value, silent=TRUE)
                if(inherits(intHaz, 'try-error')){
                  intHaz <- integrate(indRateFct, lower=0, upper=ran, stop.on.error = FALSE, rel.tol = 0.01)$value
                }
              }
              # If transformed random variate exceeds max. value of integr. hazard, we will not find a finite random waiting time.
              if(u<=intHaz){
                invHazFct <- function(x){
                  #cat('x: ',x,'\n')
                  try.res <- try(integrate(indRateFct, lower=0, upper=x)$value-u, silent=TRUE)
                  #print(try.res)
                  if(inherits(try.res, 'try-error')){
                    #cat('Seemingly, divergent intergral for ID ',id,
                    # ' in state ',currState,' at age ',currAge,' at time ',calTime, ' to state ',destState,
                    #  ' for random number: ',u,'\n')
                    try.res <- integrate(indRateFct, lower=0, upper=x, stop.on.error = FALSE, rel.tol = 0.01)$value-u
                  }
                  #cat('res: ',try.res,'\n-----\n')
                  return(try.res)
                }
                # Find random waiting time.
                timeToNext <- uniroot(invHazFct,interval=c(0,ran))$root
              } else {
                timeToNext <- Inf
              }
            }

              }
            }
        nextEventMatrix[i,1] <- destState
        nextEventMatrix[i,2] <- (timeToNext+lagToWaitingTime)*365.25    # time to next event in days
      }
      #print(nextEventMatrix)
      nE <- nextEventMatrix[which(nextEventMatrix[,2]==min(as.numeric(nextEventMatrix[,2]))),,drop=F]
      if(dim(nE)[1]>1)
        nE <- nE[1,,drop=F]
      if(nE[1,2]!=Inf){
        # Cal. time of next event of individual. (If there is one.)
        tt <- calTime + as.numeric(nE[1,2])
        #print(tt)
        #cat(nE[1,1],'---',tt,'\n')
        # Check whether next event implies school enrollment. If yes, adjust transition time to ensure that the individual
        # enters school at Sept. 1 in the year he/she turns seven.
        if(schoolEnrol){ # Is school enrollment considered in this simulation model? Yes, then continue; otherwise skip this part
          if(isSchoolEnrolment(currState,nE[1,1])){
            enYear <- getYear(tt)
            if(getMonth(tt) <= monthSchoolEnrol) {
              enDate <- getInDays_my(enYear, monthSchoolEnrol)
            } else {
              enYear <- enYear+1
              enDate <- getInDays_my(enYear, monthSchoolEnrol)
            }
            diffToEn <- as.numeric(enDate-tt)
            nE[1,2] <- as.numeric(nE[1,2]) + diffToEn
          }
        }
        # Enqueue new event (if there is one).
        queue <<- rbind(queue, c(id, t.clock, currState, currAge - lagToWaitingTime*365.25, nE[1,1], nE[1,2], birthTime, initState, mothID, fathID, partID))
      }
      #cat('\n----------\n')
      return(nE)
    }
    else
      for(i in 1:length(possTr)){
        tr <- possTr[i]
        destState <-  as.numeric(names(tr))
        cS <- codingScheme[codingScheme[,2] %in% currState, 2+c(1:nSubStates)]
        dS <- codingScheme[codingScheme[,2] %in% destState, 2+c(1:nSubStates)]

        tageInYears <- trunc(ageInYears)
        tCalTime <- trunc(1970.001+calTime/365.25)
        indRateFctDET <- function(x){
          res <- eval(do.call(tr,
                              args=list(age=tageInYears+x,calTime=tCalTime+x)))
          return(res)
        }
        ranAccuracyInDays <- (0:(trunc(ran*365.25)+0.99))/365.25
        detE <- indRateFctDET(ranAccuracyInDays)
        daysToTrInYears <- (which(detE == Inf)[1] - 1)/365.25
        if (Inf %in% detE) {
          timeToNext <- daysToTrInYears
        } else {
          u <- -log(1-runif(1))
          #cat('It: ',i,'--u: ',u,'\n')
          # Extract individual transition rate (depending on age and calendar time)
          indRateFct <- function(x){
            ageIn <- ageInYears+x
            calIn <- 1970.001+calTime/365.25+x
            res <- eval(do.call(tr, args=list(age=ageIn,calTime= calIn)))
            if(TRUE %in% (res<0))
              stop('I have found negative rate value/s for transition: ',tr,'\n
                 This is implausible. Please check this. Simulation has been stopped.\n')
            #cat('x: ',x,' -- res', res,'\n')
            #cat('\n---\n')
            return(res)
          }
          if(sum(indRateFct(0:ran))==0){ # Rate function contains only zeros.
            intHaz <- 0
          } else {
            # Integrated hazard at max. value
            intHaz <- try(integrate(indRateFct, lower=0, upper=ran)$value, silent=TRUE)
            if(inherits(intHaz, 'try-error')){
              intHaz <- integrate(indRateFct, lower=0, upper=ran, stop.on.error = FALSE, rel.tol = 0.01)$value
            }
          }
          # If transformed random variate exceeds max. value of integr. hazard, we will not find a finite random waiting time.
          if(u<=intHaz){
            invHazFct <- function(x){
              #cat('x: ',x,'\n')
              try.res <- try(integrate(indRateFct, lower=0, upper=x)$value-u, silent=TRUE)
              #print(try.res)
              if(inherits(try.res, 'try-error')){
                #cat('Seemingly, divergent intergral for ID ',id,
                # ' in state ',currState,' at age ',currAge,' at time ',calTime, ' to state ',destState,
                #  ' for random number: ',u,'\n')
                try.res <- integrate(indRateFct, lower=0, upper=x, stop.on.error = FALSE, rel.tol = 0.01)$value-u
              }
              #cat('res: ',try.res,'\n-----\n')
              return(try.res)
            }
            # Find random waiting time.
            timeToNext <- uniroot(invHazFct,interval=c(0,ran))$root
          } else {
            timeToNext <- Inf
          }
        }
        nextEventMatrix[i,1] <- destState
        nextEventMatrix[i,2] <- (timeToNext+lagToWaitingTime)*365.25    # time to next event in days
      }
    #print(nextEventMatrix)
    nE <- nextEventMatrix[which(nextEventMatrix[,2]==min(as.numeric(nextEventMatrix[,2]))),,drop=F]
    if(dim(nE)[1]>1)
      nE <- nE[1,,drop=F]
    if(nE[1,2]!=Inf){
      # Cal. time of next event of individual. (If there is one.)
      tt <- calTime + as.numeric(nE[1,2])
      #print(tt)
      #cat(nE[1,1],'---',tt,'\n')
      # Check whether next event implies school enrollment. If yes, adjust transition time to ensure that the individual
      # enters school at Sept. 1 in the year he/she turns seven.
      if(schoolEnrol){ # Is school enrollment considered in this simulation model? Yes, then continue; otherwise skip this part
        if(isSchoolEnrolment(currState,nE[1,1])){
          enYear <- getYear(tt)
          if(getMonth(tt) <= monthSchoolEnrol) {
            enDate <- getInDays_my(enYear, monthSchoolEnrol)
          } else {
            enYear <- enYear+1
            enDate <- getInDays_my(enYear, monthSchoolEnrol)
          }
          diffToEn <- as.numeric(enDate-tt)
          nE[1,2] <- as.numeric(nE[1,2]) + diffToEn
        }
      }
      # Enqueue new event (if there is one).
      queue <<- rbind(queue, c(id, t.clock, currState, currAge - lagToWaitingTime*365.25, nE[1,1], nE[1,2], birthTime, initState, mothID, fathID, partID))
    }
    #cat('\n----------\n')
    return(nE)
  }

  # ----------------------------------------------------------------------------------------------------------------------
  # ----------------------------------------------------------------------------------------------------------------------
  # E. INITIALIZATION
  # ----------------------------------------------------------------------------------------------------------------------
  # ----------------------------------------------------------------------------------------------------------------------
  # Compute next events for members of starting population
  cat('Initialization ... \n')
  time_init_start = Sys.time()
  print(paste("Starting initialization at: ", time_init_start))

  if(length(fertTr)>0){
    fertTrExpanded <- buildFertTrExpanded()
  }
  if(schoolEnrol){
    eduTrExpanded <- buildEduTrExpanded()
  }
  if(length(partTr)>0){
    partTrExpanded <- buildPartTrExpanded()
    allMateSearches <- 0
    unsuccMateCounter <- 0
    if(!(rule %in% c(0:2))) rule <- 0 # not defined partner choice rule leads to random draw
    if(rule %in% 1){# consider age difference for choosing partner
      if(is.null(ageDiffDistr)){
        ageDiffDistr <- function(ageDiff){
          return(ifelse(abs(ageDiff)<=5, 1,0))
        }
      } else { # test whether given ageDiffDistr function gives feasible results
        if(anyNA(ageDiffDistr(seq(from=-c((maxAge-18)), to=(maxAge-18), length=1000)))){
          cat("There are NA values in the function defining the success probabilities for mating along age differences.")
          indNAAgeDiff <- which(is.na(ageDiffDistr(seq(from=-c((maxAge-18)), to=(maxAge-18), length=1000))), arr.ind=TRUE)
          cat("NA values occur for the following age differences (age of potential male spouse - age of potential female spouse): ",seq(from=-c((maxAge-18)), to=(maxAge-18), length=1000)[indNAAgeDiff], "\n")
          stop("Found NAs in success probability function for mating along age differences.")
        }
      }
    }
  }
  if(length(sepTr)>0){
    sepTrExpanded <- buildSepTrExpanded()
  }
  if(length(absPartTr)>0){
    absPartTrExpanded <- buildAbsPartTrExpanded()
  }

  motherID <- rep(NA, nrow(initPop))
  fatherID <- rep(NA, nrow(initPop))
  partnerID <- rep(NA, nrow(initPop))
  initPop <- cbind.data.frame(initPop, motherID, fatherID, partnerID)

  birthTimeInDays <- getInDays(initPop[,'birthDate'])
  IN <- matrix(c(initPop[,'ID'], # ID
                 rep(99,nrow(initPop)), # currState (= initState)
                 simStartInDays-birthTimeInDays, # age
                 rep(simStartInDays,nrow(initPop)), # calenderTime
                 birthTimeInDays, # birth time in days
                 rep(99,nrow(initPop)),  # initState
                 initPop[,'motherID'], # motherID
                 initPop[,'fatherID'], # fatherID
                 initPop[,'partnerID']), # partnerID
            ncol=9, nrow=nrow(initPop))
  IN[,2] <- IN[,6] <- codingScheme[match(initPop[,'initState'], codingScheme[,1]),2]

  if(TRUE %in% (IN[,3]<0)) {
    cat("There are persons born later than simulation starting date in the initial population. Related IDs are: ")
    negAge <- IN[,1][IN[,3]<0]
    for(i in 1:length(negAge)){
      cat(negAge[i]," ")
    }
    stop("Error: Negative age in initial population.")
  }
  if(TRUE %in% (IN[,3]/365.25>maxAge)) {
    cat("In the initial population, there are persons older than maxAge at simulation starting date. Related IDs are: ")
    invalAge <- IN[,1][IN[,3]/365.25>maxAge]
    for(i in 1:length(invalAge)){
      cat(invalAge[i]," ")
    }
    stop("Error: Older than max. age in initial population.")
  }

  if(length(partTr)>0){ # create initial pool of potential mates
    if(sum(is.na(IN[,9]))>0){
      timeAvailable <- IN[,4]
      matePool <- cbind(IN[is.na(IN[,9,drop=FALSE]),c(1,2,5:8)], timeAvailable) # id, currentState, birthday in days, initState, motherID, fatherID, time when available
      matePool <- matePool[!(matePool[,2] %in% partTrExpanded[,2]),] # exclude persons with partner states
    }
    if(sum(is.na(IN[,9]))==0){
      matePool <- NULL
    }
    allS <- do.call(rbind,sapply(rownames(transitionMatrix),strsplit, split="/"))
    femStates <- rownames(transitionMatrix)[which(allS %in% "f", arr.ind=T)]
    maleStates <- rownames(transitionMatrix)[which(allS %in% "m", arr.ind=T)]
    femStatesNum <- codingScheme[match(femStates, codingScheme[,1]),2]
    maleStatesNum <- codingScheme[match(maleStates, codingScheme[,1]),2]
  }

  maxId <- max(IN[,1])

  init <- apply(IN, 1, getNextStep)

  # If immigrants enter the population, compute next events for them.
  if(!is.null(immigrPop)){

    # Check whether migrants are already born when they migrate
    if(TRUE %in% (immigrPop$immigrDate<immigrPop$birthDate)){
      cat("In the immigration population, there are persons who are not yet born when they immigrate. Related IDs are: ")
      invalImDate <- immigrPop$ID[immigrPop$immigrDate<immigrPop$birthDate]
      for(i in 1:length(invalImDate)){
        cat(invalImDate[i]," ")
      }
      stop("Error: Not yet born at immigration date.")
    }
    # Check whether all migrants migrate after simulation starting date
    if(TRUE %in% (immigrPop$immigrDate<simHorizon[1])){
      cat("In the immigration population, there are persons who are specified to migrate into the virtual population before simulation starting time. That's against MicSim's concept of migration. Related IDs are: ")
      invalImDate <- immigrPop$ID[immigrPop$immigrDate<simHorizon[1]]
      for(i in 1:length(invalImDate)){
        cat(invalImDate[i]," ")
      }
      stop("Error: Migrate before simulation starting date.")
    }
    # Check whether all migrants migrate before simulation stopping date
    if(TRUE %in% (immigrPop$immigrDate>simHorizon[2])){
      cat("In the immigration population, there are persons who are specified to migrate into the virtual population after simulation ending time. That's a bit meaningless. Related IDs are: ")
      invalImDate <- immigrPop$ID[immigrPop$immigrDate>simHorizon[2]]
      for(i in 1:length(invalImDate)){
        cat(invalImDate[i]," ")
      }
      stop("Error: Migrate after simulation ending date.")
    }

    motherID <- rep(NA, nrow(immigrPop))
    fatherID <- rep(NA, nrow(immigrPop))
    partnerID <- rep(NA, nrow(immigrPop))
    immigrPop <- cbind(immigrPop, motherID, fatherID, partnerID)

    immigrTimeInDays <- getInDays(immigrPop$immigrDate)
    birthTimeInDays <- getInDays(immigrPop$birthDate)
    IM <- matrix(c(immigrPop[,'ID'], # ID
                   rep(99,nrow(immigrPop)), # currState (= initState)
                   immigrTimeInDays - birthTimeInDays, #age
                   immigrTimeInDays, # calenderTime
                   birthTimeInDays, # birth time in days
                   rep(99,nrow(immigrPop)),  # initState
                   immigrPop[, 'motherID'], # motherID
                   immigrPop[, 'fatherID'], # fatherID
                   immigrPop[, 'partnerID']), # partnerID
                 ncol=9, nrow=nrow(immigrPop))
    IM[,2] <- IM[,6] <- codingScheme[match(immigrPop[,'immigrInitState'], codingScheme[,1]),2]

    maxId <- max(IN[,1], IM[,1])

    # Check whether all migrants are younger than maxAge when they migrate
    ageIm <- IM[,3]/365.25 # age at immigration in years
    if(TRUE %in% (ageIm>maxAge)){
      cat("In the immigration population, there are persons who older than `maxAge' when they migrate. That's a bit meaningless. Related IDs are: ")
      invalImAge <- immigrPop$ID[ageIm>maxAge]
      for(i in 1:length(invalImAge)){
        cat(invalImAge[i]," ")
      }
      stop("Error: Migrants in the input data are older than `maxAge'.")
    }

    if(length(partTr)>0){ # add to mate pool migrants
      if(sum(is.na(IM[,9]))>0){
        timeAvailable <- IM[,4]
        matePool_IM <- cbind(IM[is.na(IM[,9]),c(1,2,5:8)], timeAvailable)
        matePool_IM <- matePool_IM[!(matePool_IM[,2] %in% partTrExpanded[,2]),] # exclude persons with partner states
        matePool <- rbind(matePool, matePool_IM)
      }
    }


    imit <- apply(IM, 1, getNextStep, isIMInitEvent=T)


    immigrInitPop <- immigrPop[,c('ID','birthDate','immigrInitState','motherID', 'fatherID', 'partnerID')]
    colnames(immigrInitPop)[3] <- 'initState'
    initPop <- rbind(initPop, immigrInitPop)

  }

  time_init_end = Sys.time()
  print(paste("Ending initialization at:   ", time_init_end))
  #p_time = (time_init_end - time_init_start)
  #print(paste("#Time needed for initialization: ", p_time))

  # ----------------------------------------------------------------------------------------------------------------------
  # ----------------------------------------------------------------------------------------------------------------------
  # F. SIMULATION
  # ----------------------------------------------------------------------------------------------------------------------
  # ----------------------------------------------------------------------------------------------------------------------
  # Run simulation either until queue is empty or until simulation horizon has been reached.
  cat('Simulation is running ... \n')
  time_sim_start = Sys.time()
  print(paste("Starting simulation at: ", time_sim_start))

  currYear <- trunc(simHorizon[1]/10000)
  cat('Year: ',currYear,'\n')
  while(nrow(queue)>0 & t.clock <= simStopInDays){

    #cat('\n-----------\n')

    # Sort queue according to soonest event to happen.
    queue <- queue[order(queue[,2] + queue[,6]),,drop=F] # columns: 'ID','currTime','currState','currAge','nextState','timeToNextState','birthtime', 'initState', 'motherID', 'fatherID', 'partnerID'; in queue currTime in days since 01-01-1970
    #print_t.clock <- paste(c(trunc(getDay(t.clock)), getMonth(t.clock), getYear(t.clock)), collapse="/")
    #cat("Time: ", print_t.clock, "\n")
    #print(dim(queue)[1])
    # Enqueue individual who has the soonest event to happen.
    indS <- queue[1,]
    #cat("Event for: ",indS[1],"\n")
    #print(indS)
    #cat("\n")
    # Remove he/she from queue.
    queue <- queue[-1,,drop=F]
    # Set the global clock.
    t.clock <- indS[2] + indS[6] # in days since 01-01-1970
    cY <- getYear(t.clock) # transform days since 01-01-1970 back to years
    # If the global clock exceeds the simulation horizon, stop simulation.
    if(t.clock > simStopInDays)
      break
    if(cY>currYear){
      cat('Year: ',cY,'\n')
      currYear <- cY
    }
    # Age at current transition
    age <- indS[4] + indS[6] # in days

    # If current state is not an absorbent one, check whether newborn, partnership onset or dissolution is caused, trigger events of linked individuals
    if(!indS[5] %in% absStatesNum){

      # Current transition causes a newborn? If yes, add one to simulation population.
      if(length(fertTr)>0){
        if(isBirthEvent(indS[3],indS[5])){
          #cat("Newborn: Mother: ",indS[1], " - Father: ",indS[11], "\n")
            addNewNewborn_dur(birthTime=t.clock, motherID=indS[1], motherState=indS[5], motherbirth=indS[7], fatherID=indS[11]) #this could add mothers to relMatrix -> this is fertility for mother
        }
      }

      # Current transition cause the onset of a partnership
      if(length(partTr)>0){
        if(isPartEvent(indS[3],indS[5])){
          #print("Partner search event \n")
          allMateSearches <- allMateSearches + 1
          #print(indS)
          resPart <- createPartLink(partTime=t.clock, egoID=indS[1], egoState=indS[5], egoBirth=indS[7], egoMothID=indS[9], egoFathID=indS[10], rule=rule)
              # returned partner infos: id, currentState, birthday in days, initState, motherID, fatherID, time when available
          partID <- resPart[1]
          if(is.na(partID)){
            unsuccMateCounter <- unsuccMateCounter + 1
            #cat("No partner found for: id: ", indS[1], " for transition ", indS[3] , " -> ", indS[5], " at time ", t.clock,".\n")
          }
          if(!is.na(partID)){
            #cat("Partner found: Ego: ",indS[1], " -- Partner: ", partID, "\n")
            indS[11] <- partID
            queue <- queue[!(queue[,1] %in% partID),] # if partner is in the queue, remove her/him, register her/him, and compute new event for him
            #cat("Partner: ", resPart,"\n")
            lS_part <- resPart[2]
            cS_part <- partTrExpanded[partTrExpanded[,1] %in% resPart[2],2] # derived state for partner due to matching
            age_part <- t.clock - resPart[3]
            transitions <- rbind(transitions, c(partID, lS_part, cS_part, t.clock, age_part, resPart[5:6], indS[1])) # register partnership event for mate
             # compute next event for partner
              resP <- getNextStep(c(partID, cS_part, age_part, t.clock, resPart[c(3:6)], indS[1])) # ID, currState, age, calTime, birthtime, initState, motherID, fatherID, partID
          }
        }
      }

#      # Current transition causes the dissolution of a partnership, delete partnership ID
#      if(length(sepTr)>0){
#        partID <- indS[11]
#      }

      # Current transition causes the dissolution of a partnership also for partner; compute new event for the partner after dissolution
      if(length(sepTr)>0){
        if(isSepEvent(indS[3],indS[5]) & !is.na(indS[11])){
          #print("Separation event \n")
          #print(indS)

          partID <- indS[11] # make separation event for partner as well
          indS[11] <- NA # delete partnerID for next event
          indS_part <- queue[queue[,1] %in% partID,]
          if(length(indS_part)==0) {# partner is not in queue (since for him/her no further event is scheduled), take info from transition register
            partTrInfo <- transitions[transitions[,1] %in% partID,,drop=F]
            indS_part <- c(partTrInfo[nrow(partTrInfo),1],NA,
                           partTrInfo[nrow(partTrInfo),3],NA,NA,NA,
                           getInDays(initPop[initPop[,1] %in% partID,2]),
                           partTrInfo[1,2],partTrInfo[nrow(partTrInfo),6:7],NA)
          } else { # otherwise: partner is in queue, de-queue her/him and take info on partner from queue
            queue <- queue[!(queue[,1] %in% partID),]
          }
          #print(indS_part)
          indS_part[11] <- NA
          list_nextEvent <- sepTrExpanded[sepTrExpanded[,1] %in% indS_part[3],2] # next state

          #if( length(list_nextEvent) != length(sepTr)){
          #  cat("Inconsistency: length(list_nextSepEvent): ",length(list_nextEvent), " is not length(sepTr): ",length(sepTr),"\n")
          #}

          indS_part[5] <- list_nextEvent[sample(x=1:length(list_nextEvent),size=1,prob=probSepTr)]
          age_part <- t.clock - indS_part[7]

          transitions <- rbind(transitions, c(indS_part[c(1,3,5)], t.clock, age_part, indS_part[c(9:11)])) # register separation event of partner

          resP <- getNextStep(c(indS_part[c(1,5)], age_part, t.clock, indS_part[c(7:10)], NA)) # ID, currState, age, calTime, birthtime, initState, motherID, fatherID, partID

          matePool <- rbind(matePool, c(indS[c(1,5,7:10)],t.clock))
          matePool <- rbind(matePool, c(indS_part[c(1,5,7:10)],t.clock))
        }
      }
    }

    # Current transition cause widowhood or similar
    if(length(absPartTr)>0){
      if((indS[5] %in% absStatesNum) & !is.na(indS[11])){
        #print("Widowhooed event \n")
        partID <- indS[11] # partner experiences a widowhood or similar event
        indS_part <- queue[queue[,1] %in% partID,]
        if(length(indS_part)==0) {# partner is not in queue (since for him/her no further event is scheduled), take info from transition register
          partTrInfo <- transitions[transitions[,1] %in% partID,,drop=F]
          indS_part <- c(partTrInfo[nrow(partTrInfo),1],NA,
                         partTrInfo[nrow(partTrInfo),3],NA,NA,NA,
                         getInDays(initPop[initPop[,1] %in% partID,2]),
                         partTrInfo[1,2],partTrInfo[nrow(partTrInfo),6:7],NA)
        } else { # otherwise: partner is in queue, de-queue her/him and take info on partner from queue
          queue <- queue[!(queue[,1] %in% partID),]
        }
        indS_part[11] <- NA
        indS_part[5] <- absPartTrExpanded[absPartTrExpanded[,1] %in% indS_part[3],2] # next state
        age_part <- t.clock - indS_part[7]
        transitions <- rbind(transitions, c(indS_part[c(1,3,5)], t.clock, age_part, indS_part[c(9:11)]))

        resP <- getNextStep(c(indS_part[c(1,5)], age_part, t.clock, indS_part[c(7:10)], NA)) # ID, currState, age, calTime, birthtime, initState, mothID, fathID, partID

        matePool <- rbind(matePool, c(indS_part[c(1,5,7:10)],t.clock))
        indS[11] <- partID
      }


    # Remove people experiencing a transition to an absorbing event from mate pool
    if(indS[5] %in% absStatesNum){
      matePool <- matePool[!(matePool[,1] %in% indS[1]),]
    }
    }
    # Compute next event for EGO if current event is not an absorbent one
    if(!indS[5] %in% absStatesNum){

    res <- getNextStep(c(indS[c(1,5)], age, t.clock, indS[c(7:11)])) # ID, currState, age, calTime, birthtime, initState, isMig, mothID, fathID, partnID
    }
    #print(res)

    # Register transition.
    transitions <- rbind(transitions, c(indS[c(1,3,5)], t.clock, age, indS[c(9:11)]))

    #print(transitions)
    #cat('\n-----------\n')
  }
  transitions <- transitions[order(transitions[,1]),,drop=F] # columns: ID, From, To, transitionTime, transitionAge, motherID, partnerID

  if (nrow(transitions) == 0){

    transitionsOut <- data.frame(ID=as.numeric(initPop[,'ID']), From= rep(NA,nrow(initPop)),
                                 To=rep(NA,nrow(initPop)), transitionTime = rep(NA,nrow(initPop)),
                                 transitionAge = rep(NA,nrow(initPop)),
                                 motherID = rep(NA,nrow(initPop)),
                                 fatherID = rep(NA,nrow(initPop)),
                                 partnerID = rep(NA,nrow(initPop)),
                                 stringsAsFactors = FALSE)
    cat('Simulation has finished.\n')
    time_sim_end = Sys.time()
    print(paste("Ending simulation at: ", time_sim_end))
    cat('Beware that along the simulation horizon the individual/s considered do/es not experience any transition/s.\n')
    cat('------------------\n')

  } else {

    cat('Simulation has finished.\n')
    time_sim_end = Sys.time()
    print(paste("Ending simulation at: ", time_sim_end))
    if(length(partTr) > 0){
    cat('Number of unsuccessful partner search events: ', unsuccMateCounter, 'of', allMateSearches, '(', round(unsuccMateCounter/allMateSearches*100,2),'% )\n')
    cat('\n------------------\n')
    }

    # ----------------------------------------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------------
    # G. GENERATE OUTPUT
    # ----------------------------------------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------------

    indCodesTo <- match(transitions[,2], codingScheme[,2]) # transform numerical state codes to strings according to codingScheme
    transTo <- codingScheme[indCodesTo,1]
    indCodesFrom <- match(transitions[,3], codingScheme[,2])
    transFrom <- codingScheme[indCodesFrom,1]
    transitionsOut <- data.frame(ID=transitions[,1], From=transTo , To=transFrom,
                                 transitionTime = getInDateFormat(transitions[,4]),
                                 transitionAge = round(transitions[,5]/365.25,2),
                                 motherID = transitions[,6],
                                 fatherID = transitions[,7],
                                 partnerID = transitions[,8],
                                 stringsAsFactors = FALSE)
  }

  pop <- merge(initPop, transitionsOut, all=T, by='ID')
  pop <- pop[order(as.numeric(pop[,1])),]
  pop$motherID <- ifelse(!is.na(pop$motherID.x), pop$motherID.x, pop$motherID.y)
  pop$fatherID <- ifelse(!is.na(pop$fatherID.x), pop$fatherID.x, pop$fatherID.y)
  pop$partnerID <- ifelse(!is.na(pop$partnerID.y), pop$partnerID.y, pop$partnerID.x)
  pop <- pop[, !(colnames(pop) %in% c("motherID.x", "motherID.y", "fatherID.x", "fatherID.y", "partnerID.x", "partnerID.y"))]
  return(pop)
}

