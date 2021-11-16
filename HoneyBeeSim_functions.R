#For function storage/ developement
#TODOs

#1) Prevent mother-son mating: probably some mother-son mating in nature as well - so that's fine

# Functions
#3) Revise Pull individuals from the caste (update the cast)
#4) Create a function that create multiple virgin queens (return colony with multiple virgin queens): Jernej
#11) Synchonize the print of the Colony with the AlphaSimR Pop: Jana
#13) Replace nDrones/nWorkers with just n / nInd: Jana
#14) Write a function to remove workers/drones from the colony
#15) Put SimParam in the createColony: Jana
#17) Create a function to extract the YOB
#18) Create a function to compute the age of the queen
#9) Create a function to sample the location for the swarm: later
#10) Create a function to sample locations within a radius: later
#16) Write a setPheno function for the colony: later
#12) Create a function for the Loss of the colony


# Text
#5) THink about providing informative messages for the functions: Laura
#6) Think of a good names for the swarmed colony (the one that stay)
#"20) Write descriptions for all the functions

# Think
#7) Think about removing workers and drones in "instantaneous" functions (opposite to adding them)
#19) THink about replacing phenotypes in the swarm/supersede/split


#TODO for script
#S1) Distribute supersedure/swarming events throughout the year/seasons
#S2) Distribute colony losses events throughout the year/seasons

#' #-----------------------------------------------------------------------
#' # Class colony list
#' #-----------------------------------------------------------------------
#' #' @title ColonyList
#' #' #'
#' #' @description
#' #' The ColonyList represents the list of the colonies.
#' #' It is designed to behave like a list of colonies.
#' #'
#' #' @param x a 'ColonyList' object
#' #' @param i index of populations or colonies
#' #'
#' #' @slot colony list of \code{\link{Colony-class}} and/or
#' #' \code{ColonyList-class}
#' #'
#' #'
#' #' @export
#' setClass("ColonyList",
#'          slots=c(colonies="list"))
#' 
#' 
#' #' @describeIn ColonyList Extract ColonyList by index
#' setMethod("[",
#'           signature(x = "ColonyList"),
#'           function(x, i){
#'             x@colonies = x@colonies[i]
#'             return(x)
#'           }
#' )
#' 
#' #' @describeIn ColonyList Extract Colony by index
#' setMethod("[[",
#'           signature(x = "ColonyList"),
#'           function (x, i){
#'             return(x@colonies[[i]])
#'           }
#' )
#' 
#' #' @describeIn ColonyLIst Combine multiple ColonyLists
#' setMethod("c",
#'           signature(x = "ColonyList"),
#'           function (x, ...){
#'             for(y in list(...)){
#'               if(class(y)=="NULL"){
#'                 # Do nothing
#'               }else{
#'                 if(class(y)=="Colony"){
#'                   x@colonies = c(x@colonies, y)
#'                 }else{
#'                   stopifnot(class(y)=="MegaPop")
#'                   x@colonies = c(x@colonies, y@colonies)
#'                 }
#'               }
#'             }
#'             return(x)
#'           }
#' )
#' 
#' 
#' #' @title Create new ColonyList
#' #'
#' #' @description
#' #' Creates a new \code{\link{ColonyList-class}} from one or more
#' #' \code{\link{Colony-class}} and/or \code{\link{ColonyList-class}}
#' #' objects.
#' #'
#' #' @param ... one or more \code{\link{Colony-class}} and/or
#' #' \code{\link{ColonyList-class}} objects.
#' #'
#' #' @return Returns an empty object of \code{\link{ColonyList-class}}
#' #'
#' #' @examples
#' #' @export
#' 
#' createColonyList = function(...){
#'   input = list(...)
#'   class = sapply(input, "class")
#'   stopifnot(all(class=="Colony" | class=="ColonyList") | all(class=="NULL"))
#'   output = new("ColonyList", colonies=input)
#'   return(output)
#' }
#' 
#' #' @title Add colony to the ColonyList
#' #'
#' #' @description
#' #'
#' #' @param 
#' #'
#' #' @return 
#' #'
#' #' @examples
#' #' @export
#' 
#' addColonyToTheColonyList= function(colony, colonyList){
#'   if (class(colony) != "Colony") {
#'     message("The colony parameter is not a Colony object.")
#'   }
#'   colonyList@colonies = append(colonyList@colonies, list(colony))
#'   return(colonyList)
#' }


#-----------------------------------------------------------------------
# Create colony
#-----------------------------------------------------------------------
#' @rdname createColony
#' @method createColony
#' @title Creates a honeybee colony
#' @usage \method{createColony}(id, location, queen, drones, workers, virgin_queen, fathers, pheno, last_event)
#' @description List. Creates a honeybee colony as a list with the following elements:
#'   \id location, queen, drones, workers, virgin_queens, pheno, fathers, and last event
#'   \All elements of the list, expect for \code{last_even}, are assumed NULL if not specified
#'   \otherwise. \code{last_event} is set to "new_colony".
#'   # TODO: we will likely need queen age too - but that should go into colony@queen@misc slot!
     # TODO: can also look at hive "strength" based on number of colony workers
#' @param id ID of the colony.
#' @param location Location of the colony.
#' @param queen AlphaSimR individual object to become the queen of the colony.
#' @param drones AlphaSimR population object to become the drones of the colony.
#' @param workers AlphaSimR population object to become the workers of the colony.
#' @param virgin_queens AlphaSimR individual or population object to become the virgin queen(s) of the colony.
#' @param pheno
#' @param fathers AlphaSimR population object of the fathers of the colony - i.e. the drones the queen mated with (semen in spermatheca)
#' @param last_event Last event of the colony. Default is set to "new_colony", other possible values are TODO?????
#'
#' @example inst/examples/examples_createColony.R
#' @return Returns AlphaSimR class "Colony" object.
#' @export

setClassUnion("PopOrNULL", c("Pop", "NULL"))
setClassUnion("characterOrNULL", c("character", "NULL"))
setClassUnion("numericOrNULL", c("numeric", "NULL"))
setClassUnion("listOrNULL", c("list", "NULL"))
setClass("Colony",
         slots=c(id="character",
                 location="numericOrNULL",
                 queen="PopOrNULL",
                 drones="PopOrNULL",
                 workers="PopOrNULL",
                 virgin_queens="PopOrNULL",
                 pheno="matrix",
                 swarm="logical",
                 split="logical",
                 supersedure="logical",
                 collapse="logical",
                 #rob="logical",
                 production="logical",
                 last_event="character",
                 misc="listOrNULL"
                 ))

#' @describeIn Colony Show colony summary
setMethod("show",
          signature(object = "Colony"),
          function (object){
            cat("An object of class", 
                classLabel(class(object)), "\n")
            cat("Id:", ifelse(!is.null(object@id), object@id, 0),"\n")
            cat("Location:", ifelse(!is.null(object@location), object@location, 0),"\n")
            cat("Queens:", ifelse(!is.null(object@queen), object@queen@nInd, 0),"\n")
            cat("Virgin queens:", ifelse(!is.null(object@virgin_queens), object@virgin_queens@nInd, 0),"\n")
            cat("Drones:", ifelse(!is.null(object@drones), object@drones@nInd, 0),"\n")
            cat("Workers:", ifelse(!is.null(object@workers), object@workers@nInd, 0), "\n")
            cat("Fathers:", ifelse(!is.null(object@queen@misc$fathers), object@queen@misc$fathers@nInd, 0), "\n")
            cat("Events:", paste(if(object@swarm) "swarm", if(object@split) "split", if(object@supersedure) "supersede", if(object@collapse) "collapse"), "\n")
            invisible()
          }
)

#' @title Create new Colony
#' 
#' @description
#' Creates a new \code{\link{Colony}} 
#'
#' @param 
#'
#' @return Returns an object of \code{\link{Colony}}
#' 
#' @examples 
#' 
#' @export


#TODO: do we want to check for any conditions when creating it??????
createColony = function(id = NULL, location = NULL, queen = NULL, drones = NULL, 
                        workers = NULL, virgin_queens = NULL, fathers = NULL, 
                        pheno = NULL, swarm = FALSE, split = FALSE, supersedure =FALSE,
                        collapse = FALSE, #rob = FALSE,
                        production = FALSE,
                        last_event = NULL, yearOfBirth = NULL, misc = NULL,
                        simParam=NULL) { 
  
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  
  if(is.null(id)){
    if(!is.null(queen)){
      id = queen@id
    }
  } 
  

  if (!is.null(queen) & !is.null(fathers)) {
    virgin_queens <- randCross2(females = queen,
                                males = fathers,
                                nCrosses = 1)
  }
  
  if (!is.null(queen)) {
    queen@misc = list(fathers = fathers, yearOfBirth = yearOfBirth)
  }
  
  output = new("Colony",
              id=as.character(id),
              location=location,
              queen=queen,
              drones=drones,
              workers=workers,
              virgin_queens=virgin_queens,
              pheno=matrix(),
                          #ncol=simParam@nTraits),
              swarm=swarm,
              split=split,
              supersedure=supersedure,
              collapse=collapse,
              #rob=rob,
              production=production,
              last_event="new_colony",
              misc=list())


  return(output)
}

#=======================================================================
# createWorkers
# =======================================================================
#' @rdname createWorkers
#' @method createWorkers
#' @title Creates workers of the colony
#' @usage \method{createWorkers}(colony, nWorkers)
#' @description Creates the specified number of workers in the colony
#'       \by mating the current queen and the fathers and adds them in
#'       \ the \code{colony@workers} slot.
#' @param colony Colony class. AlphaSimR Colony object from the \code{createColony(...)} call
#' @param nDrones Numeric. Number of drones to create
#'
#' @example inst/examples/examples_createDrones.R
#'
#' @export

createWorkers = function(colony, nInd){
    if (is.null(colony@queen)) {
      stop("Missing queen!") 
    }
    if (is.null(colony@queen@misc$fathers)) {
      stop("Missing fathers!")
    }

    workerPop = randCross2(females = colony@queen,
                           males = colony@queen@misc$fathers,
                           nCrosses = nInd)
    return(workerPop)
}

#=======================================================================
# createDrones
# =======================================================================
#' @rdname createDrones
#' @method createDrones
#' @title Creates drones of the colony as double haploids
#' @usage \method{createDrones}(colony, nDrones)
#' @description Creates the specified number of drones in the colony
#'       \as double haploids from the current queen and adds them in
#'       \ the \code{colony@drones} slot.
#' @param colony Colony class. AlphaSimR Colony object from the \code{createColony(...)} call
#' @param nDrones Numeric. Number of drones to create
#'
#' @example inst/examples/examples_createDrones.R
#'
#' @export

createDrones = function(colony, nInd){
  if (is.null(colony@queen)) {
    stop("Missing queen!") 
  }
  dronePop = makeDH(pop = colony@queen, nDH = nInd)
  return(dronePop)
}


#=======================================================================
# Create DCA
# =======================================================================
#' @rdname createDCA
#' @method createDCA
#' @title Creates a drone congregation area (DCA) from the list of colonies
#' @usage \method{createDCA}(colony_list, colonyIDs)
#' @description Creates a drone congregation area (DCA) from selected colonies.
#' The function takes the list of all the colonies and a vector of IDs of the selected ones.
#' The function returns a combined population of drones.
#' @seealso \code{\link[??????]{select_colonies}}
#' @param colonies 
#' @param colonyIDs
#'
#' @example inst/examples/examples_createDrones.R
#'
#' @return Single AlphaSim population object of drones from all the selected colonies.
#' @export

createDCA = function(colonies) {
  #dcaColonyList = selectColonies(colonyList, colonyIDs)
  DCA = lapply(X = colonies, FUN = function(z) z@drones)  
  DCA = mergePops(popList = DCA)
  print(paste0("Created a DCA with ", DCA@nInd, " drones."))
  return(DCA)
}


#=======================================================================
# Select Fathers
#NOT SURE WHETHER WE NEED THIS BUT JUST TO NOTE THIS HAS TO BE MADE SOMEWHERE IN THE SIMULATION!
# =======================================================================
#' @rdname SelectFathersfromDCA
#' @method SelectFathersfromDCA
#' @title Selects and stores fathers from DCA to mate with queen 
#' @usage \method{SelectFathersfromDCA}(DCA, nFathers)
#' @description Selects a number individuals from DCA at random to become fathers 
#'              # These fathers will go on to mate with a virgin queen and create a new colony 
#'              # By storing the fathers you can keep track of the pedigree of the colonies
#' @param DCA Poplist. AlphaSimR DCA object from the \code{createDCA(...)} call
#' @param nFathers Numeric. Number of fathers to create 
#'
#' @example inst/examples/examples_SelectFathersfromDCA.R
#' 
#' @export

pullDronesFromDCA = function(DCA, nInd) {
  selectedDronesID = sample(DCA@id, size = nInd, replace = FALSE)
  sel = DCA@id %in% selectedDronesID
  selectedDrones = DCA[sel]
  updatedDCA = DCA[!sel]
  return(list(selectedDrones = selectedDrones, DCA = updatedDCA))
}

#=======================================================================
# Extract the year of birth of the queen
# =======================================================================
extractQueenYOB <- function(colony) {
  return(colony@queen@misc$yearOfBirth)
}

#=======================================================================
# Compute the age of the queen
# =======================================================================
computeQueenAge <- function(colony, currentYear) {
  return(currentYear - colony@queen@misc$yearOfBirth)
}

#=======================================================================
# addWorkers
# =======================================================================
#' @rdname addWorkers
#' @method addWorkers
#' @title 
#' @usage 
#' @description 
#'    THIS DOES NOT REPLACE THE EXISTING ANIMALS!!!
#' @param colony Colony class. AlphaSimR Colony object from the \code{createColony(...)} call
#' @param nWorkersAdd
#'
#' @example inst/examples/examples_addWorkers.R
#' @export
#' 
addWorkers = function(colony, nInd) {
  newWorkers = createWorkers(colony, nInd)
  if (!is.null(colony@workers)) {
    colony@workers = mergePops(list(colony@workers, newWorkers))
  } else {
    colony@workers = newWorkers
  }
  return(colony)
}  

#=======================================================================
# addDrones
# =======================================================================
#' @rdname addDrones
#' @method addDrones
#' @title 
#' @usage 
#' @description 
#'    THIS DOES NOT REPLACE THE EXISTING ANIMALS!!!
#' @param colony Colony class. AlphaSimR Colony object from the \code{createColony(...)} call
#' @param nDronesAdd
#'
#' @example inst/examples/examples_addDrones.R
#' @export
#'a
addDrones <- function(colony, nInd) {
  newDrones = createDrones(colony, nInd)
  if (!is.null(colony@drones)) {
    colony@drones = mergePops(list(colony@drones, newDrones))
  } else {
    colony@drones = newDrones
  }
  return(colony)
}  


#=======================================================================
# Replace workers - MAYBE WE DON?T NEED THIS SINCE WE HAVE createWorker
# =======================================================================
#' @rdname replaceWorkers
#' @method replaceWorkers
#' @title Replaces workers with new workers with new genetic information 
#' @usage \method{replaceWorkers}(colony, nWorkers)
#' @description If new genetic information is added to the colony (in the case of supersedure or a swarm),
#'              workers will be replaced to account for the new queen and fathers
#' @param colony Colony class. AlphaSimR Colony object from the \code{createColony(...)} call
#' @param nWorkers Numeric. OPTIONAL! Default set to current number of workers in the colony.
#'
#' @example inst/examples/examples_replaceWorkers.R
#' @export

replaceWorkers = function(colony, p = 1) {
  nWorkers = colony@workers@nInd
  nWorkersReplaced = round(nWorkers * p)
  if (nWorkersReplaced < nWorkers) {
    nWorkersStay <- nWorkers - nWorkersReplaced
    colony@workers <- c(selectInd(colony@workers, nInd = nWorkersStay, use = "rand"),
                        createWorkers(colony, nWorkers = nWorkersReplaced))
  } else (
    colony@workers = createWorkers(colony, nWorkersReplaced)
  )
  
  return(colony)
}

#=======================================================================
# Replace drones - MAYBE WE DON?T NEED THIS SINCE WE HAVE createDrone
# =======================================================================
#' @rdname replaceDrones
#' @method replaceDrones
#' @title Replaces drone with new drone with new genetic information 
#' @usage \method{replaceDones}(colony, nDrones)
#' @description If new genetic information is added to the colony (in the case of supersedure or a swarm),
#'              drones will be replaced to account for the new queen and fathers
#' @param colony Colony class. AlphaSimR Colony object from the \code{createColony(...)} call
#' @param nDrones Numeric. OPTIONAL! Default set to current number of drones in the colony.
#'
#' @example inst/examples/examples_replaceDrones.R
#' @export

replaceDrones = function(colony, p=1) {
  nDrones = colony@drones@nInd
  nDronesReplaced = round(nDrones * p)
  if (nDronesReplaced < nDrones) {
    nDronesStay <- nDrones - nDronesReplaced
    colony@drones <- c(selectInd(colony@drones, nInd = nDronesStay, use = "rand"),
                        createDrones(colony, nDrones = nDronesReplaced))
  } else (
    colony@drones = createDrones(colony, nDronesReplaced)
  )
  
  return(colony)
}

#=======================================================================
# Extract Individuals from a cast - this could be used for example to extract virgin queens
# NOT SURE WHETHER WE NEED THIS BUT COULD BE USEFUL!!!
# =======================================================================
pullIndFromCaste = function(colony, caste, nInd) {
  if (nInd > slot(colony, caste)@nInd) {
    stop(paste0("Not enough individuals in ", caste, " ! " ,
                nInd, " required, but ", slot(colony, caste)@nInd, " available."))
  }
  pullId = sample(slot(colony, caste)@id, nInd, replace = F)
  pullMatch = slot(colony, caste)@id %in% pullId
  stayMatch = !slot(colony, caste)@id %in% pullId
  
  indPull = slot(colony, caste)[pullMatch]
  indStay = slot(colony, caste)[stayMatch]
  
  slot(colony, caste) = indStay
  return(list(colony = colony, pulledInd = indPull))
}

#=======================================================================
# Cross colony
# =======================================================================
#' @rdname crossColony
#' @method crossColony
#' @title Crosses a colony with a virgin queen to the population of drones.
#' @usage \method{crossColony}(colony, fathers, nWorkers, nDrones)
#' @description Crosses a colony with a virgin queen to the population of drones,
#' \creates workers, drones and a new virgin queen and write them to the corresponding
#' \slots of the colony object.
#' #IF the colony is queenless - select a queen from the virgin queen - if not, mate the current queen!!!
#' @seealso \code{\link[??????]{createColony}}
#' @param colony Colony class. AlphaSimR Colony object from the \code{createColony(...)} call :
#'               INPUT SHOULD BE A COLONY WITH A VIRGIN QUEEN!!!!! 
#' @param fathers Pop Class. Father group taken from Colony object from the \code{createColony(...)} call 
#' @param nWorkers Numeric.Number of workers to create
#' @param nDrones Numeric. Number of drones to create

#'
#' @example inst/examples/examples_crossColony.R
#'
#' @export

crossColony = function(colony, fathers=NULL, nWorkers=0, nDrones=0) {
  if (is.null(colony@virgin_queens)) {
    stop("No virgin queen!")
  }
  
  if (!is.null(colony@queen)) {
    stop("Mated queen present!")
  }
  
  if(all(!is.null(colony@queen@misc$fathers), !is.null(fathers))) {
    stop("Queen already mated!")
  } else if (all(is.null(colony@queen@misc$fathers), is.null(fathers))) {
    stop("Missing fathers!")
  } else if (all(is.null(colony@queen@misc$fathers), !is.null(fathers))) {
    colony@queen@misc$fathers = fathers
  }
  
  colony@queen = selectInd(colony@virgin_queens, nInd = 1, use = "rand")
  
  if (nWorkers != 0) {
    colony@workers = createWorkers(colony, nWorkers)
  }
  if (nDrones != 0) {
    colony@drones = createDrones(colony, nDrones)
  }
  
  colony@virgin_queens = selectInd(colony@workers, nInd = 1, use = "rand")
  
  return(colony)
}


#=======================================================================
# Collapse of the colony
# =======================================================================
collapseColony <- function(colony) {
  colony@collapse <- TRUE
  return(colony)
}

#=======================================================================
# Swarm
# =======================================================================
#' @rdname swarm
#' @method swarm
#' @title Replicates the swarming process and produces two colonies.
#' @usage \method{createSwarm}(colony, perSwarm)
#' @description List. Replicates the swarming of the colony - the process in which
#' a part of the workers leave with the old queen and creates a new colony,
#' while a part of the workers stay with a new queen and the old drones.
#' 1. Compute the number of workers that will leave with the swarm 
#' 2. Compute the number of workers that will stay in the original colony
#' 3. Identify which workers will swarm and which will stay with the original colony
#'       - returns TRUE if workers swarm and FALSE if workers stay
#'       # TODO: give ids of workers that will leave with the swarm
#' 4. Create a new colony entity that represents the colony that stays!!! The swarm is the old colony. 
#' Set the virgin queens for the new colony. All the drones stay at the original location, in the new colony
#' 5. Set new set of workers to the original colony (which is just a subset of the original set)
#'    - The new queen becomes the virgin queen. The new virgin queen will be set by crossColony function 
#'       #colony@queen = selectInd(colony@virgin_queens,  nInd = 1, use = "rand")
#'  6. Change the status of the colony 
#'      #TODO: better names for this but we have to know which one stayed and which one left due to the drones
#'      #TODO: do we want to have some information about the link (i.e. mother_colony=?")

#' @seealso \code{\link[??????]{createColony}}
#' @param colony Colony class. AlphaSimR Colony object from the \code{createColony(...)} call
#' @param perSwarm Numeric. Percentage of colony that will swarm
#'
#' @example inst/examples/examples_swarm.R
#' @return Two colonies, one with the new queen and proportion of workers and
#' one with the old queen and proportion of workers.
#' @export

swarmColony = function(colony, pSwarm, crossVirginQueen = FALSE, fathers = NULL, nWorkers = 0, nDrones = 0, swarmLocation = NULL) {
  if (is.null(colony@virgin_queens)) {
    stop("Virgin queen not present in the colony, cannot swarm")
  }
  if (is.null(colony@drones)) {
    warning("No drones present in the colony!")
  }
  
  nWorkersSwarm = round(colony@workers@nInd * pSwarm, 0)
  nWorkersStay = colony@workers@nInd - nWorkersSwarm
  workersSwarmId = sample(x = colony@workers@id, size = nWorkersSwarm, replace = FALSE) 
  workersStayId = colony@workers@id[!colony@workers@id %in% workersSwarmId]

  newColony = createColony()
  newColony@virgin_queens = selectInd(colony@virgin_queens, 1, use = "rand")
  newColony@workers = colony@workers[workersStayId]
  newColony@drones = colony@drones
  newColony@id = newColony@virgin_queens@id
  newColony@location = colony@location
  
  if (crossVirginQueen) {
    if (is.null(fathers)) {
      stop("No fathers provided, cannot mate the queen!")
    }
    newColony@queen <- newColony@virgin_queens
    newColony@queen@misc$fathers <- fathers
    newColony@workers <- addWorkers(newColony, nWorkers)
    newColony@drones <- addDrones(newColony, nDrones)
    newColony@virgin_queen <- createWorkers(newColony, 1)
  }
  
  swarm = colony
  swarm@workers = colony@workers[workersSwarmId]
  swarm@virgin_queens = selectInd(swarm@workers, nInd = 1, use = "rand")
  swarm@drones = NULL
  swarm@location = swarmLocation


  newColony@last_event = "swarmStay" 
  swarm@last_event = "swarm" 
  
  newColony@swarm = TRUE
  swarm@swarm = TRUE
  newColony@production = FALSE
  swarm@production = FALSE
  
  
  message("Created two colonies.")
  return(list(swarmStay = newColony, swarm = swarm))
}


#=======================================================================
# Supersede
# =======================================================================
#' @rdname supersedeColony
#' @method supersedeColony
#' @title Replicates a supersedure of the colony and replaces the queen with a virgin queen.
#' @usage \method{supersedure}(colony)
#' @description Replicates the process of supersedure, where the
#' queen is replaced by a new virgin queen. The workers and the drones stay
#' in the colony.
#'      # Jana: I DON?T LIKE THIS - now we have a queen that in not mated --> but that should be virgin queen!
#' @seealso \code{\link[??????]{supersedure}}
#' @param colony Colony class. AlphaSimR Colony object from the \code{createColony(...)} call
#'
#' @example inst/examples/examples_supersedure.R
#' @export

supersedeColony = function(colony, crossVirginQueen = FALSE, fathers = NULL, nWorkers = 0, nDrones = 0) {
  colony@queen <- NULL
  colony@virgin_queens = selectInd(colony@virgin_queens,  nInd = 1, use = "rand")
  
  colony@id = colony@virgin_queens@id
  
  if (crossVirginQueen) {
    if (is.null(fathers)) {
      stop("No fathers provided, cannot mate the queen!")
    }
    colony@queen <- colony@virgin_queens
    colony@queen@misc$fathers <- fathers
    colony@workers <- addWorkers(colony, nWorkers)
    colony@drones <- addDrones(colony, nDrones)
    colony@virgin_queens <- createWorkers(colony, 1)
  }
  
  colony@last_event <- "superseded"
  colony@supersedure = TRUE
  colony@production = TRUE
  
  return(colony)
}


#=======================================================================
# Split colony
# =======================================================================
#' @rdname splitColony
#' @method splitColony
#' @title Split the colony in two colonies.
#' @usage \method{splitColony}(colony, per_split)
#' @description Split the colony in two colonies - one with the old queen and
#' a part of workers and drones, and one a part of workers and NO QUEEN!!!!
#' @seealso \code{\link[??????]{splitColony}}
#' @param colony Colony class. AlphaSimR Colony object from the \code{createColony(...)} call
#' @param per_split
#'
#' @example inst/examples/examples_splitColony.R
#' @export
splitColony = function(colony, pSplit, newQueen = NULL, crossVirginQueen = FALSE, fathers = NULL, nWorkers = 0, nDrones = 0, splitLocation = NULL) {
  nWorkersSplit = round(colony@workers@nInd * pSplit, 0)
  noWorkersStay = colony@workers@nInd - nWorkersSplit
  workersSplitId = sample(x = colony@workers@id, size = nWorkersSplit, replace = FALSE)
  workersStayId = colony@workers@id[!colony@workers@id %in% workersSplitId] 
  splitColony = createColony()
  splitColony@workers = colony@workers[workersSplitId]
  splitColony@location = splitLocation

  colony@workers = colony@workers[workersStayId]

  if (!is.null(newQueen)) {
    if (is.null(newQueen@misc$fathers)) {
      splitColony@virgin_queens <- newQueen
      splitColony@id <- splitColony@virgin_queens@id
    }
    if (!is.null(newQueen@misc$fathers)) {
      splitColony@queen <- newQueen
      splitColony@id <- splitColony@queen@id
    }
    
    if (crossVirginQueen) {
      if (is.null(fathers)) {
        stop("No fathers provided, cannot mate the queen!")
      }
      if (!is.null(splitColony@queen)) {
        stop("Queen already mated!")
      }
      splitColony@queen <- splitColony@virgin_queens
      splitColony@queen@misc$fathers <- fathers
      splitColony@workers <- addWorkers(splitColony, nWorkers)
      splitColony@drones <- addDrones(splitColony, nDrones)
      splitColony@virgin_queens <- createWorkers(splitColony, 1)
    }
  }
  
  #Change the status of the colony
  colony@last_event = "splitStay" 
  splitColony@last_event = "split" 
  
  colony@split = TRUE
  splitColony@split = TRUE
  
  colony@production = TRUE
  splitColony@production = FALSE
  
  message("Created two colonies.")
  return(list(splitStay = colony, split = splitColony))
}

#=======================================================================
# Select colonies
# =======================================================================
#' @rdname selectColonies
#' @method selectColonies
#' @title Select the colonies from the colony list based on IDs.
#' @usage \method{selectColonies}(colony_list, colony_ids)
#' @description Select the colonies from the list of all colonies based
#' on colony IDs and return a list of selected colonies.
#' @param colony_list 
#' @param colony_ids
#'
#' @example inst/examples/examples_selectColonies.R
#' @return A list of selected colonies.
#' @export
#' 
selectColonies <- function(colonyList, colonyIDs) {
  selColonyList = colonyList[sapply(colonyList@colonies, FUN = function(x) x@id %in% colonyIDs)]
  return(selColonyList)
}

# setPheno = function() {} # keep this one commented for now, we need some object-oriented magic for this to work on our colony and not to clash with AlphaSimR:::setPheno()





#Things to add in the future/ check list 
#genetic model --> they have either infinitesimal or finite (we don't have to worry about this)
 #simulate queen and worker effects (yet to do)
   #simulate demographic history (in progress)
  #population --> specify the size of the population (#colonies/y)
    #simulate passive population (not actively selected but exchanged queens/drones with the selected)
    #specify the rates of exchange between selected and non-selected (we can control this within the breeding program
    #specify min/max age of queen to produce queens/drones
    #culling ages
  #mating --> specify how queen of the breding and passive populations mate - free mating, AI, breeding stations
    #decide how many drones are involved
    #For mating stations specify - drone producting queens on mating stations, whether they are related, and how secure is the mating station (probability of a drone in the DCa coming from the mating stations hives)
  #selection - rates of queen selection
    #how many sisters from the same sister group (within/across family selection)
    #perform: phenotypic, genotypic, random selection on selection based on EBVs
