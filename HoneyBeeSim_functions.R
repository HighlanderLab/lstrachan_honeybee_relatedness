#For function storage/ developement

#-----------------------------------------------------------------------
# Class colony list
#-----------------------------------------------------------------------
#' @title ColonyList
#' #'
#' @description
#' The ColonyList represents the list of the colonies.
#' It is designed to behave like a list of colonies.
#'
#' @param x a 'ColonyList' object
#' @param i index of populations or colonies
#'
#' @slot colony list of \code{\link{Colony-class}} and/or
#' \code{ColonyList-class}
#'
#'
#' @export
setClass("ColonyList",
         slots=c(colonies="list"))


#' @describeIn ColonyList Extract ColonyList by index
setMethod("[",
          signature(x = "ColonyList"),
          function(x, i){
            x@colonies = x@colonies[i]
            return(x)
          }
)

#' @describeIn ColonyList Extract Colony by index
setMethod("[[",
          signature(x = "ColonyList"),
          function (x, i){
            return(x@colonies[[i]])
          }
)

#' @describeIn ColonyLIst Combine multiple ColonyLists
setMethod("c",
          signature(x = "ColonyList"),
          function (x, ...){
            for(y in list(...)){
              if(class(y)=="NULL"){
                # Do nothing
              }else{
                if(class(y)=="Colony"){
                  x@colonies = c(x@colonies, y)
                }else{
                  stopifnot(class(y)=="MegaPop")
                  x@colonies = c(x@colonies, y@colonies)
                }
              }
            }
            return(x)
          }
)


#' @title Create new ColonyList
#'
#' @description
#' Creates a new \code{\link{ColonyList-class}} from one or more
#' \code{\link{Colony-class}} and/or \code{\link{ColonyList-class}}
#' objects.
#'
#' @param ... one or more \code{\link{Colony-class}} and/or
#' \code{\link{ColonyList-class}} objects.
#'
#' @return Returns an empty object of \code{\link{ColonyList-class}}
#'
#' @examples
#' @export

createColonyList = function(...){
  input = list(...)
  class = sapply(input, "class")
  stopifnot(all(class=="Colony" | class=="ColonyList") | all(class=="NULL"))
  output = new("ColonyList", colonies=input)
  return(output)
}

#' @title Add colony to the ColonyList
#'
#' @description
#'
#' @param 
#'
#' @return 
#'
#' @examples
#' @export

addColonyToTheColonyList= function(colony, colonyList){
  if (class(colony) != "Colony") {
    message("The colony parameter is not a Colony object.")
  }
  colonyList@colonies = append(colonyList@colonies, list(colony))
  return(colonyList)
}


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
setClass("Colony",
         slots=c(id="character",
                 location="characterOrNULL",
                 queen="PopOrNULL",
                 drones="PopOrNULL",
                 workers="PopOrNULL",
                 virgin_queens="PopOrNULL",
                 fathers="PopOrNULL",
                 pheno="matrix",
                 last_event="character"
                 ))

#' @describeIn Colony Show colony summary
setMethod("show",
          signature(object = "Colony"),
          function (object){
            cat("An object of class", 
                classLabel(class(object)), "\n")
            cat("Id:", ifelse(!is.null(object@id), object@id, 0),"\n")
            cat("Location:", ifelse(!is.null(object@location), object@location, 0),"\n")
            cat("No queens:", ifelse(!is.null(object@queen), object@queen@nInd, 0),"\n")
            cat("No virgin queens:", ifelse(!is.null(object@virgin_queens), object@virgin_queens@nInd, 0),"\n")
            cat("No drones:", ifelse(!is.null(object@drones), object@drones@nInd, 0),"\n")
            cat("No workers:", ifelse(!is.null(object@workers), object@workers@nInd, 0), "\n")
            cat("No fathers:", ifelse(!is.null(object@fathers), object@fathers@nInd, 0), "\n")
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
                        pheno = NULL, last_event = NULL) { #, simParam = NULL, ...
  
  # if(is.null(simParam)){
  #   simParam = get("SP",envir=.GlobalEnv)
  # }
  
  if(is.null(id)){
    if(!is.null(queen)){
      id = queen@id
    }
  } 
  
  if(is.null(location)){
    if(!is.null(id)){
      location = id
    }
  }  
  
 output = new("Colony",
              id=as.character(id),
              location=as.character(location),
              queen=queen,
              drones=drones,
              workers=workers,
              virgin_queens=virgin_queens,
              fathers=fathers,
              pheno=matrix(),
                          #nrow=,
                          #ncol=simParam@nTraits),
              last_event="new_colony") # misc=vector("list",rawPop@nInd)


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

createWorkers = function(colony, nWorkers){
    if (is.null(colony@queen)) {
      stop("Missing queen!") 
    }
    if (is.null(colony@fathers)) {
      stop("Missing fathers!")
    }

    workerPop = randCross2(females = colony@queen,
                           males = colony@fathers,
                           nCrosses = nWorkers)
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

createDrones = function(colony, nDrones){
  if (is.null(colony@queen)) {
    stop("Missing queen!") 
  }
  dronePop = makeDH(pop = colony@queen, nDH = nDrones)
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
#' @param colony_list
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

selectFathersFromDCA = function(DCA, nFathers) {
  return(selectInd(DCA, nInd = nFathers, use = "rand"))
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
  if (all(is.null(colony@virgin_queens), is.null(colony@queen))) {
    stop("No queen or virgin queen!")
  }
  
  if(all(!is.null(colony@fathers), !is.null(fathers))) {
    stop("Fathers already present in the colony!")
  } else if (all(is.null(colony@fathers), is.null(fathers))) {
    stop("Missing fathers!")
  } else if (all(is.null(colony@fathers), !is.null(fathers))) {
    colony@fathers = fathers
  }
  
  if (is.null(colony@queen)) {
    colony@queen = selectInd(colony@virgin_queens, nInd = 1, use = "rand")
   }

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
# Extract Individuals from a cast - this could be used for example to extract virgin queens
#NOT SURE WHETHER WE NEED THIS BUT COULD BE USEFUL!!!
# =======================================================================
extractIndFromCast = function(colony, cast, nInd) {
  if (nInd > slot(colony, cast)@nInd) {
    stop(paste0("Not enough individuals in ", cast, " ! " ,
                nInd, " required, but ", slot(colony, cast)@nInd, " available."))
  }
  
  return(selectInd(slot(colony, cast), nInd = nInd, use = "rand"))
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

swarmColony = function(colony, perSwarm) {
  noWorkersSwarm = round(colony@workers@nInd * perSwarm, 0)
  noWorkersStay = colony@workers@nInd - noWorkersSwarm
  workersSwarmId = sample(x = colony@workers@id, size = noWorkersSwarm, replace = FALSE) 
  workersStayId = colony@workers@id[!colony@workers@id %in% workersSwarmId]

  newColony = createColony()
  newColony@virgin_queens = selectInd(colony@virgin_queens, 1, use = "rand")
  newColony@workers = colony@workers[workersStayId]
  newColony@drones = colony@drones
  newColony@id = newColony@virgin_queens@id
  newColony@location = colony@location
  
  swarm = colony
  swarm@workers = colony@workers[workersSwarmId]
  swarm@virgin_queens = selectInd(swarm@workers, nInd = 1, use = "rand")
  swarm@drones = NULL
  swarm@location = NULL



  newColony@last_event = "swarmStay" 
  swarm@last_event = "swarmLeave" 
  
  message("Created two colonies. The location for the swarm not set. ")
  return(list(newColony = newColony, swarm = swarm))
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

supersedeColony = function(colony) {
  colony@queen = selectInd(colony@virgin_queens,  nInd = 1, use = "rand")
  colony@virgin_queens = selectInd(swarm@workers, nInd = 1, use = "rand")
  colony@id = colony@queen@id
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
splitColony = function(colony, perSplit) {
  # Compute the number of workers that will leave with the swarm (per_swarm is the % of workers that will swarm)
  noWorkersSplit = round(colony@workers@nInd * perSplit, 0)
  # Compute the number of workers that will stay in the original colony (better to do it this way due to rounding issues)
  noWorkersStay = colony@workers@nInd - noWorkersSplit
  # Which workers are taken
  workersSplitId = sample(x = colony@workers@id, size = noWorkersSplit, replace = FALSE) # TODO gives ids of workers that will leave with the swarm
  # Which workers stay
  workersStayId = colony@workers@id[!colony@workers@id %in% workersSplitId] # tells which workers will leave (TRUE) and which won't (FALSE)
  # Create a new colony entity that represents the swarm and set its queen and workers, all the drones stay in the original colony
  splitColony = createColony()
  splitColony@workers = colony@workers[workersSplitId]

  # Set new set of workers to the original colony (which is just a subset of the original set)
  colony@workers = colony@workers[workersStayId]

  #Change the status of the colony
  colony@last_event = "splitStay" #TODO: better names for this but we have to know which one stayed and which one left due to the drones
  splitColony@last_event = "splitLeave" #TODO: do we want to have some information about the link (i.e. mother_colony=?")
  message("Created two colonies. The new, splitColony, is without a queen or a virgin queen.")
  return(list(colony = colony, splitColony = splitColony))
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
addWorkers = function(colony, nWorkersAdd) {
  newWorkers = createWorkers(colony, nWorkersAdd)
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
addDrones <- function(colony, nDronesAdd) {
  newDrones = createDrones(colony, nDronesAdd)
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

replaceWorkers = function(colony, nWorkers=colony@workers@nInd) {
  colony@workers = createWorkers(colony, nWorkers)
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

replaceDrones = function(colony, nDrones=colony@drones@nInd) {
  colony@drones = createDrones(colony, nDrones)
  return(colony)
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
