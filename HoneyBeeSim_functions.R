#For function storage/ developement
#TODOs

#1) Prevent mother-son mating: probably some mother-son mating in nature as well - so that's fine

# Functions
#1) Write function reset events
#2) Create a function to sample the location for the swarm: later
#3) Create a function to sample locations within a radius: later
#4) Write a setPheno function for the colony: later
#5) Revise last event slot (whether we use/need it)
#7) THink about how to handle the colour
#8) Call setPheno after swarm, split and supersede, crossColony, createColony (follow the AlphaSimR logic of SP parameter vs setPheno)
#9) Assign csd locus - look Laura's code
#10) Make csd function in terms of inbreeding
#11) Add names to colonies in the Colonies (name = colony id)
#12) Consider adding a combine colony function (puts workers from weak into the strong colony)
#16) Think about use = "rand/something" in the level 3 functions (no need for ID then)
#14) All functions should test the class and throw a stop (colony, colonies)
#17) Only set the id when we have the queen!
#18) All add functions should call create functions
#19) Do we just remove provided qeens and fathers in colony (level2) functions? 

# Text
#1) Think about providing informative messages for the functions: Laura
#2) Think of a good names for the swarmed colony (the one that stay)

# Think
#1) Think about removing workers and drones in "instantaneous" functions (opposite to adding them)
#2) Think about replacing phenotypes in the swarm/supersede/split


#TODO for script
#S1) Distribute supersedure/swarming events throughout the year/seasons
#S2) Distribute colony losses events throughout the year/seasons

###############################################################################
# Class Colonies----

#' @title Colonies
#' #'
#' @description
#' The Colonies represents the list of the colonies.
#' It is designed to behave like a list of colonies.
#'
#' @param x a 'Colonies' object
#' @param i index of populations or colonies
#'
#' @slot colonies of \code{\link{Colony-class}} and/or
#' \code{Colonies-class}
#'
#'
#' @export
setClass("Colonies",
         slots=c(colonies="list"))


#' @describeIn Colonies Extract Colonies by index
setMethod("[",
          signature(x = "Colonies"),
          function(x, i){
            x@colonies = x@colonies[i]
            return(x)
          }
)

#' @describeIn Colonies Extract Colony by index
setMethod("[[",
          signature(x = "Colonies"),
          function (x, i){
            return(x@colonies[[i]])
          }
)

#' @describeIn Colonies Extract Colony by index
# setMethod("[[",
#           signature(x = "Colonies", i = "character"),
#           function (x, i){
#             tmp = selectColonies(x, ID = i)
#             return(tmp@colonies)
#           }
# )

#' @describeIn Colonies Combine multiple Coloniess
setMethod("c",
          signature(x = "Colonies"),
          function (x, ...){
            for(y in list(...)){
              if(class(y)=="NULL"){
                # Do nothing
              } else {
                if(class(y)=="Colony"){
                  x@colonies = c(x@colonies, y)
                }else{
                  stopifnot(class(y)=="Colonies")
                  x@colonies = c(x@colonies, y@colonies)
                }
              }
            }
            return(x)
          }
)

#' @describeIn 
setMethod("show",
          signature(object = "Colonies"),
          function (object){
            cat("An object of class",
                classLabel(class(object)), "\n")
            cat("Number of colonies:", nColonies(object), "\n")
            invisible()
          }
)




# Class Colony----

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
#' @slot id ID of the colony.
#' @slot location Numeric, location of the colony (x, y).
#' @slot queen AlphaSimR population object to become the queen of the colony.
#' @slot drones AlphaSimR population object to become the drones of the colony.
#' @slot workers AlphaSimR population object to become the workers of the colony.
#' @slot virgin_queens AlphaSimR individual or population object to become the virgin queen(s) of the colony.
#' @slot pheno A matrix of the phenotypes of the colony
#' @slot swarm Logical, whether the colony has swarmed
#' @slot split Logical, whether the colony has split
#' @slot supersedure Logical, whether the colony has superseded
#' @slot collapse Logical, whether the colony has collapsed
#' @slot production Logical, whether the colony produces hive products
#' @slot last_event Character, the last event of the colony #TODO: WE probably don't need this
#' @slot misc A list, normally empty and exists solely as an open slot available for uses to store extra information about individuals.
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
                 #remnant="",
                 supersedure="logical",
                 collapse="logical",
                 #rob="logical",
                 production="logical",
                 last_event="character",
                 misc="listOrNULL"
                 ))

#' @describeIn Colony Show colony summary
setMethod("c",
          signature(x = "Colony"),
          function (x, ...){
            Colonies = createColonies(x, ...)
            return(Colonies)
          }
)

#' @describeIn 
# setMethod("list",
#           signature(x = "Colony"),
#           function (x, ...){
#             Colonies = createColonies(x, ...)
#             return(Colonies)
#           }
# )

#' @describeIn 
setMethod("show",
          signature(object = "Colony"),
          function (object){
            cat("An object of class",
                classLabel(class(object)), "\n")
            cat("Id:", ifelse(!is.null(object@id), object@id, 0),"\n")
            cat("Location:", ifelse(!is.null(object@location), object@location, 0),"\n")
            cat("Queens:", ifelse(!is.null(object@queen), object@queen@nInd, 0),"\n")
            cat("Virgin queens:", nVirginQueens(object),"\n")
            cat("Drones:", nDrones(object),"\n")
            cat("Workers:", nWorkers(object), "\n")
            cat("Fathers:", nFathers(object), "\n")
            cat("Events:", paste(if(object@swarm) "swarm", if(object@split) "split", 
                                 if(object@supersedure) "supersede", if(object@collapse) "collapse"), "\n")
            cat("Production:", object@production, "\n")
            invisible()
          }
)

###############################################################################
#Level 0 Auxilary Fuctions----

# Number of colonies----
nColonies <- function(colonies) {
  if (!"Colonies" %in% class(colonies)) {
    stop("colonies has to be a class Colonies")
  } 
  return(length(colonies@colonies))
}

# nWorkers----
nWorkers <- function(colony) {
  n = ifelse(!is.null(colony@workers), colony@workers@nInd, 0)
  return(n)
}


# nDrones----
nDrones <- function(colony) {
  n = ifelse(!is.null(colony@drones), colony@drones@nInd, 0)
  return(n)
}


# nVirginQueens----
nVirginQueens <- function(colony) {
  n = ifelse(!is.null(colony@virgin_queens), colony@virgin_queens@nInd, 0)
  return(n)
}

# nFathers----
nFathers <- function(colony) {
  if (is.null(colony@queen)) {
    n = 0
  } else {
    n = ifelse(!is.null(colony@queen@misc$fathers), colony@queen@misc$fathers@nInd, 0)
  }
  return(n)
}

# isQueenMated----
isQueenMated <- function(x) {
  if ("Pop" %in% class(x)) {
    return(!is.null(x@misc$fathers))
  } else if ("Colony" %in% class(x)) {
    if (!is.null(x@queen)) {
      return(!is.null(x@queen@misc$fathers))
    } else {
      return(FALSE)
    }
  }
}

# Extract the year of birth of the queen----
#' @rdname extractQueenYOB
#' @method extractQueenYOB
#' @title Extract the queen's year of birth
#' @usage \method{extractQueenYOB}(colony)
#' @description Extract the year of birth of the queen \code{colony@queen@misc$yearOfBirth} slot
#' @param colony AlphaSimR population object
#' 
#' @example 
#' #'#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=200, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' #Creates colony
#' colony1 = createColony(queen = base[1], fathers = base[2:15])
#' setQueenAge(colony, 1)
#' extractQueenYOB(colony)
#' 
#' @return Integer, the year of birth of the queen.
#' 
#' @export
#' 
extractQueenYOB <- function(colony) {
  return(colony@queen@misc$yearOfBirth)
}



# Compute the age of the queen----

#' @rdname computeQueenAge
#' @method computeQueenAge
#' @title COmputer the queen's age in years
#' @usage \method{computeQueenAge}(colony, year)
#' @description Compute the age of the queen from the \code{colony@queen@misc$yearOfBirth} slot
#' @param colony AlphaSimR population object
#' @param currentYear Integer, current year
#' 
#' @example 
#' #'#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=200, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' #Creates colony
#' colony1 = createColony(queen = base[1], fathers = base[2:15])
#' setQueenAge(colony, year = 1)
#' extractQueenYOB(colony)
#' computerQueenAge(colony, currentYear = 5)
#' 
#' @return Integer, the age of the queen.
#' 
#' @export
#' 
computeQueenAge <- function(x, currentYear) {
  if ("Pop" %in% class(x)) {
    return(currentYear - x@misc$yearOfBirth)
  } else if ("Colony" %in% class(x)) {
    if (!is.null(x@queen)) {
      return(currentYear - x@queen@misc$yearOfBirth)
    }
  }
}




# Get colony IDs from the colonies----
#' @rdname getIDs
#' @method getIDs
#' @title Get the colonies IDs from the colonies
#' @usage \method{getIDs}(colonies)
#' @description Get the colony IDs from the colonies
#' @param colonies 
#'
#' @example 
#' @return Colony IDs
#' @export
#' 
getIDs <- function(colonies, ID) {
  return(sapply(colonies@colonies, FUN = function(x) x@id))
}

###############################################################################
#Level 1 Population Fuctions----

# createFounderDrones----

createFounderDrones <- function(queenPop, nDronesPerQueen) {
  return(makeDH(queenPop, nDH = nDronesPerQueen))
}


# createWorkers----

#' @rdname createWorkers
#' @method createWorkers
#' @title Creates workers of the colony
#' @usage \method{createWorkers}(colony, nInd)
#' @description Creates the specified number of workers in the colony
#'       \by mating the current queen and the fathers in the \code{colony@queen@misc$fathers} slot.
#' @param colony AlphaSimRBee Colony object from the \code{createColony(...)} call
#' @param nInd Integer, number of workers to create
#'
#' @example 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=200, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' #Creates colony
#' colony1 = createColony(queen = base[1], fathers = base[2:15])
#' colony1@workers = createWorkers(colony1, nInd = 1000)
#' 
#' @return AlphaSim population object of created workers.
#' 
#' @export

createWorkers = function(colony, nInd){
  if (is.null(colony@queen)) {
    stop("Missing queen!") 
  }
  if (!isQueenMated(colony)) {
    stop("Missing fathers!")
  }
  
  workerPop = randCross2(females = colony@queen,
                         males = colony@queen@misc$fathers,
                         nCrosses = nInd)
  return(workerPop)
}


# createDrones----

#' @rdname createDrones
#' @method createDrones
#' @title Creates drones of the colony as double haploids
#' @usage \method{createDrones}(colony, nInd)
#' @description Creates the specified number of drones in the colony
#'       \as double haploids from the current queen.
#' @param colony AlphaSimRBee Colony object from the \code{createColony(...)} call
#' @param nInd Integer, the umber of drones to create.
#'
#' @example
#' #' colony1 <- createColony(queen = base[1], fatehrs )#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=200, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' #Creates colony
#' colony1 = createColony(queen = base[1], fathers = base[2:15])
#' colony1@workers = createWorkers(colony1, nInd = 1000)
#' colony1@drones = createDrones(colony, nInd = 200)
#' 
#' @return AlphaSim population object of created drones.
#' @export

createDrones = function(colony, nInd){
  if (is.null(colony@queen)) {
    stop("Missing queen!") 
  }
  dronePop = makeDH(pop = colony@queen, nDH = nInd)
  return(dronePop)
}



# createVirginQueens----

#' @rdname createVirginQueens
#' @method createVirginQueens
#' @title Creates virgin queen of the colony as double haploids
#' @usage \method{createVirginQueens}(colony, nInd)
#' @description Creates the specified number of virgin queen in the colony
#'       \as double haploids from the current queen.
#' @param colony AlphaSimRBee Colony object from the \code{createColony(...)} call
#' @param nInd Integer, the number of virgin queens to create.
#'
#' @example
#' #' colony1 <- createColony(queen = base[1], fatehrs )#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=200, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' #Creates colony
#' colony1 = createColony(queen = base[1], fathers = base[2:15])
#' colony1@workers = createWorkers(colony1, nInd = 1000)
#' colony1@drones = createDrones(colony, nInd = 200)
#' 
#' @return AlphaSim population object of created drones.
#' @export

createVirginQueens = function(colony, nInd){
  return(createWorkers(colony, nInd = nInd))
}



# Create DCA----

#' @rdname createDCA
#' @method createDCA
#' @title Creates a drone congregation area (DCA) from the list of colonies
#' @usage \method{createDCA}(list(colonies))
#' @description Creates a drone congregation area (DCA) from selected colonies.
#' The function takes a vector of the colonies and returns a combined population of drones.
#' @seealso \code{\link[??????]{select_colonies}}
#' @param colonies A list of colonies, each of AlphaSimRBee Colony object
#'
#' @example
#'#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=200, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' #Creates colony
#' colony1 = createColony(queen = base[1], fathers = base[2:15])
#' colony2 = createColony(virgin_queens = base[16])
#' 
#' DCA = createDCA(c(colony1, colony2))
#' 
#' @return Single AlphaSim population object of drones from the provided colonies.
#' @export

createDCA = function(colonies) {
  if ("Colonies" %in% class(colonies)) {
    DCA = lapply(X = colonies@colonies, FUN = function(z) z@drones) 
    DCA = mergePops(popList = DCA)
  } else if ("Colony" %in% class(colonies)) {
    DCA = colonies@drones
  } else {
    stop("Argument colonies must be of class Colonies or Colony")
  }
  
  print(paste0("Created a DCA with ", DCA@nInd, " drones."))
  return(DCA)
}



# pullDronesFromDCA----

#' @rdname pullDronesFromDCA
#' @method pullDronesFromDCA
#' @title Pulls the drones from the DCA
#' @usage \method{pullDronesFromDCA}(DCA, nInd)
#' @description  Pulls a specified number of drones from the DCA and updates the DCA
#' @param DCA AlphaSimR population object created with \code{createDCA(...)} call
#' @param nInd Integer, the number of drones to pull from the DCA
#'
#' @example 
#' #'#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=200, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' #Creates colony
#' colony1 = createColony(queen = base[1], fathers = base[2:15])
#' colony2 = createColony(virgin_queens = base[16])
#' 
#' DCA = createDCA(c(colony1, colony2))
#' fathers = pullDronesTheDCA(DCA, n = 14)
#' 
#' @return A list with two elements. The first element is the AlphaSimR population object of
#' selected drones. The second element is the updated DCA with the selected drones removed.
#' 
#' @export

pullDronesFromDCA = function(DCA, nInd) {
  selectedDronesID = sample(DCA@id, size = nInd, replace = FALSE)
  sel = DCA@id %in% selectedDronesID
  selectedDrones = DCA[sel]
  updatedDCA = DCA[!sel]
  message(paste0("Selected ", nInd, " fathers from DCA"))
  return(list(selectedDrones = selectedDrones, DCA = updatedDCA))
}



# Pull Drone Packages from DCA---- 

pullDronePackagesFromDCA <- function(DCA, n, nAvgFathers) {
  nFathers = rpois(n = n, lambda = nAvgFathers)
  if (sum(nFathers) > DCA@nInd) {
    stop("Not enough drones in the DCA!")
  }
  ret = vector(mode = "list", length = n)
  for (package in 1:n) {
    DCAresult = pullDronesFromDCA(DCA, nInd = nFathers[package])
    DCA = DCAresult$DCA
    ret[[package]] = DCAresult$selectedDrones
  }
  return(ret)
}

# pull Individuals from the caste----

#' @rdname pullIndFromCaste
#' @method pullIndFromCaste
#' @title Pulls a number of individuals from any caste group 
#' @usage \method{pullIndFromCaste}(colony, caste, nInd)
#' @description Pulls and separates a random number of individuals from any caste group. 
#' Two list groups are created, the group of pulled individuals and the colony.
#' 
#'@seealso \code{\link[??????]{pullIndFromCaste}}
#'@param colony Colony class. AlphaSimRBee Colony object from the \code{createColony(...)} call
#'@param caste Character. Replicating the caste class structure present in the hive (queen, drones, workers etc)
#'@nInd Integer. Number of individuals to be pulled from the caste 
#'
#'@example inst/examples/examples_pullIndFromCaste.R
#'@return Two AlphaSim population objects of the colony and the group of pulled individuals.
#'@export 
#'
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

# Cross the virgin queen----

#' @rdname crossVirginQueen
#' @method crossVirginQueen
#' @title Crosses a virgin queen to a group drones
#' @usage \method{crossVirginQueen}(virginQueen, fathers)
#' @description Crosses a virgin queen to a group of drones
#' @param virginQueen AlphaSimR population object
#' @param fathers AlphaSimR population class. 
#' 
#' @example
#' @return AlphaSim population object of a mated colony
#' @export

crossVirginQueen = function(virginQueen, fathers) {
  if (isQueenMated(virginQueen)) {
    stop("The queen is mated already!")
  }
  
  if (is.null(fathers)) {
    stop("Missing fathers!")
  }
  
  if (virginQueen@nInd > 1) {
    stop("#TODO: A function to mate multiple virgin queens at once")
  }
  
  virginQueen@misc$fathers = fathers
  
  return(virginQueen)
}

###############################################################################
#Level 2 Colony Fuctions----


#Create new Colony ----
#' @title Create new Colony
#' 
#' @description
#' Creates a new \code{\link{Colony}}.
#' The function is intended for creating initial colonies from 
#' 'FOUNDERPOP' created by \code{\link{runMacs}}.
#'
#' @param id Character, the ID of the colony, which equals the ID of the queen of not stated otherwise.
#' @param location Numeric, location of the colony (x, y).
#' @param queen AlphaSimR population object to become the queen of the colony.
#' @param drones AlphaSimR population object to become the drones of the colony.
#' @param workers AlphaSimR population object to become the workers of the colony.
#' @param virgin_queens AlphaSimR individual or population object to become the virgin queen(s) of the colony.
#' @param pheno A matrix of the phenotypes of the colony
#' @param swarm Logical, whether the colony has swarmed
#' @param split Logical, whether the colony has split
#' @param supersedure Logical, whether the colony has superseded
#' @param collapse Logical, whether the colony has collapsed
#' @param production Logical, whether the colony produces hive products
#' @param last_event Character, the last event of the colony #TODO: WE probably don't need this
#' @param misc A list, normally empty and exists solely as an open slot available for uses to store extra information about individuals.
#'
#'
#' @return Returns an object of \code{\link{Colony}}
#' 
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=200, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' #Creates colony
#' colony1 = createColony(queen = base[1], fathers = base[2:15])
#' colony2 = createColony(virgin_queens = base[16])
#' 
#' @return AlphaSim Colony object.
#' 
#' @export

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

# Set the queen's Year of Birth----

#' @rdname setQueensYOB
#' @method setQueensYOB
#' @title Set the queen's year of birth
#' @usage \method{setQueenYOB}(colony)
#' @description Set the year of birth of the queen in the \code{colony@queen@misc$yearOfBirth} slot
#' @param colony AlphaSimR population object
#' @param year Integer, the year of the birth of the queen
#' 
#' @example 
#' #'#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=200, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' #Creates colony
#' colony1 = createColony(queen = base[1], fathers = base[2:15])
#' setQueenAge(colony, 1)
#' 
#' @return AlphaSimRBee Colony object
#' 
#' @export

setQueenYOB <- function(x, year) {
  if ("Pop" %in% class(x)) {
    x@misc$yearOfBirth <- year
    return(x)
  } else if ("Colony" %in% class(x)) {
    if (!is.null(x@queen)) {
      x@queen@misc$yearOfBirth <- year
      return(x)
    }
  }
}

# addWorkers----

#' @rdname addWorkers
#' @method addWorkers
#' @title Add workers to the colony
#' @usage \method{addWorkers}(colony, nInd)
#' @description Create workers and store them in the \code{colony@workers} slot. If there is
#' already some workers present in the hive, the function will not overwrite them but instead 
#' combine the newly created and the existing workers. The function returns the updated colony.
#' @param colony AlphaSimRBee Colony object from the \code{createColony(...)} call
#' @param nInd Integer, number of workers to add.
#'
#' @example 
#' Create founder haplotypes
#' founderPop = quickHaplo(nInd=200, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' #Creates colony
#' colony1 = createColony(queen = base[1], fathers = base[2:15])
#' colony1 = addWorkers(colony1, nInd = 2000)
#' 
#' @return Updated AlphaSimRBee Colony object
#' 
#' @export
#' 
addWorkers = function(colony, nInd) {
  newWorkers = createWorkers(colony, nInd)
  if (!is.null(colony@workers)) {
    colony@workers = mergePops(list(colony@workers, newWorkers))
  } else {
    colony@workers = newWorkers
  }
  print(paste0(nInd, " workers added to the colony"))
  return(colony)
}  


# addDrones----

#' @rdname addDrones
#' @method addDrones
#' @title 
#' @usage \method{addDrones}(colony, nInd)
#' @description Create drones and store them in the \code{colony@drones} slot. If there is
#' already some drones present in the hive, the function will not overwrite them but instead 
#' combine the newly created and the existing drones The function returns the updated colony.
#' @param colony AlphaSimRBee Colony object from the \code{createColony(...)} call
#' @param nInd Integer, number of drones to add.
#'
#' @example 
#' Create founder haplotypes
#' founderPop = quickHaplo(nInd=200, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' #Creates colony
#' colony1 = createColony(queen = base[1], fathers = base[2:15])
#' colony1 = addWorkers(colony1, nInd = 2000)
#' colony1 = addDrones(colony, nInd = 100)
#' 
#' @return Updated AlphaSimRBee Colony object
#' 
#'@export

addDrones <- function(colony, nInd) {
  newDrones = createDrones(colony, nInd)
  if (!is.null(colony@drones)) {
    colony@drones = mergePops(list(colony@drones, newDrones))
  } else {
    colony@drones = newDrones
  }
  print(paste0(nInd, " drones added to the colony"))
  return(colony)
}  

# 
# # reQueenColony----
# 
reQueenColony <- function(colony, queen) {
  if(!isQueenMated(queen)) {
    colony@virgin_queens = queen
  } else {
    colony@queen = queen
    colony@id = queen@id
  }
  return(colony)
}

# 
# # Build up colony (add workers and drones)----
# 
buildUpColony = function(colony, nWorkers, nDrones) {
  colony = addWorkers(colony, nInd = (nWorkers - nWorkers(colony)))
  colony = addDrones(colony, nInd = (nDrones - nDrones(colony)))
  colony@production = TRUE
  
  return(colony)
}




# addVirginQueens----

#' @rdname addVirginQueens
#' @method addVirginQueen
#' @title Create additional virgin queens
#' @usage \method{createVirginQueens}(colony, nVirginQueens)
#' @description Creates the specified number of virgin queens in the colony
#'       \by crossing the current queen and the fathers and adds them in
#'       \ the \code{colony@virgin_queens} slot.
#' @param colony AlphaSimRBee Colony object
#' @param nVirginQueens Numeric. Number of virgin queens to create
#'
#' @example 
#'
#' @export

addVirginQueens = function(colony, nVirginQueens){
  if (is.null(colony@queen)) {
    stop("Missing queen!") 
  }
  if (is.null(colony@queen@misc$fathers)) {
    stop("Missing fathers!")
  }
  
  virginQueenPop = randCross2(females = colony@queen,
                              males = colony@queen@misc$fathers,
                              nCrosses = nVirginQueens)
  
  colony@virgin_queens = virginQueenPop
  
  return(colony)
}


# Replace workers ----

#' @rdname replaceWorkers
#' @method replaceWorkers
#' @title Replaces a proportion workers with new workers with new genetic information 
#' @usage \method{replaceWorkers}(colony, p)
#' @description Replace a proportion of workers in the new with new workers from the same queen and same fathers.
#' A user would want to replace a proportion (or all) of workers after swarming and supersedure or
#' due to the short-life span of the workers.#' 
#' @param colony AlphaSimRBee Colony object.
#' @param p Numeric, proportion of workers to be replaced with new ones.
#'
#' @example
#' #' Create founder haplotypes
#' founderPop = quickHaplo(nInd=200, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' #Creates colony
#' colony1 = createColony(queen = base[1], fathers = base[2:15])
#' colony1 = addWorkers(colony1, nInd = 2000)
#' colony1 = replaceWorkers(colony1, p = 0.2)
#' 
#' @return Updated AlphaSimRBee Colony object
#'  
#' @export

replaceWorkers = function(colony, p = 1) {
  nWorkers = colony@workers@nInd
  nWorkersReplaced = round(nWorkers * p)
  if (nWorkersReplaced < nWorkers) {
    nWorkersStay <- nWorkers - nWorkersReplaced
    colony@workers <- c(selectInd(colony@workers, nInd = nWorkersStay, use = "rand"),
                        createWorkers(colony, nInd = nWorkersReplaced))
  } else (
    colony@workers = createWorkers(colony, nWorkersReplaced)
  )
  
  return(colony)
}


# Replace drones ----

#' @rdname replaceDrones
#' @method replaceDrones
#' @title Replaces drone with new drone with new genetic information 
#' @usage \method{replaceWorkers}(colony, p)
#' @description Replace a proportion of drones in the new with new drones from the same queen.
#' A user would want to replace a proportion (or all) of drones after swarming and supersedure or
#' due to the short-life span of the drones.
#' @param colony AlphaSimRBee Colony object.
#' @param p Numeric, proportion of drones to be replaced with new ones.
#'
#' @example
#' #' Create founder haplotypes
#' founderPop = quickHaplo(nInd=200, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' #Creates colony
#' colony1 = createColony(queen = base[1], fathers = base[2:15])
#' colony1 = addDrones(colony1, nInd = 2000)
#' colony1 = replaceDrones(colony1, p = 0.2)
#' 
#' @return Updated AlphaSimRBee Colony object
#'  
#' @export

replaceDrones = function(colony, p=1) {
  nDrones = colony@drones@nInd
  nDronesReplaced = round(nDrones * p)
  if (nDronesReplaced < nDrones) {
    nDronesStay <- nDrones - nDronesReplaced
    colony@drones <- c(selectInd(colony@drones, nInd = nDronesStay, use = "rand"),
                       createDrones(colony, nInd = nDronesReplaced))
  } else (
    colony@drones = createDrones(colony, nDronesReplaced)
  )
  
  return(colony)
}



# removeWorkers----

#' @rdname removeWorkers
#' @method removeWorkers
#' @title Remove selected percentage of workers
#' @usage \method{removeWorkers}(colony, p)
#' @description To decrese the number of workers for example in winter 
#' @param colony Colony class. AlphaSimR Colony object from the \code{createColony(...)} call
#' @param p Numeric. 0<=p>=1 .
#'
#' @example inst/examples/examples_removeWorkers.R
#' @export

removeWorkers = function(colony, p) {
  if ( p > 1) {
    stop("p can not be higher than 1" )
  } else if (p < 0) {
    stop("p can not be less than 0")
  } else if (p == 1) {
    colony@workers = NULL
    warning("All workers removed!")
  } else {
    nWorkers = colony@workers@nInd
    nWorkesNew = round(nWorkers * (1 - p))
    colony@workers = selectInd(colony@workers, nInd = nWorkersNew, use = "rand")
  }
  
  return(colony)
}



# removeDrones----

#' @rdname removeDrones
#' @method removeDrones
#' @title Remove selected percentage of drones
#' @usage \method{removeWorkers}(colony, nWorkers)
#' @description To decrese the number of drones for example in winter 
#' @param colony Colony class. AlphaSimR Colony object from the \code{createColony(...)} call
#' @param p Numeric. 0<=p>=1 .
#'
#' @example inst/examples/examples_removeDrones.R
#' @export

removeDrones = function(colony, p) {
  if ( p > 1) {
    stop("p can not be higher than 1" )
  } else if (p < 0) {
    stop("p can not be less than 0")
  } else if (p == 1) {
    colony@drones = NULL
    warning("All drones removed!")
  } else {
    nDrones = colony@drones@nInd
    nDronesNew = round(nDrones * (1 - p))
    colony@drones = selectInd(colony@drones, nInd = nDronesNew, use = "rand")
  }
  
  return(colony)
}


# Reset events----

#' @rdname resetEvents
#' @method resetEvents
#' @title Reset the swarm, split, supersedure events
#' @usage \method{resetEvents}(colony)
#' @description Reset the slots swarm, split and supersedure to FALSE
#'  
#'@param colony AlphaSimRBee Colony object.
#'
#'@example
#'@return An updated AlphaSimRBee Colony object
#'
#'@export 
#'
resetEvents <- function(colony) {
  colony@swarm = FALSE
  colony@split = FALSE
  colony@supersedure = FALSE
  return(colony)
}

# Cross colony----

#' @rdname crossColony
#' @method crossColony
#' @title Crosses a colony with a virgin queen to a group of fathers pulled from the DCA.
#' @usage \method{crossColony}(colony, fathers, nWorkers, nDrones)
#' @description Crosses a colony with a virgin queen to a group of fathers pulled from the DCA
#' \creates workers, drones and a new virgin queen and write them to the corresponding
#' \slots of the colony object.
#' #IF the colony is queenless - select a queen from the virgin queen - if not, mate the current queen!!!
#' @seealso \code{\link[??????]{createColony}}
#' @param colony AlphaSimRBee Colony object with a non-mated virgin queen
#' @param fathers Pop Class. Father group pulled from the DCA. 
#' @param nWorkers Integer.Number of workers to create
#' @param nDrones Integer. Number of drones to create

#'
#' @example inst/examples/examples_crossColony.R
#' @return Single AlphaSim population object of a mated colony
#' @export

crossColony = function(colony, fathers=NULL, nWorkers=0, nDrones=0) {
  if (is.null(colony@virgin_queens)) {
    stop("No virgin queen!")
  }
  
  if (!is.null(colony@queen)) {
    stop("Mated queen present!")
  }
  
  if (is.null(fathers)) {
    stop("Missing fathers!")
  }
  
  colony@queen = selectInd(colony@virgin_queens, nInd = 1, use = "rand")
  colony@id = colony@queen@id
  colony@queen@misc$fathers = fathers
  
  if (nWorkers != 0) {
    colony@workers = createWorkers(colony, nWorkers)
  }
  if (nDrones != 0) {
    colony@drones = createDrones(colony, nDrones)
  }
  
  colony@virgin_queens = selectInd(colony@workers, nInd = 1, use = "rand")
  
  return(colony)
}



# Collapse of the colony ----

#' @rdname collapseColony 
#' @method collapseColony 
#' @title Replicates colony collapse
#' @usage \method{collapseColony}(colony)
#' @description Replicates the collapse of a colony. This can be due to winter losses, disease or other factors.
#'  
#' @seealso \code{\link[??????]{collapseColony}}
#' @param colony Colony class. AlphaSimRBee Colony object from the \code{createColony(...)} call
#' 
#' @example inst/examples/examples_collapseColony.R
#' @return Single AlphaSim population object of collapsed colony
#' @export
#' 
collapseColony <- function(colony) {
  colony@collapse <- TRUE
  return(colony)
}


# Swarm colony ----

#' @rdname swarmColony 
#' @method swarmColony 
#' @title Replicates the swarming process and produces two colonies.
#' @usage \method{swarmColony}(colony, pSwarm, crossVirginQueen. fathers, nWorkers, nDrones, swarmLocation)
#' @description List. Replicates the swarming of the colony - the process in which
#' a part of the workers leave with the old queen and creates a new colony (the swarm),
#' while a part of the workers stay with a new queen and the old drones.
#' The swarming colony contains the old mated queen,
#'  a percentage (pSwarm) of the original colonies workers, no drones and a virgin queen is created from the worker population. 
#'  A new location must be given to the new swarm colony. 
#'  The colony that stays contains the remaining workers and drones. A virgin queen is selected from the workers and mated if fathers are present. 

#' @seealso \code{\link[??????]{createColony}}
#' @param colony Colony class. AlphaSimRBee Colony object from the \code{createColony(...)} call
#' @param pSwarm Integer. Percentage of colony that will swarm
#' @param crossVirginQueen Logical. Whether a virgin queen is to be mated 
#' @param fathers AlphaSimR population object. Number of fathers pulled from the DCA
#' @param pWorkers Numeric, proportion of workers that are replaced with the workers from the new queen in the remnant colony
#' @param pDrones Numeric, proportion of drones that are replaced with the drones from the new queen in the remnant colony
#' @param swarm Location Integer. X,Y coordinates of newly made swarmed hive
#'
#' @example inst/examples/examples_swarm.R
#' @return Two colonies, one with the new queen and proportion of workers and
#' one with the old queen and proportion of workers.
#' @export

swarmColony = function(colony, pSwarm = 0.5, crossVirginQueen = FALSE, fathers = NULL, 
                       pWorkers = 1, pDrones = 1, swarmLocation = NULL) {
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
  
  remnantColony = createColony()
  remnantColony@virgin_queens = selectInd(colony@virgin_queens, 1, use = "rand")
  remnantColony@workers = colony@workers[workersStayId]
  remnantColony@drones = colony@drones
  remnantColony@location = colony@location
  
  if (crossVirginQueen) {
    if (is.null(fathers)) {
      stop("No fathers provided, cannot mate the queen!")
    }
    remnantColony@queen <- remnantColony@virgin_queens
    remnantColony@id = remnantColony@queen@id
    remnantColony@queen@misc$fathers <- fathers
    remnantColony <- replaceWorkers(remnantColony, pWorkers)
    remnantColony <- replaceDrones(remnantColony, pDrones)
    remnantColony@virgin_queens <- createWorkers(remnantColony, 1)
  }
  
  swarm = colony
  swarm@workers = colony@workers[workersSwarmId]
  swarm@virgin_queens = selectInd(swarm@workers, nInd = 1, use = "rand")
  swarm@drones = NULL
  swarm@location = swarmLocation
  
  
  remnantColony@last_event = "remnant" 
  swarm@last_event = "swarm" 
  
  remnantColony@swarm = TRUE
  swarm@swarm = TRUE
  remnantColony@production = FALSE
  swarm@production = FALSE
  
  
  message("Created two colonies.")
  
  return(list(remnant = remnantColony, swarm = swarm))
}



# Supersede colony ----

#' @rdname supersedeColony
#' @method supersedeColony
#' @title Replicates a supersedure of the colony and replaces the queen with a virgin queen.
#' @usage \method{supersedureColony}(colony, crossVirginQueen, fathers, nWorkers, nDrones)
#' @description Replicates the process of supersedure, where the
#' queen is replaced by a new virgin queen. The workers and the drones stay
#' in the colony. If no fathers are present, mating of the virgin queen does not occur. 
#' @seealso \code{\link[??????]{supersedure}}
#' @param colony Colony class. AlphaSimRBee Colony object from the \code{createColony(...)} call
#' @param crossVirginQueen Logical. Whether a virgin queen is to be mated 
#' @param fathers AlphaSimR population object. Number of fathers pulled from the DCA
#' @param pWorkers Numeric, proportion of workers that are replaced with the workers from the new queen
#' @param pDrones Numeric, proportion of drones that are replaced with the drones from the new queen
#'
#' @example inst/examples/examples_supersedeColony.R
#' @return Single AlphaSim population object of superseded colony 
#' @export

supersedeColony = function(colony, crossVirginQueen = FALSE, fathers = NULL, 
                           pWorkers = 1, pDrones = 1) {
  if (is.null(colony@queen)) {
    stop("No queen present in the colony!")
  }
  colony@queen <- NULL
  colony@virgin_queens = selectInd(colony@virgin_queens,  nInd = 1, use = "rand")
  
  if (crossVirginQueen) {
    if (is.null(fathers)) {
      stop("No fathers provided, cannot mate the queen!")
    }
    colony@queen <- colony@virgin_queens
    colony@queen@misc$fathers <- fathers
    colony <- replaceWorkers(colony, pWorkers)
    colony <- replaceDrones(colony, pDrones)
    colony@virgin_queens <- createWorkers(colony, 1)
  }
  
  colony@last_event <- "superseded"
  colony@supersedure = TRUE
  colony@production = TRUE
  
  return(colony)
}



# Split colony----

#' @rdname splitColony
#' @method splitColony
#' @title Split the colony in two colonies.
#' @usage \method{splitColony}(colony, pSplit, newQueen, crossVirginQueen, fathers, nWorkers, nDrones, splitLocation)
#' @description Spit the colony into two new colonies to prevent swarming (in managed populations) 
#' - one colony is with the old queen and a part of the workers and drones (this is the remaining colony)
#' - the split colony is taken to a new location with part of the workers. 
#'  A new mated queen can be introduced to the split colony. 
#'  If no new queen is introduced, a virgin queen must be present to mate with fathers from DCA and continue colony  
#'  #TODO: Split does not produce drones
#' @seealso \code{\link[??????]{splitColony}}
#' @param colony Colony class. AlphaSimRBee Colony object from the \code{createColony(...)} call
#' @param pSplit Integer. Percentage of hive to split 
#' @param newQueen AlphaSimR population object. A new mated queen is brought into the colony from other source 
#' @param crossVirginQueen Logical. If no mated queen is introduced, a virgin queen must be present to mate and continue colony 
#' @param fathers AlphaSimR population object. Number of fathers pulled from the DCA
#' @param pWorkers Numeric, proportion of workers that are replaced with the workers from the new queen in the split
#' @param splitLocation Integer. X,Y coordinates of newly made split hive 
#'
#' @example inst/examples/examples_splitColony.R
#' @return Two AlphaSim population objects of the split colony and the remaining colony 
#' @export
splitColony = function(colony, pSplit = 0.30, newQueen = NULL, crossVirginQueen = FALSE, fathers = NULL, 
                       pWorkers = 1, splitLocation = NULL) {
  nWorkersSplit = round(colony@workers@nInd * pSplit, 0)
  noWorkersStay = colony@workers@nInd - nWorkersSplit
  workersSplitId = sample(x = colony@workers@id, size = nWorkersSplit, replace = FALSE)
  workersStayId = colony@workers@id[!colony@workers@id %in% workersSplitId] 
  splitColony = createColony()
  splitColony@workers = colony@workers[workersSplitId]
  splitColony@location = splitLocation
  
  colony@workers = colony@workers[workersStayId]
  
  if (!is.null(newQueen)) {
    if (!isQueenMated(newQueen)) {
      splitColony@virgin_queens <- newQueen
    }
    if (isQueenMated(newQueen)) {
      splitColony@queen <- newQueen
      splitColony@id <- splitColony@queen@id
    }
    
    if (crossVirginQueen) {
      if (is.null(fathers)) {
        stop("No fathers provided, cannot mate the queen!")
      }
      splitColony@queen <- splitColony@virgin_queens
      if (isQueenMated(splitColony)) {
        stop("Queen already mated!")
      }
      splitColony@queen@misc$fathers <- fathers
    }
  }
  
  if (!is.null(splitColony@queen)) {
    splitColony <- replaceWorkers(splitColony, pWorkers)
    splitColony@virgin_queens <- createWorkers(splitColony, 1)
    
  }
  
  #Change the status of the colony
  colony@last_event = "remnant" 
  splitColony@last_event = "split" 
  
  colony@split = TRUE
  splitColony@split = TRUE
  
  colony@production = TRUE
  splitColony@production = FALSE
  
  message("Created two colonies.")
  return(list(remnant = colony, split = splitColony))
}

###############################################################################
#Level 3 Colonies Functions----


#Create  Colonies ----
#' @title Create Colonies
#'
#' @description
#' Creates a new \code{\link{Colonies-class}} from one or more
#' \code{\link{Colony-class}} and/or \code{\link{Colonies-class}}
#' objects.
#'
#' @param ... one or more \code{\link{Colony-class}} and/or
#' \code{\link{Colonies-class}} objects.
#'
#' @return Returns an empty object of \code{\link{Colonies-class}}
#'
#' @examples
#' @export

createColonies = function(..., n = NULL){
  if (is.null(n)) {
    input = list(...)
    class = sapply(input, "class")
    stopifnot(all(class=="Colony" | class=="Colonies") | all(class=="NULL"))
    output = new("Colonies", colonies=input)    
  } else {
    output = new("Colonies", colonies=vector(mode = "list", length = n))
  }
  
  return(output)
}

#Add colony to Colonies----
#' @title Add colony to the Colonies
#'
#' @description
#'
#' @param
#'
#' @return
#'
#' @examples
#' @export

addColonyToTheColonies= function(colony, Colonies){
  if (class(colony) != "Colony") {
    message("The colony parameter is not a Colony object.")
  }
  Colonies@colonies = append(Colonies@colonies, list(colony))
  return(Colonies)
}

# Select colonies ----

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
selectColonies <- function(colonies, ID) {
  return(colonies[sapply(colonies@colonies, FUN = function(x) x@id %in% ID)])
}



# Pull colonies ----

#' @rdname pullColonies
#' @method pullColonies
#' @title Pull the colonies from the colony list based on IDs.
#' @usage \method{pullColonies}(colonies, colonyIds)
#' @description Pull the colonies from the list of all colonies based
#' on colony IDs and return two lists: a list of selected colonies and 
#' updated original colonies
#' @param colonies 
#' @param colonyIDs
#'
#' @example 
#' @return Two lists: a list of selected colonies and an updated inpute colonies
#' @export
#' 
pullColonies <- function(colonies, ID) {
  pulledColonies <- selectColonies(colonies, ID)
  remainingColonies <- removeColonies(colonies, ID)
  return(list(pulledColonies = pulledColonies, remainingColonies = remainingColonies))
}


# Remove colonies ----

#' @rdname removeColonies
#' @method removeColonies
#' @title Remove the colonies from the colony list based on IDs.
#' @usage \method{removeColonies}(colony_list, colony_ids)
#' @description Remove the colonies from the list of all colonies based
#' on colony IDs and return a list of remaining colonies.
#' @param colony_list 
#' @param colony_ids
#'
#' @example 
#' @return A list of remaining colonies.
#' @export
#' 
removeColonies <- function(colonies, ID) {
  return(colonies[!sapply(colonies@colonies, FUN = function(x) x@id %in% ID)])
}


# Create multiple Virgin colonies ----

createMultipleVirginColonies = function(founderPop, nColonies) {
  ret = createColonies(n = nColonies)
  virginQueens = selectInd(founderPop, nInd = nColonies, use = "rand")
  for (colony in 1:nColonies) {
    ret@colonies[[colony]] = createColony(virgin_queens = virginQueens[colony])
  }
  return(ret)
}


# Create multiple mated colonies ----

createMultipleMatedColonies = function(founderPop, nColonies, nAvgFathers) {
  ret = createColonies(n = nColonies)
  queensID = sample(founderPop@id, size = nColonies, replace = FALSE)
  queenMatch = founderPop@id[founderPop@id %in% queensID]
  queens = founderPop[queenMatch]
  DPQMatch = founderPop@id[!founderPop@id %in% queensID]
  DPQs = founderPop[DPQMatch]
  DCA =  createFounderDrones(queenPop = DPQs, nDronesPerQueen = 10)
  fatherPackages = pullDronePackagesFromDCA(DCA, n = nColonies,nAvgFathers = nAvgFathers)
  for (colony in 1:nColonies) {
    ret@colonies[[colony]] = createColony(queen = queens[colony], 
                                          fathers = fatherPackages[[colony]])
  }
  
  return(ret)
}


# 
# # Build up colonies (add workers and drones)----
# 
buildUpColonies = function(colonies, nWorkers, nDrones) {
  nCol = nColonies(colonies)
  for (colony in 1:nCol) {
    colonies@colonies[[colony]] = buildUpColony(colony = colonies[[colony]], 
                                                nWorkers = nWorkers,
                                                nDrones = nDrones) 
  }
  return(colonies)
}

# 
# # reQueen Colonies----
# 
reQueenColonies = function(colonies, queens) {
  nCol = nColonies(colonies)
  if (nInd(queens) < nCol) {
    stop("Not enough queens!")
  }
  for (colony in 1:nCol) {
    colonies@colonies[[colony]] = reQueenColony(colony = colonies[[colony]], 
                                                queen = queens[colony]) 
  }
  return(colonies)
}





# # Collapse the colonies----
# 
# collapseColonies <- function(colonies, ID) {
#   return(removeColonies(colonies, ID))
# }




# Supersede the colonies----

supersedeColonies <- function(colonies
                              #, crossVirginQueen = FALSE, fathers = NULL, pWorkers = 1, pDrones = 1
) {
  nColonies = nColonies(colonies)
  # if (length(fathers) < nCol) {
  #   stop("Not enought fathers)
  # }
  for (colony in 1:nCol) {
    colonies@colonies[[colony]] = supersedeColony(colonies[[colony]])
  }
  return(colonies)
  
}




# Swarm the colonies----

swarmColonies <- function(colonies
                          #, crossVirginQueen = FALSE, fathers = NULL, pWorkers = 1, pDrones = 1
) {
  nCol = nColonies(colonies)
  # if (length(fathers) < nCol) {
  #   stop("Not enought fathers)
  # }
  ret = list(swarms = createColonies(n = nCol), remnants = createColonies(n = nCol))
  
  for (colony in 1:nCol) {
    tmp = swarmColony(colonies[[colony]])
    ret$swarms@colonies[[colony]] = tmp$swarm
    ret$remnants@colonies[[colony]] = tmp$remnant
  }
  return(ret)
  
}


# Split the colonies----

splitColonies <- function(colonies
                          #, crossVirginQueen = FALSE, fathers = NULL, pWorkers = 1, pDrones = 1
) {
  nCol = nColonies(colonies)
  # if (length(fathers) < nCol) {
  #   stop("Not enought fathers)
  # }
  ret = list(splits = createColonies(n = nCol), remnants = createColonies(n = nCol))
  
  for (colony in 1:nCol) {
    tmp = splitColony(colonies[[colony]])
    ret$splits@colonies[[colony]] = tmp$split
    ret$remnants@colonies[[colony]] = tmp$remnant
  }
  return(ret)
  
}

# Cross the colonies----

crossColonies <- function(colonies, DCA, nAvgFathers
                          #, crossVirginQueen = FALSE, fathers = NULL, pWorkers = 1, pDrones = 1
) {
  nCol = nColonies(colonies)
  # if (length(fathers) < nCol) {
  #   stop("Not enought fathers)
  # }
  ret = createColonies(n = nCol)
  nFathers = rpois(n = nCol, lambda = nAvgFathers)
  fatherPackages = pullDronePackagesFromDCA(DCA, n = nCol, nAvgFathers = nAvgFathers)
  
  for (colony in 1:nCol) {
    ret@colonies[[colony]] <- crossColony(colonies[[colony]],
                                          fathers = fatherPackages[[colony]]
    )
  }
  return(ret)
  
}

###############################################################################



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
