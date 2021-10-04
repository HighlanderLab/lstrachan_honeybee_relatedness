#For function storage/ developement 

#CreateColony function - to create a new colony from existing (not only base population). Arguments- list of nodes within the hive (to be populated)
createColony = function(id = NULL, location = NULL, queen = NULL, drones = NULL, workers = NULL, virgin_queens = NULL, pheno = NULL, fathers = NULL, last_event = NULL) {
  colony = vector(mode = "list",  length = 9)
  names(colony) = c("id", "location", "queen", "drones", "workers", "virgin_queens", "pheno", "fathers", "last_event")
  if (!is.null(id)) {
    colony$id = id
  }
  if (!is.null(location)) {
    colony$location = location
  }
  if (!is.null(queen)) {
    colony$queen = queen
  }
  if (!is.null(drones)) {
    colony$drones = drones
  }
  if (!is.null(workers)) {
    colony$workers = workers
  }
  if (!is.null(virgin_queens)) {
    colony$virgin_queens = virgin_queens
  }
  if (!is.null(pheno)) {
    colony$pheno = pheno
  }
  if (!is.null(fathers)) {
    colony$fathers = fathers
  }
  if (!is.null(last_event)) {
    colony$last_event = "new_colony"
  }
  return(colony)
} 
# we will likely need queen age too - but that should go into colony$queen$misc slot!
#can also look at hive "strength" based on number of colony workers 

#Drone creation (making them haplodiploid)
createDrones = function(colony, nDrones){
  colony$drones = makeDH(pop = colony$queen, nDH = nDrones)
}

#Drone congregation area of the base population made in colony_list 
createDCA = function(colony_list = NULL, colonyIDs = NULL) { # Would it be better
  dca_colony_list = select_colonies(colony_list, colonyIDs)
  DCA = lapply(X = dca_colony_list, FUN = function(z) z$drones)  
  DCA = mergePops(popList = DCA)
  
  return(popList = DCA)
}

#JANA: INPUT SHOULD BE A COLONY WITH A VIRGIN QUEEN!!!!!
crossColony = function(colony, drone_pop, nFathers, nWorkers_created, nDrones_created) {
  #TODO: do we mate the virgin queen or set the new before this???????
  #The drone pop can be created with the createDCA - but that is a separate step
  #1) Write the fathers into the colony$fathers slot
  colony$fathers = drone_pop[sample(drone_pop@id, nFathers, replace = FALSE)]
  #2) Create the workers and write them in the worker slot
  colony$queen = selectInd(colony$virgin_queens,  nInd = 1, use = "rand")
  colony$workers = randCross2(females = colony$queen, males = drone_pop, nCrosses = nWorkers_created)
  
  colony$drones = createDrones(colony, nDrones_created)
  colony$virgin_queen = selectInd(colony$workers, nInd = 1, use = "rand") #Jana: edited - simpler this way
  
}


# Generate a swarm and modify the existing colony  # colony - list, colony
# little honey
# younger queens swarm less
# build-up after swarm depends on season - should this be an extra function / or will this be in country scripts?
createSwarm = function(colony, per_swarm) {
  # Compute the number of workers that will leave with the swarm (per_swarm is the % of workers that will swarm)
  noWorkersSwarm = round(colony$workers$nInd * per_swarm, 0)
  # Compute the number of workers that will stay in the original colony (better to do it this way due to rounding issues)
  noWorkersStay = colony$workers$nInd - noWorkersSwarm
  # Which workers swarm
  workersSwarmId = sample(x = colony$workers@id, size = noWorkersSwarm, replace = FALSE) # TODO gives ids of workers that will leave with the swarm
  # Which workers stay
  workersStayId = colony$workers@id %in% sel_workers # tells which workers will leave (TRUE) and which won't (FALSE)
  # Create a new colony entity that represents the swarm and set its queen and workers, all the drones stay in the original colony
  swarm = createColony()
  swarm$queen = colony$queen 
  swarm$workers = colony$workers[workersSwarmId]
  swarm$virgin_queen = selectInd(swarm$workers, nInd = 1, use = "rand")
    # there won't be any drones this season

  # Set a new queen 
    #Jana: the new queen becomes the virgin_queen --> the new virgin queen will be set by crossColony
  colony$queen = selectInd(colony$virgin_queens,  nInd = 1, use = "rand")
  # Set new set of workers to the original colony (which is just a subset of the original set)
  colony$workers = colony$workers[!sel_workers]

  #Jana: This is not ok, since this has to come from the new fathers!!!
  #colony$virgin_queen = colony[sample(x = colony$workers,size = 1, repace = FALSE)] 
  
  # drones stay from the previous queen
  # possibly more code here - do we do mating of the new virgin queen here in this function? 
  #Change the status of the colony
  colony$last_even = "swarmStay" #TODO: better names for this but we have to know which one stayed and which one left due to the drones
  swarm$last_even = "swarmLeave" #TODO: do we want to have some information about the link (i.e. mother_colony=?")
  
  return(list(colony = colony, swarm = swarm))
}



#Splitting of the hive during supersedure and new queen made              
supersedure = function(colony) {
  colony$queen = selectInd(colony$virgin_queens,  nInd = 1, use = "rand")
}



splitColony = function(colony, per_split) {
  # Compute the number of workers that will leave with the swarm (per_swarm is the % of workers that will swarm)
  noWorkersSplit = round(colony$workers$nInd * per_split, 0)
  # Compute the number of workers that will stay in the original colony (better to do it this way due to rounding issues)
  noWorkersStay = colony$workers$nInd - noWorkersSplit
  # Which workers are taken
  workersSplitId = sample(x = colony$workers@id, size = noWorkersSplit, replace = FALSE) # TODO gives ids of workers that will leave with the swarm
  # Which workers stay
  workersStayId = colony$workers@id %in% sel_workers # tells which workers will leave (TRUE) and which won't (FALSE)
  # Create a new colony entity that represents the swarm and set its queen and workers, all the drones stay in the original colony
  splitColony = colony()
  splitColony$queen = colony$queen 
  splitColony$workers = colony$workers[workersSplitId]
  splitColony$virgin_queen = selectInd(colony$workers, nInd = 1, use = "rand")
  # there won't be any drones this season
  
  # Set a new queen
  colony$queen = selectInd(colony$virgin_queens, nInd = 1, use = "rand")
  # Set new set of workers to the original colony (which is just a subset of the original set)
  colony$workers = colony$workers[!sel_workers]

  #Jana: This again is not true - the new virgin queen has to come from tnew queen and drones
  #colony$virgin_queen = colony[sample(x = colony$workers,size = 1, replace = FALSE)] 
  #Change the status of the colony
  colony$last_even = "splitStay" #TODO: better names for this but we have to know which one stayed and which one left due to the drones
  splitColony$last_even = "splitLeave" #TODO: do we want to have some information about the link (i.e. mother_colony=?")
  
  return(list(colony = colony, splitColony = splitColony))
}


buildUpColony = function(colony, per_workers_increase, nDrones) {
  #builtColony = colony() Jana: you don't need this! COlony is an input parameter!!!! We are just changing the existing one, not creating a new one!
  nNewWorkers = round(colony$workers@nInd * per_workers_increase, 0)
  newWorkers = randCross2(females = colony$queen, males = colony$fathers, nCrosses = nNewWorkers)
  colony$workers = mergePop(colony$workers, colony$newWorkers)
  colony$drones = createDrones(colony, nDrones) #TODO: variable drone pop number 
}  

# TODO: this isn't working yet!
select_colonies <- function(colony_list, colony_ids) {
  sel_colony_list = colony_list[sapply(colony_list, FUN = function(x) x$id %in% colony_ids)]
  #sel_colony_list = sel_colony_list[sapply(sel_colony, FUN = function(x) !is.null(x))]
  return(sel_colony_list)
}

# setPheno = function() {} # keep this one commented for now, we need some object-oriented magic for this to work on our colony and not to clash with AlphaSimR:::setPheno()





#Things to add in the future/ check list 
#genetic model --> they have either infinitesimal or finite (we don't have to worry about this)
--> #simulate queen and worker effects (yet to do)
  --> #simulate demographic history (in progress)
  #population --> specify the size of the population (#colonies/y)
  --> #simulate passive population (not actively selected but exchanged queens/drones with the selected)
  --> #specify the rates of exchange between selected and non-selected (we can control this within the breeding program
  --> #specify min/max age of queen to produce queens/drones
  --> #culling ages
  #mating --> specify how queen of the breding and passive populations mate - free mating, AI, breeding stations
  --> #decide how many drones are involved
  --> #For mating stations specify - drone producting queens on mating stations, whether they are related, and how secure is the mating station (probability of a drone in the DCa coming from the mating stations hives)
  #selection - rates of queen selection
  - #how many sisters from the same sister group (within/across family selection)
  - #perform: phenotypic, genotypic, random selection on selection based on EBVs
