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
  if (is.null(workers)) {
    colony$workers = workers
  }
  if (is.null(virgin_queens)) {
    colony$virgin_queens = virgin_queens
  }
  if (is.null(pheno)) {
    colony$pheno = pheno
  }
  if (is.null(fathers)) {
    colony$fathers = fathers
    }
  if (is.null(last_event)) {
    colony$last_event = "new_colony"
    }
return(colony)
} 
# we will likely need queen age too - but that should go into colony$queen$misc slot!
#can also look at hive "strength" based on number of colony workers 



 # Generate a swarm and modify the existing colony  # colony - list, colony
# little honey
# younger queens swarm less
# build-up after swarm depends on season - should this be an extra function / or will this be in country scripts?
createSwarm = function(colony, per_swarm) {
  # Compute the number of workers that will leave with the swarm (per_swarm is the % of workers that will swarm)
  noWorkersSwarm = round(length(colony$workers$id) * per_swarm, 0)
  # Compute the number of workers that will stay in the original colony (better to do it this way due to rounding issues)
  noWorkersStay = length(colony$workers$id) - noWorkersSwarm
  # Which workers swarm
  workersSwarmId = sample(x = colony$workers@id, size = noWorkersSwarm, replace = FALSE) # TODO gives ids of workers that will leave with the swarm
  # Which workers stay
  workersStayId = colony$workers@id %in% sel_workers # tells which workers will leave (TRUE) and which won't (FALSE)
  # Create a new colony entity that represents the swarm and set its queen and workers, all the drones stay in the original colony
  swarm = colony()
  swarm$queen = colony$queen 
  swarm$workers = colony$workers[workersSwarmId]
  swarm$virgin_queen = #TODO: sample a virgin queen from the workers
  # there won't be any drones this season
  # there won't be any virgin queens either, unless supersedure happens
  # Set new set of workers to the original colony (which is just a subset of the original set)
  colony$workers = colony$workers[!sel_workers]
  # Set a new queen
  colony$queen = colony[sample(x = colony$virgin_queen$id, size = 1)]
  colony$virgin_queen = NA #TODO: select a new virgin_queen
  # drones stay from the previous queen
  # possibly more code here - do we do mating of the new virgin queen here in this function? 
  #Change the status of the colony
  colony$last_even = "swarmStay" #TODO: better names for this but we have to know which one stayed and which one left due to the drones
  swarm$last_even = "swarmLeave" #TODO: do we want to have some information about the link (i.e. mother_colony=?")
  
  return(list(colony = colony, swarm = swarm))
}

#Drone congregation area of the base population made in colony_list 
createDCA = function(colony_list = NULL, colonyIDs = NULL) { # Would it be better
  dca_colony_list = select_colonies(colony_list, colonyIDs)
  DCA = lapply(X = dca_colony_list, FUN = function(z) z$drones)  
  DCA = mergePops(popList = DCA)
  
  return(popList = DCA)
}

# TODO: this isn't working yet!
select_colonies <- function(colony_list, colony_ids) {
  sel_colony_list = lapply(colony_list, FUN = function(x) if (x$id %in% colony_ids) {return(x)})
}

#Splitting of the hive during supersedure and new queen made              
supersedure = function(create_colony) {
  nWorkersPerDrone = nBeesPerColony / nMatingDrones[create_colony]
  
  supersedure = create_colony()
  supersedure$queen = supersedure$virgin_queens
  supersedure$drones = makeDH(pop = supersedure$queen, nDH = 50) #TODO: variable number of drones  
  
  n = nInd(supersedure$workers)
  supersedure$workers = c(supersedure$workers[sample.int(n = n, size = round(n/2))], 
                          randCross2(females = supersedure$queen, males = DCA, #To do: artificial insemination here? / local DCA 
                                     nCrosses = nMatingDrones[colony], nProgeny = round(nWorkersPerDrone/2)))
  supersedure$virgin_queens = randCross2(females =  supersedure$queen, 
                                         males = base_pop, nCrosses = 1, nProgeny = 1) #to do: variable number of drones
  
  return(list(colony = supersedure))
}

#another way to make the supersedure workers ?
sel_workers = sample(x = supersedure$workers@id, size = ..., replace = FALSE) #TODO: complete size
sel_workers = supersedure$workers@id %in% sel_workers
supersedure$workers = c(supersedure$workers[sel_workers], 
                        randCross2(females = supersedure$queen, males = DCA, #To do: artificial insemination here? / local DCA 
                                   nCrosses = nMatingDrones[colony], nProgeny = round(nWorkersPerDrone/2)))



beeCross = function(colony, drone_pop, nBees_created) {
  #TODO: do we mate the virgin queen or set the new before this???????
  #The drone pop can be created with the createDCA - but that is a separate step
  #1) Write the fathers into the colony$fathers slot
  colony$fathers = drone_pop
  #2) Create the workers and write them in the worker slot
  colony$workers = randCross2(females = colony$virgin_queen, males = drone_pop, nProgeny = nBees_created)
  colony$drones = #TODO: CreateDrones functions
  colony$queen = colony$virgin_queen
  colony$virgin_queen = #TODO: sample from workers

  # Jana: I don't think this function returns anything
}

splitColony = function(colony, split_percentage) { #This is a colony split artificially by keeper 

  # TODO: Jana: I would do this the same as the swarm
  sel_workers = sample(x = colony$workers@id, size = split_percentage, replace = FALSE) #gives ids of workers that will leave with the swarm
  sel_workers = colony$workers@id %in% sel_workers # tells which workers will leave (TRUE) and which won't (FALSE)
  
  splitColony = colony()
  splitColony$queen = colony$queen 
  splitColony$workers = colony$workers[sel_workers]
  
  colony$workers = colony$workers[!sel_workers]
  sel_virgin_queen = sample(x = splitColony$virgin_queen, size = 1)
  colony$queen = colony$virgin_queen[sel_virgin_queen]
  colony$virgin_queen = NA # queen kills all the virgin queens
  # drones stay from the previous queen
  
  return(list(colony = colony, splitColony = splitColony))
}
    
 # buildUpColony = function(colony)  #still being written 
#TODO: add in workers AND drones

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
