#For function storage/ developement 

#CreateColony function - to create a new colony from existing (not only base population). Arguments- list of nodes within the hive (to be populated)
createColony = function(id = NULL, location = NULL, queen = NULL, drones = NULL, workers = NULL, virgin_queens = NULL, pheno = NULL, fathers = NULL) {
  colony = vector(mode = "list",  length = 8)
  names(colony) = c("id", "location", "queen", "drones", "workers", "virgin_queens", "pheno", "fathers")
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
return(colony)
} 
# we will likely need queen age too - but that should go into colony$queen$misc slot!
#can also look at hive "strength" based on number of colony workers 



 # Generate a swarm and modify the existing colony  # colony - list, colony
# little honey
# younger queens swarm less
# build-up after swarm depends on season - should this be an extra function / or will this be in country scripts?
createSwarm = function(colony) {
  
  sel_workers = sample(x = colony$workers@id, size = ..., replace = FALSE) # TODO gives ids of workers that will leave with the swarm
  sel_workers = colony$workers@id %in% sel_workers # tells which workers will leave (TRUE) and which won't (FALSE)
  swarm = colony()
  swarm$queen = colony$queen 
  swarm$workers = colony$workers[sel_workers]
  # there won't be any drones this season
  # there won't be any virgin queens either, unless supersedure happens
  
  colony$workers = colony$workers[!sel_workers]
  sel_virgin_queen = sample(x = swarm$workers@id, size = 1)
  colony$queen = colony$virgin_queen[sel_virgin_queen]
  colony$virgin_queen = NA # queen kills all the virgin queens
  # drones stay from the previous queen
  # possibly more code here - do we do mating of the new virgin queen here in this function? 
  
  return(list(colony = colony, swarm = swarm))
}



createDCA = function() {}

supersedure = function () {} # how do we make this one into a verb?

beeCross = function() {}

splitColony = function() {}

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
