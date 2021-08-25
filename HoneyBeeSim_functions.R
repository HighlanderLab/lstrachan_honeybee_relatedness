#For function storage/ developement 

createColony = function(id = NULL, location = NULL) {
  colony = vector(mode = "list",  length = 7)
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
} # remember to add fathers and id, we will likely need queen age too - but that should go into colony$queen$misc slot!

createSwarm = function() {} 
# no new drones from existing queen
# little honey
# younger queens swarm less
# build-up after swarm depends on season - should this be an extra function / or will this be in country scripts?

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
