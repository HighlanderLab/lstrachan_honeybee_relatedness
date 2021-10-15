library(AlphaSimR)
source("HoneyBeeSim_functions.R")

#Create founders
founder = quickHaplo(n = 1000,
                     nChr = 1,
                     segSites = 10)

#Create base population
SP = SimParam$new(founder)
base = newPop(founder)


# ---Create a colony---#
# Create colonies 1 and 2
col1 = createColony(queen = base[1], fathers = base[2:50])
col1
col2 = createColony(id = 2, queen = base[51])
col2

#########################################
#Test the functions on the Colony object

#--- createDrones ---#
# Create drones for both colonies
col1@drones <- createDrones(col1, 100)
col1@drones
col2@drones <- createDrones(col2, 50)
col2@drones

#--- createDCA ---#
# Create DCA from the drones of colonies 1 and 2
DCA <- createDCA(list(col1, col2))

#--- selectFathersFromDCA ---#
# Select fathers for both colonies
fathers1 <- selectFathersFromDCA(DCA, 14)
fathers2 <- selectFathersFromDCA(DCA, 14)

# This could also be done with a more general --- extractIndFromCast ---#
fathers1 <- extractIndFromCast(col1, "drones", 14)
fathers2 <- extractIndFromCast(col2, "drones", 14)
fathers1
class(fathers1)
fathers2

#--- crossColony ---#
# Cross colony 1 with fathers1 and col 2 with fathers 2
col1 = crossColony(col1, fathers1, nWorkers = 1000)
col1 = crossColony(col1, nWorkers = 1000)
col1
col2 = crossColony(col2, fathers2, nWorkers = 2000)
col2

#--- swarmColony---#
# Swarm colony 1 and produce colony 3
swarmed = swarmColony(col1, perSwarm=0.4)
col1 = swarmed$swarm
#ID of the colony is the queen
col1
col3 = swarmed$newColony
col3
# The new colony is not crossed!
col3 <- crossColony(col3, fathers = DCA, nWorkers = 700, nDrones = 30)
col3


#--- supersedeColony  ---#
#Supersede colony 2
col2
col2 = supersedeColony(col2)
col2

#--- splitColony ---#
# Split colony 2
col2
splitted = splitColony(col2, 0.3)
col2 = splitted$colony
col2
col4 = splitted$splitColony
col4


#--- addWorkers, addDrones ---#
# Build up colonies 2 and 4 (after splitting)
col2
col2 = addWorkers(col2, nWorkersAdd = 1000)
col2
col2 = addDrones(col2, 200)
col4
col4 = addWorkers(col4, nWorkersAdd = 2400)
col4 = addDrones(col4, nDronesAdd = 100)


#--- replaceWorkers, replaceDrones ---#
#Replace workers and drones in col2
col2
max(col2@workers@id)
col2 = replaceWorkers(col2)
col2
max(col2@workers@id)

max(col2@drones@id)
col2 = replaceDrones(col2)
max(col2@drones@id)



# #Create ColonyList
# breedingProgram = createColonyList()
# # Add to colonylist
# breedingProgram <- addColonyToTheColonyList(colony = col1, colonyList = breedingProgram)
# breedingProgram <- addColonyToTheColonyList(colony = col2, colonyList = breedingProgram)
# str(breedingProgram@colonies)
# length(breedingProgram@colonies)

# 
# extractIDs <- function(colonyList) {
#   return(sapply(colonyList@colonies, FUN = function(x) x@id))
# }
# 
# extractIDs(breedingProgram)
# 
# # Test the iteration
# for (colony in 1:1) {
#   breedingProgram[[colony]]@drones <- createDrones(breedingProgram[[colony]], 10)
# }
# breedingProgram[[1]]@drones


