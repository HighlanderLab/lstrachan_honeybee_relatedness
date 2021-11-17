library(AlphaSimR)
source("HoneyBeeSim_functions.R")

# Founder population - colonies
#Create founders
founder = quickHaplo(n = 1000,
                     nChr = 1,
                     segSites = 10)

#Create base population
SP = SimParam$new(founder)
base = newPop(founder)


####
#Parameters
noFounderColonies = 10
colonyFullSize = 1000
p1swarm = 0.05
p1supersede = 0.05
p1collapse = 0.10
p2swarm = 0.01
p2supersede = p1supersede
p2collapse = 0.10
p3collapse = 0.35


col <- createColony(queen = base[1],
                    fathers = base[2:16],
                    last_event = "period3")
col
col1 <- createColony(queen = base[17],
                    fathers = base[18:28],
                    last_event = "period3")
col1
# Period 1:
# Create short-living workers 10k-70k
col <- addWorkers(col, nInd = 400)
col
col1 <- addWorkers(col1, nInd = 300)
col1


# Create drones
col <- addDrones(col, nInd = 20)
col
col1 <- addDrones(col1, nInd = 30)
col1

# Splitting occurs
splitResult <- splitColony(col, pSplit = 0.3, newQueen = col@virgin_queens,
                          crossVirginQueen = TRUE, fathers = createDrones(col1, 14))
col <- splitResult$remnant
col
col2 <- splitResult$split
col2

# Swarming and supersedure occurs (supersedure is successful when drones are present)
swarmResult = swarmColony(col, pSwarm = 0.1, crossVirginQueen = TRUE, fathers = createDrones(col2, 12), pWorkers = 0.1, pDrones = 0.2)
col = swarmResult$swarm
col
col4 = swarmResult$remnant
col4
col1
col1 <- supersedeColony(col1, crossVirginQueen = TRUE, fathers = createDrones(col2, 12), pWorkers = 0.5, pDrones = 0.5)
col1

# First opportunity for mating
# TODO: Queen mated
# TODO: Colony loss

# Colony losses (p)
#TODO: Later

# Period 2:
# Replace ½ / 2/3 / all of spring workers and drones
col <- replaceWorkers(col)
col <- replaceDrones(col)
col


# Swarming and supersedure occurs at a lower percentage (supersedure is successful when drones are present)
#TODO: Later

# Production period (all except development)
#TODO: Later

# Colony losses (p)
#TODO: Later

# Period 3:
# Remove drones
col@drones <- NULL

# Reduce number of workers
# Replace workers with long-winter workers
col@workers  <- NULL


# Highest percentage of colony losses
#TODO: COlony loss

# Reset the events for all the colonies
col <- resetEvents(col)
col
