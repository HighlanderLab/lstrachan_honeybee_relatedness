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


col <- createColony(queen = base[1],
                    fathers = base[2:16],
                    last_event = "period3")
col

# Period 1:
# Create short-living workers 10k-70k
col <- addWorkers(col, nWorkers = 400)
col

# Create drones
col <- addDrones(col, nDrones = 20)
col

# Splitting occurs
#TODO

# Swarming and supersedure occurs (supersedure is successful when drones are present)
#TODO: Later

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