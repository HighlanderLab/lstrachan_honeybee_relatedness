library(AlphaSimR)
source("/home/jana/lstrachan_honeybee_sim/HoneyBeeSim_functions.R")

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
#TODO: Change into a vector of probabilites (age 1,2,3 of the queen)
################################################################################

#Period1
#Create 10 mated colonies from the base population
apiary1 = createMultipleMatedColonies(base, nColonies = 10, nAvgFathers = 15)

#Build-up the colonies
apiary1 = buildUpColonies(apiary1, nWorkers = colonyFullSize, nDrones = colonyFullSize * 0.1)

#Split all the colonies
tmp <- splitColonies(apiary1)
apiary1 <- tmp$remnants
apiary0 <- tmp$splits

#Create 10 virgin queens
virginQueens = createVirginQueens(apiary1[[3]], 10)

# Requeen the splits
apiary0 <- reQueenColonies(apiary0, queens = virginQueens)

#Mate the split colonies
DCA = createDCA(apiary1)
apiary0 = crossColonies(apiary0, DCA, nAvgFathers = 15)

#Build-up the splits
apiary0 = buildUpColonies(apiary0, nWorkers = colonyFullSize, nDrones = colonyFullSize * 0.1)

# Swarm a percentage of apiary1 colonies


################################################################################
################################################################################
################################################################################
DCA = createFounderDrones(queenPop = base[1:10], nDronesPerQueen = 10)
DCA




###########################################################
#Period 1
###########################################################





#Create 3 colonies
DCAresult = pullDronesFromDCA(DCA, 14)
DCA = DCAresult$DCA
col1 = createColony(queen = base[11], fathers = DCAresult$selectedDrones)
DCAresult = pullDronesFromDCA(DCA, 12)
DCA = DCAresult$DCA
col2 = createColony(queen = base[12], fathers = DCAresult$selectedDrones)
DCAresult = pullDronesFromDCA(DCA, 15)
DCA = DCAresult$DCA
col3 = createColony(queen = base[13], fathers = DCAresult$selectedDrones)
DCA


# Add workers and drones
col1 = addWorkers(col1, nInd = colonyFullSize)
col2 = addWorkers(col2, nInd = colonyFullSize)
col3 = addWorkers(col3, nInd = colonyFullSize)
col1 = addDrones(col1, nInd = colonyFullSize*0.1)
col2 = addDrones(col2, nInd = colonyFullSize*0.1)
col3 = addDrones(col3, nInd = colonyFullSize*0.1)


# Create virgin queens for the splits
col3 = addVirginQueens(col3, nVirginQueens = 3)
tmp = pullIndFromCaste(col3, "virgin_queens", 3)
col3 = tmp$colony
virgin_queens = tmp$pulledInd
virgin_queens

# Create DCA from colonies 1 - 3
DCA = createDCA(list(col1, col2, col3))
DCA

# Split all the colonies
DCAresult = pullDronesFromDCA(DCA, nInd = 11)
DCA = DCAresult$DCA
col1split <- splitColony(col1,
                         pSplit = 0.3,
                         newQueen = virgin_queens[1],
                         crossVirginQueen = TRUE,
                         fathers = DCAresult$selectedDrones)
col1 = col1split$remnant
col4 = col1split$split
col1
col4

DCAresult = pullDronesFromDCA(DCA, nInd = 13)
DCA = DCAresult$DCA
col2split <- splitColony(col2,
                         pSplit = 0.2,
                         newQueen = virgin_queens[2],
                         crossVirginQueen = TRUE,
                         fathers = DCAresult$selectedDrones)
col2 = col2split$remnant
col5 = col2split$split
col2
col5

DCAresult = pullDronesFromDCA(DCA, nInd = 16)
DCA = DCAresult$DCA
col3split <- splitColony(col3,
                         pSplit = 0.25,
                         newQueen = virgin_queens[3],
                         crossVirginQueen = TRUE,
                         fathers = DCAresult$selectedDrones)
col3 = col3split$remnant
col6 = col3split$split
col3
col6

# 

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
