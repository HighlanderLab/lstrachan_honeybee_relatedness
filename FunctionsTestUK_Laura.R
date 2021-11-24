library(AlphaSimR)
source("~/Desktop/GitHub/lstrachan_honeybee_sim/Functions_L0_auxilary.R")
source("~/Desktop/GitHub/lstrachan_honeybee_sim/Functions_L1_Pop.R")
source("~/Desktop/GitHub/lstrachan_honeybee_sim/Functions_L2_Colony.R")
source("~/Desktop/GitHub/lstrachan_honeybee_sim/Functions_L3_Colonies.R")
source("~/Desktop/GitHub/lstrachan_honeybee_sim/Classes.R")

########################       Global parameters        ########################

nFounderColonies = 10
nWorkersFull = 1000
nAvgFathers = 17

#Events percentages (based on Period)
#Period1
p1swarm = 0.05
p1supersede = 0.05
p1collapse = 0.10
#Period2
p2swarm = 0.01
p2supersede = p1supersede
p2collapse = 0.10
#Period3
p3collapse = 0.35
#TODO: Change into a vector of probabilities (age 1,2,3 of the queen)


##################     Initial simulation set-up         #######################

#Create founder population
founder_population = quickHaplo(nInd = 1000,
                                nChr = 16,
                                ploidy = 2L,
                                inbred = FALSE, 
                                segSites = 1000)

#Add simulation parameters 
SP = SimParam$new(founder_population) 

#Add traits: honey yield
# Queen and average worker effect covariance 
covA = matrix(data = c( 1.0, -0.5,
                        -0.5,  2.0),
              nrow = 2)
# Queen and individual worker effect
covA = matrix(data = c(covA[1,1], covA[1,2],
                       covA[2,1], covA[2,2]*nWorkersFull),
              nrow = 2)
covA

SP$addTraitA(nQtlPerChr = 10, 
             mean = c(0, 0), 
             var = diag(covA), 
             corA = cov2cor(covA))
covE = covA
covE[1,2] <- covE[2,1] <- 0

#Create base population 
base = newPop(founder_population)

####################             Period 1- Year 1        #######################

#Create 10 mated colonies from the base population
age1= createMultipleMatedColonies(base, nColonies = nFounderColonies, nAvgFathers = nAvgFathers)

#Build-up the colonies
age1= buildUpColonies(age1, nWorkers = nWorkersFull, nDrones = nWorkersFull * 0.1)

#Split all of the colonies----
tmp <- splitColonies(age1)
age1 <- tmp$remnants
age0p1 <- tmp$splits

#Create 10 virgin queens (Vqueens taken from hive 3 in age1)
virginQueens = createVirginQueens(age1[[3]], nColonies(age0p1))

# Requeen the splits
age0p1 <- reQueenColonies(age0p1, queens = virginQueens)

# Swarm a percentage of age1 colonies----
tmp = pullColonies(age1, p = p1swarm)
age1 = tmp$remainingColonies
tmp = swarmColonies(tmp$pulledColonies)
age0p1 = c(age0p1, tmp$remnants)
age1 = c(age1, tmp$swarms)

#Supersede a percentage of age1 colonies----
tmp = pullColonies(age1, p = p1supersede)
age1 = tmp$remainingColonies
tmp = supersedeColonies(tmp$pulledColonies)
age0p1 = c(age0p1, tmp)

#Mate the age0 virgin queens----
DCA = createDCA(age1)
age0p1 = crossColonies(age0p1, DCA, nAvgFathers = nAvgFathers)

#Collapse a percentage of age1 colonies----
age1 = selectColonies(age1, p = 1 - p1collapse)

####################             Period 2- Year 1        #######################




