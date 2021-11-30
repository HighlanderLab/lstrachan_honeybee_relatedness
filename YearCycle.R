library(AlphaSimR)
source("~/JANA/github/AlphaSimRBee/AlphaSimRBee/R/Classes.R")
source("~/JANA/github/AlphaSimRBee/AlphaSimRBee/R/Functions_L0_auxilary.R")
source("~/JANA/github/AlphaSimRBee/AlphaSimRBee/R/Functions_L1_Pop.R")
source("~/JANA/github/AlphaSimRBee/AlphaSimRBee/R/Functions_L2_Colony.R")
source("~/JANA/github/AlphaSimRBee/AlphaSimRBee/R/Functions_L3_Colonies.R")


####
#Parameters
nFounderColonies = 20
nWorkersFull = 1000
nAvgFathers = 15


#Period1
p1swarm = 0.05
p1supersede = 0.05
p1collapse = 0.10
#Period2
p2swarm = 0.01
p2supersede = p1supersede
p2collapse = 0.10
#Period3
p3collapseAge0 = 0.10
p3collapseAge1 = 0.25
#TODO: Change into a vector of probabilites (age 1,2,3 of the queen)
################################################################################

# Founder population - colonies
#Create founders
founder = quickHaplo(n = 1000,
                     nChr = 1,
                     segSites = 10)

#Create base population
SP = SimParam$new(founder)
#Add traits: honey yield
# Queen and average worker effect cov
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
base = newPop(founder)


####Year 0
#Period1
#Create 10 mated colonies from the base population

for (year in 1:2) {
  if (year == 1) {
    age1 = createMultipleMatedColonies(base, nColonies = nFounderColonies, nAvgFathers = nAvgFathers)
  } else {
    age2 = age1
    age1 = age0
    age0 = NULL
    age0p1 = NULL
    age0p2 = NULL
  }

  #########################################################################
  #Period1
  #########################################################################
  #Build-up the colonies
  age1 = buildUpColonies(age1, nWorkers = nWorkersFull, nDrones = nWorkersFull * 0.1)
  if (year > 1) {
    age2 <- buildUpColonies(age2, nWorkers = nWorkersFull, nDrones = nWorkersFull * 0.1)
  }
  
  #Split all age1 colonies
  tmp <- splitColonies(age1)
  age1 <- tmp$remnants
  # The queens of the splits are 0 years old
  age0p1 <- tmp$splits
  
  if (year > 1) {
    #Split all age2 colonies
    tmp <- splitColonies(age2)
    age2 <- tmp$remnants
    # The queens of the splits are 0 years old
    age0p1 <- c(age0p1, tmp$splits)
  }
    
  #Create 10 virgin queens
  #Sample colony for the virgin queens
  virginDonor <- sample(1:nColonies(age1), size = 1)
  virginQueens = createVirginQueens(age1[[virginDonor]], nColonies(age0p1))
  
  # Requeen the splits --> queens are now 0 years old
  age0p1 <- reQueenColonies(age0p1, queens = virginQueens)
  
  # Swarm a percentage of age1 colonies
  tmp = pullColonies(age1, p = p1swarm)
  age1 = tmp$remainingColonies
  tmp = swarmColonies(tmp$pulledColonies)
  age0p1 = c(age0p1, tmp$remnants)
  age1 = c(age1, tmp$swarms)
  
  if (year > 1) {
    # Swarm a percentage of age2 colonies
    tmp = pullColonies(age2, p = p1swarm)
    age2 = tmp$remainingColonies
    tmp = swarmColonies(tmp$pulledColonies)
    age0p1 = c(age0p1, tmp$remnants)
    age2 = c(age2, tmp$swarms)
  }
  
  #Supersede age 1 collonies
  tmp = pullColonies(age1, p = p1supersede)
  age1 = tmp$remainingColonies
  tmp = supersedeColonies(tmp$pulledColonies)
  age0p1 = c(age0p1, tmp)
  
  if (year > 1) {
    #Supersede age 1 collonies
    tmp = pullColonies(age2, p = p1supersede)
    age2 = tmp$remainingColonies
    tmp = supersedeColonies(tmp$pulledColonies)
    age0p1 = c(age0p1, tmp)
  }
  
  #Mate the split colonies
  if (year == 1) {
    DCA = createDCA(age1)
  } else {
    DCA = createDCA(c(age1, age2))
  }
  age0p1 = crossColonies(age0p1, DCA, nAvgFathers = nAvgFathers)
  
  #Collapse
  age1 = selectColonies(age1, p = 1 - p1collapse)
  if (year > 1) {
    age2 = selectColonies(age2, p = 1 - p1collapse)
  }
  
  
  #########################################################################
  #Period2
  #########################################################################
  # Swarm a percentage of age1 colonies
  tmp = pullColonies(age1, p = p2swarm)
  age1 = tmp$remainingColonies
  tmp = swarmColonies(tmp$pulledColonies)
  # The queens of the remnant colonies are of age 0
  age0p2 = tmp$remnants
  age1 = c(age1, tmp$swarms)
  
  if (year > 1) {
    # Swarm a percentage of age2 colonies
    tmp = pullColonies(age2, p = p2swarm)
    age2 = tmp$remainingColonies
    tmp = swarmColonies(tmp$pulledColonies)
    # The queens of the remnant colonies are of age 0
    age0p2 = c(age0p2, tmp$remnants)
    age2 = c(age2, tmp$swarms)
  }
  
  #Supersede a part of age1 colonies
  tmp = pullColonies(age1, p = p2supersede)
  age1 = tmp$remainingColonies
  tmp = supersedeColonies(tmp$pulledColonies)
  # The queens of superseded colonies are of age 0
  age0p2 = c(age0p2, tmp)
  
  if (year > 1) {
    #Supersede a part of age1 colonies
    tmp = pullColonies(age2, p = p2supersede)
    age2 = tmp$remainingColonies
    tmp = supersedeColonies(tmp$pulledColonies)
    # The queens of superseded colonies are of age 0
    age0p2 = c(age0p2, tmp)
  }
  
  #Mate the split colonies
  # Replace all the drones
  age1 <- replaceDronesColonies(age1)
  if (year > 1) {
    age2 <- replaceDronesColonies(age2)
  }
  
  if (year == 1) {
    DCA = createDCA(age1)
  } else {
    DCA = createDCA(c(age1, age2))
  }
  # Cross age 0 period 2 swarms and splits
  age0p2 = crossColonies(age0p2, DCA, nAvgFathers = nAvgFathers)
  
  #Collapse
  age1 = selectColonies(age1, p = 1 - p2collapse)
  if (year > 1) {
    age2 = selectColonies(age2, p = 1 - p2collapse)
  }
  
  # Merge all age 0 colonies (from both periods)
  age0 <- c(age0p1, age0p2)
  
  #########################################################################
  #Period3
  #########################################################################
  # Collapse age0 queens
  age0 <- selectColonies(age0, p = (1 - p3collapseAge0))
  age1 <- selectColonies(age1, p = (1 - p3collapseAge1))
  age2 = NULL #We don't need this but just to show the workflow!!!

}



#Phenotype
# age1 = setPhenoColonies(age1, varE = covE)
# 
# phenoColony = function(x) {
#   x@queen@pheno[, 1] + colMeans(x@workers@pheno[, 2, drop = FALSE])
# }
# age1pheno = setPhenoColonies(age1, FUN = phenoColony, varE = covE)

# phenoColony = function(x) {
#   cbind(x@queen@pheno[, "honey_yield_Q"] + colMeans(x@workers@pheno[, "honey_yield_W", drop = FALSE]),
#         x@queen@pheno[, "swarm_Q"]       + colMeans(x@workers@pheno[, "swarm_W",       drop = FALSE]))
# }