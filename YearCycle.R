library(AlphaSimR)
source("/home/jana/lstrachan_honeybee_sim/Classes.R")
source("/home/jana/lstrachan_honeybee_sim/Functions_L0_auxilary.R")
source("/home/jana/lstrachan_honeybee_sim/Functions_L1_Pop.R")
source("/home/jana/lstrachan_honeybee_sim/Functions_L2_Colony.R")
source("/home/jana/lstrachan_honeybee_sim/Functions_L3_Colonies.R")


####
#Parameters
nFounderColonies = 10
nWorkersFull = 100
nAvgFathers = 2


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



#Period1
#Create 10 mated colonies from the base population
age1 = createMultipleMatedColonies(base, nColonies = nFounderColonies, nAvgFathers = nAvgFathers)

#Build-up the colonies
age1 = buildUpColonies(age1, nWorkers = nWorkersFull, nDrones = nWorkersFull * 0.1)

#Split all the colonies
tmp <- splitColonies(age1)
age1 <- tmp$remnants
age0p1 <- tmp$splits

#Create 10 virgin queens
virginQueens = createVirginQueens(age1[[3]], nColonies(age0p1))

# Requeen the splits
age0p1 <- reQueenColonies(age0p1, queens = virginQueens)

# Swarm a percentage of age1 colonies
tmp = pullColonies(age1, p = p1swarm)
age1 = tmp$remainingColonies
tmp = swarmColonies(tmp$pulledColonies)
age0p1 = c(age0p1, tmp$remnants)
age1 = c(age1, tmp$swarms)

#Supersede
tmp = pullColonies(age1, p = p1supersede)
age1 = tmp$remainingColonies
tmp = supersedeColonies(tmp$pulledColonies)
age0p1 = c(age0p1, tmp)

#Mate the split colonies
DCA = createDCA(age1)
age0p1 = crossColonies(age0p1, DCA, nAvgFathers = nAvgFathers)

#Collapse
age1 = selectColonies(age1, p = 1 - p1collapse)

#Period2
# Swarm a percentage of age1 colonies
tmp = pullColonies(age1, p = p2swarm)
age1 = tmp$remainingColonies
tmp = swarmColonies(tmp$pulledColonies)
age0p2 = tmp$remnants
age1 = c(age1, tmp$swarms)

#Supersede
tmp = pullColonies(age1, p = p2supersede)
age1 = tmp$remainingColonies
tmp = supersedeColonies(tmp$pulledColonies)
age0p2 = c(age0p2, tmp)

#Mate the split colonies
DCA = createDCA(age1)
age0p2 = crossColonies(age0p2, DCA, nAvgFathers = nAvgFathers)

#Phenotype
age1 = setPhenoColonies(age1, varE = covE)

phenoColony = function(x) {
  x@queen@pheno[, 1] + colMeans(x@workers@pheno[, 2, drop = FALSE])
}
age1pheno = setPhenoColonies(age1, FUN = phenoColony, varE = covE)
 
# phenoColony = function(x) {
#   cbind(x@queen@pheno[, "honey_yield_Q"] + colMeans(x@workers@pheno[, "honey_yield_W", drop = FALSE]),
#         x@queen@pheno[, "swarm_Q"]       + colMeans(x@workers@pheno[, "swarm_W",       drop = FALSE]))
# }

#Collapse
age1 = selectColonies(age1, p = 1 - p2collapse)
