createMatingStationDCA <- function(colony, nDPQs = 30, nDronePerDPQ = 100) {
  DPQs <- createVirginQueens(colony, nInd = nDPQs)
  drones <- createDrones(DPQs, nInd = nDronePerDPQ)
  return(drones)
}



founderGenomes <- quickHaplo(nInd = 22, nChr = 1, segSites = 100)
SP <- SimParamBee$new(founderGenomes)
basePop <- createVirginQueens(founderGenomes)

apiary <- createColonies(x = basePop[1:10], n = 10)
openDCA <- createDrones(basePop[11:20], nInd = 100)
fatherGroups <- pullDroneGroupsFromDCA(openDCA, n = 10, nFathers = nFathersPoisson)
# Create two mating stations with basePop[21] and basePop[22] queens as the sire
# 1) Mate the queens on the DCA
sireColony1 <- createColony(basePop[21])
sireColony2 <- createColony(basePop[22])
sireColony1 <- crossColony(colony = sireColony1, drones = openDCA, nFathers = 15)
sireColony2 <- crossColony(colony = sireColony2, drones = openDCA, nFathers = 15)
matingStation1 <- createMatingStation(sireColony1, nDPQ = 30, nDronePerDPQ = 1000)
matingStation2 <- createMatingStation(sireColony2, nDPQ = 30, nDronePerDPQ = 1000)
drone1 <- createDrones(sireColony1, nInd = 1)

# Try to create a breeding plan - some of the queens will get mated to the open DCA and some to specific queen
# DO we also want to allow for a single queen/colony as the input (as sire / AI)?
matingPlan <- data.frame(VirginQueens = getId(mergePops(getVirginQueens(apiary))), # we really need colony IDs
                         Father = c(rep("openDCA", 3), rep("matingStation1", 3), rep("matingStation2", 3), "drone1"),
                         nFathers = rep("nFathersPoisson"))


# Check whether the length of the mating plan matches the length of the input
colonies <- apiary
nCol <- nColonies(colonies)
ret <- createColonies()
for (colony in seq_len(nCol)) {
  ret[[colony]] <- crossColony(
    colony = colonies[[colony]],
    drones = get(matingPlan$Father[matingPlan$VirginQueens == getVirginQueens(colonies[[1]])@id]),
    nFathers = nFathersPoisson,
    removeFathers = T,
    checkMating = T,
    simParamBee = SP
  )
}






