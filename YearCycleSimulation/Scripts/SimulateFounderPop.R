# Load packages
library(AlphaSimR)
library(SIMplyBee)

# Founder opulation parameters -------------------------------------------------------------------
nMelN = 800
nCar = 400
nChr = 16
nSegSites = 1000

#Founder population ---------------------------------------------------------
#Create a founder population of A. m. mellifera and A. m. carnica bees
print("Simulating founders")
founderGenomes <- simulateHoneyBeeGenomes(nMelN = nMelN,
                                          nCar = nCar,
                                          nChr = nChr,
                                          nSegSites = nSegSites,
					  nThreads = 16)

save(founderGenomes, file="FounderGenomes_ThreePop_16chr.RData")
