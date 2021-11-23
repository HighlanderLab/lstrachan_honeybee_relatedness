library(AlphaSimR)

#Create founder haplotypes
 founderPop = quickHaplo(nInd=200, nChr=1, segSites=10)
 
 #Set simulation parameters
 SP = SimParam$new(founderPop)
 
 #Create population
 pop = newPop(founderPop, simParam=SP)


