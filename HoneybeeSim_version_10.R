
#Packages used ------------------------------------------------------------
library(AlphaSimR)
p_supersedure = 0.5
p_swarming = 0.5
nBeesPerColony = 500

#Setting up founder population ----------------------------------------------
founder_population = quickHaplo(nInd = 10, nChr = 16, ploidy = 2L,
                                inbred = FALSE, segSites = 1000)

#TO DO: DEFINE GENETIC PARAMETERS (STDPOPSIM)

#Setting simulation parameters-----------------------------------------------
SP = SimParam$new(founder_population) 
SP$addTraitA(
  nQtlPerChr = 1000,
  mean = 0, 
  var = 1, 
  corA = NULL, 
  gamma = FALSE,
  shape = 1,
  force = FALSE) 
#No of QTLs assumed to be ~1000 to be stable 

#Adding phenotypes----
setPheno(
  founder_population,
  h2 = NULL,            # a vector of desired narrow-sense heritabilities for each trait
  H2 = NULL,            # a vector of desired broad-sense heritabilities for each trait.
  varE = NULL,          # error (co)variances for traits
  reps = 1,               
  fixEff = 1L,          # fixed effect to assign to the population
  p = NULL,             # p value for environmental covariate used by GxE traits, if NULL a value is sampled at random 
  onlyPheno = FALSE,    # if onlyPHENO= FALSE returns object of pop-class, if TRUE matrix is returned 
  simParam = NULL
)
#There are three arguments for setting the error variance of a phenotype: 
#h2, H2, and varE.
#The user should only use one of these arguments. 
#If the user supplies values for more than one, only one will be used according to order in which they are listed above.




#Define Csd locus = assuming 64 haplotypes (2^6 = 64) ---------------------------
GenMap = SP$genMap
GenMap[[3]][1:6] = 0
SP$switchGenMap(GenMap) 


#Create a base population --------------------------------------------------
base_pop = newPop(rawPop = founder_population)

template_colony = vector(mode = "list", length = 7)
names(template_colony) = c("name","id", "location", "queen", "drones", "workers", "virgin_queens", "pheno")
#We will expand the above list with further data, say (location,pheno, spatial(vector of x,y))

#Number of colonies controller--------------- 
nColonies <- 10 # to do: variable number of Colonies 

#Create a colony list -----------------------------------------------------------
colony_list = vector(mode = "list", length = nColonies)
for (colony in 1:nColonies) {
  # colony = 1 
  colony_list[[colony]] = template_colony
  
  #populating the queens group
  colony_list[[colony]]$queen = randCross(pop = base_pop, nCrosses = 1, nProgeny = 1)
  
  
  #populating the drones group and making the drones Haploid 
  colony_list[[colony]]$drones = makeDH(pop = colony_list[[colony]]$queen, nDH = 50) #to do: variable number of drones
  
  
  #populating the workers group
  colony_list[[colony]]$workers = randCross2(females =  colony_list[[colony]]$queen, 
                                         males = base_pop, nCrosses = 5, nProgeny = nBeesPerColony) #to do: variable number of drones
  
  #csd of workers (gets rid of unviable diploids)
  Temp = pullSegSiteGeno(colony_list[[colony]]$workers)[3, 1:6]
  sel = apply(Temp, MARGIN = 1, FUN = function(z) any(z == 1))
  colony_list[[colony]]$workers = colony_list[[colony]]$workers[sel]
  
  #populating virgin queens   
  #n_virgin_queens = rpois(n = 1, lambda = 25)  #lambda = number of v_q produced per year (STILL TBD)
  colony_list[[colony]]$virgin_queens = randCross2(females =  colony_list[[colony]]$queen, 
                                               males = base_pop, nCrosses = 1, nProgeny = 1) #to do: variable number of drones
  
}

#Merge the pops / create DCA and virgin queen group  --------
DCA = lapply(X = colony_list, FUN = function(z) z$drones)  
DCA = mergePops(popList = DCA)

virgin_queens_group = lapply(X = colony_list, FUN = function(z) z$virgin_queens)
virgin_queens_group[1] #Reminder = first queen in the queen group 

#Queens age (variable )
colony$queen@misc$age = 1
#so e.g if  colony$queen@misc$age = >1 then the probability of her survival goes down 

#Create new colony with 1 queen and DCA -----------------------------------------------
nMatingDrones = rpois(n = nColonies, lambda = 17)  # sampling random number of drones for this queen
for (colony in 1:nColonies) {
  #colony=1 
  change = sample(x = c("no change", "supersedure", "swarming", "collapse"), n = 1 ) # to do : add the probabilities
  if (change == "supersedure") {
    nWorkersPerDrone = nBeesPerColony / nMatingDrones[colony] # to do: variable number of workers 
    
    colony_list[[colony]]$queen = colony_list[[colony]]$virgin_queens #TO DO : THINK ABOUT HOW TO DO SELECTION AMONG VIRGIN QUEENS AND REPLACEMENT 
    colony_list[[colony]]$drones = makeDH(pop= colony_list[[colony]]$queen, nDH = 50) #to do: variable number of drones
    n = nInd(colony_list[[colony]]$workers)
    colony_list[[colony]]$workers = c(colony_list[[colony]]$workers[sample.int(n = n, size = round(n/2))], 
                                  randCross2(females = colony_list[[colony]]$queen, males = DCA, #To do: artificial insemination here? / local DCA 
                                             nCrosses = nMatingDrones[colony], nProgeny = round(nWorkersPerDrone/2)))
    colony_list[[colony]]$virgin_queens = randCross2(females =  colony_list[[colony]]$queen, 
                                                 males = base_pop, nCrosses = 1, nProgeny = 1) #to do: variable number of drones
  }
  if (change == "swarming") {
    swarm = colony_list[[colony]]
    n = nInd(swarm$workers)  # TODO: make n variable
    swarm$workers = swarm$workers[sample.int(n = n, size = round(n/2))]
    colony_list = c(colony_list, swarm) #TODO: change location of the swarm 
    #TODO : write function to split colony  input is 1 colony, output is 2 colonys 1 colony is queenless 
  }
  if (change == "collapse") {
    
  }
}



















