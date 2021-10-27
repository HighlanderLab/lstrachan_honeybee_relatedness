library(AlphaSimR)
source("~/Desktop/GitHub/lstrachan_honeybee_sim/HoneyBeeSim_functions.R")

#Create founder population ----

founder_population = quickHaplo(nInd = 1000, nChr = 16, ploidy = 2L,
                                inbred = FALSE, segSites = 1000)

#Create base population
SP = SimParam$new(founder_population) 
SP$addTraitA(nQtlPerChr = 1000)
base = newPop(founder_population)

#Create colonies -made 5 colonies 
col1 = createColony(id = 1, queen = base[1])
col2 = createColony(id = 2, queen = base[200])
col3 = createColony(id = 3, queen = base [20])
col4 = createColony(id = 4, queen = base [53])
col5 = createColony(id = 5, queen = base [643])
#check colonies
#col?
#can't have fathers at this stage if you need to crossColony with DCA fathers


#populate drones - range from 200-300 drones per colony 
col1@drones = createDrones(colony = col1, nDrones = 200)
col2@drones = createDrones(colony = col2, nDrones = 300)
col3@drones = createDrones(colony = col3, nDrones = 250)
col4@drones = createDrones(colony = col4, nDrones = 225)
col5@drones = createDrones(colony = col5, nDrones= 275)
#check drones
#col?@drones


#createDCA
DCA <- createDCA(list(col1, col2, col3, col4, col5))
DCA

#select fathers from the DCA - range from 16-20 fathers 
fathers1 = selectFathersFromDCA(DCA, nFathers = 16)
fathers2 = selectFathersFromDCA(DCA, nFathers = 17)
fathers3 = selectFathersFromDCA(DCA, nFathers = 18)
fathers4 = selectFathersFromDCA(DCA, nFathers = 19)
fathers5 = selectFathersFromDCA(DCA, nFathers = 20)
#check fathers
#fathers?
  
#cross founder colonies with corresponding fathers , workers range 1000-1500, drones range 50-60
col1 = crossColony(colony = col1, fathers = fathers1, nWorkers = 1000, nDrones = 50)
col2 = crossColony(colony = col2, fathers = fathers2, nWorkers = 1250, nDrones = 60)
col3 = crossColony(colony = col3, fathers = fathers3, nWorkers = 1500, nDrones = 50)  
col4 = crossColony(colony = col4, fathers = fathers4, nWorkers = 1300, nDrones = 55)
col5 = crossColony(colony = col5, fathers = fathers5, nWorkers = 1400, nDrones = 58)
#check colony
#col? 

#First year ----
#colonies 1 +2 swarm, colony 3 supersedure,colony 4 nothing happens, colony 5 dies
#this year chance of swarm is 0.6 
swarm1 = swarmColony(col1, perSwarm=0.6)
col1 = swarm1$swarm #same ID, no location 
col6 = swarm1$newColony  #has new Id but location 1 
swarm2 = swarmColony(col2, perSwarm=0.6)
col2 = swarm2$swarm
col7 = swarm2$newColony


#col6 and col7 must both be crossed using crossColony 
fathers6 = selectFathersFromDCA(DCA, nFathers = 21)
col6 = crossColony(colony = col6, fathers = fathers6, nWorkers = 1200, nDrones =50)
fathers7 = selectFathersFromDCA(DCA, nFathers = 21)
col7 = crossColony(colony = col7, fathers = fathers7, nWorkers = 1400, nDrones =55)

#build up swarmed colonies (col1, col2)
col1 = addWorkers(col1, nWorkersAdd = 500)
col1 = addDrones(col1, nDrones = 50)
col2 = addWorkers(col2, nWorkersAdd = 500)
col2 = addDrones(col2, nDrones = 50)

#col3 supersedure 
col3 = supersedeColony(col3)

#TODO: create wintercolonylosses fuction - removes from colonyList 