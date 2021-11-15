
library(AlphaSimR)
source('C:/Users/jernejb/Desktop/BeeSim/lstrachan_honeybee_sim-main/HoneyBeeSim_functions.R') 

### nardim genome ###

founderGenomes = quickHaplo(nInd = 1000,
                             nChr = 16, 
                             segSites = 40,)

SP = SimParam$new(founderGenomes)

cebele = newPop(founderGenomes)
        
### nardim družine ### ce ardis 4 druzine ti jih po difolti da na 4 lokacije mors na roke jih nastvt na 1, morde bi blu smiselnu poenostavt kreiranje druzin da ne rabis na roke
# bo difolti bi blu dobru da je lokacija null
col1 = createColony(id = 1, 
                    queen = cebele[1],
                    fathers = cebele[2:18],
                    location = 1) # funkcija crate koloni mora avtomatsko nastavit virgin queen, èe nastaviš matico in fotre kot njihovo potomko

col2 = createColony(id = 2, 
                    queen = cebele[19],
                    fathers = cebele[19:35],
                    location = 1)

col3 = createColony(id = 3, 
                    queen = cebele[36],
                    fathers = cebele[37:50],
                    location = 1)

col4 = createColony(id = 4, 
                    queen = cebele[51],
                    fathers = cebele[52:70],
                    location = 1)

col5 = createColony(id = 5)

### add drones to colonies

col1@drones = createDrones(col1, nDrones = 1000)
col2@drones = createDrones(col2, nDrones = 1500)
col3@drones = createDrones(col3, nDrones = 1200)
col4@drones = createDrones(col4, nDrones = 1300)

### create DCA

dca = createDCA(list(col1, col2, col3, col4))

### adding workers to colones

col1@workers = createWorkers(col1, nWorkers = 10000)
col2@workers = createWorkers(col2, nWorkers = 15000)
col3@workers = createWorkers(col3, nWorkers = 12000)
col4@workers = createWorkers(col4, nWorkers = 13000)

### selecting fathers for further crossings of the colonies

fotr1 = selectFathersFromDCA(DCA = dca, nFathers = sample(10:25, size = 1))
fotr2 = selectFathersFromDCA(DCA = dca, nFathers = sample(10:25, size = 1))
fotr3 = selectFathersFromDCA(DCA = dca, nFathers = sample(10:25, size = 1))

### creating workers

col1@virgin_queens = createWorkers(col1, nWorkers = 1)
col2@virgin_queens = createWorkers(col2, nWorkers = 1)
col3@virgin_queens = createWorkers(col3, nWorkers = 1)
col4@virgin_queens = createWorkers(col4, nWorkers = 1)

### swarming, catching a swarm and crossing swarmed colony

produkt_rojenja = swarmColony(col1, perSwarm = sample((0.4:0.75), size = 1 )) # narest je treba èek, èe so prisotne virgin queens, èe jih ni te ne spusti èez
col1 = produkt_rojenja$newColony
col6 = produkt_rojenja$swarm # ne rabiš create colony kr je itak produkt colony

### adding workers and drones to swarm and swarmed colony

col1 = crossColony(col1, fathers = fotr1)
col1 = addWorkers(col1, nWorkersAdd = 5000)
col1 = addDrones(col1, nDronesAdd = 3000)                                   

col6 = addWorkers(col6, nWorkersAdd = 5000)
col6 = addDrones(col6, nDronesAdd = 700)

### supersedure

col2 = supersedeColony(col2) #ali dodati virgin queen into queen slot or and send a message the queen is not mated yet, ali samo spraznit queen slot...
                             #popraviti namesto swarm v drugi vrstici na colony
### splitting the colony, adding virgin to split and crss it

narejenc = splitColony(col3, perSplit = 0.3)
col3 = narejenc$colony
col7 = narejenc$splitColony #### rabmo funkcijo add virgin queen
col7@virgin_queens = col4@virgin_queens
col7 = crossColony(col7, fathers = selectFathersFromDCA(DCA = dca, nFathers = sample(10:25, size = 1)))### poglej kaj je gregor nalimu lauri v kanal

### killed the queen, empty her spermatheca

col7@queen = NULL## we need function kill the queen
col7@fathers = NULL

### add virgin queen and cross it

col7@queen = col4@virgin_queens
col7 = crossColony(col7, fathers = selectFathersFromDCA(DCA = dca, nFathers = sample(10:25, size = 1)))
col7


apiary <- list(col1, col2, col3, col4, col5, col6, col7)
class(apiary)

apiary[[1]]@queen

for(colony in apiary){
 if (colony@workers = 0){
   print('no workers')
 }else {
  nuc = splitColony(colony, perSplit = 0.3)
  
 # col7 = nuc$splitColony#### kle je treba še pogruntat kaku kreirat nove kolonije +1
  colony = nuc$colony ### bi mouglu bit ok, naredi, še da preskoèi prazne družine
  #apiary[length(apiary) + 1] = colony
}}



