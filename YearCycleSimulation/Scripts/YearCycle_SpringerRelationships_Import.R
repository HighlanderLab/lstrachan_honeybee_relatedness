#setwd("~/Desktop/GitHub/lstrachan_honeybee_sim/YearCycleSimulation")
# Clean workspace
rm(list = ls())

# Define functions
computeRelationship_genomic <- function(x, csd = TRUE,  ColonyAF = FALSE, AF0.5 = FALSE,
                                        SingleBaseAF = NULL, MultiBaseAF = NULL,
                                        SingleBaseAFcsdLocus = NULL, MultiBaseAFcsdLocus = NULL,
                                        SingleBaseAFcsdChr = NULL, MultiBaseAFcsdChr = NULL) {
  if (isColony(x)) {
    # Build the colony up to 1,000 workers and 200 drones
    colony <- buildUp(x = x,
                      nWorkers = 1000,
                      nDrones = 200,
                      exact = TRUE)
    # Extract the genotypes of all the colony members
    geno <- getSegSiteGeno(colony, caste = "all", collapse = TRUE)
    # Collect the sex of all the colony members in the order than genotypes
    sex <- getCasteSex(colony, caste = "all", collapse = TRUE)
    #get colony IDs
    id <- getCasteId(colony, caste = "all")

    #Determine sister type:
    CarWorkerFathers <- colony@workers@father #collect the fathers information for the workers
    FathersPop <- colony@queen@misc[[1]][["fathers"]] # Extract father Pop object
    WorkersFatherTable <- data.frame(workers= colony@workers@id, fathers = colony@workers@father)
    WorkersFatherTable$DPQ <-  sapply(WorkersFatherTable$father, FUN = function(x) FathersPop[x]@mother)

  } else if (isPop(x)) {
    geno <- getSegSiteGeno(x)
    sex <- x@sex
    id <- x@id
    # Queens populations- making a null object to simplify for now
    CarWorkerFathers = NULL
    WorkersFatherTable = NULL
    }

  ##################################################
  # IBS
  ##################################################
  # TODO: Combine generations 1 and 10 - so that they have the same reference population (regarding allele frequencies)
  print("IBS")
  print(Sys.time())

  ### --- WHOLE GENOME --- ###
  # ALLELE FREQUENCY 0.5
  if (AF0.5) {
    ibsAF0.5 <- calcBeeGRMIbs(x = geno,
                              sex = sex,
                              alleleFreq = rep(0.5, ncol(geno))) #modified allele freq
  } else {
    ibsAF0.5 <- NULL
  }

  #COLONY ALLELE FREQUENCY
  if (ColonyAF){
    ColonyAF_wg = calcBeeAlleleFreq(x = geno, sex = sex)
    ibsColonyAF <- calcBeeGRMIbs(x = geno,
                                 sex = sex,
                                 alleleFreq = ColonyAF_wg)
  } else {
    ibsColonyAF <- NULL
  }

  #SINGLE POPULATION BASE QUEENS ALLELE FREQUENCY
  if (is.null(SingleBaseAF)){
    ibsSingleBaseAF <- NULL
  } else {
    ibsSingleBaseAF <- calcBeeGRMIbs(x = geno,
                                     sex = sex,
                                     alleleFreq = SingleBaseAF)
  }

  #MULIPLE POPULATION BASE QUEENS ALLELE FREQUENCY
  if (is.null(MultiBaseAF)){
    ibsMultiBaseAF <- NULL
  } else {
    ibsMultiBaseAF <- calcBeeGRMIbs(x = geno,
                                    sex = sex,
                                    alleleFreq = MultiBaseAF)
  }

  if (csd) {
    ### --- CSD CHROMOSOME --- ###
    csdChr <- SP$csdChr

    if (isColony(x)) {
    csdChrGeno <- pullMarkerGeno(getCastePop(colony, caste = "all", collapse = TRUE),
                                 markers = names(SP$genMap[[SP$csdChr]]))
    sex <- getCasteSex(colony, caste = "all", collapse = T)
    } else {
      csdChrGeno <- pullMarkerGeno(x, markers = names(SP$genMap[[SP$csdChr]]))
      sex <- getCasteSex(x, caste = "all", collapse = T)
    }

    # ALLELE FREQUENCY 0.5
    if (AF0.5) {
      ibsAF0.5_csdChr <- calcBeeGRMIbs(x = csdChrGeno,
                                       sex = sex,
                                       alleleFreq = rep(0.5, length(SP$genMap[[SP$csdChr]])))
    } else {
      ibsAF0.5_csdChr <- NULL
    }

    #COLONY ALLELE FREQUENCY
    if (ColonyAF){
      ColonyAF_csdChr = calcBeeAlleleFreq(x = csdChrGeno, sex = sex)
      ibsColonyAF_csdChr <- calcBeeGRMIbs(x = csdChrGeno,
                                          sex = sex,
                                          alleleFreq = ColonyAF_csdChr)
    } else {
      ibsColonyAF_csdChr <- NULL
    }

    #SINGLE POPULATION BASE QUEENS ALLELE FREQUENCY
    if (is.null(SingleBaseAFcsdChr)) {
      ibsSingle_csdChr <- NULL
    } else {
      ibsSingle_csdChr <- calcBeeGRMIbs(x =  csdChrGeno,
                                        sex = sex,
                                        alleleFreq = SingleBaseAFcsdChr)
    }

    #MULIPLE POPULATION BASE QUEENS ALLELE FREQUENCY
    if (is.null(MultiBaseAFcsdChr)) {
      ibsMulti_csdChr <- NULL
    } else {
      ibsMulti_csdChr <- calcBeeGRMIbs(x =  csdChrGeno,
                                       sex = sex,
                                       alleleFreq = MultiBaseAFcsdChr)
    }


    #### --- CSD LOCUS ---###
    if (isColony(x)) {
    csdLocusGeno <- getCsdGeno(colony, caste = "all", collapse = TRUE)
    } else {
      csdLocusGeno <- getCsdGeno(x, caste = "all", collapse = TRUE)
    }

    # ALLELE FREQUENCY 0.5
    if (AF0.5) {
      ibsAF0.5_csdLocus <- calcBeeGRMIbs(x = csdLocusGeno,
                                         sex = sex,
                                         alleleFreq = rep(0.5, length(SP$csdPosStart:SP$csdPosStop)))
    } else {
      ibsAF0.5_csdLocus <- NULL
    }

    #COLONY ALLELE FREQUENCY
    if (ColonyAF){
      ColonyAF_csdLocus = calcBeeAlleleFreq(x = csdLocusGeno, sex = sex)
      ibsColonyAF_csdLocus <- calcBeeGRMIbs(x = csdLocusGeno,
                                            sex = sex,
                                            alleleFreq = ColonyAF_csdLocus)
    } else {
      ibsColonyAF_csdLocus <- NULL
    }

    #MULIPLE POPULATION BASE QUEENS ALLELE FREQUENCY
    if (is.null(MultiBaseAFcsdLocus)){
      ibsMultiBaseCsdLocus <- NULL
    } else {
      ibsMultiBaseCsdLocus <- calcBeeGRMIbs(x = csdLocusGeno,
                                            sex = sex,
                                            alleleFreq = MultiBaseAFcsdLocus)
    }

    #SINGLE POPULATION BASE QUEENS ALLELE FREQUENCY
    if (is.null(SingleBaseAFcsdLocus)) {
      ibsSingleBaseCsdLocus <- NULL
    } else {
      ibsSingleBaseCsdLocus <- calcBeeGRMIbs(x = csdLocusGeno,
                                             sex = sex,
                                             alleleFreq = SingleBaseAFcsdLocus)
    }
  }
  ##################################################
  # IBD realized
  ##################################################
  # Compute the IBD relationship matrices
  print("IBD")
  print(Sys.time())
  ### --- WHOLE GENOME --- ###
  if (isColony(x)) {
    haplo <- getIbdHaplo(colony, caste = "all", collapse = TRUE)
  } else if (isPop(x)) {
    haplo <- getIbdHaplo(x)
  }
  ibd <- calcBeeGRMIbd(x = haplo)
  ibd <- ibd$indiv

  if (csd) {
    ### --- CSD CHROMOSOME --- ###
    print("csd IBD")
    print(Sys.time())
    if (isColony(x)) {
    csdChrHaplo <- getIbdHaplo(x = colony, caste = "all", collapse = TRUE, chr = SP$csdChr)
    } else {
      csdChrHaplo <- getIbdHaplo(x = x, caste = "all", collapse = TRUE, chr = SP$csdChr)
    }
    ibd_csdChr <- calcBeeGRMIbd(x = csdChrHaplo)
    ibd_csdChr <- ibd_csdChr$indiv

    ### --- CSD LOCUS --- ###
    ibd_csdLocus <- calcBeeGRMIbd(x = haplo[, paste(SP$csdChr,
                                                    SP$csdPosStart:SP$csdPosStop,
                                                    sep = "_")])
    ibd_csdLocus <- ibd_csdLocus$indiv
  }

  if (csd){
    return(list(IBSAF0.5 = ibsAF0.5, IBScolonyAF = ibsColonyAF, IBSsingleAF = ibsSingleBaseAF, IBSmultiAF = ibsMultiBaseAF, IBD = ibd,
                IBSAF0.5CSDChr = ibsAF0.5_csdChr, IBScolonyAF_CSDChr = ibsColonyAF_csdChr, IBSsingleCSDChr = ibsSingle_csdChr, IBSmultiCSDChr = ibsMulti_csdChr, IBDcsdChr = ibd_csdChr,
                IBSAF0.5CSDLocus = ibsAF0.5_csdLocus, colonyAF_CSDLocus = ibsColonyAF_csdLocus, IBSsingleAFcsdLocus = ibsSingleBaseCsdLocus, IBSmultiAFCsdLocus = ibsMultiBaseCsdLocus, IBDCsdLocus = ibd_csdLocus,
                ID = id, CarWorkerFathers = CarWorkerFathers, WorkersFatherTable = WorkersFatherTable))
  } else {
    return(list(IBSAF0.5 = ibsAF0.5, IBScolonyAF = ibsColonyAF, IBSsingleAF = ibsSingleBaseAF, IBSmultiAF = ibsMultiBaseAF, IBD = ibd,
                IBSAF0.5CSDChr = NULL, IBScolonyAF_CSDChr = NULL, IBSsingleCSDChr = NULL, IBSmultiCSDChr = NULL, IBDcsdChr = NULL,
                IBSAF0.5CSDLocus = NULL, IBScolonyAF_CSDChr = NULL, IBSsingleAFcsdLocus = NULL, IBSmultiAFCsdLocus = NULL, IBDCsdLocus = NULL,
                ID = id, CarWorkerFathers = CarWorkerFathers, WorkersFatherTable = WorkersFatherTable))
  }
}

computeRelationship_pedigree <- function(pedigree) {
  pedigree <- as.data.frame(pedigree)
  # nadiv needs missing as NA
  pedigree$mother[pedigree$mother == 0] <- NA
  pedigree$father[pedigree$father == 0] <- NA
  pedigree$ID <- rownames(pedigree)

  # Females are 1, males are 0
  colnames(pedigree) <- c("Dam", "Sire", "Sex", "ID")
  # The order for nadiv should be ID, Dam, Sire, Sex
  pedigree <- pedigree[, c("ID", "Dam", "Sire", "Sex")]
  pedigree$Sire[pedigree$Sex == 1] <- NA

  tmp <- makeS(pedigree = pedigree, heterogametic = "1", returnS = FALSE)
  IBDe <- tmp #$S
  #  dimnames(IBDe) <- list(rownames(pedigree), rownames(pedigree))
  return(IBDe)
}

maintainApiarySize <- function(age0 = NULL, age1 = NULL) {
  if ((nColonies(age0) + nColonies(age1)) > apiarySize) { # check if the sum of all colonies is greater than apiary size
    IDsplits <- getId(age0)[hasSplit(age0)] # get the IDs of age 0 that are splits
    splits0 <- pullColonies(age0, ID = IDsplits) # pull the splits out of age 0
    age0split <- splits0$pulled # create an object for age 0 splits
    age0swarm <- splits0$remnant # create an object for swarms and superseded colonies
    age0needed <- apiarySize - nColonies(age1) # calculate the number of age 0 colonies that are needed to fill up the apiary
    splitsNeeded <- age0needed - nColonies(age0swarm) # calculate the number of splits needed
    if (age0needed <= nColonies(age0swarm)) { # check if the number of age 0 colonies needed is lower or equal to age 0 swarms
      swarmID <- sample(getId(age0swarm), age0needed) # if yes, select the ids of swarms that will stay in apiary
      swarmTMP <- pullColonies(age0swarm, ID = swarmID) # pull out those selected age0 swarms
      age0 <- swarmTMP$pulled # put selected swarms to age 0 object
    } else if (age0needed > nColonies(age0swarm)) { # in case when age 0 needed is grater than number of swarm select splits
      nSplitsNeeded <- age0needed - nColonies(age0swarm) # calculate the number of splits needed
      splitId <- sample(getId(age0split), nSplitsNeeded) # select ids of splits
      splitTmp <- pullColonies(age0split, ID = splitId) # pull the splits
      splits <- splitTmp$pulled # select pulled splits
      age0 <- c(age0swarm, splits) # combine splits and swarms in age 0 object
    }
    return(age0)
  }
}
# Load packages
library(AlphaSimR)
library(ggplot2)
library(tictoc)
library(R6)
library(nadiv)
library(Matrix)
library(SIMplyBee)
library(dplyr)
library(tidyr)
# TODO: replace with devtools installation from Github once the package is operational
# Source the development version of AlphaSimR


# Founder population parameters -------------------------------------------------------------------
nMelN = 800
nCar = 400
nChr = 16
nDronesPerQueen = 50
nSegSites = 1000
# Population parameters -------------------------------------------------------------------
# Number of repeats
nRep <- 1
# Number of years
nYear <- 10
# Number of colonies in the apiary
apiarySize <- 300
# Number of workers in a full colony
nWorkers <- 10 # TODO: change to 20K
# Number of drones in a full colony
nDrones <- 50 #nWorkers * 0.2
# Number of drones the queen mates with (could also be a function)
pFathers <- nFathersPoisson
# Number of created virgin queens
nVirginQueens <- 1

# Period parameters -------------------------------------------------------------------
# Period 1 (spring)
# Percentage of colonies that swarm in period 1
p1swarm <- 0.05
# Percentage of colonies that supersede in period 1
p1supersede <- 0.05
# Percentage of colonies that collapse in period 1
p1collapse <- 0.10

# Period2 (summer)
# Percentage of colonies that swarm in period 2
p2swarm <- 0.01
# Percentage of colonies that supersede in period 2
p2supersede <- p1supersede
# Percentage of colonies that collapse in period 2
p2collapse <- p1collapse

# Period3 (winter)
# Percentage of age 0 colonies that collapse in period 3
p3collapseAge0 <- 0.25
# Percentage of age 2 colonies that collapse in period 3
p3collapseAge1 <- 0.3

# Percentage import from carnica to mellifera
pImport <- 0.5

# Create data frames for recording the number of age0 and age1 colonies, csd variability and for recording cpu time
loopTime <- data.frame(Rep = NA, tic = NA, toc = NA, msg = NA, time = NA)

# Prepare recording function
data_rec <- function(datafile, colonies, year, population) {
  queens = mergePops(getQueen(colonies))
  datafile = rbind(datafile,
                   data.frame(colonies             = deparse(substitute(colonies)),
                              population           = population,
                              year                 = year,
                              Id                   = queens@id,
                              MId                  = queens@mother,
                              FId                  = queens@father,
                              nFathers             = nFathers(queens),
                              nDPQ                 = sapply(getFathers(queens), function(x) length(unique(x@mother))),
                              nCsdAlColony         = sapply(colonies@colonies, function(x) nCsdAlleles(x, collapse = TRUE)),
                              nCsdApiary           = rep(nCsdAlleles(colonies, collapse = TRUE), queens@nInd),
                              pHomBrood            = calcQueensPHomBrood(queens),
                              gvQueens_QueenTrait  = sapply(getGv(colonies, caste = "queen"), function(x) x[1,1]),
                              gvQueens_WorkerTrait = sapply(getGv(colonies, caste = "queen"), function(x) x[1,2])
                   ))}
colonyRecords = NULL

# Start of the rep-loop ---------------------------------------------------------------------
for (Rep in 1:nRep) {
  # Rep <- 1
  cat(paste0("Rep: ", Rep, "/", nRep, "\n"))
  # Measure cpu time
  tic(paste0(nYear, 'y loop'))
  # Start profiling
  Rprof()


  # Founder population ---------------------------------------------------------
  # Create a founder population of A. m. mellifera and A. m. carnica bees
  #founderGenomes <- simulateHoneyBeeGenomes(nMelN = nMelN,
  #                                           nCar = nCar,
  #                                           nChr = nChr,
  #                                           nSegSites = nSegSites)
  #
  #save(founderGenomes, file="founderGenomes_ThreePop.RData")
  print("Loading in the founderData")
  load("FounderGenomes_ThreePop_16chr.RData")
  #load("~/Desktop/GitHub/lstrachan_honeybee_sim/YearCycleSimulation/PlottingData/FounderGenomes_ThreePop_16chr.RData")
  # Create SP object and write in the global simulation/population parameters
  SP <- SimParamBee$new(founderGenomes, csdChr = ifelse(nChr >= 3, 3, 1), nCsdAlleles = 128)
  SP$nWorkers <- nWorkers
  SP$nDrones <- nDrones
  SP$nFathers <- pFathers
  SP$nVirginQueens <- nVirginQueens
  SP$swarmP <- 0.5
  SP$splitP <- 0.3
  # Track the pedigree
  SP$setTrackPed(TRUE)
  # Track the recombination
  SP$setTrackRec(TRUE)
  # Add a SNP chip with 3 SNPs per chromosome
  SP$addSnpChip(nSnpPerChr = 3)
  # define csd chromomsome
  csdChr <- SP$csdChr
  # Add traits - taken from the QuantGen vignettte
  mean <- c(20, 0)
  varA <- c(1, 1 / SP$nWorkers)
  corA <- matrix(data = c( 1.0, -0.5,
                           -0.5,  1.0), nrow = 2, byrow = TRUE)
  SP$addTraitA(nQtlPerChr = 100, mean = mean, var = varA, corA = corA,
               name = c("queenTrait", "workersTrait"))

  varE <- c(3, 3 / SP$nWorkers)
  # TODO: what is a reasonable environmental correlation between queen and worker effects?
  corE <- matrix(data = c(1.0, 0.3,
                          0.3, 1.0), nrow = 2, byrow = TRUE)
  SP$setVarE(varE = varE, corE = corE)


  # Create a base population for A. m. mellifera, A. m. mellifera cross, and A. m. carnica
  virginQueens <- list(Mel = createVirginQueens(x = founderGenomes[1:(nMelN/2)]),
                       MelCross = createVirginQueens(x = founderGenomes[(nMelN/2 + 1):nMelN]),
                       Car = createVirginQueens(x = founderGenomes[(nMelN +1):(nMelN + nCar)]))
  # Create drones for A. m. mellifera, A. m. mellifera cross, and A. m. carnica
  drones <- list(Mel = createDrones(x = virginQueens$Mel[(apiarySize+1):(nMelN/2)], nInd = nDronesPerQueen),
                 MelCross = createDrones(x = virginQueens$MelCross[(apiarySize+1):(nMelN/2)], nInd = nDronesPerQueen),
                 Car = createDrones(x = virginQueens$Car[(apiarySize+1):nCar], nInd = nDronesPerQueen))
  # Get fathers for Mel, MelCross and Car
  fathersMel <- pullDroneGroupsFromDCA(drones$Mel, n = nInd(virginQueens$Mel[1:apiarySize]), nDrones = nFathersPoisson)
  fathersMelCross <- pullDroneGroupsFromDCA(drones$MelCross, n = nInd(virginQueens$MelCross[1:apiarySize]), nDrones = nFathersPoisson)
  fathersCar <- pullDroneGroupsFromDCA(drones$Car, n = nInd(virginQueens$Car[1:apiarySize]), nDrones = nFathersPoisson)
  # Force two drones for the first colony to be brothers (add two drones from the same mother)
  fathersMel[[1]] <- c(fathersMel[[1]], createDrones(virginQueens$Mel[fathersMel[[1]]@mother[1]], 2))
  fathersMelCross[[1]] <- c(fathersMelCross[[1]], createDrones(virginQueens$MelCross[fathersMelCross[[1]]@mother[1]], 2))
  fathersCar[[1]] <- c(fathersCar[[1]], createDrones(virginQueens$Car[fathersCar[[1]]@mother[1]], 2))

  # Mate virgin queens with fathers
  queens <- list(Mel = SIMplyBee::cross(x = virginQueens$Mel[1:apiarySize], drones = fathersMel),
                 MelCross = SIMplyBee::cross(x = virginQueens$MelCross[1:apiarySize], drones = fathersMelCross),
                 Car = SIMplyBee::cross(x = virginQueens$Car[1:apiarySize], drones = fathersCar))

  #Set allele frequency for queens
  tmp <- c(virginQueens$Mel, virginQueens$Car)
  alleleFreqBaseQueens <- calcBeeAlleleFreq(x = getSegSiteGeno(tmp),
                                            sex = tmp@sex)

  alleleFreqBaseQueensCar <- calcBeeAlleleFreq(x = getSegSiteGeno(virginQueens$Car),
                                               sex = virginQueens$Car@sex)

  alleleFreqBaseQueensMel <- calcBeeAlleleFreq(x = getSegSiteGeno(virginQueens$Mel),
                                               sex = virginQueens$Mel@sex)
  #Get allele freq for csd locus
  csdLocus <- paste0(SP$csdChr, "_", SP$csdPosStart:SP$csdPosStop)
  alleleFreqCsdLocusBaseQueens <- alleleFreqBaseQueens[csdLocus]
  alleleFreqCsdLocusBaseCar <- alleleFreqBaseQueensCar[csdLocus]
  alleleFreqCsdLocusBaseMel <- alleleFreqBaseQueensMel[csdLocus]

  #Get allele freq for csd Chromosome - this pulls out only the 3rd chromosome
  alleleFreqCsdChrBaseQueens <- t(as.data.frame(alleleFreqBaseQueens))[, grepl(pattern = paste0("^", csdChr, "_"), x = colnames(t(as.data.frame(alleleFreqBaseQueens))))] %>% t()
  alleleFreqCsdChrBaseCar <- t(as.data.frame(alleleFreqBaseQueensCar))[, grepl(pattern = paste0("^", csdChr, "_"), x = colnames(t(as.data.frame(alleleFreqBaseQueensCar))))] %>% t()
  alleleFreqCsdChrBaseMel <- t(as.data.frame(alleleFreqBaseQueensMel))[, grepl(pattern = paste0("^", csdChr, "_"), x = colnames(t(as.data.frame(alleleFreqBaseQueensMel))))] %>% t()


  # Start the year-loop ------------------------------------------------------------------
  for (year in 1:nYear) {
    print("Starting the cycle")
    #year <- 1
    #year <- year + 1
    cat(paste0("Year: ", year, "/", nYear, "\n"))
    # If this is the first year, create some colonies to start with
    # If not, promote the age0 to age1, age1 to age2 and remove age2 colonies
    if (year == 1) {
      print("Creating initial colonies")
      age1 <- list(Mel = createMultiColony(x = queens$Mel, n = apiarySize),
                   MelCross = createMultiColony(x = queens$MelCross, n = apiarySize),
                   Car = createMultiColony(x = queens$Car, n = apiarySize))
      print("Record initial colonies")
      colonyRecords <- data_rec(datafile = colonyRecords, colonies = age1$Mel, year = year, population = "Mel")
      colonyRecords <- data_rec(datafile = colonyRecords, colonies = age1$MelCross, year = year, population = "MelCross")
      colonyRecords <- data_rec(datafile = colonyRecords, colonies = age1$Car, year = year, population = "Car")

    } else {
      age2 <- list(Mel = age1$Mel, MelCross = age1$MelCross, Car = age1$Car)
      age1 <- list(Mel = age0$Mel, MelCross = age0$MelCross, Car = age0$Car)
      age0 <- list(Mel = NULL, MelCross = NULL, Car = NULL)
      age0p1 <- list(Mel = NULL, MelCross = NULL, Car = NULL)
      age0p2 <- list(Mel = NULL, MelCross = NULL, Car = NULL)
    }

    #In year 1, inspect the relationship in one of the colonies
    if (year == 1) {
      print("Computing carnica generation1 colony relationship")
      print(Sys.time())
      springerColony1_Car <- computeRelationship_genomic(x = age1$Car[[1]],
                                                         csd = isCsdActive(SP),
                                                         ColonyAF = TRUE,
                                                         AF0.5 = TRUE,
                                                         SingleBaseAF = alleleFreqBaseQueensCar,
                                                         MultiBaseAF = alleleFreqBaseQueens,
                                                         SingleBaseAFcsdLocus = alleleFreqCsdLocusBaseCar,
                                                         MultiBaseAFcsdLocus  = alleleFreqCsdLocusBaseQueens,
                                                         SingleBaseAFcsdChr = alleleFreqCsdChrBaseCar,
                                                         MultiBaseAFcsdChr = alleleFreqCsdChrBaseQueens)
       print("Computing queens generation1 colony relationship")
      print(Sys.time())
      springerQueens1_MelAF <- computeRelationship_genomic(x = c(queens$Mel, queens$MelCross, queens$Car),
                                                     csd = isCsdActive(SP),
                                                     ColonyAF = TRUE,
                                                     AF0.5 = TRUE,
                                                     SingleBaseAF = alleleFreqBaseQueensMel,
                                                     MultiBaseAF = alleleFreqBaseQueens,
                                                     SingleBaseAFcsdLocus = alleleFreqCsdLocusBaseMel,
                                                     MultiBaseAFcsdLocus  = alleleFreqCsdLocusBaseQueens,
                                                     SingleBaseAFcsdChr = alleleFreqCsdChrBaseMel,
                                                     MultiBaseAFcsdChr = alleleFreqCsdChrBaseQueens)

      springerQueens1_CarAF <- computeRelationship_genomic(x = c(queens$Mel, queens$MelCross, queens$Car),
                                                           csd = isCsdActive(SP),
                                                           ColonyAF = TRUE,
                                                           AF0.5 = TRUE,
                                                           SingleBaseAF = alleleFreqBaseQueensCar,
                                                           MultiBaseAF = alleleFreqBaseQueens,
                                                           SingleBaseAFcsdLocus = alleleFreqCsdLocusBaseCar,
                                                           MultiBaseAFcsdLocus  = alleleFreqCsdLocusBaseQueens,
                                                           SingleBaseAFcsdChr = alleleFreqCsdChrBaseCar,
                                                           MultiBaseAFcsdChr = alleleFreqCsdChrBaseQueens)

      springerQueensPop1 <- rbind(data.frame(ID = sapply(getQueen(age1$Mel), FUN = function(x) x@id), Pop = "Mel"),
                                  data.frame(ID = sapply(getQueen(age1$MelCross), FUN = function(x) x@id), Pop = "MelCross"),
                                  data.frame(ID = sapply(getQueen(age1$Car), FUN = function(x) x@id), Pop = "Car"))

    }
    print("Done computing initial relationships")

    # Period1 ------------------------------------------------------------------
    # Build-up the colonies
    print(paste0("Building up the colonies to ", nWorkers, " and ", nDrones))
    print(Sys.time())
    age1 <- list(Mel = buildUp(age1$Mel),
                 MelCross = buildUp(age1$MelCross),
                 Car = buildUp(age1$Car))
    if (year > 1) {
      age2 <- list(Mel = buildUp(age2$Mel),
                   MelCross = buildUp(age2$MelCross),
                   Car = buildUp(age2$Car))
    }

    # Split all age1 colonies
    print("Splitting the colonies")
    print(Sys.time())
    tmp <- list(Mel = split(age1$Mel),
                MelCross = split(age1$MelCross),
                Car = split(age1$Car))
    age1 <- list(Mel = tmp$Mel$remnant,
                 MelCross = tmp$MelCross$remnant,
                 Car = tmp$Car$remnant)
    # The queens of the splits are 0 years old
    age0p1 <- list(Mel = tmp$Mel$split, MelCross = tmp$MelCross$split, Car = tmp$Car$split)

    if (year > 1) {
      # Split all age2 colonies
      tmp <- list(Mel = split(age2$Mel),
                  MelCross = split(age2$MelCross),
                  Car = split(age2$Car))
      age2 <- list(Mel = tmp$Mel$remnant,
                   MelCross = tmp$MelCross$remnant,
                   Car = tmp$Car$remnant)
      # The queens of the splits are 0 years old
      age0p1 <- list(Mel = c(age0p1$Mel, tmp$Mel$split),
                     MelCross = c(age0p1$MelCross, tmp$MelCross$split),
                     Car = c(age0p1$Car, tmp$Car$split))
    }

    # Create virgin queens
    # Sample colony for the virgin queens
    # TODO: this is likely an inbreeding bottleneck - just one colony producing all virgin queens!
    #       maybe we should use selectColonies() and then create virginQueens from them, but we
    #       need createVirginQueens() to work with multiple colonies too, and equally for createWorkers()
    #       and createDrones()
    print("Create virgin queens, period 1")
    print(Sys.time())
    virginDonor <- list(Mel = sample.int(n = nColonies(age1$Mel), size = 1),
                        MelCross = sample.int(n = nColonies(age1$MelCross), size = 1),
                        Car = sample.int(n = nColonies(age1$Car), size = 1))
    # Virgin queens for splits!
    virginQueens <- list(Mel = createVirginQueens(age1$Mel[[virginDonor$Mel]], nInd = nColonies(age0p1$Mel)),
                         MelCross = createVirginQueens(age1$MelCross[[virginDonor$MelCross]], nInd = nColonies(age0p1$MelCross)),
                         Car = createVirginQueens(age1$Car[[virginDonor$Car]], nInd = nColonies(age0p1$Car)))

    # Requeen the splits --> queens are now 0 years old
    age0p1 <- list(Mel = reQueen(age0p1$Mel, queen = virginQueens$Mel),
                   MelCross = reQueen(age0p1$MelCross, queen = virginQueens$MelCross),
                   Car = reQueen(age0p1$Car, queen = virginQueens$Car))

    # Swarm a percentage of age1 colonies
    print("Swarm colonies, P1")
    print(Sys.time())
    tmp <- list(Mel = pullColonies(age1$Mel, p = p1swarm),
                MelCross = pullColonies(age1$MelCross, p = p1swarm),
                Car = pullColonies(age1$Car, p = p1swarm))
    age1 <- list(Mel = tmp$Mel$remnant,
                 MelCross = tmp$MelCross$remnant,
                 Car = tmp$Car$remnant)
    tmp <- list(Mel = swarm(tmp$Mel$pulled),
                MelCross = swarm(tmp$MelCross$pulled),
                Car = swarm(tmp$Car$pulled))
    age0p1 <- list(Mel = c(age0p1$Mel, tmp$Mel$remnant),
                   MelCross = c(age0p1$MelCross, tmp$MelCross$remnant),
                   Car = c(age0p1$Car, tmp$Car$remnant))
    age1 <- list(Mel = c(age1$Mel, tmp$Mel$swarm),
                 MelCross = c(age1$MelCross, tmp$MelCross$swarm),
                 Car = c(age1$Car, tmp$Car$swarm))


    if (year > 1) {
      # Swarm a percentage of age2 colonies
      tmp <- list(Mel = pullColonies(age2$Mel, p = p1swarm),
                  MelCross = pullColonies(age2$MelCross, p = p1swarm),
                  Car = pullColonies(age2$Car, p = p1swarm))
      age2 <- list(Mel = tmp$Mel$remainingColonies,
                   MelCross = tmp$MelCross$remnant,
                   Car = tmp$Car$remnant)
      tmp <- list(Mel = swarm(tmp$Mel$pulled),
                  MelCross = swarm(tmp$MelCross$pulled),
                  Car = swarm(tmp$Car$pulled))
      age0p1 <- list(Mel = c(age0p1$Mel, tmp$Mel$remnant),
                     MelCross = c(age0p1$MelCross, tmp$MelCross$remnant),
                     Car = c(age0p1$Car, tmp$Car$remnant))
      age2 <- list(Mel = c(age2$Mel, tmp$Mel$swarm),
                   MelCross = c(age2$MelCross, tmp$MelCross$swarm),
                   Car = c(age2$Car, tmp$Car$swarm))
    }

    # Supersede age1 colonies
    print("Supersede colonies, P1")
    print(Sys.time())
    tmp <- list(Mel = pullColonies(age1$Mel, p = p1supersede),
                MelCross = pullColonies(age1$MelCross, p = p1supersede),
                Car = pullColonies(age1$Car, p = p1supersede))
    age1 <- list(Mel = tmp$Mel$remnant,
                 MelCross = tmp$MelCross$remnant,
                 Car = tmp$Car$remnant)
    tmp <- list(Mel = supersede(tmp$Mel$pulled),
                MelCross = supersede(tmp$MelCross$pulled),
                Car = supersede(tmp$Car$pulled))
    age0p1 <- list(Mel = c(age0p1$Mel, tmp$Mel),
                   MelCross = c(age0p1$MelCross, tmp$MelCross),
                   Car = c(age0p1$Car, tmp$Car))


    if (year > 1) {
      # Supersede age2 colonies
      tmp <- list(Mel = pullColonies(age2$Mel, p = p1supersede),
                  MelCross = pullColonies(age2$MelCross, p = p1supersede),
                  Car = pullColonies(age2$Car, p = p1supersede))
      age2 <- list(Mel = tmp$Mel$remnant,
                   MelCross = tmp$MelCross$remnant,
                   Car = tmp$Car$remnant)
      tmp <- list(Mel = supersede(tmp$Mel$pulled),
                  MelCross = supersede(tmp$MelCross$pulled),
                  Car = supersede(tmp$Car$pulled))
      age0p1 <- list(Mel = c(age0p1$Mel, tmp$Mel),
                     MelCross = c(age0p1$MelCross, tmp$MelCross),
                     Car = c(age0p1$Car, tmp$Car))
    }

    # Mate the split colonies
    print("Mate split colonies, P1")
    print(Sys.time())
    if (year == 1) {
      DCAMel <- createDCA(age1$Mel)
      age0p1$Mel <- cross(age0p1$Mel, drones = pullDroneGroupsFromDCA(DCA = DCAMel, n = nColonies(age0p1$Mel), nDrones = nFathersPoisson))
      DCAMelCross <- createDCA(c(age1$MelCross,
                                 selectColonies(age1$Car, n = round(nColonies(age1$MelCross) * pImport, 0)))) #Is this 0 the P argument?
      age0p1$MelCross <- cross(age0p1$MelCross, drones = pullDroneGroupsFromDCA(DCA = DCAMelCross, n = nColonies(age0p1$MelCross), nDrones = nFathersPoisson))
      DCACar <- createDCA(age1$Car)
      age0p1$Car <- cross(age0p1$Car, drones = pullDroneGroupsFromDCA(DCA = DCACar, n = nColonies(age0p1$Car), nDrones = nFathersPoisson))
    } else {
      #Mel
      DCAMel <- createDCA(c(age1$Mel, age2$Mel))
      fathersMel <-  pullDroneGroupsFromDCA(DCA = DCAMel, n = nColonies(age0p1$Mel), nDrones = nFathersPoisson)
      for (x in 1:length(fathersMel)){
        fathersMel[[x]] <- c(fathersMel[[x]], createDrones(age1$Mel[[1]], nInd = 2))
      }
      age0p1$Mel <- cross(age0p1$Mel, drones = fathersMel)
      
      #MelCross
      DCAMelCross <- createDCA(c(age1$MelCross,
                                 selectColonies(age1$Car, n = round(nColonies(age1$MelCross) * pImport, 0)),
                                 age2$MelCross,
                                 selectColonies(age2$Car, n = round(nColonies(age2$MelCross) * pImport, 0))))
      fathersMelCross <-  pullDroneGroupsFromDCA(DCA = DCAMelCross, n = nColonies(age0p1$MelCross), nDrones = nFathersPoisson)
      for (x in 1:length(fathersMelCross)){
        fathersMelCross[[x]] <- c(fathersMelCross[[x]], createDrones(age1$MelCross[[1]], nInd = 2))
      }
      age0p1$MelCross <- cross(age0p1$MelCross, drones = fathersMelCross)
      
      #Car
      DCACar <- createDCA(c(age1$Car, age2$Car))
      fathersCar <-  pullDroneGroupsFromDCA(DCA = DCACar, n = nColonies(age0p1$Car), nDrones = nFathersPoisson)
      for (x in 1:length(fathersCar)){
        fathersCar[[x]] <- c(fathersCar[[x]], createDrones(age1$Car[[1]], nInd = 2))
      }
      age0p1$Car <- cross(age0p1$Car, drones = fathersCar)
    }
    
    # Collapse
    print("Collapse colonies, P1")
    print(Sys.time())
    age1 <- list(Mel = selectColonies(age1$Mel, p = 1 - p1collapse),
                 MelCross = selectColonies(age1$MelCross, p = 1 - p1collapse),
                 Car = selectColonies(age1$Car, p = 1 - p1collapse))
    if (year > 1) {
      age2 <- list(Mel = selectColonies(age2$Mel, p = 1 - p1collapse),
                   MelCross = selectColonies(age2$MelCross, p = 1 - p1collapse),
                   Car = selectColonies(age2$Car, p = 1 - p1collapse))
    }

    # Period2 ------------------------------------------------------------------
    print("PERIOD 2")
    # Swarm a percentage of age1 colonies
    # Mellifera
    print("Swarm colonies, P2")
    print(Sys.time())
    tmp <- list(Mel = pullColonies(age1$Mel, p = p2swarm),
                MelCross = pullColonies(age1$MelCross, p = p2swarm),
                Car = pullColonies(age1$Car, p = p2swarm))
    age1 <- list(Mel = tmp$Mel$remnant,
                 MelCross = tmp$MelCross$remnant,
                 Car = tmp$Car$remnant)
    tmp <- list(Mel = swarm(tmp$Mel$pulled),
                MelCross = swarm(tmp$MelCross$pulled),
                Car = swarm(tmp$Car$pulled))
    # The queens of the remnant colonies are of age 0
    age0p2 <- list(Mel = tmp$Mel$remnant,
                   MelCross = tmp$MelCross$remnant,
                   Car = tmp$Car$remnant)
    age1 <- list(Mel = c(age1$Mel, tmp$Mel$swarm),
                 MelCross = c(age1$MelCross, tmp$MelCross$swarm),
                 Car = c(age1$Car, tmp$Car$swarm))

    if (year > 1) {
      # Swarm a percentage of age2 colonies
      tmp <- list(Mel = pullColonies(age2$Mel, p = p2swarm),
                  MelCross = pullColonies(age2$MelCross, p = p2swarm),
                  Car = pullColonies(age2$Car, p = p2swarm))
      age2 <- list(Mel = tmp$Mel$remnant,
                   MelCross = tmp$MelCross$remnant,
                   Car = tmp$Car$remnant)
      tmp <- list(Mel = swarm(tmp$Mel$pulled),
                  MelCross = swarm(tmp$MelCross$pulled),
                  Car = swarm(tmp$Car$pulled))
      # The queens of the remnant colonies are of age 0
      age0p2 <- list(Mel = tmp$Mel$remnant,
                     MelCross = tmp$MelCross$remnant,
                     Car = tmp$Car$remnant)
      age2 <- list(Mel = c(age2$Mel, tmp$Mel$swarm),
                   MelCross = c(age2$MelCross, tmp$MelCross$swarm),
                   Car = c(age2$Car, tmp$Car$swarm))
    }

    # Supersede a part of age1 colonies
    print("Supersede colonies, P2")
    print(Sys.time())

    tmp <- list(Mel = pullColonies(age1$Mel, p = p2supersede),
                MelCross = pullColonies(age1$MelCross, p = p2supersede),
                Car = pullColonies(age1$Car, p = p2supersede))
    age1 <- list(Mel = tmp$Mel$remnant,
                 MelCross = tmp$MelCross$remnant,
                 Car = tmp$Car$remnant)
    tmp <- list(Mel = supersede(tmp$Mel$pulled),
                MelCross = supersede(tmp$MelCross$pulled),
                Car = supersede(tmp$Car$pulled))
    # The queens of superseded colonies are of age 0
    age0p2 <- list(Mel = c(age0p2$Mel, tmp$Mel),
                   MelCross = c(age0p2$MelCross, tmp$MelCross),
                   Car = c(age0p2$Car, tmp$Car))

    if (year > 1) {
      # Supersede a part of age2 colonies
      tmp <- list(Mel = pullColonies(age2$Mel, p = p2supersede),
                  MelCross = pullColonies(age2$MelCross, p = p2supersede),
                  Car = pullColonies(age2$Car, p = p2supersede))
      age2 <- list(Mel = tmp$Mel$remnant,
                   MelCross = tmp$MelCross$remnant,
                   Car = tmp$Car$remnant)
      tmp <- list(Mel = supersede(tmp$Mel$pulled),
                  MelCross = supersede(tmp$MelCross$pulled),
                  Car = supersede(tmp$Car$pulled))
      # The queens of superseded colonies are of age 0
      age0p2 <- list(Mel = c(age0p2$Mel, tmp$Mel),
                     MelCross = c(age0p2$MelCross, tmp$MelCross),
                     Car = c(age0p2$Car, tmp$Car))
    }

    # Replace all the drones
    print("Replace Drones, P2")
    print(Sys.time())

    age1$Mel <- replaceDrones(age1$Mel)
    age1$MelCross <- replaceDrones(age1$MelCross)
    age1$Car <- replaceDrones(age1$Car)
    if (year > 1) {
      age2$Mel <- replaceDrones(age2$Mel)
      age2$MelCross <- replaceDrones(age2$MelCross)
      age2$Car <- replaceDrones(age2$Car)
    }

    # Mate the colonies
    # Import p percentage of carnica colonies into mellifera DCA
    print("Mate colonies, P2")
    print(Sys.time())


    if (year == 1) {
      DCAMel <- createDCA(age1$Mel)
      age0p2$Mel <- cross(age0p2$Mel, drones = pullDroneGroupsFromDCA(DCA = DCAMel, n = nColonies(age0p2$Mel), nDrones = nFathersPoisson))
      DCAMelCross <- createDCA(c(age1$MelCross,
                                 selectColonies(age1$Car, n = round(nColonies(age1$MelCross) * pImport, 0))))
      age0p2$MelCross <- cross(age0p2$MelCross, drones = pullDroneGroupsFromDCA(DCA = DCAMelCross, n = nColonies(age0p2$MelCross), nDrones = nFathersPoisson))
      DCACar <- createDCA(age1$Car)
      age0p2$Car <- cross(age0p2$Car, drones = pullDroneGroupsFromDCA(DCA = DCACar, n = nColonies(age0p2$Car), nDrones = nFathersPoisson))
    } else {
      #Mel
      DCAMel <- createDCA(c(age1$Mel, age2$Mel))
      fathersMel <- pullDroneGroupsFromDCA(DCA = DCAMel, n = nColonies(age0p2$Mel), nDrones = nFathersPoisson)
      for (x in 1:length(fathersMel)){
        fathersMel[[x]] <- c(fathersMel[[x]], createDrones(age1$Mel[[1]], nInd = 2))
      }
      age0p2$Mel <- cross(age0p2$Mel, drones = fathersMel)
      
      #MelCross
      DCAMelCross <- createDCA(c(age1$MelCross,
                                 selectColonies(age1$Car, n = round(nColonies(age1$MelCross) * pImport, 0)),
                                 age2$MelCross,
                                 selectColonies(age2$Car, n = round(nColonies(age2$MelCross) * pImport, 0))))
      fathersMelCross <- pullDroneGroupsFromDCA(DCA = DCAMelCross, n = nColonies(age0p2$MelCross), nDrones = nFathersPoisson)
      for (x in 1:length(fathersMelCross)){
        fathersMelCross[[x]] <- c(fathersMelCross[[x]], createDrones(age1$MelCross[[1]], nInd = 2))
      }
      age0p2$MelCross <- cross(age0p2$MelCross, drones = fathersMelCross)
      
      #Car
      DCACar <- createDCA(c(age1$Car, age2$Car))
      fathersCar <-  pullDroneGroupsFromDCA(DCA = DCACar, n = nColonies(age0p2$Car), nDrones = nFathersPoisson)
      for (x in 1:length(fathersCar)){
        fathersCar[[x]] <- c(fathersCar[[x]], createDrones(age1$Car[[1]], nInd = 2))
      }
      age0p2$Car <- cross(age0p2$Car, drones = fathersCar)
    }

    # Collapse
    age1 <- list(Mel = selectColonies(age1$Mel, p = 1 - p2collapse),
                 MelCross = selectColonies(age1$MelCross, p = 1 - p2collapse),
                 Car = selectColonies(age1$Car, p = 1 - p2collapse))
    if (year > 1) {
      age2 <- list(Mel = selectColonies(age2$Mel, p = 1 - p2collapse),
                   MelCross = selectColonies(age2$MelCross, p = 1 - p2collapse),
                   Car = selectColonies(age2$Car, p = 1 - p2collapse))
    }

    # Merge all age 0 colonies (from both periods)
    age0 <- list(Mel = c(age0p1$Mel, age0p2$Mel),
                 MelCross = c(age0p1$MelCross, age0p2$MelCross),
                 Car = c(age0p1$Car, age0p2$Car))
    colonyRecords <- data_rec(datafile = colonyRecords, colonies = age0$Mel, year = year, population = "Mel")
    colonyRecords <- data_rec(datafile = colonyRecords, colonies = age0$MelCross, year = year, population = "MelCross")
    colonyRecords <- data_rec(datafile = colonyRecords, colonies = age0$Car, year = year, population = "Car")

    # Period3 ------------------------------------------------------------------
    # Collapse age0 queens
    print("PERIOD 3")
    print("Collapse colonies, P3")
    print(Sys.time())

    age0 <- list(Mel = selectColonies(age0$Mel, p = (1 - p3collapseAge0)),
                 MelCross = selectColonies(age0$MelCross, p = (1 - p3collapseAge0)),
                 Car = selectColonies(age0$Car, p = (1 - p3collapseAge0)))
    age1 <- list(Mel = selectColonies(age1$Mel, p = (1 - p3collapseAge1)),
                 MelCross = selectColonies(age1$MelCross, p = (1 - p3collapseAge1)),
                 Car = selectColonies(age1$Car, p = (1 - p3collapseAge1)))
    age2 <- list(Mel = NULL, MelCross=NULL, Car = NULL) #We don't need this but just to show the workflow!!!


    # Maintain the number of colonies ------------------------------------------
    # Keep all of age1, age0 swarmed so we build it up with some splits, while we remove (sell) the other splits
    print("Maintain the number, P2")
    print(Sys.time())

    age0$Mel <- maintainApiarySize(age0 = age0$Mel, age1 = age1$Mel)
    age0$MelCross <- maintainApiarySize(age0 = age0$MelCross, age1 = age1$MelCross)
    age0$Car <- maintainApiarySize(age0 = age0$Car, age1 = age1$Car)

    for (subspecies in c("Mel", "MelCross", "Car")) {
      if ((nColonies(age0[[subspecies]]) + nColonies(age1[[subspecies]])) != apiarySize) {
        stop(paste0("The number of colonies for ", subspecies, " does not match the apiary size!"))
      }
    }

  } # Year-loop

  a <- toc()
  loopTime <- rbind(loopTime, c(Rep, a$tic, a$toc, a$msg, (a$toc - a$tic)))

  print("Computing carnica year 10 relationships")
  print(Sys.time())
  pHomBroodYear10 <- sample(names(which(pHomBrood(age1$Car) > 0.2)), 1)
  springerColony10_Car <- computeRelationship_genomic(x = age1$Car[[pHomBroodYear10]],
                                                      csd = isCsdActive(SP),
                                                      ColonyAF = TRUE,
                                                      AF0.5 = TRUE,
                                                      SingleBaseAF = alleleFreqBaseQueensCar,
                                                      MultiBaseAF = alleleFreqBaseQueens,
                                                      SingleBaseAFcsdLocus = alleleFreqCsdLocusBaseCar,
                                                      MultiBaseAFcsdLocus = alleleFreqCsdLocusBaseQueens,
                                                      SingleBaseAFcsdChr = alleleFreqCsdChrBaseCar,
                                                      MultiBaseAFcsdChr = alleleFreqCsdChrBaseQueens)
  print("Computing queens relationships")
  print(Sys.time())
  springerQueens10_MelAF <- computeRelationship_genomic(x = c(mergePops(getQueen(age1$Mel)),
                                                        mergePops(getQueen(age1$MelCross)),
                                                        mergePops(getQueen(age1$Car))),
                                                  csd = isCsdActive(SP),
                                                  ColonyAF = TRUE,
                                                  AF0.5 = TRUE,
                                                  MultiBaseAF = alleleFreqBaseQueens,
                                                  SingleBaseAF = alleleFreqBaseQueensMel,
                                                  MultiBaseAFcsdLocus = alleleFreqCsdLocusBaseQueens,
                                                  SingleBaseAFcsdLocus = alleleFreqCsdLocusBaseMel,
                                                  MultiBaseAFcsdChr = alleleFreqCsdChrBaseQueens,
                                                  SingleBaseAFcsdChr = alleleFreqCsdChrBaseMel)

  springerQueens10_CarAF <- computeRelationship_genomic(x = c(mergePops(getQueen(age1$Mel)),
                                                        mergePops(getQueen(age1$MelCross)),
                                                        mergePops(getQueen(age1$Car))),
                                                  csd = isCsdActive(SP),
                                                  ColonyAF = TRUE,
                                                  AF0.5 = TRUE,
                                                  MultiBaseAF = alleleFreqBaseQueens,
                                                  SingleBaseAF = alleleFreqBaseQueensCar,
                                                  MultiBaseAFcsdLocus = alleleFreqCsdLocusBaseQueens,
                                                  SingleBaseAFcsdLocus = alleleFreqCsdLocusBaseCar,
                                                  MultiBaseAFcsdChr = alleleFreqCsdChrBaseQueens,
                                                  SingleBaseAFcsdChr = alleleFreqCsdChrBaseCar)

  springerQueensPop10 <- rbind(data.frame(ID = sapply(getQueen(age1$Mel), FUN = function(x) x@id), Pop = "Mel"),
                               data.frame(ID = sapply(getQueen(age1$MelCross), FUN = function(x) x@id), Pop = "MelCross"),
                               data.frame(ID = sapply(getQueen(age1$Car), FUN = function(x) x@id), Pop = "Car"))

  # Compute the pedigree relationship matrix
  print("Computing pedigree relationships")
  print(Sys.time())
  pedigree <- SP$pedigree
  caste <- SP$caste

  save(pedigree, caste,
       springerColony1_Car, springerColony10_Car,
       springerQueens1_MelAF, springerQueens10_CarAF, springerQueensPop1,
       springerQueens10_MelAF, springerQueens10_CarAF, springerQueensPop10,
       colonyRecords, file = "SpringerSimulation_import_objects.RData")

  IBDe <- computeRelationship_pedigree(SP$pedigree)
  save(IBDe, file = "IBDe_SpringerSimulation.RData")

} # Rep-loop

print("Saving image data")
save.image("SpringerSimulation_import.RData")








