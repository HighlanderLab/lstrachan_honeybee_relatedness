#setwd("~/Desktop/GitHub/lstrachan_honeybee_sim/YearCycleSimulation")
# Clean workspace
rm(list = ls())

# Define functions
computeRelationship_genomic <- function(x, alleleFreq = NULL, useOwnAlleleFreq = FALSE, alleleFreqCsd = NULL, csd = TRUE) {
  if (isColony(x)) {
    # Build the colony up to 1,000 workers and 200 drones
    colony <- buildUp(x = x,
                      nWorkers = 1000,
                      nDrones = 200)
    # Extract the genotypes of all the colony members
    geno <- rbind(getSegSiteGeno(colony, caste = "queen"),
                  getSegSiteGeno(colony, caste = "workers"),
                  getSegSiteGeno(colony, caste = "drones"),
                  getSegSiteGeno(colony, caste = "fathers"))
    # Collect the sex of all the colony members in the order than genotypes
    sex <- c(rep("F", nQueens(colony)),
             rep("F", nWorkers(colony)),
             rep("M", nDrones(colony)),
             rep("M", nFathers(colony)))
  } else if (SIMplyBee:::isPop(x)) {
    geno <- getSegSiteGeno(x)
    sex <- x@sex
  }
  # Compute the IBS relationship matrix
  # TODO: Combine generations 1 and 10 - so that they have the same reference population (regarding allele frequencies)
  print("IBS")
  print(Sys.time())
  if (is.null(alleleFreq) & isFALSE(useOwnAlleleFreq)) { #no allele freq info given
    ibs <- calcBeeGRMIbs(x = geno,
                         sex = sex,
                         alleleFreq = rep(0.5, ncol(geno))) #modified allele freq 
    ibsOwnFreq <- NULL
  } else if (isTRUE(useOwnAlleleFreq)) {
    ownAlleleFreq = calcBeeAlleleFreq(x = geno, sex = sex) #use function defaults (geno must be assigned to x) 
    if (is.null(alleleFreq)) {
      ibsOwnFreq <- calcBeeGRMIbs(x = geno,
                                  sex = sex,
                                  alleleFreq = ownAlleleFreq)
      ibs = NULL
    } else {
      ibs <- calcBeeGRMIbs(x = geno,
                           sex = sex,
                           alleleFreq = alleleFreq)
      ibsOwnFreq <- calcBeeGRMIbs(x = geno,
                                  sex = sex,
                                  alleleFreq = ownAlleleFreq)
    }
  }


  #Csd
  if (is.null(alleleFreqCsd)) { #DO THIS
    alleleFreqCsd <- rep(0.5, length(SP$csdPosStart:SP$csdPosStop))
  }
  # Collect the IBD haplotypes for all the colony members
  print("IBD")
  print(Sys.time())
  if (isColony(x)) {
    haplo <- rbind(getIbdHaplo(colony, caste = "queen"),
                   getIbdHaplo(colony, caste = "workers"),
                   getIbdHaplo(colony, caste = "drones"),
                   getIbdHaplo(colony, caste = "fathers"))
  } else if (SIMplyBee:::isPop(x)) {
    haplo <- getIbdHaplo(x)
  }
  # Compute the IBD relationship matrix
  ibd <- calcBeeGRMIbd(x = haplo)
  ibd <- ibd$indiv

  print("csd IBD")
  print(Sys.time())
  if (csd) {
    # Only chromosome 3
    csdChr <- SP$csdChr
    csdChrHaplo <- haplo[, grepl(pattern = paste0(csdChr, "_"), x = colnames(haplo))]
    ibd_csdChr <- calcBeeGRMIbd(x = csdChrHaplo)
    ibd_csdChr <- ibd_csdChr$indiv
    csdChrGeno <- geno[, grepl(pattern = paste0(csdChr, "_"), x = colnames(geno))]
    ibs_csdChr <- calcBeeGRMIbs(x =  csdChrGeno,
                                sex = sex,
                                alleleFreq = NULL) #THIS SHOULD BE NULL! 

    # Only csd locus
    ibd_csd <- calcBeeGRMIbd(x = haplo[, paste(SP$csdChr,
                                               SP$csdPosStart:SP$csdPosStop,
                                               sep = "_")])
    ibd_csd <- ibd_csd$indiv
    ibs_csd <- calcBeeGRMIbs(x = geno[, paste(SP$csdChr,
                                              SP$csdPosStart:SP$csdPosStop,
                                              sep = "_")],
                             sex = sex,
                             alleleFreq = alleleFreqCsd)
  } else {
    ibd_csdChr = ibd_csdChr = ibd_csd = ibs_csd = NULL
  }

  if (isColony(x)) {
    id <- getCasteId(colony, caste = "all")
  } else if (SIMplyBee:::isPop(x)) {
    id <- x@id
  }



  return(list(IBS = ibs, IBSOwnFreq = ibsOwnFreq, IBD = ibd,
              IBScsdChr = ibs_csdChr, IBDcsdChr = ibd_csdChr,
              IBSCsd = ibs_csd, IBDCsd = ibd_csd, ID = id))
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


# TODO: replace with devtools installation from Github once the package is operational
# Source the development version of AlphaSimR


library(SIMplyBee)

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
                   )
  )
}
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
  queens <- list(Mel = cross(x = virginQueens$Mel[1:apiarySize], drones = fathersMel),
                 MelCross = cross(x = virginQueens$MelCross[1:apiarySize], drones = fathersMelCross),
                 Car = cross(x = virginQueens$Car[1:apiarySize], drones = fathersCar))

  #Set allele frequency for queens
  tmp <- c(virginQueens$Mel, virginQueens$Car)
  alleleFreqBaseQueens <- calcBeeAlleleFreq(x = getSegSiteGeno(tmp),
                                            sex = tmp@sex)
  #Get allele freq for csd chromomsome
  csdLoci <- paste0(SP$csdChr, "_", SP$csdPosStart:SP$csdPosStop)
  alleleFreqCsdBaseQueens <- alleleFreqBaseQueens[csdLoci]

  # Start the year-loop ------------------------------------------------------------------
  for (year in 1:nYear) {
    print("Starting the cycle")
    year <- 1
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
       print("Computing initial relationships")
       # Choose the first colony of age 1 to inspect relationship in the base population
       # print("Computing mellifera generation1 colony relationship")
       # print(Sys.time())
      #  springerColony1_Mel <- computeRelationship_genomic(x = age1$Mel[[1]],
      #                                                     csd = isCsdActive(SP),
      #                                                     alleleFreq = alleleFreqBaseQueens,
      #                                                     alleleFreqCsd = alleleFreqCsdBaseQueens,
      #                                                     useOwnAlleleFreq = TRUE)
      # print("Computing mellifera cross generation1 colony relationship")
      #  print(Sys.time())
      #  springerColony1_MelCross <- computeRelationship_genomic(x = age1$MelCross[[1]],
      #                                                     csd = isCsdActive(SP),
      #                                                     alleleFreq = alleleFreqBaseQueens,
      #                                                     alleleFreqCsd = alleleFreqCsdBaseQueens,
      #                                                     useOwnAlleleFreq = TRUE)
      print("Computing carnica generation1 colony relationship")
      print(Sys.time())
      springerColony1_Car <- computeRelationship_genomic(x = age1$Car[[1]],
                                                         csd = isCsdActive(SP),
                                                         alleleFreq = alleleFreqBaseQueens,
                                                         alleleFreqCsd = alleleFreqCsdBaseQueens,
                                                         useOwnAlleleFreq = TRUE)
      print("Computing queens generation1 colony relationship")
      print(Sys.time())
      springerQueens1 <- computeRelationship_genomic(x = c(queens$Mel, queens$MelCross, queens$Car),
                                                     csd = isCsdActive(SP),
                                                     alleleFreq = alleleFreqBaseQueens,
                                                     alleleFreqCsd = alleleFreqCsdBaseQueens,
                                                     useOwnAlleleFreq = TRUE)
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
      DCAMel <- createDCA(c(age1$Mel, age2$Mel))
      age0p1$Mel <- cross(age0p1$Mel, drones = pullDroneGroupsFromDCA(DCA = DCAMel, n = nColonies(age0p1$Mel), nDrones = nFathersPoisson))
      DCAMelCross <- createDCA(c(age1$MelCross,
                                 selectColonies(age1$Car, n = round(nColonies(age1$MelCross) * pImport, 0)),
                                 age2$MelCross,
                                 selectColonies(age2$Car, n = round(nColonies(age2$MelCross) * pImport, 0))))
      age0p1$MelCross <- cross(age0p1$MelCross, drones = pullDroneGroupsFromDCA(DCA = DCAMelCross, n = nColonies(age0p1$MelCross), nDrones = nFathersPoisson))
      DCACar <- createDCA(c(age1$Car, age2$Car))
      age0p1$Car <- cross(age0p1$Car, drones = pullDroneGroupsFromDCA(DCA = DCACar, n = nColonies(age0p1$Car), nDrones = nFathersPoisson))
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
      DCAMel <- createDCA(c(age1$Mel, age2$Mel))
      age0p2$Mel <- cross(age0p2$Mel, drones = pullDroneGroupsFromDCA(DCA = DCAMel, n = nColonies(age0p2$Mel), nDrones = nFathersPoisson))
      DCAMelCross <- createDCA(c(age1$MelCross,
                                 selectColonies(age1$Car, n = round(nColonies(age1$MelCross) * pImport, 0)),
                                 age2$MelCross,
                                 selectColonies(age2$Car, n = round(nColonies(age2$MelCross) * pImport, 0))))
      age0p2$MelCross <- cross(age0p2$MelCross, drones = pullDroneGroupsFromDCA(DCA = DCAMelCross, n = nColonies(age0p2$MelCross), nDrones = nFathersPoisson))
      DCACar <- createDCA(c(age1$Car, age2$Car))
      age0p2$Car <- cross(age0p2$Car, drones = pullDroneGroupsFromDCA(DCA = DCACar, n = nColonies(age0p2$Car), nDrones = nFathersPoisson))
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

  print("Computing final relationships")
  print(Sys.time())

  # Take the first colony of age1 and extract genotypes for relationship computation
  # print("Computing mellifera relationships")
  # print(Sys.time())
  # springerColony10_Mel <- computeRelationship_genomic(x = age1$Mel[[1]],
  #                                                     csd = isCsdActive(SP),
  #                                                     alleleFreq = alleleFreqBaseQueens,
  #                                                     alleleFreqCsd = alleleFreqCsdBaseQueens,
  #                                                     useOwnAlleleFreq = TRUE)
  # print("Computing mellifera cross relationships")
  # print(Sys.time())
  # springerColony10_MelCross <- computeRelationship_genomic(x = age1$MelCross[[1]],
  #                                                          csd = isCsdActive(SP),
  #                                                          alleleFreq = alleleFreqBaseQueens,
  #                                                          alleleFreqCsd = alleleFreqCsdBaseQueens,
  #                                                          useOwnAlleleFreq = TRUE)
  print("Computing carnica relationships")
  print(Sys.time())
  springerColony10_Car <- computeRelationship_genomic(x = age1$Car[[1]],
                                                      csd = isCsdActive(SP),
                                                      alleleFreq = alleleFreqBaseQueens,
                                                      alleleFreqCsd = alleleFreqCsdBaseQueens,
                                                      useOwnAlleleFreq = TRUE)
  print("Computing queens relationships")
  print(Sys.time())
  springerQueens10 <- computeRelationship_genomic(x = c(mergePops(getQueen(age1$Mel)),
                                                        mergePops(getQueen(age1$MelCross)),
                                                        mergePops(getQueen(age1$Car))),
                                                  csd = isCsdActive(SP),
                                                  alleleFreq = alleleFreqBaseQueens,
                                                  alleleFreqCsd = alleleFreqCsdBaseQueens,
                                                  useOwnAlleleFreq = TRUE)
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
       springerQueens1, springerQueensPop1,
       springerQueens10, springerQueensPop10,
       colonyRecords, file = "SpringerSimulation_import_objects.RData")

  IBDe <- computeRelationship_pedigree(SP$pedigree)
  save(IBDe, file = "IBDe_SpringerSimulation.RData")

} # Rep-loop

print("Saving image data")
save.image("SpringerSimulation_import.RData")








