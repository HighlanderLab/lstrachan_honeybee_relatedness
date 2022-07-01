#setwd("~/Desktop/GitHub/lstrachan_honeybee_sim/YearCycleSimulation")
# Clean workspace
rm(list = ls())

# Define functions
computeRelationship_genomic <- function(x, alleleFreq = NULL, alleleFreqCsd = NULL, csd = TRUE) {
  if (isColony(x)) {
    # Build the colony up to 1,000 workers and 200 drones
    colony <- buildUpColony(colony = x,
                            nWorkers = 1000,
                            nDrones = 200)
    # Extract the genotypes of all the colony members
    geno <- rbind(getCasteSegSiteGeno(colony, caste = "queen"),
                  getCasteSegSiteGeno(colony, caste = "workers"),
                  getCasteSegSiteGeno(colony, caste = "drones"),
                  getCasteSegSiteGeno(colony, caste = "fathers"))
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
  if (is.null(alleleFreq)) {
    alleleFreq <- rep(0.5, ncol(geno))
  }
  if (is.null(alleleFreqCsd)) {
    alleleFreqCsd <- rep(0.5, length(SP$csdPosStart:SP$csdPosStop))
  }
  ibs <- calcBeeGRMIbs(x = geno,
                       sex = sex,
                       alleleFreq = alleleFreq)
  # Collect the IBD haplotypes for all the colony members
  print("IBD")
  print(Sys.time())
  if (isColony(x)) {
    haplo <- rbind(getCasteIbdHaplo(colony, caste = "queen"),
                   getCasteIbdHaplo(colony, caste = "workers"),
                   getCasteIbdHaplo(colony, caste = "drones"),
                   getCasteIbdHaplo(colony, caste = "fathers"))
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
                                alleleFreq = alleleFreq)

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
  
  

  return(list(IBS = ibs, IBD = ibd,
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

getCsdInfo <- function (colonies, subspecies = NULL) {
  totalCsd <- nCsdAlleles(colonies, collapse = TRUE)
  csdVariability <- data.frame()
  pDiploidDrones <- data.frame()
  for (n in 1:nColonies(colonies)){
    csdVariability <- rbind(csdVariability,
                            c(Rep = Rep, year = year, id = colonies[[n]]@id,
                              nCSD = nCsdAlleles(colonies[[n]], collapse = TRUE),
                              totalCSD = totalCsd, 
                              subspecies = subspecies))
    pDiploidDrones <- rbind(pDiploidDrones,
                            c(Rep = Rep, year = year, id = colonies[[n]]@id,
                              pQueenHomBrood = calcQueensPHomBrood(colonies[[n]]), subspecies = subspecies))
  }
  colnames(csdVariability) <- c("Rep", "year", "id", "nCSD", "totalCSD", "subspecies")
  colnames(pDiploidDrones) <- c("Rep", "year", "id", "pQueenHomBrood", "subspecies")
  return(list(csdVariability = csdVariability, pDiploidDrones = pDiploidDrones))
}


maintainApiarySize <- function(age0 = NULL, age1 = NULL) {
  if ((nColonies(age0) + nColonies(age1)) > apiarySize) { # check if the sum of all colonies is greater than apiary size
    IDsplits <- getId(age0)[hasSplit(age0)] # get the IDs of age 0 that are splits
    splits0 <- pullColonies(age0, ID = IDsplits) # pull the splits out of age 0
    age0split <- splits0$pulledColonies # create an object for age 0 splits
    age0swarm <- splits0$remainingColonies # create an object for swarms and superseded colonies
    age0needed <- apiarySize - nColonies(age1) # calculate the number of age 0 colonies that are needed to fill up the apiary
    splitsNeeded <- age0needed - nColonies(age0swarm) # calculate the number of splits needed
    if (age0needed <= nColonies(age0swarm)) { # check if the number of age 0 colonies needed is lower or equal to age 0 swarms
      swarmID <- sample(getId(age0swarm), age0needed) # if yes, select the ids of swarms that will stay in apiary
      swarmTMP <- pullColonies(age0swarm, ID = swarmID) # pull out those selected age0 swarms
      age0 <- swarmTMP$pulledColonies # put selected swarms to age 0 object
    } else if (age0needed > nColonies(age0swarm)) { # in case when age 0 needed is grater than number of swarm select splits
      nSplitsNeeded <- age0needed - nColonies(age0swarm) # calculate the number of splits needed
      splitId <- sample(getId(age0split), nSplitsNeeded) # select ids of splits
      splitTmp <- pullColonies(age0split, ID = splitId) # pull the splits
      splits <- splitTmp$pulledColonies # select pulled splits
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
nChr = 1
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
nWorkers <- 2 # TODO: change to 20K
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
csdVariability <- data.frame(Rep = NA, year = NA, id = NA, nCSD = NA, totalCSD = NA, subspecies = NA)
pDiploidDrones <- data.frame(Rep = NA, year = NA, id = NA, pQueenHomBrood = NA, subspecies = NA)

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
  #save(founderGenomes, file="FounderGenomes_ThreePop.RData")
  print("Loading in the founderData")
  load("FounderGenomes_ThreePop.RData")
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

  # Create a base population for A. m. mellifea, A. m. mellifera cross, and A. m. carnica
  virginQueens <- list(Mel = createVirginQueens(x = founderGenomes[1:(nMelN/2)]),
                       MelCross = createVirginQueens(x = founderGenomes[(nMelN/2 + 1):nMelN]),
                       Car = createVirginQueens(x = founderGenomes[(nMelN +1):(nMelN + nCar)]))
  # Create drones for A. m. mellifea, A. m. mellifera cross, and A. m. carnica
  drones <- list(Mel = createDrones(x = virginQueens$Mel[(apiarySize+1):(nMelN/2)], nInd = nDronesPerQueen),
                 MelCross = createDrones(x = virginQueens$MelCross[(apiarySize+1):(nMelN/2)], nInd = nDronesPerQueen),
                 Car = createDrones(x = virginQueens$Car[(apiarySize+1):nCar], nInd = nDronesPerQueen))
  # Mate A. m. mellifera queens with A. m. mellifera drones
  queens <- list(Mel = crossVirginQueen(pop = virginQueens$Mel[1:apiarySize], drones = drones$Mel),
                 MelCross = crossVirginQueen(pop = virginQueens$MelCross[1:apiarySize], drones = drones$MelCross),
                 Car = crossVirginQueen(pop = virginQueens$Car[1:apiarySize], drones = drones$Car))

  tmp <- c(queens$Mel, queens$Car)
  alleleFreqBaseQueens <- calcBeeAlleleFreq(x = getSegSiteGeno(tmp),
                                            sex = tmp@sex)
  #TODO: get allele freq for csd chromomsome 
  csdLoci <- paste0(SP$csdChr, "_", SP$csdPosStart:SP$csdPosStop)
  alleleFreqCsdBaseQueens <- alleleFreqBaseQueens[csdLoci]
  
  # Start the year-loop ------------------------------------------------------------------
  for (year in 1:nYear) {
    print("Starting the cycle")
    # year <- 1
    # year <- year + 1
    cat(paste0("Year: ", year, "/", nYear, "\n"))
    # If this is the first year, create some colonies to start with
    # If not, promote the age0 to age1, age1 to age2 and remove age2 colonies
    if (year == 1) {
      print("Creating initial colonies")
      age1 <- list(Mel = createColonies(x = queens$Mel, n = apiarySize),
                   MelCross = createColonies(x = queens$MelCross, n = apiarySize),
                   Car = createColonies(x = queens$Car, n = apiarySize))
    } else {
      age2 <- list(Mel = age1$Mel, MelCross = age1$MelCross, Car = age1$Car)
      age1 <- list(Mel = age0$Mel, MelCross = age0$MelCross, Car = age0$Car)
      age0 <- list(Mel = NULL, MelCross = NULL, Car = NULL)
      age0p1 <- list(Mel = NULL, MelCross = NULL, Car = NULL)
      age0p2 <- list(Mel = NULL, MelCross = NULL, Car = NULL)
    }

    # In year 1, inspect the relationship in one of the colonies
    if (year == 1) {
      print("Computing initial relationships")
      # Choose the first colony of age 1 to inspect relationship in the base population
      print("Computing mellifera generation1 colony relationship")
      print(Sys.time())
      springerColony1_Mel <- computeRelationship_genomic(x = age1$Mel[[1]],
                                                         csd = isCsdActive(SP),
                                                         alleleFreq = alleleFreqBaseQueens,
                                                         alleleFreqCsd = alleleFreqCsdBaseQueens)
      print("Computing mellifera cross generation1 colony relationship")
      print(Sys.time())
      springerColony1_MelCross <- computeRelationship_genomic(x = age1$MelCross[[1]],
                                                              csd = isCsdActive(SP),
                                                              alleleFreq = alleleFreqBaseQueens,
                                                              alleleFreqCsd = alleleFreqCsdBaseQueens)
      print("Computing carnica generation1 colony relationship")
      print(Sys.time())
      springerColony1_Car <- computeRelationship_genomic(x = age1$Car[[1]],
                                                         csd = isCsdActive(SP),
                                                         alleleFreq = alleleFreqBaseQueens,
                                                         alleleFreqCsd = alleleFreqCsdBaseQueens)
      print("Computing queens generation1 colony relationship")
      print(Sys.time())
      springerQueens1 <- computeRelationship_genomic(x = c(queens$Mel, queens$MelCross, queens$Car),
                                                     csd = isCsdActive(SP),
                                                     alleleFreq = alleleFreqBaseQueens,
                                                     alleleFreqCsd = alleleFreqCsdBaseQueens)
      springerQueensPop1 <- rbind(data.frame(ID = sapply(getQueen(age1$Mel), FUN = function(x) x@id), Pop = "Mel"),
                                  data.frame(ID = sapply(getQueen(age1$MelCross), FUN = function(x) x@id), Pop = "MelCross"),
                                  data.frame(ID = sapply(getQueen(age1$Car), FUN = function(x) x@id), Pop = "Car"))

    }
    print("Done computing initial relationships")

    # Period1 ------------------------------------------------------------------
    # Build-up the colonies
    print(paste0("Building up the colonies to ", nWorkers, " and ", nDrones))
    print(Sys.time())
    age1 <- list(Mel = buildUpColonies(age1$Mel),
                 MelCross = buildUpColonies(age1$MelCross),
                 Car = buildUpColonies(age1$Car))
    if (year > 1) {
      age2 <- list(Mel = buildUpColonies(age2$Mel),
                   MelCross = buildUpColonies(age2$MelCross),
                   Car = buildUpColonies(age2$Car))
    }

    # Split all age1 colonies
    print("Splitting the colonies")
    print(Sys.time())
    tmp <- list(Mel = splitColonies(age1$Mel),
                MelCross = splitColonies(age1$MelCross),
                Car = splitColonies(age1$Car))
    age1 <- list(Mel = tmp$Mel$remnants,
                 MelCross = tmp$MelCross$remnants,
                 Car = tmp$Car$remnants)
    # The queens of the splits are 0 years old
    age0p1 <- list(Mel = tmp$Mel$splits, MelCross = tmp$MelCross$splits, Car = tmp$Car$splits)

    if (year > 1) {
      # Split all age2 colonies
      tmp <- list(Mel = splitColonies(age2$Mel),
                  MelCross = splitColonies(age2$MelCross),
                  Car = splitColonies(age2$Car))
      age2 <- list(Mel = tmp$Mel$remnants,
                   MelCross = tmp$MelCross$remnants,
                   Car = tmp$Car$remnants)
      # The queens of the splits are 0 years old
      age0p1 <- list(Mel = c(age0p1$Mel, tmp$Mel$splits),
                     MelCross = c(age0p1$MelCross, tmp$MelCross$splits),
                     Car = c(age0p1$Car, tmp$Car$splits))
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
    age0p1 <- list(Mel = reQueenColonies(age0p1$Mel, queens = virginQueens$Mel),
                   MelCross = reQueenColonies(age0p1$MelCross, queens = virginQueens$MelCross),
                   Car = reQueenColonies(age0p1$Car, queens = virginQueens$Car))

    # Swarm a percentage of age1 colonies
    print("Swarm colonies, P1")
    print(Sys.time())
    tmp <- list(Mel = pullColonies(age1$Mel, p = p1swarm),
                MelCross = pullColonies(age1$MelCross, p = p1swarm),
                Car = pullColonies(age1$Car, p = p1swarm))
    age1 <- list(Mel = tmp$Mel$remainingColonies,
                 MelCross = tmp$MelCross$remainingColonies,
                 Car = tmp$Car$remainingColonies)
    tmp <- list(Mel = swarmColonies(tmp$Mel$pulledColonies),
                MelCross = swarmColonies(tmp$MelCross$pulledColonies),
                Car = swarmColonies(tmp$Car$pulledColonies))
    age0p1 <- list(Mel = c(age0p1$Mel, tmp$Mel$remnants),
                   MelCross = c(age0p1$MelCross, tmp$MelCross$remnants),
                   Car = c(age0p1$Car, tmp$Car$remnants))
    age1 <- list(Mel = c(age1$Mel, tmp$Mel$swarms),
                 MelCross = c(age1$MelCross, tmp$MelCross$swarms),
                 Car = c(age1$Car, tmp$Car$swarms))


    if (year > 1) {
      # Swarm a percentage of age2 colonies
      tmp <- list(Mel = pullColonies(age2$Mel, p = p1swarm),
                  MelCross = pullColonies(age2$MelCross, p = p1swarm),
                  Car = pullColonies(age2$Car, p = p1swarm))
      age2 <- list(Mel = tmp$Mel$remainingColonies,
                   MelCross = tmp$MelCross$remainingColonies,
                   Car = tmp$Car$remainingColonies)
      tmp <- list(Mel = swarmColonies(tmp$Mel$pulledColonies),
                  MelCross = swarmColonies(tmp$MelCross$pulledColonies),
                  Car = swarmColonies(tmp$Car$pulledColonies))
      age0p1 <- list(Mel = c(age0p1$Mel, tmp$Mel$remnants),
                     MelCross = c(age0p1$MelCross, tmp$MelCross$remnants),
                     Car = c(age0p1$Car, tmp$Car$remnants))
      age2 <- list(Mel = c(age2$Mel, tmp$Mel$swarms),
                   MelCross = c(age2$MelCross, tmp$MelCross$swarms),
                   Car = c(age2$Car, tmp$Car$swarms))
    }

    # Supersede age1 colonies
    print("Supersede colonies, P1")
    print(Sys.time())
    tmp <- list(Mel = pullColonies(age1$Mel, p = p1supersede),
                MelCross = pullColonies(age1$MelCross, p = p1supersede),
                Car = pullColonies(age1$Car, p = p1supersede))
    age1 <- list(Mel = tmp$Mel$remainingColonies,
                 MelCross = tmp$MelCross$remainingColonies,
                 Car = tmp$Car$remainingColonies)
    tmp <- list(Mel = supersedeColonies(tmp$Mel$pulledColonies),
                MelCross = supersedeColonies(tmp$MelCross$pulledColonies),
                Car = supersedeColonies(tmp$Car$pulledColonies))
    age0p1 <- list(Mel = c(age0p1$Mel, tmp$Mel),
                   MelCross = c(age0p1$MelCross, tmp$MelCross),
                   Car = c(age0p1$Car, tmp$Car))


    if (year > 1) {
      # Supersede age2 colonies
      tmp <- list(Mel = pullColonies(age2$Mel, p = p1supersede),
                  MelCross = pullColonies(age2$MelCross, p = p1supersede),
                  Car = pullColonies(age2$Car, p = p1supersede))
      age2 <- list(Mel = tmp$Mel$remainingColonies,
                   MelCross = tmp$MelCross$remainingColonies,
                   Car = tmp$Car$remainingColonies)
      tmp <- list(Mel = supersedeColonies(tmp$Mel$pulledColonies),
                  MelCross = supersedeColonies(tmp$MelCross$pulledColonies),
                  Car = supersedeColonies(tmp$Car$pulledColonies))
      age0p1 <- list(Mel = c(age0p1$Mel, tmp$Mel),
                     MelCross = c(age0p1$MelCross, tmp$MelCross),
                     Car = c(age0p1$Car, tmp$Car))
    }

    # Mate the split colonies
    print("Mate split colonies, P1")
    print(Sys.time())
    if (year == 1) {
      DCAMel <- createDCA(age1$Mel)
      age0p1$Mel <- crossColonies(age0p1$Mel, drones = DCAMel)
      DCAMelCross <- createDCA(c(age1$MelCross,
                                 selectColonies(age1$Car, n = round(nColonies(age1$MelCross) * pImport, 0))))
      age0p1$MelCross <- crossColonies(age0p1$MelCross, drones = DCAMelCross)
      DCACar <- createDCA(age1$Car)
      age0p1$Car <- crossColonies(age0p1$Car, drones = DCACar)
    } else {
      DCAMel <- createDCA(c(age1$Mel, age2$Mel))
      age0p1$Mel <- crossColonies(age0p1$Mel, drones = DCAMel)
      DCAMelCross <- createDCA(c(age1$MelCross,
                                 selectColonies(age1$Car, n = round(nColonies(age1$MelCross) * pImport, 0)),
                                 age2$MelCross,
                                 selectColonies(age2$Car, n = round(nColonies(age2$MelCross) * pImport, 0))))
      age0p1$MelCross <- crossColonies(age0p1$MelCross, drones = DCAMelCross)
      DCACar <- createDCA(c(age1$Car, age2$Car))
      age0p1$Car <- crossColonies(age0p1$Car, drones = DCACar)
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
    age1 <- list(Mel = tmp$Mel$remainingColonies,
                 MelCross = tmp$MelCross$remainingColonies,
                 Car = tmp$Car$remainingColonies)
    tmp <- list(Mel = swarmColonies(tmp$Mel$pulledColonies),
                MelCross = swarmColonies(tmp$MelCross$pulledColonies),
                Car = swarmColonies(tmp$Car$pulledColonies))
    # The queens of the remnant colonies are of age 0
    age0p2 <- list(Mel = tmp$Mel$remnants,
                   MelCross = tmp$MelCross$remnants,
                   Car = tmp$Car$remnants)
    age1 <- list(Mel = c(age1$Mel, tmp$Mel$swarms),
                 MelCross = c(age1$MelCross, tmp$MelCross$swarms),
                 Car = c(age1$Car, tmp$Car$swarms))

    if (year > 1) {
      # Swarm a percentage of age2 colonies
      tmp <- list(Mel = pullColonies(age2$Mel, p = p2swarm),
                  MelCross = pullColonies(age2$MelCross, p = p2swarm),
                  Car = pullColonies(age2$Car, p = p2swarm))
      age2 <- list(Mel = tmp$Mel$remainingColonies,
                   MelCross = tmp$MelCross$remainingColonies,
                   Car = tmp$Car$remainingColonies)
      tmp <- list(Mel = swarmColonies(tmp$Mel$pulledColonies),
                  MelCross = swarmColonies(tmp$MelCross$pulledColonies),
                  Car = swarmColonies(tmp$Car$pulledColonies))
      # The queens of the remnant colonies are of age 0
      age0p2 <- list(Mel = tmp$Mel$remnants,
                     MelCross = tmp$MelCross$remnants,
                     Car = tmp$Car$remnants)
      age2 <- list(Mel = c(age2$Mel, tmp$Mel$swarms),
                   MelCross = c(age2$MelCross, tmp$MelCross$swarms),
                   Car = c(age2$Car, tmp$Car$swarms))
    }

    # Supersede a part of age1 colonies
    print("Supersede colonies, P2")
    print(Sys.time())

    tmp <- list(Mel = pullColonies(age1$Mel, p = p2supersede),
                MelCross = pullColonies(age1$MelCross, p = p2supersede),
                Car = pullColonies(age1$Car, p = p2supersede))
    age1 <- list(Mel = tmp$Mel$remainingColonies,
                 MelCross = tmp$MelCross$remainingColonies,
                 Car = tmp$Car$remainingColonies)
    tmp <- list(Mel = supersedeColonies(tmp$Mel$pulledColonies),
                MelCross = supersedeColonies(tmp$MelCross$pulledColonies),
                Car = supersedeColonies(tmp$Car$pulledColonies))
    # The queens of superseded colonies are of age 0
    age0p2 <- list(Mel = c(age0p2$Mel, tmp$Mel),
                   MelCross = c(age0p2$MelCross, tmp$MelCross),
                   Car = c(age0p2$Car, tmp$Car))

    if (year > 1) {
      # Supersede a part of age2 colonies
      tmp <- list(Mel = pullColonies(age2$Mel, p = p2supersede),
                  MelCross = pullColonies(age2$MelCross, p = p2supersede),
                  Car = pullColonies(age2$Car, p = p2supersede))
      age2 <- list(Mel = tmp$Mel$remainingColonies,
                   MelCross = tmp$MelCross$remainingColonies,
                   Car = tmp$Car$remainingColonies)
      tmp <- list(Mel = supersedeColonies(tmp$Mel$pulledColonies),
                  MelCross = supersedeColonies(tmp$MelCross$pulledColonies),
                  Car = supersedeColonies(tmp$Car$pulledColonies))
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
      age0p2$Mel <- crossColonies(age0p2$Mel, drones = DCAMel)
      DCAMelCross <- createDCA(c(age1$MelCross,
                                 selectColonies(age1$Car, n = round(nColonies(age1$MelCross) * pImport, 0))))
      age0p2$MelCross <- crossColonies(age0p2$MelCross, drones = DCAMelCross)
      DCACar <- createDCA(age1$Car)
      age0p2$Car <- crossColonies(age0p2$Car, drones = DCACar)
    } else {
      DCAMel <- createDCA(c(age1$Mel, age2$Mel))
      age0p2$Mel <- crossColonies(age0p2$Mel, drones = DCAMel)
      DCAMelCross <- createDCA(c(age1$MelCross,
                                 selectColonies(age1$Car, n = round(nColonies(age1$MelCross) * pImport, 0)),
                                 age2$MelCross,
                                 selectColonies(age2$Car, n = round(nColonies(age2$MelCross) * pImport, 0))))
      age0p2$MelCross <- crossColonies(age0p2$MelCross, drones = DCAMelCross)
      DCACar <- createDCA(c(age1$Car, age2$Car))
      age0p2$Car <- crossColonies(age0p2$Car, drones = DCACar)
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

    #Collect CSD info
    print("Collect csd info")
    print(Sys.info())
    csdInfoMel <- getCsdInfo(age0$Mel, subspecies = "Mel")
    csdInfoMelCross <- getCsdInfo(age0$MelCross, subspecies = "MelCross")
    csdInfoCar <- getCsdInfo(age0$Car, subspecies = "Car")
    csdVariability <- rbind(csdVariability, csdInfoMel$csdVariability)
    csdVariability <- rbind(csdVariability, csdInfoMelCross$csdVariability)
    csdVariability <- rbind(csdVariability, csdInfoCar$csdVariability)
    pDiploidDrones <-rbind(pDiploidDrones, csdInfoMel$pDiploidDrones)
    pDiploidDrones <-rbind(pDiploidDrones, csdInfoMelCross$pDiploidDrones)
    pDiploidDrones <-rbind(pDiploidDrones, csdInfoCar$DiploidDrones)


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
  print("Computing mellifera relationships")
  print(Sys.time())
  springerColony10_Mel <- computeRelationship_genomic(x = age1$Mel[[1]],
                                                      csd = isCsdActive(SP),
                                                      alleleFreq = alleleFreqBaseQueens,
                                                      alleleFreqCsd = alleleFreqCsdBaseQueens)
  print("Computing mellifera cross relationships")
  print(Sys.time())
  springerColony10_MelCross <- computeRelationship_genomic(x = age1$MelCross[[1]],
                                                           csd = isCsdActive(SP),
                                                           alleleFreq = alleleFreqBaseQueens,
                                                           alleleFreqCsd = alleleFreqCsdBaseQueens)
  print("Computing carnica relationships")
  print(Sys.time())
  springerColony10_Car <- computeRelationship_genomic(x = age1$Car[[1]],
                                                      csd = isCsdActive(SP),
                                                      alleleFreq = alleleFreqBaseQueens,
                                                      alleleFreqCsd = alleleFreqCsdBaseQueens)
  print("Computing queens relationships")
  print(Sys.time())
  springerQueens10 <- computeRelationship_genomic(x = c(mergePops(getQueen(age1$Mel)),
                                                        mergePops(getQueen(age1$MelCross)),
                                                        mergePops(getQueen(age1$Car))),
                                                  csd = isCsdActive(SP),
                                                  alleleFreq = alleleFreqBaseQueens,
                                                  alleleFreqCsd = alleleFreqCsdBaseQueens)
  springerQueensPop10 <- rbind(data.frame(ID = sapply(getQueen(age1$Mel), FUN = function(x) x@id), Pop = "Mel"),
                               data.frame(ID = sapply(getQueen(age1$MelCross), FUN = function(x) x@id), Pop = "MelCross"),
                               data.frame(ID = sapply(getQueen(age1$Car), FUN = function(x) x@id), Pop = "Car"))

  # Compute the pedigree relationship matrix
  print("Computing pedigree relationships")
  print(Sys.time())
  pedigree <- SP$pedigree
  caste <- SP$caste

  save(pedigree, caste,
       springerColony1_Mel, springerColony10_Mel,
       springerColony1_MelCross, springerColony10_MelCross,
       springerColony1_Car, springerColony10_Car,
       springerQueens1, springerQueensPop1,
       springerQueens10, springerQueensPop10,
       csdVariability, pDiploidDrones, file = "SpringerSimulation_import_objects.RData")

  IBDe <- computeRelationship_pedigree(SP$pedigree)
  save(IBDe, file = "IBDe_SpringerSimulation.RData")

} # Rep-loop

print("Saving image data")
save.image("SpringerSimulation_import.RData")








