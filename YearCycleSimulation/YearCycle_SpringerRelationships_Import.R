setwd("~/github/lstrachan_honeybee_sim/YearCycleSimulation/")
# Clean workspace
rm(list = ls())


# Define functions
computeRelationship_genomic <- function(colony) {
  # Build the colony up to 10,000 workers and 200 drones
  colony <- buildUpColony(colony = colony,
                                   nWorkers = 1000,
                                   nDrones = 200)
  # Extract the genotypes of all the colony members
  colony_geno <- rbind(getCasteSegSiteGeno(colony, caste = "queen"),
                       getCasteSegSiteGeno(colony, caste = "workers"),
                       getCasteSegSiteGeno(colony, caste = "drones"),
                       getCasteSegSiteGeno(colony, caste = "fathers"))
  # Collect the sex of all the colony members in the order than genotypes
  sex_colony <- c(rep("F", nQueens(colony)),
                  rep("F", nWorkers(colony)),
                  rep("M", nDrones(colony)),
                  rep("M", nFathers(colony)))
  # Compute the IBS relationship matrix
  # TODO: Combine generations 1 and 10 - so that they have the same reference population (regarding allele frequencies)
  ibs_colony <- calcBeeGRMIbs(x = colony_geno,
                              sex = sex_colony)
  # Collect the IBD haplotypes for all the colony members
  colony_haplo <- rbind(getCasteIbdHaplo(colony, caste = "queen"),
                        getCasteIbdHaplo(colony, caste = "workers"),
                        getCasteIbdHaplo(colony, caste = "drones"),
                        getCasteIbdHaplo(colony, caste = "fathers"))
  # Compute the IBD relationship matrix
  ibd_colony <- calcBeeGRMIbd(x = colony_haplo)
  ibd_colony <- ibd_colony$indiv

  # Only chromosome 3
  ibd_colonychr3 <- calcBeeGRMIbd(x = colony_haplo[, grepl(pattern = "3_",
                                                           x = colnames(colony_haplo))])
  # Only csd locus
  ibd_colonycsd <- calcBeeGRMIbd(x = colony_haplo[, paste(SP$csdChr,
                                                          SP$csdPosStart:SP$csdPosStop,
                                                          sep = "_")])

  colony_id <- getCasteId(colony, caste = "all")

  return(list(IBS = ibs_colony, IBD = ibd_colony, IBDChr = ibd_colonychr3, IBDCsd = ibd_colonycsd, ID = colony_id))
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

  tmp <- makeS(pedigree = pedigree, heterogametic = "1", returnS = TRUE)
  IBDe <- tmp$S
  dimnames(IBDe) <- list(rownames(pedigree), rownames(pedigree))
  return(IBDe)
}

getCsdInfo <- function (colonies, subspecies = NULL) {
  totalCsd <- nCsdAlleles(colonies, collapse = TRUE)
  for (n in 1:nColonies(colonies)){
    csdVariability <- c(Rep = Rep, year = year, id = colonies[[n]]@id,
                        nCSD = nCsdAlleles(colonies[[n]], collapse = TRUE), totalCSD = totalCsd, subspecies = subspecies)
    pDiploidDrones <- c(Rep = Rep, year = year, id = colonies[[n]]@id,
                        pQueenHomBrood = computeQueensPHomBrood(colonies), subspecies = subspecies)
    return(list(csdVariability = csdVariability, pDiploidDrones = pDiploidDrones))
  }
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

# Population parameters -------------------------------------------------------------------
# Number of repeats
nRep <- 1
# Number of years
nYear <- 10
# Number of colonies in the apiary
apiarySize <- 10
# Number of workers in a full colony
nWorkers <- 0 # TODO: change to 20K
# Number of drones in a full colony
nDrones <- 100 #nWorkers * 0.2
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
pImport <- 0.2

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
  # Create a founder population of A. m. mellifera [1:30] and carnica [31:60] bees
  # founderGenomes <- simulateHoneyBeeGenomes(nMelN = 30,
  #                                           nCar = 30,
  #                                           nChr = 1,
  #                                           nSegSites = 100)
  load("~/Documents/TwoPopFounderPop_Bee.Rdata")
  # Create SP object and write in the global simulation/population parameters
  SP <- SimParamBee$new(founderGenomes, csdChr = 3, nCsdAlleles = 128)
  SP$nWorkers <- nWorkers
  SP$nDrones <- nDrones
  SP$nFathers <- pFathers
  SP$nVirginQueens <- nVirginQueens
  SP$pSwarm <- 0.5
  SP$pSplit <- 0.3
  # Track the pedigree
  SP$setTrackPed(TRUE)
  # Track the recombination
  SP$setTrackRec(TRUE)
  # Add a SNP chip with 3 SNPs per chromosome
  SP$addSnpChip(nSnpPerChr = 3)

  # Create a base population for A. m. mellifera
  melVirginQueens <- createVirginQueens(x = founderGenomes[1:30])
  # Create a base population for A. m. carnica
  carVirginQueens <- createVirginQueens(x = founderGenomes[30:60])
  # Create A. m. mellifera drones
  melDrones <- createDrones(x = melVirginQueens[11:30], nInd = 20)
  # Create A. m. carnica drones
  carDrones <- createDrones(x = carVirginQueens[11:30], nInd = 20)
  # Mate A. m. mellifera queens with A. m. mellifera drones
  melQueens <- crossVirginQueen(pop = melVirginQueens[1:10], drones = melDrones)
  # Mate A. m. carnica queens with A. m. carnica drones
  carQueens <- crossVirginQueen(pop = carVirginQueens[1:10], drones = carDrones)

  # Start the year-loop ------------------------------------------------------------------
  for (year in 1:nYear) {
    # year <- 1
    # year <- year + 1
    cat(paste0("Year: ", year, "/", nYear, "\n"))
    # If this is the first year, create some colonies to start with
    # If not, promote the age0 to age1, age1 to age2 and remove age2 colonies
    if (year == 1) {
      age1 <- list(Mel = createColonies(x = melQueens, n = apiarySize,
                                        simParamBee = SP),
                   Car = createColonies(x = carQueens, n = apiarySize,
                                        simParamBee = SP))
    } else {
      age2 <- list(Mel = age1$Mel, Car = age1$Car)
      age1 <- list(Mel = age0$Mel, Car = age0$Car)
      age0 <- list(Mel = NULL, Car = NULL)
      age0p1 <- list(Mel = NULL, Car = NULL)
      age0p2 <- list(Mel = NULL, Car = NULL)
    }

    # In year 1, inspect the relationship in one of the colonies
    if (year == 1) {
      # Choose the first colony of age 1 to inspect relationship in the base population
      springerColony1_Mel <- computeRelationship_genomic(colony = age1$Mel[[1]])
      springerColony1_Car <- computeRelationship_genomic(colony = age1$Car[[1]])
    }

    # Period1 ------------------------------------------------------------------
    # Build-up the colonies
    age1 <- list(Mel = buildUpColonies(age1$Mel),
                 Car = buildUpColonies(age1$Car))
    if (year > 1) {
      age2 <- list(Mel = buildUpColonies(age2$Mel),
                   Car = buildUpColonies(age2$Car))
    }

    # Split all age1 colonies
    tmpMel <- splitColonies(age1$Mel)
    tmpCar <- splitColonies(age1$Car)
    age1 <- list(Mel = tmpMel$remnants,
                 Car = tmpCar$remnants)
    # The queens of the splits are 0 years old
    age0p1 <- list(Mel = tmpMel$splits, Car = tmpCar$splits)

    if (year > 1) {
      # Split all age2 colonies
      tmpMel <- splitColonies(age2$Mel)
      tmpCar <- splitColonies(age2$Car)
      age2 <- list(Mel = tmpMel$remnants,
                   Car = tmpCar$remnants)
      # The queens of the splits are 0 years old
      age0p1 <- list(Mel = c(age0p1$Mel, tmpMel$splits), Car = c(age0p1$Car, tmpCar$splits))
    }

    # Create virgin queens
    # Sample colony for the virgin queens
    # TODO: this is likely an inbreeding bottleneck - just one colony producing all virgin queens!
    #       maybe we should use selectColonies() and then create virginQueens from them, but we
    #       need createVirginQueens() to work with multiple colonies too, and equally for createWorkers()
    #       and createDrones()
    virginDonor <- list(Mel = sample.int(n = nColonies(age1$Mel), size = 1),
                        Car = sample.int(n = nColonies(age1$Car), size = 1))
    # Virgin queens for splits!
    virginQueens <- list(Mel = createVirginQueens(age1$Mel[[virginDonor$Mel]], nInd = nColonies(age0p1$Mel)),
                         Car = createVirginQueens(age1$Car[[virginDonor$Car]], nInd = nColonies(age0p1$Car)))

    # Requeen the splits --> queens are now 0 years old
    age0p1 <- list(Mel = reQueenColonies(age0p1$Mel, queens = virginQueens$Mel),
                   Car = reQueenColonies(age0p1$Car, queens = virginQueens$Car))

    # Swarm a percentage of age1 colonies
    tmpMel <- pullColonies(age1$Mel, p = p1swarm)
    tmpCar <- pullColonies(age1$Car, p = p1swarm)
    age1 <- list(Mel = tmpMel$remainingColonies,
                 Car = tmpCar$remainingColonies)
    tmpMel <- swarmColonies(tmpMel$pulledColonies)
    tmpCar <- swarmColonies(tmpCar$pulledColonies)
    age0p1 <- list(Mel = c(age0p1$Mel, tmpMel$remnants),
                   Car = c(age0p1$Car, tmpCar$remnants))
    age1 <- list(Mel = c(age1$Mel, tmpMel$swarms),
                 Car = c(age1$Car, tmpCar$swarms))


    if (year > 1) {
      # Swarm a percentage of age2 colonies
      tmpMel <- pullColonies(age2$Mel, p = p1swarm)
      tmpCar <- pullColonies(age2$Car, p = p1swarm)
      age2 <- list(Mel = tmpMel$remainingColonies,
                   Car = tmpCar$remainingColonies)
      tmpMel <- swarmColonies(tmpMel$pulledColonies)
      tmpCar <- swarmColonies(tmpCar$pulledColonies)
      age0p1 <- list(Mel = c(age0p1$Mel, tmpMel$remnants),
                     Car = c(age0p1$Car, tmpCar$remnants))
      age2 <- list(Mel = c(age2$Mel, tmpMel$swarms),
                   Car = c(age2$Car, tmpCar$swarms))
    }

    # Supersede age1 colonies
    tmpMel <- pullColonies(age1$Mel, p = p1supersede)
    tmpCar <- pullColonies(age1$Car, p = p1supersede)
    age1 <- list(Mel = tmpMel$remainingColonies,
                 Car = tmpCar$remainingColonies)
    tmpMel <- supersedeColonies(tmpMel$pulledColonies)
    tmpCar <- supersedeColonies(tmpCar$pulledColonies)
    age0p1 <- list(Mel = c(age0p1$Mel, tmpMel),
                   Car = c(age0p1$Car, tmpCar))


    if (year > 1) {
      # Supersede age2 colonies
      tmpMel <- pullColonies(age2$Mel, p = p1supersede)
      tmpCar <- pullColonies(age2$Car, p = p1supersede)
      age2 <- list(Mel = tmpMel$remainingColonies,
                   Car = tmpCar$remainingColonies)
      tmpMel <- supersedeColonies(tmpMel$pulledColonies)
      tmpCar <- supersedeColonies(tmpCar$pulledColonies)
      age0p1 <- list(Mel = c(age0p1$Mel, tmpMel),
                     Car = c(age0p1$Car, tmpCar))
    }

    # Mate the split colonies
    if (year == 1) {
      DCAMel <- createDCA(c(age1$Mel,
                            selectColonies(age1$Car, n = round(nColonies(age1$Mel) * pImport, 0))))
      age0p1$Mel <- crossColonies(age0p1$Mel, drones = DCAMel, nFathers = SP$nFathers) #TODO: Remove this
      DCACar <- createDCA(age1$Car)
      age0p1$Car <- crossColonies(age0p1$Car, drones = DCACar, nFathers = SP$nFathers) #TODO: Remove this
    } else {
      DCAMel <- createDCA(c(age1$Mel,
                            selectColonies(age1$Car, n = round(nColonies(age1$Mel) * pImport, 0)),
                            age2$Mel,
                            selectColonies(age2$Car, n = round(nColonies(age2$Mel) * pImport, 0))))
      age0p1$Mel <- crossColonies(age0p1$Mel, drones = DCAMel, nFathers = SP$nFathers) #TODO: Remove this
      DCACar <- createDCA(c(age1$Car, age2$Car))
      age0p1$Car <- crossColonies(age0p1$Car, drones = DCACar, nFathers = SP$nFathers) #TODO: Remove this
    }

    # Collapse
    age1 <- list(Mel = selectColonies(age1$Mel, p = 1 - p1collapse),
                 Car = selectColonies(age1$Car, p = 1 - p1collapse))
    if (year > 1) {
      age2 <- list(Mel = selectColonies(age2$Mel, p = 1 - p1collapse),
                   Car = selectColonies(age2$Car, p = 1 - p1collapse))
    }

    # Period2 ------------------------------------------------------------------

    # Swarm a percentage of age1 colonies
    # Mellifera
    tmpMel <- pullColonies(age1$Mel, p = p2swarm)
    tmpCar <- pullColonies(age1$Car, p = p2swarm)
    age1 <- list(Mel = tmpMel$remainingColonies,
                 Car = tmpCar$remainingColonies)
    tmpMel <- swarmColonies(tmpMel$pulledColonies)
    tmpCar <- swarmColonies(tmpCar$pulledColonies)
    # The queens of the remnant colonies are of age 0
    age0p2 <- list(Mel = tmpMel$remnants,
                   Car = tmpCar$remnants)
    age1 <- list(Mel = c(age1$Mel, tmpMel$swarms),
                 Car = c(age1$Car, tmpCar$swarms))

    if (year > 1) {
      # Swarm a percentage of age2 colonies
      tmpMel <- pullColonies(age2$Mel, p = p2swarm)
      tmpCar <- pullColonies(age2$Car, p = p2swarm)
      age2 <- list(Mel = tmpMel$remainingColonies,
                   Car = tmpCar$remainingColonies)
      tmpMel <- swarmColonies(tmpMel$pulledColonies)
      tmpCar <- swarmColonies(tmpCar$pulledColonies)
      # The queens of the remnant colonies are of age 0
      age0p2 <- list(Mel = tmpMel$remnants,
                     Car = tmpCar$remnants)
      age2 <- list(Mel = c(age2$Mel, tmpMel$swarms),
                   Car = c(age2$Car, tmpCar$swarms))
    }

    # Supersede a part of age1 colonies
    tmpMel <- pullColonies(age1$Mel, p = p2supersede)
    tmpCar <- pullColonies(age1$Car, p = p2supersede)
    age1 <- list(Mel = tmpMel$remainingColonies,
                 Car = tmpCar$remainingColonie)
    tmpMel <- supersedeColonies(tmpMel$pulledColonies)
    tmpCar <- supersedeColonies(tmpCar$pulledColonies)
    # The queens of superseded colonies are of age 0
    age0p2 <- list(Mel = c(age0p2$Mel, tmpMel),
                   Car = c(age0p2$Car, tmpCar))

    if (year > 1) {
      # Supersede a part of age2 colonies
      tmpMel <- pullColonies(age2$Mel, p = p2supersede)
      tmpCar <- pullColonies(age2$Car, p = p2supersede)
      age2 <- list(Mel = tmpMel$remainingColonies,
                   Car = tmpCar$remainingColonies)
      tmpMel <- supersedeColonies(tmpMel$pulledColonies)
      tmpCar <- supersedeColonies(tmpCar$pulledColonies)
      # The queens of superseded colonies are of age 0
      age0p2 <- list(Mel = c(age0p2$Mel, tmpMel),
                     Car = c(age0p2$Car, tmpCar))
    }

    # Replace all the drones
    age1$Mel <- replaceDrones(age1$Mel)
    age1$Car <- replaceDrones(age1$Car)
    if (year > 1) {
      age2$Mel <- replaceDrones(age2$Mel)
      age2$Car <- replaceDrones(age2$Car)
    }

    # Mate the colonies
    # Import p percentage of carnica colonies into mellifera DCA
    if (year == 1) {
      DCAMel <- createDCA(c(age1$Mel,
                            selectColonies(age1$Car, n = round(nColonies(age1$Mel) * pImport, 0))))
      DCACar <- createDCA(age1$Car)
    } else {
      DCAMel <- createDCA(c(age1$Mel,
                            selectColonies(age1$Car, n = round(nColonies(age1$Mel) * pImport, 0)),
                            age2$Mel,
                            selectColonies(age2$Car, n = round(nColonies(age2$Mel) * pImport, 0))))
      DCACar <- createDCA(c(age1$Car, age2$Car))
    }

    # Cross age 0 period 2 swarms and splits
    age0p2$Mel <- crossColonies(age0p2$Mel, drones = DCAMel, nFathers = nFathers) #TODO: REMOVE
    age0p2$Car <- crossColonies(age0p2$Car, drones = DCACar, nFathers = nFathers) #TODO: REMOVE

    # Collapse
    age1 <- list(Mel= selectColonies(age1$Mel, p = 1 - p2collapse),
                 Car = selectColonies(age1$Car, p = 1 - p2collapse))
    if (year > 1) {
      age2 <- list(Mel = selectColonies(age2$Mel, p = 1 - p2collapse),
                   Car = selectColonies(age2$Car, p = 1 - p2collapse))
    }

    # Merge all age 0 colonies (from both periods)
    age0 <- list(Mel = c(age0p1$Mel, age0p2$Mel),
                 Car = c(age0p1$Car, age0p2$Car))

    # Period3 ------------------------------------------------------------------
    # Collapse age0 queens
    age0 <- list(Mel = selectColonies(age0$Mel, p = (1 - p3collapseAge0)),
                 Car = selectColonies(age0$Car, p = (1 - p3collapseAge0)))
    age1 <- list(Mel = selectColonies(age1$Mel, p = (1 - p3collapseAge1)),
                 Car = selectColonies(age1$Car, p = (1 - p3collapseAge1)))
    age2 <- list(Mel = NULL, Car = NULL) #We don't need this but just to show the workflow!!!

    #Collect CSD info
    csdInfoMel <- getCsdInfo(age0$Mel, subspecies = "Mel")
    csdInfoCar <- getCsdInfo(age0$Car, subspecies = "Car")
    csdVariability <- rbind(csdVariability, csdInfoMel$csdVariability)
    csdVariability <- rbind(csdVariability, csdInfoCar$csdVariability)
    pDiploidDrones <-rbind(pDiploidDrones, csdInfoMel$pDiploidDrones)
    pDiploidDrones <-rbind(pDiploidDrones, csdInfoCar$DiploidDrones)


    # Maintain the number of colonies ------------------------------------------
    # Keep all of age1, age0 swarmed so we build it up with some splits, while we remove (sell) the other splits
    age0$Mel <- maintainApiarySize(age0 = age0$Mel, age1 = age1$Mel)
    age0$Car <- maintainApiarySize(age0 = age0$Car, age1 = age1$Car)

  } # Year-loop

  data.frame(Rep = NA, Year = NA, Age0 = NA, Age1 = NA, sum = NA, subspecies = NA)
  for (subspecies in c("Mel", "Car")) {
    if ((nColonies(age0[[subspecies]]) + nColonies(age1[[subspecies]])) != apiarySize) {
      stop(paste0("The number of colonies for ", subspecies, " does not match the apiary size!"))
    }
  }

  a <- toc()
  loopTime <- rbind(loopTime, c(Rep, a$tic, a$toc, a$msg, (a$toc - a$tic)))

  # Take the first colony of age1 and extract genotypes for relationship computation
  springerColony10_Mel <- computeRelationship_genomic(colony = age1$Mel[[1]])
  springerColony10_Car <- computeRelationship_genomic(colony = age1$Car[[1]])

  # Compute the pedigree relationship matrix
  IBDe <- computeRelationship_pedigree(SP$pedigree)

} # Rep-loop


save.image(file = "SpringerSimulation_import.Rdata")






