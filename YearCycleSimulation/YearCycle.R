# Clean workspace
rm(list = ls())

# Load packages
library(AlphaSimR)
library(ggplot2)
library(tictoc)
library(R6)

# TODO: replace with devtools installation from Github once the package is operational
# Source the development version of AlphaSimR
# Laura
#AlphaSimRBeeFolder <- "~/Desktop/GitHub/Fork/AlphaSimRBee"
# Gregor
#AlphaSimRBeeFolder <- "~/Documents/5_Storages/GitBox/AlphaSimRBee/AlphaSimRBee"
# Jernej
#AlphaSimRBeeFolder <- "C:/Users/jernejb/Desktop/git/AlphaSimRBee/SIMplyBee"

#source(paste0(AlphaSimRBeeFolder, "/R/Class-SimParamBee.R"))
#source(paste0(AlphaSimRBeeFolder, "/R/Class-Colony.R"))
#source(paste0(AlphaSimRBeeFolder, "/R/Class-Colonies.R"))
#source(paste0(AlphaSimRBeeFolder, "/R/Functions_L0_auxilary.R"))
#source(paste0(AlphaSimRBeeFolder, "/R/Functions_L1_Pop.R"))
#source(paste0(AlphaSimRBeeFolder, "/R/Functions_L2_Colony.R"))
#source(paste0(AlphaSimRBeeFolder, "/R/Functions_L3_Colonies.R"))

library(SIMplyBee)

# Parameters -------------------------------------------------------------------

apiarySize <- 20
nWorkersFull <- 20 # TODO: change to 20K
nAvgFathers <- 15
nDronesFull <- nWorkersFull * 0.2
if (nDronesFull < (nAvgFathers * 2)) {
  nDronesFull <-nAvgFathers * 2
}

nSNPChr <- 3
nChromo <- 16
# Period1
p1swarm <- 0.05
p1supersede <- 0.05
p1collapse <- 0.10

# Period2
p2swarm <- 0.01
p2supersede <- p1supersede
p2collapse <- p1collapse

# Period3
p3collapseAge0 <- 0.25
p3collapseAge1 <- 0.3
# TODO: Change into a vector of probabilities (age 1,2,3 of the queen)

# Create df for recording the number of age0 and age1 colonies, csd variability and for recording cpu time
loopTime <- data.frame(Rep = NA, tic = NA, toc = NA, msg = NA, time = NA)
noQueens <- data.frame(Rep = NA, Year = NA, Age0 = NA, Age1 = NA, sum = NA)
csdVariability <- data.frame(Rep = NA, year = NA, id = NA, nCSDage0 = NA, totalCSDage0 = NA)
pDiploidDrones <- data.frame(Rep = NA, year = NA, id = NA, pDidrA0 = NA, pDidrA1 = NA)
ped <- data.frame(ID = NA, mother=NA, father=NA, isDH=NA)


# Rep-loop ---------------------------------------------------------------------

nRep <- 1 
for (Rep in 1:nRep) {
  # Rep <- 1
  cat(paste0("Rep: ", Rep, "/", nRep, "\n"))
  # Measure cpu time
  tic('20y loop')
  
  # prepare matrix for genotypes
  
  SNPgenoAllMat <- matrix(ncol = (nSNPChr*nChromo))
  # Founder population ---------------------------------------------------------

  founderGenomes <- quickHaplo(nInd = 1000,
                               nChr = 16,
                               segSites = 3)
  SP <- SimParamBee$new(founderGenomes, csdChr = NULL)
  SP$setTrackPed(TRUE)
  SP$addSnpChip(nSnpPerChr = 3)
  base <- asVirginQueen(newPop(founderGenomes))

  # Year-loop ------------------------------------------------------------------

  nYear <- 3
  for (year in 1:nYear) {
    # year <- 1
    # year <- year + 1
    cat(paste0("Year: ", year, "/", nYear, "\n"))
    if (year == 1) {

      #age1 <- createColonies(pop = selectInd(base, nCol = apiarySize * 2, use = "rand"),
                            # n = apiarySize, nAvgFathers = nAvgFathers)
      age1 <- createColonies(pop = base, nCol = apiarySize, mated = TRUE,
                                         nAvgFathers = 15, nDronesPerQueen = 100,
                                         simParamBee = SP)

    } else {
      age2 <- age1
      age1 <- age0
      age0 <- NULL
      age0p1 <- NULL
      age0p2 <- NULL
    }

    # Period1 ------------------------------------------------------------------

    # Build-up the colonies
    age1 <- buildUpColonies(age1, nWorkers = nWorkersFull, nDrones = nDronesFull, exact = FALSE)
    if (year > 1) {
      age2 <- buildUpColonies(age2, nWorkers = nWorkersFull, nDrones = nDronesFull, exact = FALSE)
    }

    # Split all age1 colonies
    tmp <- splitColonies(age1)
    age1 <- tmp$remnants
    # The queens of the splits are 0 years old
    age0p1 <- tmp$splits

    if (year > 1) {
      # Split all age2 colonies
      tmp <- splitColonies(age2)
      age2 <- tmp$remnants
      # The queens of the splits are 0 years old
      age0p1 <- c(age0p1, tmp$splits)
    }

    # Create virgin queens
    # Sample colony for the virgin queens
    # TODO: this is likely an inbreeding bottleneck - just one colony producing all virgin queens!
    #       maybe we should use selectColonies() and then create virginQueens from them, but we
    #       need createVirginQueens() to work with multiple colonies too, and equally for createWorkers()
    #       and createDrones()
    virginDonor <- sample.int(n = nColonies(age1), size = 1)
    virginQueens <- createVirginQueens(age1[[virginDonor]], nInd = nColonies(age0p1) * 10)

    # Requeen the splits --> queens are now 0 years old
    age0p1 <- reQueenColonies(age0p1, queens = virginQueens$virginQueens)

    # Swarm a percentage of age1 colonies

    TMP <- lapply(age1@colonies, FUN = addVirginQueens, nInd = 50)
    age1@colonies <- TMP

    if (year > 1) {
    TMP <- lapply(age2@colonies, FUN = addVirginQueens, nInd = 50)
    age2@colonies <- TMP
    }

    tmp <- pullColonies(age1, p = p1swarm)
    age1 <- tmp$remainingColonies
    tmp <- swarmColonies(tmp$pulledColonies)
    age0p1 <- c(age0p1, tmp$remnants)
    age1 <- c(age1, tmp$swarms)

    if (year > 1) {
      # Swarm a percentage of age2 colonies
      tmp <- pullColonies(age2, p = p1swarm)
      age2 <- tmp$remainingColonies
      tmp <- swarmColonies(tmp$pulledColonies)
      age0p1 <- c(age0p1, tmp$remnants)
      age2 <- c(age2, tmp$swarms)
    }

    # Supersede age1 colonies


    TMP <- lapply(age1@colonies, FUN = addVirginQueens, nInd = 50)
    age1@colonies <- TMP

    if (year > 1) {
    TMP <- lapply(age2@colonies, FUN = addVirginQueens, nInd = 50)
    age2@colonies <- TMP
    }

    tmp <- pullColonies(age1, p = p1supersede)
    age1 <- tmp$remainingColonies
    tmp <- supersedeColonies(tmp$pulledColonies)
    age0p1 <- c(age0p1, tmp)

    if (year > 1) {
      # Supersede age2 colonies
      tmp <- pullColonies(age2, p = p1supersede)
      age2 <- tmp$remainingColonies
      tmp <- supersedeColonies(tmp$pulledColonies)
      age0p1 <- c(age0p1, tmp)
    }

    # Mate the split colonies
    if (year == 1) {
      DCA <- createDCA(age1, nInd = nDronesFull)
    } else {
      DCA <- createDCA(c(age1, age2), nInd = nDronesFull)
    }
    age0p1 <- crossColonies(age0p1, DCA = DCA, nAvgFathers = nAvgFathers)

    # Collapse
    age1 <- selectColonies(age1, p = 1 - p1collapse)
    if (year > 1) {
      age2 <- selectColonies(age2, p = 1 - p1collapse)
    }

    # Period2 ------------------------------------------------------------------

    # Swarm a percentage of age1 colonies


    TMP <- lapply(age1@colonies, FUN = addVirginQueens, nInd = 50)
    age1@colonies <- TMP

    if (year > 1) {
    TMP <- lapply(age2@colonies, FUN = addVirginQueens, nInd = 50)
    age2@colonies <- TMP
    }

    tmp <- pullColonies(age1, p = p2swarm)
    age1 <- tmp$remainingColonies
    tmp <- swarmColonies(tmp$pulledColonies)
    # The queens of the remnant colonies are of age 0
    age0p2 <- tmp$remnants
    age1 <- c(age1, tmp$swarms)

    if (year > 1) {
      # Swarm a percentage of age2 colonies
      tmp <- pullColonies(age2, p = p2swarm)
      age2 <- tmp$remainingColonies
      tmp <- swarmColonies(tmp$pulledColonies)
      # The queens of the remnant colonies are of age 0
      age0p2 <- c(age0p2, tmp$remnants)
      age2 <- c(age2, tmp$swarms)
    }

    # Supersede a part of age1 colonies


    TMP <- lapply(age1@colonies, FUN = addVirginQueens, nInd = 50)
    age1@colonies <- TMP

    if (year > 1) {
    TMP <- lapply(age2@colonies, FUN = addVirginQueens, nInd = 50)
    age2@colonies <- TMP
    }

    tmp <- pullColonies(age1, p = p2supersede)
    age1 <- tmp$remainingColonies
    tmp <- supersedeColonies(tmp$pulledColonies)
    # The queens of superseded colonies are of age 0
    age0p2 <- c(age0p2, tmp)

    if (year > 1) {
      # Supersede a part of age2 colonies
      tmp <- pullColonies(age2, p = p2supersede)
      age2 <- tmp$remainingColonies
      tmp <- supersedeColonies(tmp$pulledColonies)
      # The queens of superseded colonies are of age 0
      age0p2 <- c(age0p2, tmp)
    }

    # Replace all the drones
    age1 <- replaceDrones(age1)
    if (year > 1) {
      age2 <- replaceDrones(age2)
    }

    # Mate the colonies
    if (year == 1) {
      DCA <- createDCA(age1, nInd = nDronesFull)
    } else {
      DCA <- createDCA(c(age1, age2), nInd = nDronesFull)
    }
    # Cross age 0 period 2 swarms and splits
    age0p2 <- crossColonies(age0p2, DCA = DCA, nAvgFathers = nAvgFathers)

    # Collapse
    age1 <- selectColonies(age1, p = 1 - p2collapse)
    if (year > 1) {
      age2 <- selectColonies(age2, p = 1 - p2collapse)
    }

    # Merge all age 0 colonies (from both periods)
    age0 <- c(age0p1, age0p2)

    # Period3 ------------------------------------------------------------------

    # Get the intracolonial csd variability-------------------------------------
   # totalCsdAge0 <- nCsdAlleles(age0, collapse = TRUE)
    #for (n in 1:nColonies(age0)){
     # csdVariability <- rbind(csdVariability, c(Rep, year, age0[[n]]@id,
      #                                          nCsdAlleles(age0[[n]], collapse = TRUE), totalCsdAge0))
      #pDiploidDrones <- rbind( c(Rep = Rep, year = year, id = age0[[n]]@id, pDidrA0 = phb0, pDidrA1 = phb1))
    #}


    # Collapse age0 queens
    age0 <- selectColonies(age0, p = (1 - p3collapseAge0))
    age1 <- selectColonies(age1, p = (1 - p3collapseAge1))
    age2 <- NULL #We don't need this but just to show the workflow!!!



    # Maintain the number of colonies ------------------------------------------

    # keep all of age1, age0 swarmed so we build it up with some splits, while we remove (sell) the other splits
    if ((nColonies(age0) + nColonies(age1)) > apiarySize) { # check if the sum of all colonies is greater than apiary size
      IDsplits <- getId(age0)[hasSplit(age0)] # get the IDs of age 0 that are splits
      splits0 <- pullColonies(age0, ID = IDsplits) # pull the splits out of age 0
      age0split <- splits0$pulledColonies # create an object for age 0 splits
      age0swarm <- splits0$remainingColonies # create an object for swarms and superseded colonies
      # age0 <- NULL
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
    }
  # extract genotype ----------------------------------------------------------- 
    SNPgenoWA0 <- do.call("rbind", getCasteSnpGeno(age0, caste = "workers", nInd = 2)) #na tak nači, še za ostale
    SNPgenoDA0 <- do.call("rbind", getCasteSnpGeno(age0, caste = "drones", nInd = 2))
    SNPgenoFA0 <- do.call("rbind", getCasteSnpGeno(age0, caste = "fathers", nInd = 1))
    SNPgenoQA0 <- do.call("rbind", getCasteSnpGeno(age0, caste = "queen"))
    
    SNPgenoAllMat <- do.call("rbind", list(SNPgenoWA0, SNPgenoDA0, SNPgenoFA0, SNPgenoQA0))
    
    
  } # Year-loop
  noQueens <- rbind(noQueens, c(Rep, year, nColonies(age0), nColonies(age1), (nColonies(age0) + nColonies(age1))))
  a <- toc()
  loopTime <- rbind(loopTime, c(Rep, a$tic, a$toc, a$msg, (a$toc - a$tic)))
  
  # write getoype info
  
  #write.csv(SNPgenoAllMat, paste("SNPGeno", rep, ".csv"), quote = FALSE, row.names = FALSE)
  
  
  
} # Rep-loop

#polot the number of queens
#ggplot(noQueens, aes(Rep, sum)) +
#  geom_line(aes(Rep, sum))

# write pedigree to df
ped <- SP$pedigree

#write.csv(noQueens, paste0("NoQueens", rep, ".csv"), quote = FALSE, row.names = FALSE)
#write.csv(csdVariability, paste0("CsdVariability", rep, ".csv"), quote = FALSE, row.names = FALSE)
#write.csv(pDiploidDrones, paste0("pDiploidDrones", rep, ".csv"), quote = FALSE, row.names = FALSE)
#write.csv(ped, paste0("ped", rep, ".csv"), quote = FALSE, row.names = FALSE)


#### rabmo neko tabelo z kastami, individual id in colony ID
 caste <- getCaste(age1)
 indID <- getId(age1)
 colID <- age1

