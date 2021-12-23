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
AlphaSimRBeeFolder <- "C:/Users/jernejb/Desktop/git/AlphaSimRBee/SIMplyBee"

source(paste0(AlphaSimRBeeFolder, "/R/Class-SimParamBee.R"))
source(paste0(AlphaSimRBeeFolder, "/R/Class-Colony.R"))
source(paste0(AlphaSimRBeeFolder, "/R/Class-Colonies.R"))
source(paste0(AlphaSimRBeeFolder, "/R/Functions_L0_auxilary.R"))
source(paste0(AlphaSimRBeeFolder, "/R/Functions_L1_Pop.R"))
source(paste0(AlphaSimRBeeFolder, "/R/Functions_L2_Colony.R"))
source(paste0(AlphaSimRBeeFolder, "/R/Functions_L3_Colonies.R"))

library(SIMplyBee)

# Parameters -------------------------------------------------------------------

apiarySize <- 20
nWorkersFull <- 200 # TODO: change to 20K
nAvgFathers <- 15
nDronesFull <- nWorkersFull * 0.2
if (nDronesFull < (nAvgFathers * 2)) {
  nDronesFull <-nAvgFathers * 2
}
nVirginQueensFull <- 30

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
noQueens <- data.frame(Rep = NA, Age0 = NA, Age1 = NA, sum = NA)
csdVariability <- data.frame(Rep = NA, year = NA, avgCSDage0 = NA )

# Rep-loop ---------------------------------------------------------------------

nRep <- 1
for (Rep in 1:nRep) {
  # Rep <- 1
  cat(paste0("Rep: ", Rep, "/", nRep, "\n"))
  # Measure cpu time
  tic('20y loop')

  # Founder population ---------------------------------------------------------

  founderGenomes <- quickHaplo(nInd = apiarySize * 2,
                               nChr = 16,
                               segSites = 1000)
  SP <- SimParamBee$new(founderGenomes, csdChr = 3, nCsdAlleles = 128)
  SP <- SimParamBee$new(founderGenomes, csdChr = NULL, nCsdAlleles = 128)
  SP$nWorkers <- nWorkersFull
  SP$nDrones <- nDronesFull
  SP$nVirginQueens <- nVirginQueensFull
  base <- newPop(founderGenomes)

  # Year-loop ------------------------------------------------------------------

  nYear <- 20
  for (year in 1:nYear) {
    # year <- 1
    # year <- year + 1
    cat(paste0("Year: ", year, "/", nYear, "\n"))
    if (year == 1) {
      age1 <- createColonies(pop = base, nCol = apiarySize)
    } else {
      age2 <- age1
      age1 <- age0
      age0 <- NULL
      age0p1 <- NULL
      age0p2 <- NULL
    }

    # Period1 ------------------------------------------------------------------
    cat("Period1 -----------------------------------------------------------\n")

    cat("# Build-up the colonies\n")
    age1 <- buildUpColonies(age1)
    if (year > 1) {
      age2 <- buildUpColonies(age2)
    }

    cat("# Split all age1 colonies\n")
    tmpS1 <- splitColonies(age1)
    age1 <- tmpS1$remnants
    age0p1 <- tmpS1$splits # The queens of the splits are 0 years old

    if (year > 1) {
      cat("# Split all age2 colonies\n")
      tmpS2 <- splitColonies(age2)
      age2 <- tmpS2$remnants
      age0p1 <- c(age0p1, tmpS2$splits) # The queens of the splits are 0 years old
    }

    cat("# Create virgin queens\n")
    # Sample colony for the virgin queens
    virginDonor <- selectColonies(age1, n = 1)
    virginQueens <- getVirginQueens(virginDonor[[1]])

    cat("# Requeen the splits\n")
    age0p1 <- reQueenColonies(age0p1, queens = virginQueens) # queens are now 0 years old

    cat("# Swarm a percentage of age1 colonies\n")
    tmpPa <- pullColonies(age1, p = p1swarm)
    age1 <- tmpPa$remainingColonies
    tmpW1 <- swarmColonies(tmpPa$pulledColonies)
    age0p1 <- c(age0p1, tmpW1$remnants)
    age1 <- c(age1, tmpW1$swarms)

    if (year > 1) {
      cat("# Swarm a percentage of age2 colonies\n")
      tmpPb <- pullColonies(age2, p = p1swarm)
      age2 <- tmpPb$remainingColonies
      tmpW2 <- swarmColonies(tmpPb$pulledColonies)
      age0p1 <- c(age0p1, tmpW2$remnants)
      age2 <- c(age2, tmpW2$swarms)
    }

    cat("# Supersede age1 colonies\n")
    tmpPc <- pullColonies(age1, p = p1supersede)
    age1 <- tmpPc$remainingColonies
    tmpQ1 <- supersedeColonies(tmpPc$pulledColonies)
    age0p1 <- c(age0p1, tmpQ1)

    if (year > 1) {
      cat("# Supersede age2 colonies\n")
      tmpPd <- pullColonies(age2, p = p1supersede)
      age2 <- tmpPd$remainingColonies
      tmpQ2 <- supersedeColonies(tmpPd$pulledColonies)
      age0p1 <- c(age0p1, tmpQ2)
    }

    cat("# Mate the split colonies\n")
    if (year == 1) {
      DCA <- createDCA(age1, nInd = nDronesFull)
    } else {
      DCA <- createDCA(c(age1, age2), nInd = nDronesFull)
    }
    age0p1 <- crossColonies(age0p1, DCA = DCA, nAvgFathers = nAvgFathers)

    cat("# Collapse collonies\n")
    age1 <- selectColonies(age1, p = 1 - p1collapse)
    if (year > 1) {
      age2 <- selectColonies(age2, p = 1 - p1collapse)
    }

    # Period2 ------------------------------------------------------------------
    cat("Period2 -----------------------------------------------------------\n")

    cat("# Swarm a percentage of age1 colonies\n")
    tmpPa <- pullColonies(age1, p = p2swarm)
    age1 <- tmpPa$remainingColonies
    tmpW1 <- swarmColonies(tmpPa$pulledColonies)
    age0p2 <- tmpW1$remnants # The queens of the remnant colonies are of age 0
    age1 <- c(age1, tmpW1$swarms)

    if (year > 1) {
      cat("# Swarm a percentage of age2 colonies\n")
      tmpPb <- pullColonies(age2, p = p2swarm)
      age2 <- tmpPb$remainingColonies
      tmpW2 <- swarmColonies(tmpPb$pulledColonies)
      age0p2 <- c(age0p2, tmpW2$remnants) # The queens of the remnant colonies are of age 0
      age2 <- c(age2, tmpW2$swarms)
    }

    cat("# Supersede a part of age1 colonies\n")
    tmpPc <- pullColonies(age1, p = p2supersede)
    age1 <- tmpPc$remainingColonies
    tmpQ1 <- supersedeColonies(tmpPc$pulledColonies)
    age0p2 <- c(age0p2, tmpQ1) # The queens of superseded colonies are of age 0

    if (year > 1) {
      cat("# Supersede a part of age2 colonies\n")
      tmpPd <- pullColonies(age2, p = p2supersede)
      age2 <- tmpPd$remainingColonies
      tmpQ2 <- supersedeColonies(tmpPd$pulledColonies)
      age0p2 <- c(age0p2, tmpQ2) # The queens of superseded colonies are of age 0
    }

    cat("# Replace all the drones\n")
    age1 <- replaceDronesColonies(age1)
    if (year > 1) {
      age2 <- replaceDronesColonies(age2)
    }

    cat("# Mate age 0 period 2 swarms and splits\n")
    if (year == 1) {
      DCA <- createDCA(age1, nInd = nDronesFull)
    } else {
      DCA <- createDCA(c(age1, age2), nInd = nDronesFull)
    }
    age0p2 <- crossColonies(age0p2, DCA = DCA, nAvgFathers = nAvgFathers)

    cat("# Collapse collonies\n")
    age1 <- selectColonies(age1, p = 1 - p2collapse)
    if (year > 1) {
      age2 <- selectColonies(age2, p = 1 - p2collapse)
    }

    cat("# Merge all age 0 colonies (from both periods)\n")
    age0 <- c(age0p1, age0p2)

    # Period3 ------------------------------------------------------------------
    cat("Period3 -----------------------------------------------------------\n")

    cat("# Collapse age0 queens\n")
    age0 <- selectColonies(age0, p = (1 - p3collapseAge0))
    age1 <- selectColonies(age1, p = (1 - p3collapseAge1))
    age2 <- NULL # We don't need this but just to show the workflow!!!

    # Get the intra-colonial csd variability -----------------------------------
    cat("Get the intra-colonial csd variability ----------------------------\n")

    if (isCsdActive()) {
      tmp <- nCsdAlleles(age0)
      # TODO: we need collapse here!
      nCsdAll <- sapply(X = tmp, FUN = function(z) z$workers)
      nCSD <- sum(nCsdAll) / nColonies(age0)
      csdVariability <- rbind(csdVariability, c(Rep, year, nCSD))
    }

    # Maintain the number of colonies ------------------------------------------
    cat("Maintain the number of colonies -----------------------------------\n")

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

  } # Year-loop
  noQueens <- rbind(noQueens, c(Rep, nColonies(age0), nColonies(age1), (nColonies(age0) + nColonies(age1))))
  a <- toc()
  loopTime <- rbind(loopTime, c(Rep, a$tic, a$toc, a$msg, (a$toc - a$tic)))

  nVirginQueens(age2)
  nVirginQueens(age1)
  nVirginQueens(age0p1)

} # Rep-loop

ggplot(noQueens, aes(Rep, sum)) +
  geom_line(aes(Rep, sum))

