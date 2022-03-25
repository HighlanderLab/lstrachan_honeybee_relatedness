setwd("~/github/lstrachan_honeybee_sim/YearCycleSimulation/")
# Clean workspace
rm(list = ls())

# Load packages
library(AlphaSimR)
library(ggplot2)
library(tictoc)
library(R6)
library(nadiv)
library(Matrix)
library(pedigreemm)
library(ggcorrplot)



# TODO: replace with devtools installation from Github once the package is operational
# Source the development version of AlphaSimR


#library(SIMplyBee)

# Parameters -------------------------------------------------------------------
nRep <- 1
nYear <- 10

apiarySize <- 10
nWorkers <- 0 # TODO: change to 20K
nFathers <- 15
nDrones <- 100 #nWorkers * 0.2
nVirginQueens <- 1

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

# Create df for recording the number of age0 and age1 colonies, csd variability and for recording cpu time
loopTime <- data.frame(Rep = NA, tic = NA, toc = NA, msg = NA, time = NA)
nQueens <- data.frame(Rep = NA, Year = NA, Age0 = NA, Age1 = NA, sum = NA)
csdVariability <- data.frame(Rep = NA, year = NA, id = NA, nCSDage0 = NA, totalCSDage0 = NA)
pDiploidDrones <- data.frame(Rep = NA, year = NA, id = NA, pQueenHomBrood_age0 = NA, pQueenHomBrood_age1 = NA)

# Rep-loop ---------------------------------------------------------------------


for (Rep in 1:nRep) {
  # Rep <- 1
  cat(paste0("Rep: ", Rep, "/", nRep, "\n"))
  # Measure cpu time
  tic(paste0(nYear, 'y loop'))
  # Start profiling
  Rprof()


  # Founder population ---------------------------------------------------------

  founderGenomes <- quickHaplo(nInd = 30,
                               nChr = 16,
                               segSites = 100)
  SP <- SimParamBee$new(founderGenomes, csdChr = 3, nCsdAlleles = 128)
  SP$nWorkers <- nWorkers
  SP$nDrones <- nDrones
  SP$nFathers <- nFathers
  SP$nVirginQueens <- nVirginQueens
  SP$pSwarm <- 0.5
  SP$pSplit <- 0.3
  SP$setTrackPed(TRUE)
  SP$setTrackRec(TRUE)
  SP$addSnpChip(nSnpPerChr = 3)

  base <- createVirginQueens(founderGenomes)

  # Year-loop ------------------------------------------------------------------


  for (year in 1:nYear) {
    # year <- 1
    # year <- year + 1
    cat(paste0("Year: ", year, "/", nYear, "\n"))
    if (year == 1) {

      age1 <- createColonies(pop = base, n = apiarySize, mated = TRUE,
                             nDronesPerQueen = 30,
                             simParamBee = SP)

    } else {
      age2 <- age1
      age1 <- age0
      age0 <- NULL
      age0p1 <- NULL
      age0p2 <- NULL
    }

    # Choose the first colony of age 1 to inspect relationship in the base population
    springerColony1 <- age1[[1]]
    # Building the colony up to 10,000 workers
    springerColony1 <- buildUpColony(colony = springerColony1,
                                     nWorkers = 10000,
                                     nDrones = 2000)
    # Extract the genotypes of all the colony members
    springerColony1_geno <- rbind(getCasteSegSiteGeno(springerColony1, caste = "queen"),
                                  getCasteSegSiteGeno(springerColony1, caste = "workers"),
                                  getCasteSegSiteGeno(springerColony1, caste = "drones"),
                                  getCasteSegSiteGeno(springerColony1, caste = "fathers"))
    # Collect the sex of all the colony members in the order than genotypes
    sex_springerColony1 <- c(rep("F", nQueens(springerColony1)),
                            rep("F", nWorkers(springerColony1)),
                            rep("M", nDrones(springerColony1)),
                            rep("M", nFathers(springerColony1)))
    # Compute the IBD relationship matrix
    ibs_springerColony1 <- calcBeeGRMIbs(x = springerColony1_geno,
                                        sex = sex_springerColony1)
    # Collect the IBD haplotypes for all the colony members
    springerColony1_haplo <- rbind(getCasteIbdHaplo(springerColony1, caste = "queen"),
                                   getCasteIbdHaplo(springerColony1, caste = "workers"),
                                   getCasteIbdHaplo(springerColony1, caste = "drones"),
                                   getCasteIbdHaplo(springerColony1, caste = "fathers"))
    # Compute the IBD relationship matrix
    ibd_springerColony1 <- calcBeeGRMIbd(x = springerColony1_haplo)
    ibd_springerColony1 <- ibd_springerColony1$indiv

    springerColony1_id <- getCasteId(springerColony1, caste = "all")

    # Period1 ------------------------------------------------------------------

    # Build-up the colonies
    age1 <- buildUpColonies(age1)
    if (year > 1) {
      age2 <- buildUpColonies(age2)
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
    # Virgin queens for splits!
    virginQueens <- createVirginQueens(age1[[virginDonor]], nInd = nColonies(age0p1))

    # Requeen the splits --> queens are now 0 years old
    age0p1 <- reQueenColonies(age0p1, queens = virginQueens)

    # Swarm a percentage of age1 colonies
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
      DCA <- createDCA(age1)
    } else {
      DCA <- createDCA(c(age1, age2))
    }
    age0p1 <- crossColonies(age0p1, drones = DCA, nFathers = SP$nFathers) #TODO: Remove this

    # Collapse
    age1 <- selectColonies(age1, p = 1 - p1collapse)
    if (year > 1) {
      age2 <- selectColonies(age2, p = 1 - p1collapse)
    }

    # Period2 ------------------------------------------------------------------

    # Swarm a percentage of age1 colonies
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
      DCA <- createDCA(age1)
    } else {
      DCA <- createDCA(c(age1, age2))
    }
    # Cross age 0 period 2 swarms and splits
    age0p2 <- crossColonies(age0p2, drones = DCA, nFathers = nFathers) #TODO: REMOVE

    # Collapse
    age1 <- selectColonies(age1, p = 1 - p2collapse)
    if (year > 1) {
      age2 <- selectColonies(age2, p = 1 - p2collapse)
    }

    # Merge all age 0 colonies (from both periods)
    age0 <- c(age0p1, age0p2)

    # Period3 ------------------------------------------------------------------

    # Collapse age0 queens
    age0 <- selectColonies(age0, p = (1 - p3collapseAge0))
    age1 <- selectColonies(age1, p = (1 - p3collapseAge1))
    age2 <- NULL #We don't need this but just to show the workflow!!!


    # Maintain the number of colonies ------------------------------------------
    totalCsdAge0 <- nCsdAlleles(age0, collapse = TRUE)
    for (n in 1:nColonies(age0)){
     csdVariability <- rbind(csdVariability, c(Rep, year, age0[[n]]@id,
                                              nCsdAlleles(age0[[n]], collapse = TRUE), totalCsdAge0))
     pDiploidDrones <- rbind(pDiploidDrones,
                            c(Rep = Rep, year = year, id = age0[[n]]@id,
                              pQueenHomBrood_age0 = computeQueensPHomBrood(age0),
                              pQueenHomBrood_age1 = computeQueensPHomBrood(age1)))
    }

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


  nQueens <- rbind(nQueens, c(Rep, year, nColonies(age0), nColonies(age1), (nColonies(age0) + nColonies(age1))))
  a <- toc()
  loopTime <- rbind(loopTime, c(Rep, a$tic, a$toc, a$msg, (a$toc - a$tic)))

  # Take the first colony of age1 and extract genotypes for relationship computation
  springerColony10 <- age1[[1]]
  springerColony10 <- buildUpColony(colony = springerColony10,
                                    nWorkers = 10000,
                                    nDrones = 2000)
  springerColony10_geno <- rbind(getCasteSegSiteGeno(springerColony10, caste = "queen"),
                                 getCasteSegSiteGeno(springerColony10, caste = "workers"),
                                 getCasteSegSiteGeno(springerColony10, caste = "drones"),
                                 getCasteSegSiteGeno(springerColony10, caste = "fathers"))

  sex_springerColony10 <- c(rep("F", nQueens(springerColony10)),
                            rep("F", nWorkers(springerColony10)),
                            rep("M", nDrones(springerColony10)),
                            rep("M", nFathers(springerColony10)))
  ibs_springerColony10 <- calcBeeGRMIbs(x = springerColony10_geno,
                                        sex = sex_springerColony10)

  springerColony10_haplo <- rbind(getCasteIbdHaplo(springerColony10, caste = "queen"),
                                  getCasteIbdHaplo(springerColony10, caste = "workers"),
                                  getCasteIbdHaplo(springerColony10, caste = "drones"),
                                  getCasteIbdHaplo(springerColony10, caste = "fathers"))

  ibd_springerColony10 <- calcBeeGRMIbd(x = springerColony10_haplo)
  ibd_springerColony10 <- ibd_springerColony10$indiv

  springerColony10_id <- getCasteId(springerColony10, caste = "all")

  # Compute the pedigree relationship matrix
  pedigree <- as.data.frame(SP$pedigree)
  # nadiv needs missing as NA
  pedigree$mother[pedigree$mother == 0] <- NA
  pedigree$father[pedigree$father == 0] <- NA
  pedigree$ID <- rownames(pedigree)

  # Females are 1, males are 0
  pedigree$isDH <-as.numeric(pedigree$isDH == 0)
  colnames(pedigree) <- c("Dam", "Sire", "Sex", "ID")
  # The order for nadiv should be ID, Dam, Sire, Sex
  pedigree <- pedigree[, c("ID", "Dam", "Sire", "Sex")]
  pedigree$Sire[pedigree$Sex == 0] <- NA

  A <- makeS(pedigree = pedigree, heterogametic = "1", returnS = TRUE)
  S <- A$S

} # Rep-loop


save.image(file = "~/github/lstrachan_honeybee_sim/YearCycleSimulation/SpringerSimulation.Rdata")

write.csv(nQueens, paste0("nQueens", Rep, ".csv"), quote = FALSE, row.names = FALSE)
write.csv(csdVariability, paste0("CsdVariability", Rep, ".csv"), quote = FALSE, row.names = FALSE)
write.csv(pDiploidDrones, paste0("pDiploidDrones", Rep, ".csv"), quote = FALSE, row.names = FALSE)
#write.csv(ped, paste0("ped", Rep, ".csv"), quote = FALSE, row.names = FALSE)






