library(dplyr)
library(tidyr)
library(ggplot2)
library(Matrix)
library(SIMplyBee)
library(gridExtra)

print("Reading in the data")
#Laura's laptop data
data <- load("~/Desktop/GitHub/lstrachan_honeybee_sim/YearCycleSimulation/SpringerSimulation_import.RData")
Sinv <- readMM("~/Desktop/GitHub/lstrachan_honeybee_sim/YearCycleSimulation/Sinv.mm")
#Eddie data
#data <- load("SpringerSimulation_import_objects.RData")
#Sinv <- readMM("Sinv.mm")
#save.image("~/Documents/")

# The data contains two populations - mellifera and carnica
# The carnica stays "pure" throughout the simulation
# Mellifera gets mated with a proportion of carnica drones

#Plotting help
# The colourblind palette with grey:
cbPalette <- c("#FF6DB6", "#490092", "#6DB6FF")


print("Assigning objects")
ped <- pedigree
caste <- caste
# Colony in year 1
# Carnica
colonyCar1 <- springerColony1_Car
ibsCar1 <- colonyCar1$IBS
ibsCcar1 <- colonyCar1$IBSOwnFreq
ibsCar1_csdChr <- colonyCar1$IBScsdChr
ibsCar1_csd <- colonyCar1$IBSCsd
ibdCar1 <- colonyCar1$IBD
ibdCar1_csdChr <- colonyCar1$IBDcsdChr
ibdCar1_csd <- colonyCar1$IBDCsd
idCar1 <- colonyCar1$ID
# Queens in year 1
queens1 <- springerQueens1
ibsQueens1 <- queens1$IBS
ibsCQueens1 <- queens1$IBSOwnFreq
ibsQueens1_csdChr <- queens1$IBScsdChr
ibsQueens1_csd <- queens1$IBSCsd
ibdQueens1 <- queens1$IBD
ibdQueens1_csdChr <- queens1$IBDcsdChr
ibdQueens1_csd <- queens1$IBDCsd
idQueens1 <- queens1$ID
idPopQueens1 <- springerQueensPop1

# Colony in year 10
# Carnica
colonyCar10 <- springerColony10_Car
ibsCar10 <- colonyCar10$IBS
ibsCcar10 <- colonyCar10$IBSOwnFreq
ibsCar10_csdChr <- colonyCar10$IBScsdChr
ibsCar10_csd <- colonyCar10$IBSCsd
ibdCar10 <- colonyCar10$IBD
ibdCar10_csdChr <- colonyCar10$IBDcsdChr
ibdCar10_csd <- colonyCar10$IBDCsd
idCar10 <- colonyCar10$ID
# Queens in year 10
queens10 <- springerQueens10
ibsQueens10 <- queens10$IBS
ibsCQueens10 <- queens10$IBSOwnFreq
ibsQueens10_csdChr <- queens10$IBScsdChr
ibsQueens10_csd <- queens10$IBSCsd
ibdQueens10 <- queens10$IBD
ibdQueens10_csdChr <- queens10$IBDcsdChr
ibdQueens10_csd <- queens10$IBDCsd
idQueens10 <- queens10$ID
idPopQueens10 <- springerQueensPop10


getS <- function(Sinv, ids, with = ids, diagOnly = FALSE, vector = FALSE) {
  ids <- as.numeric(ids)
  with <- as.numeric(with)
  x <- sparseMatrix(i = ids, j = 1:length(ids), dims = c(nrow(Sinv), length(ids)))
  M1 <- as(x, "dMatrix")
  Sids <- solve(Sinv, M1)[with,]

  if (diagOnly) {
    Sids <- diag(Sids)
  }
  if (vector) {
    Sids <- c(as.matrix(Sids))
  }
  return(Sids)
}

# Plotting functions
prepareDataForPlotting_Colony <- function(ibsDF = NULL, ibsCdf = NULL, ibdDF = NULL, Sinv = NULL,  idDF) {
  # workers to workers
  print("IBSb WW")
  tmp <- ibsDF[idDF$workers, idDF$workers]
  IBSb_WW1 <- c(tmp[lower.tri(tmp, diag = FALSE)])
  ret <- data.frame(Value = IBSb_WW1, Rel = "WW", Type = "IBS")

  print("IBSc WW")
  tmp <- ibsCdf[idDF$workers, idDF$workers]
  IBSc_WW1 <- c(tmp[lower.tri(tmp, diag = FALSE)])
  ret <- rbind(ret, data.frame(Value = IBSc_WW1, Rel = "WW", Type = "IBSOwnFreq"))

  print("IBDr WW")
  tmp <- ibdDF[idDF$workers, idDF$workers]
  IBDr_WW1 <- c(tmp[lower.tri(tmp, diag = FALSE)])
  ret <- rbind(ret, data.frame(Value = IBDr_WW1, Rel = "WW", Type = "IBDr"))

  if (!is.null(Sinv)) {
    print("IBDe WW")
    tmp <- getS(Sinv, ids = idDF$workers, vector = TRUE)
    IBDe_WW <- c(tmp[lower.tri(tmp, diag = FALSE)])
    ret <- rbind(ret, data.frame(Value = IBDe_WW, Rel = "WW", Type = "IBDe"))
  }

  # workers vs drones
  print("IBSb WD")
  IBSb_WD1 <- c(ibsDF[idDF$workers, idDF$drones])
  ret <- rbind(ret, data.frame(Value = IBSb_WD1, Rel = "WD", Type = "IBS"))

  print("IBSc WD")
  IBSc_WD1 <- c(ibsCdf[idDF$workers, idDF$drones])
  ret <- rbind(ret, data.frame(Value = IBSc_WD1, Rel = "WD", Type = "IBSOwnFreq"))

  print("IBDr WD")
  IBDr_WD1 <- c(ibdDF[idDF$workers, idDF$drones])
  ret <- rbind(ret, data.frame(Value = IBDr_WD1, Rel = "WD", Type = "IBDr"))

  if (!is.null(Sinv)) {
    print("IBDe WD")
    ret <- rbind(ret, data.frame(Value = getS(Sinv, ids = idDF$workers, with = idDF$drones, vector = TRUE),
                                 Rel = "WD", Type = "IBDe"))
  }

  # drones vs drones
  print("IBSb DD")
  tmp <- ibsDF[idDF$drones, idDF$drones]
  IBSb_DD1 <- c(tmp[lower.tri(tmp, diag = FALSE)])
  ret <- rbind(ret, data.frame(Value = IBSb_DD1, Rel = "DD", Type = "IBS"))

  print("IBSc DD")
  tmp <- ibsCdf[idDF$drones, idDF$drones]
  IBSc_DD1 <- c(tmp[lower.tri(tmp, diag = FALSE)])
  ret <- rbind(ret, data.frame(Value = IBSc_DD1, Rel = "DD", Type = "IBSOwnFreq"))

  print("IBDr DD")
  tmp <- ibdDF[idDF$drones, idDF$drones]
  IBDr_DD1 <- c(tmp[lower.tri(tmp, diag = FALSE)])
  ret <- rbind(ret, data.frame(Value = IBDr_DD1, Rel = "DD", Type = "IBDr"))

  if (!is.null(Sinv)) {
    print("IBDe DD")
    tmp <- getS(Sinv, ids = idDF$drones, vector = TRUE)
    IBDe_DD <- c(tmp[lower.tri(tmp, diag = FALSE)])
    ret <- rbind(ret, data.frame(Value = IBDe_DD, Rel = "DD", Type = "IBDe"))
  }

  # queen vs workers
  print("IBSb QW")
  ret <- rbind(ret, data.frame(Value = c(ibsDF[idDF$workers, idDF$queen]),
                               Rel = "QW", Type = "IBS"))

  print("IBSc QW")
  ret <- rbind(ret, data.frame(Value = c(ibsCdf[idDF$workers, idDF$queen]),
                               Rel = "QW", Type = "IBSOwnFreq"))

  print("IBDr QW")
  ret <- rbind(ret, data.frame(Value = c(ibdDF[idDF$workers, idDF$queen]),
                               Rel = "QW", Type = "IBDr"))

  if (!is.null(Sinv)) {
    print("IBDe QW")
    ret <- rbind(ret, data.frame(Value = getS(Sinv, ids = idDF$queen, with = idDF$workers, vector = TRUE),
                                 Rel = "QW", Type = "IBDe"))
  }

  # queen vs drones
  print("IBSb QD")
  ret <- rbind(ret, data.frame(Value = c(ibsDF[idDF$queen, idDF$drones]),
                               Rel = "QD", Type = "IBS"))

  print("IBSc QD")
  ret <- rbind(ret, data.frame(Value = c(ibsCdf[idDF$queen, idDF$drones]),
                               Rel = "QD", Type = "IBSOwnFreq"))

  print("IBDr QD")
  ret <- rbind(ret, data.frame(Value = c(ibdDF[idDF$queen, idDF$drones]),
                               Rel = "QD", Type = "IBDr"))

  if (!is.null(Sinv)) {
    print("IBDe QD")
    ret <- rbind(ret, data.frame(Value = getS(Sinv, ids = idDF$queen, with = idDF$drones, vector = TRUE),
                                 Rel = "QD", Type = "IBDe"))
  }

  return(ret)
}

prepareDataForPlotting_ColonyDiag <- function(ibsDF = NULL, ibsCdf = NULL, ibdDF = NULL, Sinv = NULL,  idDF) {
  # workers
  print("IBSb WW")
  tmp <- diag(ibsDF[idDF$workers, idDF$workers])
  ret <- data.frame(Value = tmp, Rel = "WW", Type = "IBS")

  print("IBSc WW")
  tmp <- diag(ibsCdf[idDF$workers, idDF$workers])
  ret <- rbind(ret, data.frame(Value = tmp, Rel = "WW", Type = "IBSOwnFreq"))

  print("IBDr WW")
  tmp <- diag(ibdDF[idDF$workers, idDF$workers])
  ret <- rbind(ret, data.frame(Value = tmp, Rel = "WW", Type = "IBDr"))

  if (!is.null(Sinv)) {
    print("IBDe WW")
    tmp <- getS(Sinv, ids = idDF$workers, vector = TRUE, diagOnly = TRUE)
    ret <- rbind(ret, data.frame(Value = tmp, Rel = "WW", Type = "IBDe"))
  }
  # drones
  print("IBSb DD")
  tmp <- diag(ibsDF[idDF$drones, idDF$drones])
  ret <- rbind(ret, data.frame(Value = tmp, Rel = "DD", Type = "IBS"))

  print("IBSc DD")
  tmp <- diag(ibsCdf[idDF$workers, idDF$workers])
  ret <- rbind(ret, data.frame(Value = tmp, Rel = "DD", Type = "IBSOwnFreq"))

  print("IBDr DD")
  tmp <- diag(ibdDF[idDF$drones, idDF$drones])
  ret <- rbind(ret, data.frame(Value = tmp, Rel = "DD", Type = "IBDr"))

  if (!is.null(Sinv)) {
    print("IBDe DD")
    tmp <- getS(Sinv, ids = idDF$drones, vector = TRUE, diagOnly = TRUE)
    ret <- rbind(ret, data.frame(Value = tmp, Rel = "DD", Type = "IBDe"))
  }
}

prepareDataForPlotting_ColonyCSD <- function(ibsDF = NULL, ibdDF = NULL, Sinv = NULL,  idDF) {
  # workers to workers
  print("IBSb WW")
  tmp <- ibsDF[idDF$workers, idDF$workers]
  IBSb_WW1 <- c(tmp[lower.tri(tmp, diag = FALSE)])
  ret <- data.frame(Value = IBSb_WW1, Rel = "WW", Type = "IBS")

  print("IBDr WW")
  tmp <- ibdDF[idDF$workers, idDF$workers]
  IBDr_WW1 <- c(tmp[lower.tri(tmp, diag = FALSE)])
  ret <- rbind(ret, data.frame(Value = IBDr_WW1, Rel = "WW", Type = "IBDr"))

  if (!is.null(Sinv)) {
    print("IBDe WW")
    tmp <- getS(Sinv, ids = idDF$workers, vector = TRUE)
    IBDe_WW <- c(tmp[lower.tri(tmp, diag = FALSE)])
    ret <- rbind(ret, data.frame(Value = IBDe_WW, Rel = "WW", Type = "IBDe"))
  }

  # workers vs drones
  print("IBSb WD")
  IBSb_WD1 <- c(ibsDF[idDF$workers, idDF$drones])
  ret <- rbind(ret, data.frame(Value = IBSb_WD1, Rel = "WD", Type = "IBS"))

  print("IBDr WD")
  IBDr_WD1 <- c(ibdDF[idDF$workers, idDF$drones])
  ret <- rbind(ret, data.frame(Value = IBDr_WD1, Rel = "WD", Type = "IBDr"))

  if (!is.null(Sinv)) {
    print("IBDe WD")
    ret <- rbind(ret, data.frame(Value = getS(Sinv, ids = idDF$workers, with = idDF$drones, vector = TRUE),
                                 Rel = "WD", Type = "IBDe"))
  }

  # drones vs drones
  print("IBSb DD")
  tmp <- ibsDF[idDF$drones, idDF$drones]
  IBSb_DD1 <- c(tmp[lower.tri(tmp, diag = FALSE)])
  ret <- rbind(ret, data.frame(Value = IBSb_DD1, Rel = "DD", Type = "IBS"))

  print("IBDr DD")
  tmp <- ibdDF[idDF$drones, idDF$drones]
  IBDr_DD1 <- c(tmp[lower.tri(tmp, diag = FALSE)])
  ret <- rbind(ret, data.frame(Value = IBDr_DD1, Rel = "DD", Type = "IBDr"))

  if (!is.null(Sinv)) {
    print("IBDe DD")
    tmp <- getS(Sinv, ids = idDF$drones, vector = TRUE)
    IBDe_DD <- c(tmp[lower.tri(tmp, diag = FALSE)])
    ret <- rbind(ret, data.frame(Value = IBDe_DD, Rel = "DD", Type = "IBDe"))
  }

  # queen vs workers
  print("IBSb QW")
  ret <- rbind(ret, data.frame(Value = c(ibsDF[idDF$workers, idDF$queen]),
                               Rel = "QW", Type = "IBS"))

  print("IBDr QW")
  ret <- rbind(ret, data.frame(Value = c(ibdDF[idDF$workers, idDF$queen]),
                               Rel = "QW", Type = "IBDr"))

  if (!is.null(Sinv)) {
    print("IBDe QW")
    ret <- rbind(ret, data.frame(Value = getS(Sinv, ids = idDF$queen, with = idDF$workers, vector = TRUE),
                                 Rel = "QW", Type = "IBDe"))
  }

  # queen vs drones
  print("IBSb QD")
  ret <- rbind(ret, data.frame(Value = c(ibsDF[idDF$queen, idDF$drones]),
                               Rel = "QD", Type = "IBS"))

  print("IBDr QD")
  ret <- rbind(ret, data.frame(Value = c(ibdDF[idDF$queen, idDF$drones]),
                               Rel = "QD", Type = "IBDr"))

  if (!is.null(Sinv)) {
    print("IBDe QD")
    ret <- rbind(ret, data.frame(Value = getS(Sinv, ids = idDF$queen, with = idDF$drones, vector = TRUE),
                                 Rel = "QD", Type = "IBDe"))
  }

  return(ret)
}

prepareDataForPlotting_Queens <- function(ibsDF = NULL, ibsCdf = NULL, ibdDF = NULL, Sinv = NULL,  idPopDF) {
  #Get populations IDs
  melID <- idPopDF$ID[idPopDF$Pop == "Mel"]
  melCrossID <- idPopDF$ID[idPopDF$Pop == "MelCross"]
  carID <- idPopDF$ID[idPopDF$Pop == "Car"]

  #Compute relationships between queens of different populations (QQ)
  #IBSb
  IBSbMelMelCross <- data.frame(Value = as.vector(list(ibsDF[melID, melCrossID])[[1]]),
                               Pops = "Mel_MelCross", Rel = "QQ")
  IBSbMelCar <- data.frame(Value = as.vector(list(ibsDF[melID, carID])[[1]]),
                          Pops = "Mel_Car", Rel = "QQ")
  IBSbMelCrossCar <- data.frame(Value = as.vector(list(ibsDF[melCrossID, carID])[[1]]),
                               Pops = "MelCross_Car", Rel = "QQ")
  IBSb <-rbind(IBSbMelMelCross, IBSbMelCar, IBSbMelCrossCar)
  IBSb$Type <- "IBS"

  #IBSc
  IBScMelMelCross <- data.frame(Value = as.vector(list(ibsCdf[melID, melCrossID])[[1]]),
                               Pops = "Mel_MelCross", Rel = "QQ")
  IBScMelCar <- data.frame(Value = as.vector(list(ibsCdf[melID, carID])[[1]]),
                          Pops = "Mel_Car", Rel = "QQ")
  IBScMelCrossCar <- data.frame(Value = as.vector(list(ibsCdf[melCrossID, carID])[[1]]),
                               Pops = "MelCross_Car", Rel = "QQ")
  IBSc <-rbind(IBScMelMelCross, IBScMelCar, IBScMelCrossCar)
  IBSc$Type <- "IBSOwnFreq"

  #IBDr
  IBDrMelMelCross <- data.frame(Value = as.vector(list(ibdDF[melID, melCrossID])[[1]]),
                                Pops = "Mel_MelCross", Rel = "QQ")
  IBDrMelCar <- data.frame(Value = as.vector(list(ibdDF[melID, carID])[[1]]),
                           Pops = "Mel_Car", Rel = "QQ")
  IBDrMelCrossCar <- data.frame(Value = as.vector(list(ibdDF[melCrossID, carID])[[1]]),
                                Pops = "MelCross_Car", Rel = "QQ")
  IBDr <- rbind(IBDrMelMelCross, IBDrMelCar, IBDrMelCrossCar)
  IBDr$Type <- "IBDr"

  #IBDe
  if (!is.null(Sinv)) {
    IBDeMelMelCross <- data.frame(Value = getS(Sinv, ids = melID, with = melCrossID, vector = TRUE),
                                  Pops = "Mel_MelCross", Rel = "QQ", Type = "IBDe")
    IBDeMelCar <- data.frame(Value = getS(Sinv, ids = melID, with = carID, vector = TRUE),
                             Pops = "Mel_Car", Rel = "QQ", Type = "IBDe")
    IBDeMelCrossCar <- data.frame(Value = getS(Sinv, ids = melCrossID, with = carID, vector = TRUE),
                                  Pops = "MelCross_Car", Rel = "QQ", Type = "IBDe")
    IBDe <- rbind(IBDeMelMelCross, IBDeMelCar, IBDeMelCrossCar)
  }

  # Compute relationships between queens of the same populations/ inbreeding without diagonal (Q)
  #IBSb
  tmp <- ibsDF[melID, melID]
  IBSbMel_Mel <- c(tmp[lower.tri(tmp, diag = FALSE)])
  IBSbMel_Mel <- data.frame(Value = as.vector(list(IBSbMel_Mel)[[1]]), Rel = "Q", Type = "IBS", Pops = "Mel")

  tmp <- ibsDF[carID, carID]
  IBSbCar_Car <- c(tmp[lower.tri(tmp, diag = FALSE)])
  IBSbCar_Car <-data.frame(Value = as.vector(list(IBSbCar_Car)[[1]]), Rel = "Q", Type = "IBS", Pops = "Car")

  tmp <- ibsDF[melCrossID, melCrossID]
  IBSbmelCross_melCross <-  c(tmp[lower.tri(tmp, diag = FALSE)])
  IBSbmelCross_melCross <- data.frame(Value = as.vector(list(IBSbmelCross_melCross)[[1]]), Rel = "Q", Type = "IBS", Pops = "MelCross")

  IBSb <- rbind(IBSb, IBSbMel_Mel, IBSbCar_Car, IBSbmelCross_melCross)

  #IBSc
  tmp <- ibsCdf[melID, melID]
  IBScMel_Mel <- c(tmp[lower.tri(tmp, diag = FALSE)])
  IBScMel_Mel <- data.frame(Value = as.vector(list(IBScMel_Mel)[[1]]), Rel = "Q", Type = "IBSOwnFreq", Pops = "Mel")

  tmp <- ibsCdf[carID, carID]
  IBScCar_Car <- c(tmp[lower.tri(tmp, diag = FALSE)])
  IBScCar_Car <-data.frame(Value = as.vector(list(IBScCar_Car)[[1]]), Rel = "Q", Type = "IBSOwnFreq", Pops = "Car")

  tmp <- ibsCdf[melCrossID, melCrossID]
  IBScmelCross_melCross <-  c(tmp[lower.tri(tmp, diag = FALSE)])
  IBScmelCross_melCross <- data.frame(Value = as.vector(list(IBScmelCross_melCross)[[1]]), Rel = "Q", Type = "IBSOwnFreq", Pops = "MelCross")

  IBSc <- rbind(IBSc, IBScMel_Mel, IBScCar_Car, IBScmelCross_melCross)

  #IBDr
  tmp <- ibdDF[melID, melID]
  IBDrMel_Mel <- c(tmp[lower.tri(tmp, diag = FALSE)])
  IBDrMel_Mel <- data.frame(Value = as.vector(list(IBDrMel_Mel)[[1]]), Rel = "Q", Type = "IBDr", Pops = "Mel")

  tmp <- ibdDF[carID, carID]
  IBDrCar_Car <- c(tmp[lower.tri(tmp, diag = FALSE)])
  IBDrCar_Car <- data.frame(Value = as.vector(list(IBDrCar_Car)[[1]]), Rel = "Q", Type = "IBDr", Pops = "Car")

  tmp <- ibdDF[melCrossID, melCrossID]
  IBDrmelCross_melCross <-  c(tmp[lower.tri(tmp, diag = FALSE)])
  IBDrmelCross_melCross <- data.frame(Value = as.vector(list(IBDrmelCross_melCross)[[1]]), Rel = "Q", Type = "IBDr", Pops = "MelCross")

  IBDr <- rbind(IBDr, IBDrMel_Mel, IBDrCar_Car, IBDrmelCross_melCross)

  #IBDe
  if (!is.null(Sinv)) {
    tmp1 <- getS(Sinv, melID, diagOnly = FALSE, vector = FALSE)
    tmp2 <- getS(Sinv, melCrossID, diagOnly = FALSE, vector = FALSE)
    tmp3 <- getS(Sinv, carID, diagOnly = FALSE, vector = FALSE)
    IBDe <- rbind(IBDe,
                    data.frame(Value = c(tmp1[lower.tri(tmp1, diag = FALSE)]),
                                Pops = "Mel", Rel = "Q", Type = "IBDe"),
                     data.frame(Value = c(tmp2[lower.tri(tmp2, diag = FALSE)]),
                                Pops = "MelCross", Rel = "Q", Type = "IBDe"),
                     data.frame(Value = c(tmp3[lower.tri(tmp3, diag = FALSE)]),
                                Pops = "Car", Rel = "Q", Type = "IBDe"))

  }


  # Inbreeding (diagonal!!!) if queens (F)
  inbIBSb <- rbind(data.frame(Value = diag(ibsDF[melID, melID]),
                             Pops = "Mel"),
                  data.frame(Value = diag(ibsDF[melCrossID, melCrossID]),
                             Pops = "MelCross"),
                  data.frame(Value = diag(ibsDF[carID, carID]),
                             Pops = "Car"))
  inbIBSb$Type = "IBS"
  inbIBSb$Rel = "F"
  IBSb <- rbind(IBSb, inbIBSb)

  inbIBSc <- rbind(data.frame(Value = diag(ibsCdf[melID, melID]),
                             Pops = "Mel"),
                  data.frame(Value = diag(ibsCdf[melCrossID, melCrossID]),
                             Pops = "MelCross"),
                  data.frame(Value = diag(ibsCdf[carID, carID]),
                             Pops = "Car"))
  inbIBSc$Type = "IBSOwnFreq"
  inbIBSc$Rel = "F"
  IBSc <- rbind(IBSc, inbIBSc)

  inbIBDr <- rbind(data.frame(Value = diag(ibdDF[melID, melID]),
                              Pops = "Mel"),
                   data.frame(Value = diag(ibdDF[melCrossID, melCrossID]),
                              Pops = "MelCross"),
                   data.frame(Value = diag(ibdDF[carID, carID]),
                              Pops = "Car"))
  inbIBDr$Type = "IBDr"
  inbIBDr$Rel = "F"
  IBDr <- rbind(IBDr, inbIBDr)

  ret <- rbind(IBSb, IBSc, IBDr)

  if (!is.null(Sinv)) { #TODO: CHANGE THIS WHEN YOU GET PEDIGREE ESTIMATES
    inbIBDe <- rbind(data.frame(Value = getS(Sinv, melID, diagOnly = TRUE),
                                Pops = "Mel"),
                     data.frame(Value = getS(Sinv, melCrossID, diagOnly = TRUE),
                                Pops = "MelCross"),
                     data.frame(Value = getS(Sinv, carID, diagOnly = TRUE),
                                Pops = "Car"))
    inbIBDe$Type = "IBDe"
    inbIBDe$Rel = "F"
    IBDe <- rbind(IBDe, inbIBDe)

    ret <- rbind(ret, IBDe)
  }

  return(ret)
}

prepareDataForPlotting_QueensCSD <- function(ibsDF = NULL, ibdDF = NULL, Sinv = NULL,  idPopDF) {
  #Get populations IDs
  melID <- idPopDF$ID[idPopDF$Pop == "Mel"]
  melCrossID <- idPopDF$ID[idPopDF$Pop == "MelCross"]
  carID <- idPopDF$ID[idPopDF$Pop == "Car"]

  #Compute relationships between queens of different populations (QQ)
  #IBSb
  IBSbMelMelCross <- data.frame(Value = as.vector(list(ibsDF[melID, melCrossID])[[1]]),
                                Pops = "Mel_MelCross", Rel = "QQ")
  IBSbMelCar <- data.frame(Value = as.vector(list(ibsDF[melID, carID])[[1]]),
                           Pops = "Mel_Car", Rel = "QQ")
  IBSbMelCrossCar <- data.frame(Value = as.vector(list(ibsDF[melCrossID, carID])[[1]]),
                                Pops = "MelCross_Car", Rel = "QQ")
  IBSb <-rbind(IBSbMelMelCross, IBSbMelCar, IBSbMelCrossCar)
  IBSb$Type <- "IBS"

  #IBSc

  #IBDr
  IBDrMelMelCross <- data.frame(Value = as.vector(list(ibdDF[melID, melCrossID])[[1]]),
                                Pops = "Mel_MelCross", Rel = "QQ")
  IBDrMelCar <- data.frame(Value = as.vector(list(ibdDF[melID, carID])[[1]]),
                           Pops = "Mel_Car", Rel = "QQ")
  IBDrMelCrossCar <- data.frame(Value = as.vector(list(ibdDF[melCrossID, carID])[[1]]),
                                Pops = "MelCross_Car", Rel = "QQ")
  IBDr <- rbind(IBDrMelMelCross, IBDrMelCar, IBDrMelCrossCar)
  IBDr$Type <- "IBDr"

  #IBDe
  if (!is.null(Sinv)) {
    IBDeMelMelCross <- data.frame(Value = getS(Sinv, ids = melID, with = melCrossID, vector = TRUE),
                                  Pops = "Mel_MelCross", Rel = "QQ", Type = "IBDe")
    IBDeMelCar <- data.frame(Value = getS(Sinv, ids = melID, with = carID, vector = TRUE),
                             Pops = "Mel_Car", Rel = "QQ", Type = "IBDe")
    IBDeMelCrossCar <- data.frame(Value = getS(Sinv, ids = melCrossID, with = carID, vector = TRUE),
                                  Pops = "MelCross_Car", Rel = "QQ", Type = "IBDe")
    IBDe <- rbind(IBDeMelMelCross, IBDeMelCar, IBDeMelCrossCar)
  }

  # Compute relationships between queens of the same populations/ inbreeding without diagonal (Q)
  #IBSb
  tmp <- ibsDF[melID, melID]
  IBSbMel_Mel <- c(tmp[lower.tri(tmp, diag = FALSE)])
  IBSbMel_Mel <- data.frame(Value = as.vector(list(IBSbMel_Mel)[[1]]), Rel = "Q", Type = "IBS", Pops = "Mel")

  tmp <- ibsDF[carID, carID]
  IBSbCar_Car <- c(tmp[lower.tri(tmp, diag = FALSE)])
  IBSbCar_Car <-data.frame(Value = as.vector(list(IBSbCar_Car)[[1]]), Rel = "Q", Type = "IBS", Pops = "Car")

  tmp <- ibsDF[melCrossID, melCrossID]
  IBSbmelCross_melCross <-  c(tmp[lower.tri(tmp, diag = FALSE)])
  IBSbmelCross_melCross <- data.frame(Value = as.vector(list(IBSbmelCross_melCross)[[1]]), Rel = "Q", Type = "IBS", Pops = "MelCross")

  IBSb <- rbind(IBSb, IBSbMel_Mel, IBSbCar_Car, IBSbmelCross_melCross)

  #IBDr
  tmp <- ibdDF[melID, melID]
  IBDrMel_Mel <- c(tmp[lower.tri(tmp, diag = FALSE)])
  IBDrMel_Mel <- data.frame(Value = as.vector(list(IBDrMel_Mel)[[1]]), Rel = "Q", Type = "IBDr", Pops = "Mel")

  tmp <- ibdDF[carID, carID]
  IBDrCar_Car <- c(tmp[lower.tri(tmp, diag = FALSE)])
  IBDrCar_Car <- data.frame(Value = as.vector(list(IBDrCar_Car)[[1]]), Rel = "Q", Type = "IBDr", Pops = "Car")

  tmp <- ibdDF[melCrossID, melCrossID]
  IBDrmelCross_melCross <-  c(tmp[lower.tri(tmp, diag = FALSE)])
  IBDrmelCross_melCross <- data.frame(Value = as.vector(list(IBDrmelCross_melCross)[[1]]), Rel = "Q", Type = "IBDr", Pops = "MelCross")

  IBDr <- rbind(IBDr, IBDrMel_Mel, IBDrCar_Car, IBDrmelCross_melCross)

  #IBDe
  if (!is.null(Sinv)) {
    tmp1 <- getS(Sinv, melID, diagOnly = FALSE, vector = FALSE)
    tmp2 <- getS(Sinv, melCrossID, diagOnly = FALSE, vector = FALSE)
    tmp3 <- getS(Sinv, carID, diagOnly = FALSE, vector = FALSE)
    IBDe <- rbind(IBDe,
                  data.frame(Value = c(tmp1[lower.tri(tmp1, diag = FALSE)]),
                             Pops = "Mel", Rel = "Q", Type = "IBDe"),
                  data.frame(Value = c(tmp2[lower.tri(tmp2, diag = FALSE)]),
                             Pops = "MelCross", Rel = "Q", Type = "IBDe"),
                  data.frame(Value = c(tmp3[lower.tri(tmp3, diag = FALSE)]),
                             Pops = "Car", Rel = "Q", Type = "IBDe"))

  }


  # Inbreeding (diagonal!!!) if queens (F)
  inbIBSb <- rbind(data.frame(Value = diag(ibsDF[melID, melID]),
                              Pops = "Mel"),
                   data.frame(Value = diag(ibsDF[melCrossID, melCrossID]),
                              Pops = "MelCross"),
                   data.frame(Value = diag(ibsDF[carID, carID]),
                              Pops = "Car"))
  inbIBSb$Type = "IBS"
  inbIBSb$Rel = "F"
  IBSb <- rbind(IBSb, inbIBSb)

  inbIBDr <- rbind(data.frame(Value = diag(ibdDF[melID, melID]),
                              Pops = "Mel"),
                   data.frame(Value = diag(ibdDF[melCrossID, melCrossID]),
                              Pops = "MelCross"),
                   data.frame(Value = diag(ibdDF[carID, carID]),
                              Pops = "Car"))
  inbIBDr$Type = "IBDr"
  inbIBDr$Rel = "F"
  IBDr <- rbind(IBDr, inbIBDr)

  ret <- rbind(IBSb, IBDr)

  if (!is.null(Sinv)) { #TODO: CHANGE THIS WHEN YOU GET PEDIGREE ESTIMATES
    inbIBDe <- rbind(data.frame(Value = getS(Sinv, melID, diagOnly = TRUE),
                                Pops = "Mel"),
                     data.frame(Value = getS(Sinv, melCrossID, diagOnly = TRUE),
                                Pops = "MelCross"),
                     data.frame(Value = getS(Sinv, carID, diagOnly = TRUE),
                                Pops = "Car"))
    inbIBDe$Type = "IBDe"
    inbIBDe$Rel = "F"
    IBDe <- rbind(IBDe, inbIBDe)

    ret <- rbind(ret, IBDe)
  }

  return(ret)
}

prepareDataForPlottingHeatMap_Queens <- function(ibsDF = NULL, ibsCdf = NULL, ibdDF = NULL, Sinv = NULL, idDF = NULL) {
  ibsDF <- as.data.frame(ibsDF)
  columns <- colnames(ibsDF)
  ibsDF$ID <- rownames(ibsDF)
  ibsDFL <- ibsDF %>% pivot_longer(cols = all_of(columns))
  ibsDFL$Method <- "IBS"

  ibsCdf <- as.data.frame(ibsCdf)
  columns <- colnames(ibsCdf)
  ibsCdf$ID <- rownames(ibsCdf)
  ibsCdfL <- ibsCdf %>% pivot_longer(cols = all_of(columns))
  ibsCdfL$Method <- "IBSOwnFreq"

  ibdrDF <- as.data.frame(ibdDF)
  columns <- colnames(ibdrDF)
  ibdrDF$ID <- rownames(ibdrDF)
  ibdrDFL <- ibdrDF %>% pivot_longer(cols = all_of(columns))
  ibdrDFL$Method <- "IBDr"
  ret <- rbind(ibsDFL, ibsCdfL, ibdrDFL)

  if (!is.null(Sinv)) {
    ibdeDF <-  as.data.frame(as.matrix(getS(Sinv, ids = idDF)))
    rownames(ibdeDF) <- idDF
    colnames(ibdeDF) <- idDF
    columns <- colnames(ibdeDF)
    ibdeDF$ID <- as.character(idDF)
    ibdeDFL <- ibdeDF %>% pivot_longer(cols = all_of(columns))
    ibdeDFL$Method <- "IBDe"
    ret <- rbind(ret, ibdeDFL)
  }
  return(ret)
}

prepareDataForPlottingHeatMap_QueensCSD <- function(ibsDF = NULL, ibdDF = NULL, Sinv = NULL, idDF = NULL) {
  ibsDF <- as.data.frame(ibsDF)
  columns <- colnames(ibsDF)
  ibsDF$ID <- rownames(ibsDF)
  ibsDFL <- ibsDF %>% pivot_longer(cols = all_of(columns))
  ibsDFL$Method <- "IBS"

  ibdrDF <- as.data.frame(ibdDF)
  columns <- colnames(ibdrDF)
  ibdrDF$ID <- rownames(ibdrDF)
  ibdrDFL <- ibdrDF %>% pivot_longer(cols = all_of(columns))
  ibdrDFL$Method <- "IBDr"
  ret <- rbind(ibsDFL, ibdrDFL)

  if (!is.null(Sinv)) {
    ibdeDF <-  as.data.frame(as.matrix(getS(Sinv, ids = idDF)))
    rownames(ibdeDF) <- idDF
    colnames(ibdeDF) <- idDF
    columns <- colnames(ibdeDF)
    ibdeDF$ID <- as.character(idDF)
    ibdeDFL <- ibdeDF %>% pivot_longer(cols = all_of(columns))
    ibdeDFL$Method <- "IBDe"
    ret <- rbind(ret, ibdeDFL)
  }
  return(ret)
}

plotColony <- function(df, rel = c("WD", "WW", "DD"), type = c("IBDr", "IBDe"), x_axis = NULL, show.legend = NULL, strip.text = NULL) {
  df$Type <- factor(df$Type, levels = c("IBDe", "IBDr", "IBS", "IBSOwnFreq"))
  df$Rel <- factor(df$Rel, levels= c("WW", "WD", "DD", "QW", "QD"))
  p <- ggplot(df[df$Rel %in% rel & df$Type %in% type, ],
                 aes(x = Value, fill = Rel)) +
    geom_histogram(binwidth = 0.01, position = "identity", show.legend = show.legend) + facet_grid(rows = vars(Type), scales = "free") + xlim(x_axis)
  plot <- p + scale_fill_manual("", values=cbPalette, aesthetics = c("colour","fill")) + theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), strip.text = strip.text)
  return(plot)
}

plotColonyF <- function(df, rel = c("WW", "DD"), type = c("IBDr", "IBDe"), x_axis = NULL, show.legend = NULL){
  df$Type <- factor(df$Type, levels = c("IBDe", "IBDr", "IBS", "IBSOwnFreq"))
  df$Rel <- factor(df$Rel, levels= c("WW", "WD", "DD", "QW", "QD"))
  p <- ggplot(df[df$Rel %in% rel & df$Type %in% type, ],
              aes(x = Value, fill = Rel)) +
    geom_histogram(binwidth = 0.01, position = "identity", show.legend = show.legend) + facet_wrap(.~Type, scales = "free_y") + xlim(x_axis)
  plot <- p + scale_fill_manual("", values=cbPalette, aesthetics = c("colour","fill")) + theme_classic()
  return(plot)
}

plotQueensQQ <- function(df, rel = c("QQ"), type = c("IBDr", "IBDe", "IBS"), x_axis = NULL, show.legend = NULL, strip.text = NULL) {
  df$Type <- factor(df$Type, levels = c("IBDe", "IBDr", "IBS", "IBSOwnFreq"))
  df$Pops <- factor(df$Pops, levels = c("Mel_MelCross", "Mel_Car", "MelCross_Car"))
  p <- ggplot(df[df$Rel %in% rel & df$Type %in% type, ],
                   aes(x = Value, fill = Pops)) +
      geom_histogram(binwidth = 0.01, position = "identity", show.legend = show.legend) +
    facet_wrap(.~Type, scales = "free_y") + xlim(x_axis) #IBS not to scale so increased scale
  plot <- p + scale_fill_manual("", values=cbPalette, aesthetics = c("colour","fill")) + theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), strip.text = strip.text)
  return(plot)
}

scatterQueens <-  function(df, rel = c("QQ"), type = c("IBDr", "IBS"), pops = c("Mel_MelCross", "Mel_Car", "MelCross_Car"), x_axis = NULL, show.legend = NULL ) {
  df$Type <- factor(df$Type, levels = c("IBDe", "IBDr", "IBS", "IBSOwnFreq"))
  df$Pops <- factor(df$Pops, levels = c("Mel_MelCross", "Mel_Car", "MelCross_Car", "Car", "Mel", "MelCross"))
  a <- pivot_wider(df[df$Rel %in% rel & df$Type %in% type & df$Pops %in% pops, ], names_from = Type, values_from = Value, values_fn = list)
  b <- unnest(a, cols = all_of(type))
  c <- ggplot(data = b, aes(x = IBDr, y = IBS)) + geom_point(aes(colour = Pops), show.legend = show.legend) + xlim(x_axis)
  plot <- c + scale_fill_manual("", values= cbPalette, aesthetics = c("colour","fill")) + theme_classic()
  return(plot)
}

plotQueensQ <- function(df, rel = c("Q", "F"), type = c("IBDr", "IBDe", "IBS"), x_axis = NULL, show.legend = NULL, strip.text = NULL) {
  df$Type <- factor(df$Type, levels = c("IBDe", "IBDr", "IBS", "IBSOwnFreq"))
  df$Rel <- factor(df$Rel, levels= c("F", "Q"))
  p <- ggplot(df[df$Rel %in% rel & df$Type %in% type, ],
                   aes(x = Value, fill = Pops)) +
      geom_histogram(binwidth = 0.01, position = "identity", show.legend = show.legend) +
      facet_wrap(.~Type, scales = "free_y") + xlim(x_axis) #xlim removed
  plot <- p + scale_fill_manual("", values=cbPalette, aesthetics = c("colour","fill")) + theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), strip.text = strip.text)
  return(plot)
}

plotQueensF <- function(df, rel = c("Q", "F"), type = c("IBDr", "IBDe", "IBS"), x_axis = NULL, show.legend = NULL) {
  df$Type <- factor(df$Type, levels = c("IBDe", "IBDr", "IBS", "IBSOwnFreq"))
  df$Rel <- factor(df$Rel, levels= c("F", "Q"))
  p <- ggplot(df[df$Rel %in% rel & df$Type %in% type, ],
              aes(x = Value - 1, fill = Pops)) +
    geom_histogram(binwidth = 0.01, position = "identity", show.legend = show.legend) +
    facet_wrap(.~Type, scales = "free_y")+ xlim(x_axis) #xlim removed
  plot <- p + scale_fill_manual("", values=cbPalette, aesthetics = c("colour","fill")) + theme_classic()
  return(plot)
}

plotQueens_heatmap <- function(df, Pop = FALSE, PopIdDF = NULL, show.legend = NULL) {
  df$ID <- as.factor(as.numeric(df$ID))
  df$name <- as.factor(as.numeric(df$name))
  df$Method <- factor(df$Method, levels = c("IBDe", "IBDr", "IBS", "IBSOwnFreq"))

  if (Pop) {
    df <- merge(df, PopIdDF, by = "ID")
    df <- merge(df, PopIdDF, by.x = "name", by.y = "ID")

    df$PopId1 <- paste0(df$Pop.x, df$ID)
    df$PopId2 <- paste0(df$Pop.y, df$name)

    n <- length(unique(df$name[df$Pop.y == "Mel"]))/2
    breaks = list(df %>% group_by(Pop.y) %>% summarise(mean = unique(name)[n]) %>% unite(PopMean, Pop.y, mean, sep=""))[[1]]$PopMean
    plot <- ggplot(data = df, aes(x = PopId1, y = PopId2, fill = value)) + geom_tile(show.legend = show.legend) + scale_fill_gradient(low = "white", high = "blue") + theme_classic() +
      scale_x_discrete(breaks = breaks, labels = c("Car", "Mel", "MelCross")) +
      scale_y_discrete(breaks = breaks, labels = c("Car", "Mel", "MelCross")) +
      xlab("") + ylab("") +
      facet_grid(rows = vars(Method))
  } else {
    plot <- ggplot(data = df, aes(x = ID, y = name, fill = value)) + geom_tile() +
      facet_grid(rows = vars(Method))
  }
  return(plot)
}

########################################################
### --- FIGURE 1&2: Pure subspecies (carnica) in years 1 and 10 ---###
########################################################
print("Plot carnica year 1")
relCar1 <- prepareDataForPlotting_Colony(ibsDF = ibsCar1, ibsCdf = ibsCcar1, ibdDF = ibdCar1, Sinv = Sinv, idDF = idCar1)
print("Plot carnica year 10")
relCar10 <- prepareDataForPlotting_Colony(ibsDF = ibsCar10, ibsCdf = ibsCcar10, ibdDF = ibdCar10, Sinv = Sinv, idDF = idCar10)

#Carnica WW/WD/DD
plotCar1 <- plotColony(relCar1, type = c("IBDe", "IBDr", "IBS", "IBSOwnFreq"), rel = c("WW", "WD", "DD"), x_axis = c(-0.55, 2.0), show.legend = FALSE, strip.text = element_blank())
plotCar1
#plotCar10 WW/WD/DD
plotCar10 <- plotColony(relCar10, type = c("IBDe", "IBDr","IBS","IBSOwnFreq"), rel = c("WW", "WD", "DD"), x_axis = c(-0.55, 2.0), show.legend = TRUE, strip.text = element_text())
plotCar10

plotCarWW <- grid.arrange(plotCar1, plotCar10, ncol = 2)

#Carnica QW/QD
plotCar1Q <- plotColony(relCar1, type = c("IBDe", "IBDr", "IBS", "IBSOwnFreq"), rel = c("QW", "QD"), x_axis = c(-0.5, 1.5), show.legend = FALSE, strip.text = element_blank() )
plotCar1Q
#plotCar10 QW/QD
plotCar10Q <- plotColony(relCar10, type = c("IBDe", "IBDr", "IBS"), rel = c("QW", "QD"), x_axis = c(-0.5, 1.5), show.legend = TRUE, strip.text = element_text())
plotCar10Q

plotCarQW <- grid.arrange(plotCar1Q, plotCar10Q, ncol = 2)

#Plot Car inbreeding
relCar1F <- prepareDataForPlotting_ColonyDiag(ibsDF = ibsCar1, ibsCdf = ibsCcar1, ibdDF = ibdCar1, Sinv = Sinv, idDF = idCar1)
relCar10F <- prepareDataForPlotting_ColonyDiag(ibsDF = ibsCar10, ibsCdf = ibsCcar10, ibdDF = ibdCar10, Sinv = Sinv, idDF = idCar10)

plotCar1F <- plotColonyF(relCar1F, rel = c("WW", "DD"), type = c("IBDe", "IBDr", "IBS", "IBSOwnFreq"), x_axis = c(0.0, 2.0), show.legend = FALSE)
plotCar1F
plotCar10F <- plotColonyF(relCar10F, rel = c("WW", "DD"), type = c("IBDe", "IBDr", "IBS", "IBSOwnFreq"), x_axis = c(0.0, 2.0), show.legend = TRUE)
plotCar10F

plotCarF <- grid.arrange(plotCar1F, plotCar10F, ncol = 2)

#Plot CAR csd Locus
#Year 1
relCar1_csd <- prepareDataForPlotting_ColonyCSD(ibsDF = ibsCar1_csd, ibdDF = ibdCar1_csd, Sinv = Sinv, idDF = idCar1)
relCar10_csd <- prepareDataForPlotting_ColonyCSD(ibsDF = ibsCar10_csd, ibdDF = ibdCar10_csd, Sinv = Sinv, idDF = idCar10)

plotCar1_csdLoc <- plotColony(relCar1_csd, type = c("IBDe", "IBDr", "IBS"), x_axis = c(-1.25, 2.0), show.legend = FALSE, strip.text = element_blank())
plotCar1_csdLoc
plotCar10_csdLoc <- plotColony(relCar10_csd, type = c("IBDe", "IBDr", "IBS"), x_axis = c(-1.25, 2.0), show.legend = TRUE, strip.text = element_text())
plotCar10_csdLoc

plotCarCsdLoc <- grid.arrange(plotCar1_csdLoc, plotCar10_csdLoc, ncol = 2)

#Plot CAR csd Chromosome
relCar1_csdChr <- prepareDataForPlotting_ColonyCSD(ibsDF = ibsCar1_csdChr, ibdDF = ibdCar1_csdChr, Sinv = Sinv, idDF = idCar1)
relCar10_csdChr <- prepareDataForPlotting_ColonyCSD(ibsDF = ibsCar10_csdChr, ibdDF = ibdCar10_csdChr, Sinv = Sinv, idDF = idCar10)


plotCar1_csdChr <- plotColony(relCar1_csdChr, type = c("IBDe", "IBDr", "IBS"), x_axis = c(-0.01, 2.0), show.legend = FALSE, strip.text = element_blank())
plotCar1_csdChr
plotCar10_csdChr <- plotColony(relCar10_csdChr, type = c("IBDe", "IBDr", "IBS"), x_axis = c(-0.01, 2.0), show.legend = TRUE, strip.text = element_text())
plotCar10_csdChr

plotCarCsdChr <- grid.arrange(plotCar1_csdChr, plotCar10_csdChr, ncol = 2)

########################################################
### --- FIGURE 4: Between queens of different subspecies (carnica vs. mellifera) ---###
########################################################
#Plot queens Year 1
print("Plot queens")
relQueens1 <- prepareDataForPlotting_Queens(ibsDF = ibsQueens1, ibsCdf = ibsCQueens1, ibdDF = ibdQueens1, Sinv = Sinv, idPopDF = idPopQueens1)
relQueens10 <- prepareDataForPlotting_Queens(ibsDF = ibsQueens10, ibsCdf = ibsCQueens10, ibdDF = ibdQueens10, Sinv = Sinv, idPopDF = idPopQueens10)

print("Between populations ")
#between populations
plotQueens1QQ <- plotQueensQQ(relQueens1, type = c("IBDr", "IBDe", "IBS"), x_axis = c(-1.0, 1.5), show.legend = FALSE, strip.text = element_text())
plotQueens1QQ
plotQueens10QQ <- plotQueensQQ(relQueens10, type = c("IBDe", "IBDr", "IBS"), x_axis = c(-1.0, 1.5), show.legend = TRUE, strip.text = element_text())
plotQueens10QQ

plotQueensQQ <- grid.arrange(plotQueens1QQ, plotQueens10QQ, nrow = 2)

#between populations as a scatter plot
plotQueensScatter1 <- scatterQueens(relQueens1, type = c("IBDr", "IBS"), pops = c("Mel_MelCross", "Mel_Car", "MelCross_Car"), x_axis = c(-0.01, 1), show.legend = FALSE)
plotQueensScatter1
plotQueensScatter10 <- scatterQueens(relQueens10, type = c("IBDr", "IBS"), pops = c("Mel_MelCross", "Mel_Car", "MelCross_Car"), x_axis = c(-0.01, 1), show.legend = TRUE)
plotQueensScatter10

plotQueensScatter <- grid.arrange(plotQueensScatter1, plotQueensScatter10, ncol = 2)

print("Within populations")
#within population - non-diagonal
plotQueens1Q <- plotQueensQ(relQueens1, rel = c("Q"), type = c("IBDr", "IBDe", "IBS"), x_axis = c(-0.5, 2.0), show.legend = FALSE, strip.text = element_text())
plotQueens1Q
plotQueens10Q <- plotQueensQ(relQueens10, rel = c("Q"), type = c("IBDr", "IBDe", "IBS"), x_axis = c(-0.5, 2.0), show.legend = TRUE, strip.text = element_text())
plotQueens10Q

plotQueensQ <- grid.arrange(plotQueens1Q, plotQueens10Q, nrow = 2)

#within population - diagonal
plotQueens1F <- plotQueensF(relQueens1, rel = c("F"), type = c("IBDr", "IBDe", "IBS"), x_axis = c(-0.1, 0.8), show.legend = FALSE)
plotQueens1F
plotQueens10F <- plotQueensF(relQueens10, rel = c("F"), type = c("IBDr", "IBDe", "IBS"), x_axis = c(-0.1, 0.8), show.legend = TRUE)
plotQueens10F

plotQueensF <- grid.arrange(plotQueens1F, plotQueens10F, nrow = 2)

print("HeatMaps")
relQueens1h <- prepareDataForPlottingHeatMap_Queens(ibsDF = ibsQueens1, ibsCdf = ibsCQueens1, ibdDF = ibdQueens1, Sinv = Sinv, idDF = idQueens1)
relQueens10h <- prepareDataForPlottingHeatMap_Queens(ibsDF = ibsQueens10, ibsCdf = ibsCQueens10, ibdDF = ibdQueens10, Sinv = Sinv, idDF = idQueens10)

plotQueens1h <- plotQueens_heatmap(relQueens1h, Pop = TRUE, PopIdDF = idPopQueens1, show.legend = FALSE)
plotQueens1h
plotQueens10h <- plotQueens_heatmap(relQueens10h, Pop = T, PopIdDF = idPopQueens10, show.legend = TRUE)
plotQueens10h

plotQueensHeatMap <- grid.arrange(plotQueens1h, plotQueens10h, ncol = 2)

#CSD STUFF
#Plot csd Loc as a histogram
relQueens1_csd_hist <- prepareDataForPlotting_QueensCSD(ibsDF = ibsQueens1_csd, ibdDF = ibdQueens1_csd, Sinv = Sinv, idPopDF = idPopQueens1)
relQueens10_csd_hist <- prepareDataForPlotting_QueensCSD(ibsDF = ibsQueens10_csd, ibdDF = ibdQueens10_csd, Sinv = Sinv, idPopDF = idPopQueens10)


plotQueens1_csdLoc<- plotQueensQQ(relQueens1_csd_hist, type = c("IBDr", "IBDe", "IBS"), x_axis = c(-1.5, 1.5), show.legend = FALSE, strip.text = element_text())
plotQueens1_csdLoc
plotQueens10_csdLoc <- plotQueensQQ(relQueens10_csd_hist, type = c("IBDr", "IBDe", "IBS"), x_axis = c(-1.5, 1.5), show.legend = TRUE, strip.text = element_text())
plotQueens10_csdLoc

plotQueensCsdLoc <- grid.arrange(plotQueens1_csdLoc, plotQueens10_csdLoc, nrow = 2)

#Plot csd loc as a heatmap
relQueens1_csd <- prepareDataForPlottingHeatMap_QueensCSD(ibsDF = ibsQueens1_csd, ibdDF = ibdQueens1_csd, Sinv = Sinv, idDF = idQueens1)
relQueens10_csd <- prepareDataForPlottingHeatMap_QueensCSD(ibsDF = ibsQueens10_csd, ibdDF = ibdQueens10_csd, Sinv = Sinv, idDF = idQueens10)


plotQueens1_csdLoc_heat <- plotQueens_heatmap(relQueens1_csd, Pop = TRUE, PopIdDF = idPopQueens1, show.legend = FALSE)
plotQueens1_csdLoc_heat
plotQueens10_csdLoc_heat <- plotQueens_heatmap(relQueens10_csd, Pop = TRUE, PopIdDF = idPopQueens10, show.legend = TRUE)
plotQueens10_csdLoc_heat

plotQueensCsdLocHeatMap <- grid.arrange(plotQueens1_csdLoc_heat, plotQueens10_csdLoc_heat, ncol = 2)

#Plot csd Chr as heatmap
relQueens1_csdChr <- prepareDataForPlottingHeatMap_QueensCSD(ibsDF = ibsQueens1_csdChr, ibdDF = ibdQueens1_csdChr, Sinv = Sinv, idDF = idQueens1)
relQueens10_csdChr <- prepareDataForPlottingHeatMap_QueensCSD(ibsDF = ibsQueens10_csdChr, ibdDF = ibdQueens10_csdChr, Sinv = Sinv, idDF = idQueens10)

plotQueens1_csdChr_heat <- plotQueens_heatmap(relQueens1_csdChr, Pop = TRUE, PopIdDF = idPopQueens1, show.legend = FALSE)
plotQueens1_csdChr_heat
plotQueens10_csdChr_heat <-plotQueens_heatmap(relQueens10_csdChr, Pop = TRUE, PopIdDF = idPopQueens10, show.legend = TRUE)
plotQueens10_csdChr_heat

plotQueensCSDChrHeapMap <- grid.arrange(plotQueens1_csdChr_heat, plotQueens10_csdChr_heat, ncol = 2)


### Csd
csdVariability <- csdVariability[!is.na(csdVariability$subspecies),]
pDiploidDrones <- pDiploidDrones[-1,]

csdVariability$nCSD <- as.numeric(csdVariability$nCSD)
csdVariability$totalCSD <- as.numeric(csdVariability$totalCSD)
csdVariability$year <- as.numeric(csdVariability$year)
csdMean <- csdVariability %>%  group_by(subspecies, year) %>%  summarize(meanCSD = mean(nCSD))
csdMean$subspecies <- factor(csdMean$subspecies, levels = c("Mel", "MelCross", "Car"))
csdVariability$subspecies <- factor(csdVariability$subspecies, levels = c("Mel", "MelCross", "Car"))
ggplot(data = csdMean, aes(x = year, y = meanCSD, colour = subspecies)) + geom_line()

ggplot(data = csdVariability, aes(x = year, y = totalCSD, colour = subspecies)) + geom_line()



pDiploidDrones$pQueenHomBrood <- as.numeric(pDiploidDrones$pQueenHomBrood)
pDiploidDrones$year <- as.numeric(pDiploidDrones$year)
pDiploidMean <- pDiploidDrones %>%  group_by(subspecies, year) %>%  summarize(meanHom = mean(pQueenHomBrood))
pDiploidMean$subspecies <- factor(pDiploidMean$subspecies, levels = c("Mel", "MelCross", "Car"))

ggplot(data = pDiploidMean, aes(x = year, y = meanHom, colour = subspecies)) + geom_line()


