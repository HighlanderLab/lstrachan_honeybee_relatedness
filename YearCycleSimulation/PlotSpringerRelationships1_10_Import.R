rm(list = ls())
setwd("~/Desktop/GitHub/lstrachan_honeybee_sim/YearCycleSimulation")

library(dplyr)
library(tidyr)
library(ggplot2)
library(Matrix)
library(SIMplyBee)
library(gridExtra)
library(ggh4x)
library(viridisLite)


print("Reading in the data")
#Laura's laptop data
data <- load("~/Desktop/GitHub/lstrachan_honeybee_sim/YearCycleSimulation/SpringerSimulation_import.RData")
Sinv <- readMM("~/Desktop/GitHub/lstrachan_honeybee_sim/YearCycleSimulation/Sinv.mm")
#Eddie data
#data <- load("SpringerSimulation_import_objects.RData")
#Sinv <- readMM("Sinv.mm")

# The data contains two populations - mellifera and carnica
# The carnica stays "pure" throughout the simulation
# Mellifera gets mated with a proportion of carnica drones

#Plotting help
# The colourblind palette with grey:
cbPaletteQ <- c("#FF6DB6", "#490092", "#6DB6FF")  #Use within the same population
cbPaletteQQ <- c("#993F00", "#FF8E32", "#B6DBFF")  #Use on different populations 
paletteViridis <- plasma(n = 3, begin = 0.3, end = 0.9) #Use on Carnica population

print("Assigning objects")
ped <- pedigree
caste <- caste

#Assign Car colony and queens 
colonyCar1 <- springerColony1_Car
colonyCar10 <- springerColony10_Car
queens1 <- springerQueens1
queens10 <- springerQueens10


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
prepareDataForPlotting_Colony <- function(ibsMultiDF = NULL, ibsSingleDF = NULL, ibsOwnDF = NULL, ibdDF = NULL, Sinv = NULL,  idDF) {
  # workers to workers
  print("IBS WW")
  tmp <- ibsMultiDF[idDF$workers, idDF$workers]
  IBS_WW1 <- c(tmp[lower.tri(tmp, diag = FALSE)])
  ret <- data.frame(Value = IBS_WW1, Rel = "WW", Type = "IBSmultiBF")
  
  print("IBSc WW")
  tmp <- ibsOwnDF[idDF$workers, idDF$workers]
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
  print("IBS WD")
  IBS_WD1 <- c(ibsMultiDF[idDF$workers, idDF$drones])
  ret <- rbind(ret, data.frame(Value = IBS_WD1, Rel = "WD", Type = "IBS"))
  
  print("IBSc WD")
  IBSc_WD1 <- c(ibsOwnDF[idDF$workers, idDF$drones])
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
  print("IBS DD")
  tmp <- ibsMultiDF[idDF$drones, idDF$drones]
  IBS_DD1 <- c(tmp[lower.tri(tmp, diag = FALSE)])
  ret <- rbind(ret, data.frame(Value = IBS_DD1, Rel = "DD", Type = "IBS"))
  
  print("IBSc DD")
  tmp <- ibsOwnDF[idDF$drones, idDF$drones]
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
  print("IBS QW")
  ret <- rbind(ret, data.frame(Value = c(ibsMultiDF[idDF$workers, idDF$queen]),
                               Rel = "QW", Type = "IBS"))
  
  print("IBSc QW")
  ret <- rbind(ret, data.frame(Value = c(ibsOwnDF[idDF$workers, idDF$queen]),
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
  print("IBS QD")
  ret <- rbind(ret, data.frame(Value = c(ibsMultiDF[idDF$queen, idDF$drones]),
                               Rel = "QD", Type = "IBS"))
  
  print("IBSc QD")
  ret <- rbind(ret, data.frame(Value = c(ibsOwnDF[idDF$queen, idDF$drones]),
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

prepareDataForPlotting_ColonyDiag <- function(ibsMultiDF = NULL, ibsOwnDF = NULL, ibdDF = NULL, Sinv = NULL,  idDF) {
  # workers
  print("IBS WW")
  tmp <- diag(ibsMultiDF[idDF$workers, idDF$workers])
  ret <- data.frame(Value = tmp, Rel = "WW", Type = "IBS")
  
  print("IBSc WW")
  tmp <- diag(ibsOwnDF[idDF$workers, idDF$workers])
  ret <- rbind(ret, data.frame(Value = tmp, Rel = "WW", Type = "IBSOwnFreq"))
  
  print("IBDr WW")
  tmp <- diag(ibdDF[idDF$workers, idDF$workers])
  ret <- rbind(ret, data.frame(Value = tmp, Rel = "WW", Type = "IBDr"))
  
  if (!is.null(Sinv)) {
    print("IBDe WW")
    tmp <- getS(Sinv, ids = idDF$workers, vector = TRUE, diagOnly = TRUE)
    ret <- rbind(ret, data.frame(Value = tmp, Rel = "WW", Type = "IBDe"))
  }
  
  #workers to drones 
  print("IBS WD")
  tmp <- diag(ibsMultiDF[idDF$workers, idDF$drones])
  ret <- rbind(ret, data.frame(Value = tmp, Rel = "WD", Type = "IBS"))
  
  print("IBSc WD")
  tmp <- diag(ibsOwnDF[idDF$workers, idDF$drones])
  ret <- rbind(ret, data.frame(Value = tmp, Rel = "WD", Type = "IBSOwnFreq"))
  
  print("IBDr WD")
  tmp <- diag(ibdDF[idDF$workers, idDF$drones])
  ret <- rbind(ret, data.frame(Value = tmp, Rel = "WD", Type = "IBDr"))
  
  if (!is.null(Sinv)) {
    print("IBDe WD")
    tmp <- getS(Sinv, ids = idDF$drones, with = idDF$workers, vector = TRUE, diagOnly = TRUE)
    ret <- rbind(ret, data.frame(Value = tmp, Rel = "WD", Type = "IBDe"))
  }

  # drones
  print("IBS DD")
  tmp <- diag(ibsMultiDF[idDF$drones, idDF$drones])
  ret <- rbind(ret, data.frame(Value = tmp, Rel = "DD", Type = "IBS"))
  
  print("IBSc DD")
  tmp <- diag(ibsOwnDF[idDF$workers, idDF$workers])
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

#Same as plot colony but without cDF 
prepareDataForPlotting_ColonyCSD <- function(ibsMultiDF = NULL, ibdDF = NULL, Sinv = NULL,  idDF) {
  # workers to workers
  print("IBS WW")
  tmp <- ibsMultiDF[idDF$workers, idDF$workers]
  IBS_WW1 <- c(tmp[lower.tri(tmp, diag = FALSE)])
  ret <- data.frame(Value = IBS_WW1, Rel = "WW", Type = "IBS")
  
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
  print("IBS WD")
  IBS_WD1 <- c(ibsMultiDF[idDF$workers, idDF$drones])
  ret <- rbind(ret, data.frame(Value = IBS_WD1, Rel = "WD", Type = "IBS"))
  
  print("IBDr WD")
  IBDr_WD1 <- c(ibdDF[idDF$workers, idDF$drones])
  ret <- rbind(ret, data.frame(Value = IBDr_WD1, Rel = "WD", Type = "IBDr"))
  
  if (!is.null(Sinv)) {
    print("IBDe WD")
    ret <- rbind(ret, data.frame(Value = getS(Sinv, ids = idDF$workers, with = idDF$drones, vector = TRUE),
                                 Rel = "WD", Type = "IBDe"))
  }
  
  # drones vs drones
  print("IBS DD")
  tmp <- ibsMultiDF[idDF$drones, idDF$drones]
  IBS_DD1 <- c(tmp[lower.tri(tmp, diag = FALSE)])
  ret <- rbind(ret, data.frame(Value = IBS_DD1, Rel = "DD", Type = "IBS"))
  
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
  print("IBS QW")
  ret <- rbind(ret, data.frame(Value = c(ibsMultiDF[idDF$workers, idDF$queen]),
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
  print("IBS QD")
  ret <- rbind(ret, data.frame(Value = c(ibsMultiDF[idDF$queen, idDF$drones]),
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

prepareDataForPlotting_Queens <- function(ibsMultiDF = NULL, ibdDF = NULL, Sinv = NULL,  idPopDF) {
  #Get populations IDs
  melID <- idPopDF$ID[idPopDF$Pop == "Mel"]
  melCrossID <- idPopDF$ID[idPopDF$Pop == "MelCross"]
  carID <- idPopDF$ID[idPopDF$Pop == "Car"]
  
  #Compute relationships between queens of different populations (QQ)
  #IBS
  IBSMelMelCross <- data.frame(Value = as.vector(list(ibsMultiDF[melID, melCrossID])[[1]]),
                                Pops = "Mel_MelCross", Rel = "QQ")
  IBSMelCar <- data.frame(Value = as.vector(list(ibsMultiDF[melID, carID])[[1]]),
                           Pops = "Mel_Car", Rel = "QQ")
  IBSMelCrossCar <- data.frame(Value = as.vector(list(ibsMultiDF[melCrossID, carID])[[1]]),
                                Pops = "MelCross_Car", Rel = "QQ")
  IBS <-rbind(IBSMelMelCross, IBSMelCar, IBSMelCrossCar)
  IBS$Type <- "IBS"
  
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
  #IBS
  tmp <- ibsMultiDF[melID, melID]
  IBSMel_Mel <- c(tmp[lower.tri(tmp, diag = FALSE)])
  IBSMel_Mel <- data.frame(Value = as.vector(list(IBSMel_Mel)[[1]]), Rel = "Q", Type = "IBS", Pops = "Mel")
  
  tmp <- ibsMultiDF[carID, carID]
  IBSCar_Car <- c(tmp[lower.tri(tmp, diag = FALSE)])
  IBSCar_Car <-data.frame(Value = as.vector(list(IBSCar_Car)[[1]]), Rel = "Q", Type = "IBS", Pops = "Car")
  
  tmp <- ibsMultiDF[melCrossID, melCrossID]
  IBSmelCross_melCross <-  c(tmp[lower.tri(tmp, diag = FALSE)])
  IBSmelCross_melCross <- data.frame(Value = as.vector(list(IBSmelCross_melCross)[[1]]), Rel = "Q", Type = "IBS", Pops = "MelCross")
  
  IBS <- rbind(IBS, IBSMel_Mel, IBSCar_Car, IBSmelCross_melCross)
  
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
  inbIBS <- rbind(data.frame(Value = diag(ibsMultiDF[melID, melID]),
                              Pops = "Mel"),
                   data.frame(Value = diag(ibsMultiDF[melCrossID, melCrossID]),
                              Pops = "MelCross"),
                   data.frame(Value = diag(ibsMultiDF[carID, carID]),
                              Pops = "Car"))
  inbIBS$Type = "IBS"
  inbIBS$Rel = "F"
  IBS <- rbind(IBS, inbIBS)
  
  inbIBDr <- rbind(data.frame(Value = diag(ibdDF[melID, melID]),
                              Pops = "Mel"),
                   data.frame(Value = diag(ibdDF[melCrossID, melCrossID]),
                              Pops = "MelCross"),
                   data.frame(Value = diag(ibdDF[carID, carID]),
                              Pops = "Car"))
  inbIBDr$Type = "IBDr"
  inbIBDr$Rel = "F"
  IBDr <- rbind(IBDr, inbIBDr)
  
  ret <- rbind(IBS, IBDr)
  
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

prepareDataForPlotting_QueensCSD <- function(ibsMultiDF = NULL, ibdDF = NULL, Sinv = NULL,  idPopDF) {
  #Get populations IDs
  melID <- idPopDF$ID[idPopDF$Pop == "Mel"]
  melCrossID <- idPopDF$ID[idPopDF$Pop == "MelCross"]
  carID <- idPopDF$ID[idPopDF$Pop == "Car"]
  
  #Compute relationships between queens of different populations (QQ)
  #IBS
  IBSMelMelCross <- data.frame(Value = as.vector(list(ibsMultiDF[melID, melCrossID])[[1]]),
                                Pops = "Mel_MelCross", Rel = "QQ")
  IBSMelCar <- data.frame(Value = as.vector(list(ibsMultiDF[melID, carID])[[1]]),
                           Pops = "Mel_Car", Rel = "QQ")
  IBSMelCrossCar <- data.frame(Value = as.vector(list(ibsMultiDF[melCrossID, carID])[[1]]),
                                Pops = "MelCross_Car", Rel = "QQ")
  IBS <-rbind(IBSMelMelCross, IBSMelCar, IBSMelCrossCar)
  IBS$Type <- "IBS"
  
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
  #IBS
  tmp <- ibsMultiDF[melID, melID]
  IBSMel_Mel <- c(tmp[lower.tri(tmp, diag = FALSE)])
  IBSMel_Mel <- data.frame(Value = as.vector(list(IBSMel_Mel)[[1]]), Rel = "Q", Type = "IBS", Pops = "Mel")
  
  tmp <- ibsMultiDF[carID, carID]
  IBSCar_Car <- c(tmp[lower.tri(tmp, diag = FALSE)])
  IBSCar_Car <-data.frame(Value = as.vector(list(IBSCar_Car)[[1]]), Rel = "Q", Type = "IBS", Pops = "Car")
  
  tmp <- ibsMultiDF[melCrossID, melCrossID]
  IBSmelCross_melCross <-  c(tmp[lower.tri(tmp, diag = FALSE)])
  IBSmelCross_melCross <- data.frame(Value = as.vector(list(IBSmelCross_melCross)[[1]]), Rel = "Q", Type = "IBS", Pops = "MelCross")
  
  IBS <- rbind(IBS, IBSMel_Mel, IBSCar_Car, IBSmelCross_melCross)
  
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
  inbIBS <- rbind(data.frame(Value = diag(ibsMultiDF[melID, melID]),
                              Pops = "Mel"),
                   data.frame(Value = diag(ibsMultiDF[melCrossID, melCrossID]),
                              Pops = "MelCross"),
                   data.frame(Value = diag(ibsMultiDF[carID, carID]),
                              Pops = "Car"))
  inbIBS$Type = "IBS"
  inbIBS$Rel = "F"
  IBS <- rbind(IBS, inbIBS)
  
  inbIBDr <- rbind(data.frame(Value = diag(ibdDF[melID, melID]),
                              Pops = "Mel"),
                   data.frame(Value = diag(ibdDF[melCrossID, melCrossID]),
                              Pops = "MelCross"),
                   data.frame(Value = diag(ibdDF[carID, carID]),
                              Pops = "Car"))
  inbIBDr$Type = "IBDr"
  inbIBDr$Rel = "F"
  IBDr <- rbind(IBDr, inbIBDr)
  
  ret <- rbind(IBS, IBDr)
  
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

prepareDataForPlottingHeatMap_Queens <- function(ibsMultiDF = NULL, ibsOwnDF = NULL, ibdDF = NULL, Sinv = NULL, idDF = NULL) {
  ibsMultiDF <- as.data.frame(ibsMultiDF)
  columns <- colnames(ibsMultiDF)
  ibsMultiDF$ID <- rownames(ibsMultiDF)
  ibsDFL <- ibsMultiDF %>% pivot_longer(cols = all_of(columns))
  ibsDFL$Method <- "IBS"
  
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

prepareDataForPlottingHeatMap_QueensCSD <- function(ibsMultiDF = NULL, ibdDF = NULL, Sinv = NULL, idDF = NULL) {
  ibsMultiDF <- as.data.frame(ibsMultiDF)
  columns <- colnames(ibsMultiDF)
  ibsMultiDF$ID <- rownames(ibsMultiDF)
  ibsDFL <- ibsMultiDF %>% pivot_longer(cols = all_of(columns))
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

plotColony <- function(df, rel = NULL, type = c("IBDe", "IBDr", "IBSmultiBF", "IBSsingleBF", "IBSOwnFreq"), years = NULL) {
  df$Type <- factor(df$Type, levels = c("IBDe", "IBDr", "IBSmultiBF", "IBSsingleBF","IBSOwnFreq"))
  df$Rel <- factor(df$Rel, levels= c("WW", "WD", "DD", "QW", "QD"))
  type_labels <- c("IBDe", "IBDr", "IBSmultiBF", "IBSsingleBF" ,"IBSOwnFreq")
  names(type_labels) <- c("IBDe", "IBDr", "IBSmultiBF", "IBSsingleBF" ,"IBSOwnFreq")
   
  df$Year <- factor(df$Year, levels = years)
  year_labels <- unlist(paste("Year", years, collapse = "_") %>% stringr::str_split("_"))
  names(year_labels) <- years
  
  p <- ggplot(df[df$Rel %in% rel & df$Type %in% type, ],
              aes(x = Value, fill = Rel)) +
    geom_histogram(binwidth = 0.01, position = "identity") +
    facet_grid2(rows = vars(Type), cols = vars(Year), scales = "free", independent = "y", labeller = labeller(Type = type_labels, Year = year_labels))
  plot <- p +
    scale_fill_manual("", values=paletteViridis, aesthetics = c("colour","fill")) +
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "top") +
    scale_x_continuous(breaks=seq(-2, 2.5, 0.25))
  return(plot)
}

plotQueensQQ <- function(df, rel = c("QQ"), type = c("IBDe", "IBDr", "IBSmultiBF", "IBSsingleBF", "IBSOwnFreq")) {
  df$Type <- factor(df$Type, levels = c("IBDe", "IBDr", "IBSmultiBF", "IBSsingleBF","IBSOwnFreq"))
  df$Pops <- factor(df$Pops, levels = c("Mel_MelCross", "Mel_Car", "MelCross_Car"))
  df$Rel  <- factor(df$Rel, levels = "QQ")
  type_labels <- c("IBDe", "IBDr", "IBSmultiBF", "IBSsingleBF","IBSOwnFreq")
  year_labels <- c("Year 1", "Year 10")
  names(type_labels) <- c("IBDe", "IBDr", "IBSmultiBF", "IBSsingleBF","IBSOwnFreq")
  names(year_labels) <- c("1", "10")
  
  plot <- ggplot(df[df$Rel %in% rel & df$Type %in% type, ],
                 aes(x = Value, fill = Pops)) +
          geom_histogram(binwidth = 0.01, position = "identity") +
          facet_grid2(rows = vars(Type), cols = vars(Year), scales = "free", independent = "y", labeller = labeller(Type = type_labels, Year = year_labels)) +
          scale_fill_manual("", values = cbPaletteQQ, aesthetics = c("colour","fill")) +
          theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "top") +
          scale_x_continuous(breaks=seq(-2, 2, 0.25))
  return(plot)
}

scatterQueens <-  function(df, rel = NULL, type = c("IBDr", "IBS"), pops = NULL, legend.position = NULL, palette = NULL, years = NULL) {
  df$Type <- factor(df$Type, levels = c("IBDe", "IBDr", "IBS", "IBSOwnFreq"))
  df$Pops <- factor(df$Pops, levels = c("Mel_MelCross", "Mel_Car", "MelCross_Car", "Car", "Mel", "MelCross"))
  df$Rel  <- factor(df$Rel, levels = c("QQ", "Q"))
  type_labels <- c("IBDe", "IBDr", "IBSbaseAF")
  names(type_labels) <- c("IBDe", "IBDr", "IBS")

  
  df$Year <- factor(df$Year, levels = years)
  year_labels <- unlist(paste("Year", years, collapse = "_") %>% stringr::str_split("_"))
  names(year_labels) <- years
  
  a <- pivot_wider(df[df$Rel %in% rel & df$Type %in% type & df$Pops %in% pops, ], names_from = Type, values_from = Value, values_fn = list)
  b <- unnest(a, cols = all_of(type))
  
  plot <- ggplot(data = b, aes(x = IBDr, y = IBS)) +
          geom_point(aes(colour = Pops)) + 
          facet_grid2(cols = vars(Year), scales = "free", independent = "y", labeller = labeller(Type = type_labels, Year = year_labels)) +
          scale_fill_manual("", values = palette, aesthetics = c("colour","fill")) + 
          theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = legend.position) +
          scale_x_continuous(breaks=seq(-2, 2, 0.25))
  return(plot)
}

plotQueensQ <- function(df, rel = NULL, type = c("IBDe", "IBDr", "IBS")) {
  df$Type <- factor(df$Type, levels = c("IBDe", "IBDr", "IBS", "IBSOwnFreq"))
  df$Rel <- factor(df$Rel, levels= c("F", "Q"))
  type_labels <- c("IBDe", "IBDr", "IBSbaseAF")
  year_labels <- c("Year 1", "Year 10")
  names(type_labels) <- c("IBDe", "IBDr", "IBS")
  names(year_labels) <- c("1", "10")
  plot <- ggplot(df[df$Rel %in% rel & df$Type %in% type, ],
                 aes(x = Value, fill = Pops)) +
       geom_histogram(binwidth = 0.01, position = "identity") +
       facet_grid2(rows = vars(Type), cols = vars(Year), scales = "free", independent = "y", labeller = labeller(Type = type_labels, Year = year_labels)) +
       scale_fill_manual("", values = cbPaletteQ, aesthetics = c("colour","fill")) + 
       theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "top") + 
       scale_x_continuous(breaks=seq(-2, 2.5, 0.25))
  return(plot)
}

plotQueens_heatmapBIND <- function(df, PopIdDF = NULL, years = NULL) {
  df$ID <- as.factor(as.numeric(df$ID))
  df$name <- as.factor(as.numeric(df$name))
  df$Method <- factor(df$Method, levels = c("IBDe", "IBDr", "IBS", "IBSOwnFreq"))
  method_labels <- c("IBDe", "IBDr", "IBSbaseAF", "IBSColonyAF")
  names(method_labels) <- c("IBDe", "IBDr", "IBS", "IBSOwnFreq")
  
  df$Year <- factor(df$Year, levels = years)
  year_labels <- unlist(paste("Year", years, collapse = "_") %>% stringr::str_split("_"))
  names(year_labels) <- years
  
  df <- merge(df, PopIdDF, by = "ID")
  df <- merge(df, PopIdDF, by.x = "name", by.y = "ID")
  df$PopId1 <- paste0(df$Pop.x, df$ID)
  df$PopId2 <- paste0(df$Pop.y, df$name)
  
  n <- length(unique(df$name[df$Pop.y == "Mel"]))/2
  breaks = df %>%
    arrange(Year, Pop.y) %>%
    dplyr::mutate(row = row_number()) %>%
    group_by(Year) %>%
    dplyr::mutate(min = min(row)) %>%
    group_by(Year, Pop.y) %>%
    summarise(breaks = max(row_number()), min = mean(min)) %>%
    ungroup() %>% group_by(Year) %>%
    dplyr::mutate(breaks2 = ifelse(test = row_number() > 1, yes = breaks/2 + breaks*(row_number()-1), no =  breaks/2)) %>%
    ungroup() %>% 
    dplyr::mutate(breaks3 = dplyr::select(., min, breaks2) %>% 
                                     rowSums(.)) %>% 
    dplyr::select(breaks = breaks3)
  
  plot <- ggplot(data = df, aes(x = PopId1, y = PopId2, fill = value)) + geom_tile() + scale_fill_gradient(low = "white", high = "blue") + 
    theme(panel.background = element_blank(), legend.position = "top") +
    scale_x_discrete(breaks = breaks$breaks, labels = rep(c("Car", "Mel", "MelCross"), length(unique(df$Year)))) +
    scale_y_discrete(labels = rep(c("Car", "Mel", "MelCross"), length(unique(df$Year)))) +
    xlab("") +
    ylab("") +
    facet_grid2(rows = vars(Method), cols = vars(Year), labeller = labeller(Method = method_labels, Year = year_labels), scales = "free", independent = "y")
  
  return(plot)
}

plotQueens_heatmapSOLO <- function(df, Pop = FALSE, PopIdDF = NULL, legend.position = NULL) {
  df$ID <- as.factor(as.numeric(df$ID))
  df$name <- as.factor(as.numeric(df$name))
  df$Method <- factor(df$Method, levels = c("IBDe", "IBDr", "IBS", "IBSOwnFreq"))
  method_labels <- c("IBDe", "IBDr", "IBSbaseAF", "IBSColonyAF")
  names(method_labels) <- c("IBDe", "IBDr", "IBS", "IBSOwnFreq")
  
  if (Pop) {
    df <- merge(df, PopIdDF, by = "ID")
    df <- merge(df, PopIdDF, by.x = "name", by.y = "ID")
    
    df$PopId1 <- paste0(df$Pop.x, df$ID)
    df$PopId2 <- paste0(df$Pop.y, df$name)
    
    n <- length(unique(df$name[df$Pop.y == "Mel"]))/2
    breaks = list(df %>% group_by(Pop.y) %>% summarise(mean = unique(name)[n]) %>% unite(PopMean, Pop.y, mean, sep=""))[[1]]$PopMean
    plot <- ggplot(data = df, aes(x = PopId1, y = PopId2, fill = value)) + geom_tile() + scale_fill_gradient(low = "white", high = "blue") + 
      theme(panel.background = element_blank(), legend.position = legend.position) +
      scale_x_discrete(breaks = breaks, labels = c("Car", "Mel", "MelCross")) +
      scale_y_discrete(breaks = breaks, labels = c("Car", "Mel", "MelCross")) +
      xlab("") + ylab("") +
      facet_grid(rows = vars(Method))
  } else {
    plot <- ggplot(data = df, aes(x = ID, y = name, fill = value)) + geom_tile() +
      facet_grid2(rows = vars(Method), labeller = labeller(Method = method_labels), scales = "free")
  }
  return(plot)
}

########################################################
### --- FIGURE 1&2: Pure subspecies (carnica) in years 1 and 10 ---###
########################################################
print("Plot carnica year 1")
relCar1 <- prepareDataForPlotting_Colony(ibsMultiDF = colonyCar1$IBSmultiBF,
                                         ibsSingleDF = colonyCar1$IBSsingleBF,
                                         ibsOwnDF = colonyCar1$IBSOwnFreq, 
                                         ibdDF = colonyCar1$IBD, 
                                         Sinv = Sinv,
                                         idDF = colonyCar1$ID,
                                         inbreeding = FALSE)

print("Plot carnica year 10")
relCar10 <- prepareDataForPlotting_Colony(ibsMultiDF = colonyCar10$IBS, ibsOwnDF = colonyCar10$IBSOwnFreq, ibdDF = colonyCar10$IBD, Sinv = Sinv, idDF = colonyCar10$ID)

relCar1$Year <- 1
relCar10$Year <- 10
dataCar <- rbind(relCar1, relCar10)

CarWWplot<- plotColony(dataCar, type = c("IBDe", "IBDr", "IBS", "IBSOwnFreq"), rel = c("WW", "WD", "DD"))
CarWWplot

CarQWplot <- plotColony(dataCar, type = c("IBDe", "IBDr", "IBS", "IBSOwnFreq"), rel = c("QW", "QD"))
CarQWplot

#Plot Car inbreeding
relCar1F <- prepareDataForPlotting_Colony(ibsMultiDF = colonyCar1$IBS,
                                          ibsOwnDF = colonyCar1$IBSOwnFreq,
                                          ibdDF = colonyCar1$IBD,
                                          Sinv = Sinv,
                                          idDF = colonyCar1$ID, 
                                          diag = TRUE)
relCar10F <- prepareDataForPlotting_ColonyDiag(ibsMultiDF = colonyCar10$IBS, ibsOwnDF = colonyCar10$IBSOwnFreq, ibdDF = colonyCar10$IBD, Sinv = Sinv, idDF = colonyCar10$ID)
relCar1F$Year <- 1
relCar10F$Year <- 10 

dataCarF <- rbind(relCar1F, relCar10F)

CarFplot <- plotColony(dataCarF, type = c("IBDe", "IBDr", "IBS", "IBSOwnFreq"), rel = c("WW", "WD", "DD"))
CarFplot
#Plot CAR csd Locus
#Year 1
relCar1_csd <- prepareDataForPlotting_ColonyCSD(ibsMultiDF = colonyCar1$IBSCsd, ibdDF = colonyCar1$IBDCsd, Sinv = Sinv, idDF = colonyCar1$ID)
relCar10_csd <- prepareDataForPlotting_ColonyCSD(ibsMultiDF = colonyCar10$IBSCsd, ibdDF = colonyCar10$IBDCsd, Sinv = Sinv, idDF = colonyCar10$ID)
relCar1_csd$Year <- 1
relCar10_csd$Year <- 10 

dataCarCSDloc <- rbind(relCar1_csd, relCar10_csd)

CarCSDlocplot <- plotColony(dataCarCSDloc, type = c("IBDe", "IBDr", "IBS"), rel = c("WW", "WD", "DD"))
CarCSDlocplot

#Plot CAR csd Chromosome
relCar1_csdChr <- prepareDataForPlotting_ColonyCSD(ibsMultiDF = colonyCar1$IBScsdChr, ibdDF = colonyCar1$IBDcsdChr, Sinv = Sinv, idDF = colonyCar1$ID)
relCar10_csdChr <- prepareDataForPlotting_ColonyCSD(ibsMultiDF = colonyCar10$IBScsdChr, ibdDF = colonyCar10$IBDcsdChr, Sinv = Sinv, idDF = colonyCar10$ID)
relCar1_csdChr$Year <- 1
relCar10_csdChr$Year <- 10 

dataCarCSDchr <- rbind(relCar1_csdChr, relCar10_csdChr)

CarCSDchr <- plotColony(dataCarCSDchr, type = c("IBDe", "IBDr", "IBS"))
CarCSDchr 
########################################################
### --- FIGURE 4: Between queens of different subspecies (carnica vs. mellifera) ---###
########################################################
#Plot queens Year 1
print("Plot queens")
relQueens1 <- prepareDataForPlotting_Queens(ibsMultiDF =  queens1$IBS, ibdDF = queens1$IBD, Sinv = Sinv, idPopDF = springerQueensPop1)
relQueens10 <- prepareDataForPlotting_Queens(ibsMultiDF = queens10$IBS, ibdDF = queens10$IBD, Sinv = Sinv, idPopDF = springerQueensPop10)
relQueens1$Year <- 1
relQueens10$Year <- 10

dataQueens <- rbind(relQueens1, relQueens10)
print("Between populations ")
#between populations
QueensQQ <- plotQueensQQ(dataQueens, rel = "QQ", type = c("IBDr", "IBDe", "IBS"))
QueensQQ

print("Within populations")
#within population - non-diagonal
QueensQ <- plotQueensQ(dataQueens, rel = "Q",  type = c("IBDr", "IBDe", "IBS"))
QueensQ
#within population - diagonal
QueensF <- plotQueensF(dataQueens, rel = "F", type = c("IBDr", "IBDe", "IBS"))
QueensF

#scatter plot looking at Q and QQ
QueensScatterQQ <- scatterQueens(relQueens10, rel = "QQ", pops = c("Mel_MelCross", "Mel_Car", "MelCross_Car"), palette = cbPaletteQQ, years = 10 , legend.position = "top")
QueensScatterQ <- scatterQueens(relQueens10, rel = "Q", pops = c("Mel", "Car", "MelCross"), palette = cbPaletteQ, years= 10, legend.position = "top")
QueensScatter <- grid.arrange(QueensScatterQ, QueensScatterQQ, ncol = 2)



print("HeatMaps")
relQueens1h <- prepareDataForPlottingHeatMap_Queens(ibsMultiDF =  queens1$IBS, ibdDF = queens1$IBD, Sinv = Sinv, idDF = queens1$ID)
relQueens10h <- prepareDataForPlottingHeatMap_Queens(ibsMultiDF = queens10$IBS, ibdDF = queens10$IBD, Sinv = Sinv, idDF = queens10$ID)
relQueens1h$Year <- 1
relQueens10h$Year <- 10
dataQh <- rbind(relQueens1h, relQueens10h)
idPopQueens <- rbind(springerQueensPop1, springerQueensPop10)


plotQueens1h <- plotQueens_heatmapSOLO(relQueens1h, Pop = TRUE, PopIdDF = springerQueensPop1, legend.position = "left")
plotQueens1h
plotQueens10h <- plotQueens_heatmapSOLO(relQueens10h, Pop = TRUE, PopIdDF = springerQueensPop10, legend.position = "right")
plotQueens10h

arrangedPlots <- grid.arrange(plotQueens1h, plotQueens10h, ncol = 2)

QueenH <- plotQueens_heatmap(dataQh, PopIdDF = idPopQueens, years = c(1,10))
QueenH

#CSD STUFF
#Plot csd Loc as a histogram (QQ)
relQueens1_csdLoc_hist <- prepareDataForPlotting_QueensCSD(ibsMultiDF = queens1$IBSCsd, ibdDF = queens1$IBDCsd, Sinv = Sinv, idPopDF = springerQueensPop1)
relQueens10_csdLoc_hist <- prepareDataForPlotting_QueensCSD(ibsMultiDF = queens10$IBSCsd, ibdDF = queens10$IBDCsd, Sinv = Sinv, idPopDF = springerQueensPop10)
relQueens1_csdLoc_hist$Year <- 1
relQueens10_csdLoc_hist$Year <- 10 
dataQueens_csdLoc <- rbind(relQueens1_csdLoc_hist, relQueens10_csdLoc_hist)

Queens_csdLoc_hist <- plotQueensQQ(dataQueens_csdLoc, rel = "QQ", type = c("IBDr", "IBDe", "IBS", "IBSOwnFreq"))
Queens_csdLoc_hist

#Plot csd loc as a heatmap TODO: 

relQueens1_csdLoc_heat <- prepareDataForPlottingHeatMap_QueensCSD(ibsMultiDF = queens1$IBSCsd, ibdDF = queens1$IBDCsd, Sinv = Sinv, idDF = queens1$ID)
relQueens10_csdLoc_heat <- prepareDataForPlottingHeatMap_QueensCSD(ibsMultiDF = queens10$IBSCsd, ibdDF = queens10$IBDCsd, Sinv = Sinv, idDF = queens10$ID)
relQueens1_csdLoc_heat$Year <- 1
relQueens10_csdLoc_heat$Year <- 10 
dataQueens_csdLoc_heat <- rbind(relQueens1_csdLoc_heat, relQueens10_csdLoc_heat)

QueensCSDloc1h <- plotQueens_heatmapSOLO(relQueens1_csdLoc_heat, Pop = TRUE, PopIdDF = springerQueensPop1, legend.position = "left")
QueensCSDloc10h <- plotQueens_heatmapSOLO(relQueens10_csdLoc_heat, Pop = TRUE, PopIdDF = springerQueensPop10, legend.position = "right")

QueensCSDLocHeat <- grid.arrange(QueensCSDloc1h, QueensCSDloc10h, ncol = 2)


#Plot csd Chr as Histogram
relQueen1_csdChr_hist <- prepareDataForPlotting_QueensCSD(ibsMultiDF = queens1$IBScsdChr, ibdDF = queens1$IBDcsdChr, Sinv = Sinv, idPopDF = springerQueensPop1)
relQueen10_csdChr_hist <- prepareDataForPlotting_QueensCSD(ibsMultiDF = queens10$IBScsdChr, ibdDF = queens10$IBDcsdChr, Sinv = Sinv, idPopDF = springerQueensPop10)
relQueen1_csdChr_hist$Year <- 1
relQueen10_csdChr_hist$Year <- 10 
dataQueens_csdChr_hist <- rbind(relQueen1_csdChr_hist, relQueen10_csdChr_hist)

Queens_csdChr_hist <- plotQueensQQ(dataQueens_csdChr_hist, rel = "QQ", type = c("IBDr", "IBDe", "IBS", "IBSOwnFreq"))


#Plot csd Chr as heatmap
relQueens1_csdChr_heat <- prepareDataForPlottingHeatMap_QueensCSD(ibsMultiDF = queens1$IBScsdChr, ibdDF = queens1$IBDcsdChr, Sinv = Sinv, idDF = queens1$ID)
relQueens10_csdChr_heat <- prepareDataForPlottingHeatMap_QueensCSD(ibsMultiDF = queens10$IBScsdChr, ibdDF = queens10$IBDcsdChr, Sinv = Sinv, idDF = queens10$ID)
relQueens1_csdChr_heat$Year <- 1
relQueens10_csdChr_heat$Year <- 10
dataQueens_csdChr_heat <- rbind(relQueens1_csdChr_heat, relQueens10_csdChr_heat)

QueensCSDchrHeat1 <- plotQueens_heatmapSOLO(relQueens1_csdChr_heat, Pop = TRUE, PopIdDF = springerQueensPop1, legend.position = "left")
QueensCSDchrHeat10 <- plotQueens_heatmapSOLO(relQueens10_csdChr_heat, Pop = TRUE, PopIdDF = springerQueensPop10, legend.position = "right")

QueenCSDchrHeat <- grid.arrange(QueensCSDchrHeat1, QueensCSDchrHeat10, ncol = 2)

Queens_CSDchr_head <- plotQueens_heatmapBIND(dataQueens_csdChr_heat, PopIdDF = idPopQueens, years= c(1,10))

### Csd
# csdVariability <- csdVariability[!is.na(csdVariability$subspecies),]
# pDiploidDrones <- pDiploidDrones[-1,]
# 
# csdVariability$nCSD <- as.numeric(csdVariability$nCSD)
# csdVariability$totalCSD <- as.numeric(csdVariability$totalCSD)
# csdVariability$year <- as.numeric(csdVariability$year)
# csdMean <- csdVariability %>%  group_by(subspecies, year) %>%  summarize(meanCSD = mean(nCSD))
# csdMean$subspecies <- factor(csdMean$subspecies, levels = c("Mel", "MelCross", "Car"))
# csdVariability$subspecies <- factor(csdVariability$subspecies, levels = c("Mel", "MelCross", "Car"))
# ggplot(data = csdMean, aes(x = year, y = meanCSD, colour = subspecies)) + geom_line()
# 
# ggplot(data = csdVariability, aes(x = year, y = totalCSD, colour = subspecies)) + geom_line()
# 
# 
# 
# pDiploidDrones$pQueenHomBrood <- as.numeric(pDiploidDrones$pQueenHomBrood)
# pDiploidDrones$year <- as.numeric(pDiploidDrones$year)
# pDiploidMean <- pDiploidDrones %>%  group_by(subspecies, year) %>%  summarize(meanHom = mean(pQueenHomBrood))
# pDiploidMean$subspecies <- factor(pDiploidMean$subspecies, levels = c("Mel", "MelCross", "Car"))
# 
# ggplot(data = pDiploidMean, aes(x = year, y = meanHom, colour = subspecies)) + geom_line()


