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
data <- load("~/Desktop/GitHub/lstrachan_honeybee_sim/YearCycleSimulation/PlottingData/SpringerSimulation_import.RData")
Sinv <- readMM("~/Desktop/GitHub/lstrachan_honeybee_sim/YearCycleSimulation/PLottingData/Sinv.mm")
#Eddie data
#data <- load("SpringerSimulation_import_objects.RData")
#Sinv <- readMM("Sinv.mm")

# The data contains two populations - mellifera and carnica
# The carnica stays "pure" throughout the simulation
# Melifera gets mated with a proportion of carnica drones

#Plotting help
# The colourblind palettes:
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
prepareDataForPlotting_Colony <- function(ibsMultiDF = NULL, ibsSingleDF = NULL, ibsColonyDF = NULL, ibdDF = NULL, Sinv = NULL, idDF, inbreeding = FALSE) {

  #Workers on workers
  print("IBSmulti WW")
  start_time <- Sys.time()
  ibsMultiWW <- data.frame(Value = ifelse(test = isTRUE(inbreeding),
                                          yes = data.frame(diag(ibsMultiDF[idDF$workers, idDF$workers])),
                                          no = data.frame(ibsMultiDF[idDF$workers, idDF$workers] %>%
                                                               c(.[lower.tri(., diag = FALSE)]))) ,
                           Rel = "WW",
                           Type = "IBSmultiBF")
  colnames(ibsMultiWW) <- c("Value", "Rel", "Type")
  Sys.time() - start_time #FOR A TIME CHECK

  print("IBSsingle WW")
  start_time <- Sys.time()
  ibsSingleWW <- data.frame(Value = data.frame(ifelse(test = isTRUE(inbreeding),
                                                      yes = data.frame(diag(ibsSingleDF[idDF$workers, idDF$workers])),
                                                      no = data.frame(ibsSingleDF[idDF$workers, idDF$workers] %>%
                                                           c(.[lower.tri(., diag = FALSE)])))),
                            Rel = "WW",
                            Type = "IBSsingleBF")
  colnames(ibsSingleWW) <- c("Value", "Rel", "Type")

  if(!is.null(ibsColonyDF)){ #if dataframe is CSD info then ibdOwnDF is NULL
    print("IBScolony WW")
    ibsColonyWW <- data.frame(Value = ifelse(test = isTRUE(inbreeding),
                                          yes = data.frame(diag(ibsColonyDF[idDF$workers, idDF$workers])),
                                          no = data.frame(ibsColonyDF[idDF$workers, idDF$workers] %>%
                                                c(.[lower.tri(., diag = FALSE)]))),
                           Rel = "WW",
                           Type = "IBScolonyF")
    colnames(ibsColonyWW) <- c("Value", "Rel", "Type")

  }

  print("IBDr WW")
  ibdrWW <- data.frame(Value = ifelse(test = isTRUE(inbreeding),
                                      yes = data.frame(diag(ibdDF[idDF$workers, idDF$workers])),
                                      no = data.frame(ibdDF[idDF$workers, idDF$workers] %>%
                                            c(.[lower.tri(., diag = FALSE)]))),
                         Rel = "WW",
                       Type = "IBDr")
  colnames(ibdrWW) <- c("Value", "Rel", "Type")

  if (!is.null(Sinv)) {
    print("IBDe WW")
    ibdeWW <- data.frame(Value = ifelse(test = isTRUE(inbreeding),
                                        yes = data.frame(getS(Sinv, ids = idDF$workers, vector = TRUE, diagOnly = TRUE)),
                                        no = data.frame(getS(Sinv, ids = idDF$workers, vector = TRUE) %>%
                                                        c(.[lower.tri(., diag = FALSE)]))),
                         Rel = "WW",
                         Type = "IBDe")
    colnames(ibdeWW) <- c("Value", "Rel", "Type")
  }
  # Workers vs drones
  if (isFALSE(inbreeding)) {
    print("IBSmulti WD")
    ibsMultiWD <- data.frame(Value = c(ibsMultiDF[idDF$workers, idDF$drones]),
                             Rel = "WD",
                             Type = "IBSmultiBF")

    print("IBSsingle WD")
    ibsSingleWD <- data.frame(Value =  c(ibsSingleDF[idDF$workers, idDF$drones]),
                              Rel = "WD",
                              Type = "IBSsingleBF")

    if (!is.null(ibsColonyDF)){
      print("IBScolony WD")
      ibsColonyWD <- data.frame(Value = c(ibsColonyDF[idDF$workers, idDF$drones]),
                             Rel = "WD",
                             Type = "IBScolonyF")
    }

    print("IBDr WD")
    ibdrWD <- data.frame(Value = c(ibdDF[idDF$workers, idDF$drones]),
                         Rel = "WD",
                         Type = "IBDr")

    if (!is.null(Sinv)) {
      print("IBDe WD")
      ibdeWD <- data.frame(Value = getS(Sinv, ids = idDF$workers, with = idDF$drones, vector = TRUE),
                           Rel = "WD",
                           Type = "IBDe")
    }} else {
      ibsMultiWD <- NULL
      ibsSingleWD <- NULL
      ibsColonyWD <- NULL
      ibdrWD <- NULL
      ibdeWD <- NULL
    }

  #Drones vs drones
  print("IBSmulti DD")
  start_time <- Sys.time()
  ibsMultiDD <- data.frame(Value = ifelse(test = isTRUE(inbreeding),
                                          yes = data.frame(diag(ibsMultiDF[idDF$drones, idDF$drones])),
                                          no = data.frame(ibsMultiDF[idDF$drones, idDF$drones] %>%
                                            c(.[lower.tri(., diag = FALSE)]))),
                           Rel = "DD",
                           Type = "IBSmultiBF")
  colnames(ibsMultiDD) <- c("Value", "Rel", "Type")
  #Sys.time() - start_time FOR A TIME CHECK

  print("IBSsingle DD")
  start_time <- Sys.time()
  ibsSingleDD <- data.frame(Value = ifelse(test = isTRUE(inbreeding),
                                           yes = data.frame(diag(ibsSingleDF[idDF$drones, idDF$drones])),
                                           no = data.frame(ibsSingleDF[idDF$drones, idDF$drones] %>%
                                             c(.[lower.tri(., diag = FALSE)]))),
                            Rel = "DD",
                            Type = "IBSsingleBF")
  colnames(ibsSingleDD) <- c("Value", "Rel", "Type")
  if(!is.null(ibsColonyDF)){ #if dataframe is CSD info then ibdCdf is NULL
    print("IBScolony DD")
    ibsColonyDD <- data.frame(Value = ifelse(test = isTRUE(inbreeding),
                                          yes = data.frame(diag(ibsColonyDF[idDF$drones, idDF$drones])),
                                          no = data.frame(ibsColonyDF[idDF$drones, idDF$drones] %>%
                                            c(.[lower.tri(., diag = FALSE)]))),
                           Rel = "DD",
                           Type = "IBScolonyF")
    colnames(ibsColonyDD) <- c("Value", "Rel", "Type")
  }

  print("IBDr DD")
  ibdrDD <- data.frame(Value = ifelse(test = isTRUE(inbreeding),
                                      yes = data.frame(diag(ibdDF[idDF$drones, idDF$drones])),
                                      no = data.frame(ibdDF[idDF$drones, idDF$drones] %>%
                                        c(.[lower.tri(., diag = FALSE)]))),
                       Rel = "DD",
                       Type = "IBDr")
  colnames(ibdrDD) <- c("Value", "Rel", "Type")

  if (!is.null(Sinv)) {
    print("IBDe DD")
    ibdeDD <- data.frame(Value = ifelse(test = isTRUE(inbreeding),
                                         yes = data.frame(getS(Sinv, ids = idDF$drones, vector = TRUE, diagOnly = TRUE)),
                                         no = data.frame(getS(Sinv, ids = idDF$drones, vector = TRUE) %>%
                                          c(.[lower.tri(., diag = FALSE)]))),
                         Rel = "DD",
                         Type = "IBDe")
    colnames(ibdeDD) <- c("Value", "Rel", "Type")
  }

  # Queens vs workers
  if (isFALSE(inbreeding)){
    print("IBSmulti QW")
    ibsMultiQW <- data.frame(Value = c(ibsMultiDF[idDF$workers, idDF$queen]),
                             Rel = "QW",
                             Type = "IBSmultiBF")

    print("IBSsingle QW")
    ibsSingleQW <- data.frame(Value = c(ibsSingleDF[idDF$workers, idDF$queen]),
                              Rel = "QW",
                              Type = "IBSsingleBF")

    if (!is.null(ibsColonyDF)){
      print("IBScolony QW")
      ibsColonyQW <- data.frame(Value = c(ibsColonyDF[idDF$workers, idDF$queen]),
                             Rel = "QW",
                             Type = "IBScolonyF")
    }

    print("IBDr QW")
    ibdrQW <- data.frame(Value = c(ibdDF[idDF$workers, idDF$queen]),
                         Rel = "QW",
                         Type = "IBDr")

    if (!is.null(Sinv)) {
      print("IBDe QW")
      ibdeQW <- data.frame(Value = getS(Sinv, ids = idDF$workers, with = idDF$queen, vector = TRUE),
                           Rel = "QW",
                           Type = "IBDe")
    }} else {
      ibsMultiQW <- NULL
      ibsSingleQW <- NULL
      ibsColonyQW <- NULL
      ibdrQW <- NULL
      ibdeQW <- NULL
    }
  # Queens vs drones
  if (isFALSE(inbreeding)){
    print("IBSmulti QD")
    ibsMultiQD <- data.frame(Value = c(ibsMultiDF[idDF$drones, idDF$queen]),
                             Rel = "QD",
                             Type = "IBSmultiBF")

    print("IBSsingle QD")
    ibsSingleQD <- data.frame(Value = c(ibsSingleDF[idDF$drones, idDF$queen]),
                              Rel = "QD",
                              Type = "IBSsingleBF")

    if (!is.null(ibsColonyDF)){
      print("IBScolony QD")
      ibsColonyQD <- data.frame(Value = c(ibsColonyDF[idDF$drones, idDF$queen]),
                             Rel = "QD",
                             Type = "IBScolonyF")
    }

    print("IBDr QD")
    ibdrQD <- data.frame(Value = c(ibdDF[idDF$drones, idDF$queen]),
                         Rel = "QD",
                         Type = "IBDr")

    if (!is.null(Sinv)) {
      print("IBDe QD")
      ibdeQD <- data.frame(Value = getS(Sinv, ids = idDF$drones, with = idDF$queen, vector = TRUE),
                           Rel = "QD",
                           Type = "IBDe")
    }} else {
      ibsMultiQD <- NULL
      ibsSingleQD <- NULL
      ibsColonyQD <- NULL
      ibdrQD <- NULL
      ibdeQD <- NULL
    }
  return(if (is.null(ibsColonyDF)){
    bind_rows(
      ibsMultiWW, ibsMultiWD, ibsMultiDD, ibsMultiQD, ibsMultiQW,
      ibsSingleWW, ibsSingleWD, ibsSingleDD, ibsSingleQD, ibsSingleQW,
      ibsColonyWW = NULL, ibsColonyWD = NULL, ibsColonyDD =  NULL, ibsColonyQD = NULL, ibsColonyQW = NULL,
      ibdrWW, ibdrWD, ibdrDD, ibdrQD, ibdrQW,
      ibdeWW, ibdeWD, ibdeDD, ibdeQD, ibdeQW)
  } else {
    bind_rows(
      ibsMultiWW, ibsMultiWD, ibsMultiDD, ibsMultiQD, ibsMultiQW,
      ibsSingleWW, ibsSingleWD, ibsSingleDD, ibsSingleQD, ibsSingleQW,
      ibsColonyWW, ibsColonyWD, ibsColonyDD, ibsColonyQD, ibsColonyQW,
      ibdrWW, ibdrWD, ibdrDD, ibdrQD, ibdrQW,
      ibdeWW, ibdeWD, ibdeDD, ibdeQD, ibdeQW)
  })
}

prepareDataForPlotting_Queens <- function(ibsMultiDF = NULL, ibsSingleDF = NULL, ibsColonyDF = NULL, ibdDF = NULL, Sinv = NULL, idPopDF) {
  #Get populations IDs
  melID <- idPopDF$ID[idPopDF$Pop == "Mel"]
  melCrossID <- idPopDF$ID[idPopDF$Pop == "MelCross"]
  carID <- idPopDF$ID[idPopDF$Pop == "Car"]

  #Compute relationships between queens of different populations (QQ)
  #IBSmulti
  print("IBSmultiQQ")
  IBSmultiQQ <- rbind(data.frame(Value = as.vector(list(ibsMultiDF[melID, melCrossID])[[1]]),
                                 Pops = "Mel_MelCross"),
                      data.frame(Value = as.vector(list(ibsMultiDF[melID, carID])[[1]]),
                                 Pops = "Mel_Car"),
                      data.frame(Value = as.vector(list(ibsMultiDF[melCrossID, carID])[[1]]),
                                 Pops = "MelCross_Car"))
  IBSmultiQQ$Rel = "QQ"
  IBSmultiQQ$Type = "IBSmultiBF"

  #IBSsingle
  print("IBSsingleQQ")
  IBSsingleQQ <- rbind(data.frame(Value = as.vector(list(ibsSingleDF[melID, melCrossID])[[1]]),
                                  Pops = "Mel_MelCross"),
                       data.frame(Value = as.vector(list(ibsSingleDF[melID, carID])[[1]]),
                                  Pops = "Mel_Car"),
                       data.frame(Value = as.vector(list(ibsSingleDF[melCrossID, carID])[[1]]),
                                  Pops = "MelCross_Car"))
  IBSsingleQQ$Rel = "QQ"
  IBSsingleQQ$Type = "IBSsingleBF"

  #IBScolony
  if(is.null(ibsColonyDF)){
    IBScolonyQQ <- NULL
  } else{
    print("IBScolonyQQ")
    IBScolonyQQ <- rbind(data.frame(Value = as.vector(list(ibsColonyDF[melID, melCrossID])[[1]]),
                                 Pops = "Mel_MelCross"),
                      data.frame(Value = as.vector(list(ibsColonyDF[melID, carID])[[1]]),
                                 Pops = "Mel_Car"),
                      data.frame(Value = as.vector(list(ibsColonyDF[melCrossID, carID])[[1]]),
                                 Pops = "MelCross_Car"))
    IBScolonyQQ$Rel = "QQ"
    IBScolonyQQ$Type = "IBScolonyF"
  }

  #IBDr
  print("IBDrQQ")
  IBDrQQ <- rbind(data.frame(Value = as.vector(list(ibdDF[melID, melCrossID])[[1]]),
                             Pops = "Mel_MelCross"),
                  data.frame(Value = as.vector(list(ibdDF[melID, carID])[[1]]),
                             Pops = "Mel_Car"),
                  data.frame(Value = as.vector(list(ibdDF[melCrossID, carID])[[1]]),
                             Pops = "MelCross_Car"))
  IBDrQQ$Rel <- "QQ"
  IBDrQQ$Type <- "IBDr"


  #IBDeQQ
  if (!is.null(Sinv)) {
    print("IBDeQQ")
    IBDeQQ <- rbind(data.frame(Value = getS(Sinv, ids = melID, with = melCrossID, vector = TRUE),
                               Pops = "Mel_MelCross"),
                    data.frame(Value = getS(Sinv, ids = melID, with = carID, vector = TRUE),
                               Pops = "Mel_Car"),
                    data.frame(Value = getS(Sinv, ids = melCrossID, with = carID, vector = TRUE),
                               Pops = "MelCross_Car"))
    IBDeQQ$Rel <- "QQ"
    IBDeQQ$Type <- "IBDe"
  }

  # Compute relationships between queens of the same populations/ without diagonal (Q)
  #IBSmultiBF
  print("IBSmultiQ")
  IBSmultiQ <- rbind(data.frame(Value = as.vector(list(c(ibsMultiDF[melID, melID][lower.tri(ibsMultiDF[melID, melID],
                                                                                            diag = FALSE)]))[[1]]),
                                Pops = "Mel"),
                     data.frame(Value = as.vector(list(c(ibsMultiDF[carID, carID][lower.tri(ibsMultiDF[carID, carID],
                                                                                            diag = FALSE)]))[[1]]),
                                Pops = "Car"),
                     data.frame(Value = as.vector(list(c(ibsMultiDF[melCrossID, melCrossID][lower.tri(ibsMultiDF[melCrossID, melCrossID],
                                                                                                      diag = FALSE)]))[[1]]),
                                Pops = "MelCross"))
  IBSmultiQ$Rel <- "Q"
  IBSmultiQ$Type <- "IBSmultiBF"

  #IBSsingleBF
  print("IBSsingleQ")
  IBSsingleQ <- rbind(data.frame(Value = as.vector(list(c(ibsSingleDF[melID, melID][lower.tri(ibsSingleDF[melID, melID],
                                                                                              diag = FALSE)]))[[1]]),
                                 Pops = "Mel"),
                      data.frame(Value = as.vector(list(c(ibsSingleDF[carID, carID][lower.tri(ibsSingleDF[carID, carID],
                                                                                              diag = FALSE)]))[[1]]),
                                 Pops = "Car"),
                      data.frame(Value = as.vector(list(c(ibsSingleDF[melCrossID, melCrossID][lower.tri(ibsSingleDF[melCrossID, melCrossID],
                                                                                                        diag = FALSE)]))[[1]]),
                                 Pops = "MelCross"))
  IBSsingleQ$Rel <- "Q"
  IBSsingleQ$Type <- "IBSsingleBF"


  #IBSsingleBF
  if (is.null(ibsColonyDF)){
    IBScolonyQ <- NULL
  } else {
    print("IBScolonyQ")
    IBScolonyQ <- rbind(data.frame(Value = as.vector(list(c(ibsColonyDF[melID, melID][lower.tri(ibsColonyDF[melID, melID], diag = FALSE)]))[[1]]),
                                Pops = "Mel"),
                     data.frame(Value = as.vector(list(c(ibsColonyDF[carID, carID][lower.tri(ibsColonyDF[carID, carID], diag = FALSE)]))[[1]]),
                                Pops = "Car"),
                     data.frame(Value = as.vector(list(c(ibsColonyDF[melCrossID, melCrossID][lower.tri(ibsColonyDF[melCrossID, melCrossID], diag = FALSE)]))[[1]]),
                                Pops = "MelCross"))
    IBScolonyQ$Rel <- "Q"
    IBScolonyQ$Type <- "IBScolonyF"
  }


  #IBDr
  print("IBDrQ")
  IBDrQ <- rbind(data.frame(Value = as.vector(list(c(ibdDF[melID, melID][lower.tri(ibdDF[melID, melID], diag = FALSE)]))[[1]]),
                            Pops = "Mel"),
                 data.frame(Value = as.vector(list(c(ibdDF[carID, carID][lower.tri(ibdDF[carID, carID], diag = FALSE)]))[[1]]),
                            Pops = "Car"),
                 data.frame(Value = as.vector(list(c(ibdDF[melCrossID, melCrossID][lower.tri(ibdDF[melCrossID, melCrossID], diag = FALSE)]))[[1]]),
                            Pops = "MelCross"))



  IBDrQ$Rel = "Q"
  IBDrQ$Type = "IBDr"

  #IBDe
  if (!is.null(Sinv)) {
    print("IBDeQ")
    tmp1 <- getS(Sinv, melID, diagOnly = FALSE, vector = FALSE)
    tmp2 <- getS(Sinv, melCrossID, diagOnly = FALSE, vector = FALSE)
    tmp3 <- getS(Sinv, carID, diagOnly = FALSE, vector = FALSE)
    IBDeQ <- rbind(data.frame(Value = c(tmp1[lower.tri(tmp1, diag = FALSE)]),
                              Pops = "Mel"),
                   data.frame(Value = c(tmp2[lower.tri(tmp2, diag = FALSE)]),
                              Pops = "MelCross"),
                   data.frame(Value = c(tmp3[lower.tri(tmp3, diag = FALSE)]),
                              Pops = "Car"))
    IBDeQ$Rel <- "Q"
    IBDeQ$Type <- "IBDe"
  }



  # Inbreeding (diagonal!!!) if queens (F)
  #IBSfMulti
  print("IBSfMulti")
  IBSfMulti <- rbind(data.frame(Value = diag(ibsMultiDF[melID, melID]),
                                Pops = "Mel"),
                     data.frame(Value = diag(ibsMultiDF[melCrossID, melCrossID]),
                                Pops = "MelCross"),
                     data.frame(Value = diag(ibsMultiDF[carID, carID]),
                                Pops = "Car"))
  IBSfMulti$Type = "IBSmultiBF"
  IBSfMulti$Rel = "F"

  #IBSfSingle
  print("IBSfSingle")
  IBSfSingle <- rbind(data.frame(Value = diag(ibsSingleDF[melID, melID]),
                                 Pops = "Mel"),
                      data.frame(Value = diag(ibsSingleDF[melCrossID, melCrossID]),
                                 Pops = "MelCross"),
                      data.frame(Value = diag(ibsSingleDF[carID, carID]),
                                 Pops = "Car"))
  IBSfSingle$Type = "IBSsingleBF"
  IBSfSingle$Rel = "F"

  #IBSfColony
  if (is.null(ibsColonyDF)){
    IBSfColony <- NULL
  } else {
    print("IBSfColony")
    IBSfColony <- data.frame(rbind(data.frame(Value = diag(ibsColonyDF[melID, melID]),
                                           Pops = "Mel"),
                                data.frame(Value = diag(ibsColonyDF[melCrossID, melCrossID]),
                                           Pops = "MelCross"),
                                data.frame(Value = diag(ibsColonyDF[carID, carID]),
                                           Pops = "Car")))
    IBSfColony$Type <- "IBScolonyF"
    IBSfColony$Rel <- "F"
  }

  #IBDrF
  print("IBDrF")
  IBDrF <- rbind(data.frame(Value = diag(ibdDF[melID, melID]),
                            Pops = "Mel"),
                 data.frame(Value = diag(ibdDF[melCrossID, melCrossID]),
                            Pops = "MelCross"),
                 data.frame(Value = diag(ibdDF[carID, carID]),
                            Pops = "Car"))
  IBDrF$Type = "IBDr"
  IBDrF$Rel = "F"

  if (!is.null(Sinv)) {
    print("IBDeF")
    IBDeF <- rbind(data.frame(Value = getS(Sinv, melID, diagOnly = TRUE),
                              Pops = "Mel"),
                   data.frame(Value = getS(Sinv, melCrossID, diagOnly = TRUE),
                              Pops = "MelCross"),
                   data.frame(Value = getS(Sinv, carID, diagOnly = TRUE),
                              Pops = "Car"))
    IBDeF$Type = "IBDe"
    IBDeF$Rel = "F"
  }

  return(
    bind_rows(
      IBSmultiQQ, IBSsingleQQ, IBScolonyQQ, IBDrQQ, IBDeQQ,
      IBSmultiQ, IBSsingleQ, IBScolonyQ, IBDrQ, IBDeQ,
      IBSfMulti, IBSfSingle, IBSfColony, IBDrF, IBDeF)
  )
}

prepareDataForPlottingHeatMap_Queens <- function(ibsMultiDF = NULL, ibsSingleDF = NULL, ibsColonyDF = NULL, ibdDF = NULL, Sinv = NULL, idDF = NULL) {
  #multi
  ibsMultiDF <- as.data.frame(ibsMultiDF)
  columns <- colnames(ibsMultiDF)
  ibsMultiDF$ID <- rownames(ibsMultiDF)
  ibsMultiDFL <- ibsMultiDF %>% pivot_longer(cols = all_of(columns))
  ibsMultiDFL$Type <- "IBSmultiBF"

  #single
  ibsSingleDF <- as.data.frame(ibsSingleDF)
  columns <- colnames(ibsSingleDF)
  ibsSingleDF$ID <- rownames(ibsSingleDF)
  ibsSingleDFL <- ibsSingleDF %>% pivot_longer(cols = all_of(columns))
  ibsSingleDFL$Type <- "IBSsingleBF"

  #own
  if (!is.null(ibsColonyDF)) {
  ibsColonyDF <- as.data.frame(ibsColonyDF)
  columns <- colnames(ibsColonyDF)
  ibsColonyDF$ID <- rownames(ibsColonyDF)
  ibsColonyDFL <- ibsColonyDF %>% pivot_longer(cols = all_of(columns))
  ibsColonyDFL$Type <- "IBScolonyF"
  }

  #IBDr
  ibdrDF <- as.data.frame(ibdDF)
  columns <- colnames(ibdrDF)
  ibdrDF$ID <- rownames(ibdrDF)
  ibdrDFL <- ibdrDF %>% pivot_longer(cols = all_of(columns))
  ibdrDFL$Type <- "IBDr"

  #IBDe
  if (!is.null(Sinv)) {
    ibdeDF <-  as.data.frame(as.matrix(getS(Sinv, ids = idDF)))
    rownames(ibdeDF) <- idDF
    colnames(ibdeDF) <- idDF
    columns <- colnames(ibdeDF)
    ibdeDF$ID <- as.character(idDF)
    ibdeDFL <- ibdeDF %>% pivot_longer(cols = all_of(columns))
    ibdeDFL$Type <- "IBDe"
  }
  return(rbind(
         ibsMultiDFL, ibsSingleDFL, ibsColonyDFL, ibdrDFL, ibdeDFL))
}

plotColony <- function(df, rel = NULL, type = NULL, years = NULL) {
  df$Type <- factor(df$Type, levels = c("IBDe", "IBDr", "IBScolonyF", "IBSsingleBF",  "IBSmultiBF"))
  df$Rel <- factor(df$Rel, levels= c("WW", "WD", "DD", "QW", "QD"))

  type_labels <- c("IBDe", "IBDr", "IBScolonyF", "IBSsingleBF",  "IBSmultiBF")
  names(type_labels) <- c("IBDe", "IBDr", "IBScolonyF", "IBSsingleBF",  "IBSmultiBF")

  df$Year <- factor(df$Year, levels = years)
  year_labels <- unlist(paste("Year", years, collapse = "_") %>% stringr::str_split("_"))
  names(year_labels) <- years

  plot <- ggplot(df[df$Rel %in% rel & df$Type %in% type, ],
              aes(x = Value, fill = Rel)) +
          geom_vline(xintercept = c(0, 0.25, 0.5, 0.75), linewidth = 0.25, colour = "grey") +
          geom_histogram(binwidth = 0.01, position = "identity") +
          facet_grid2(rows = vars(Type), cols = vars(Year), scales = "free_y",independent = "y", labeller = labeller(Type = type_labels, Year = year_labels)) +
          scale_fill_manual("", values = paletteViridis, aesthetics = c("colour","fill")) +
          theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "top") +
          scale_x_continuous(breaks=seq(-2, 10, 0.25))
  return(plot)
}

plotQueens <- function(df, rel = NULL, type = NULL, pops = NULL, years = NULL, palette = NULL) { #Use palette = cbPaletteQ with same pop and palette = cbPaletteQQ with different pops
  df$Type <- factor(df$Type, levels = c("IBDe", "IBDr", "IBScolonyF", "IBSsingleBF",  "IBSmultiBF"))
  df$Pops <- factor(df$Pops, levels = c("Mel_MelCross", "Mel_Car", "MelCross_Car", "Mel", "Car", "MelCross"))
  df$Rel  <- factor(df$Rel, levels = c("QQ", "Q", "F"))

  #TODO new type labels required?
  type_labels <- c("IBDe", "IBDr", "IBScolonyF", "IBSsingleBF",  "IBSmultiBF")
  names(type_labels) <- c("IBDe", "IBDr", "IBScolonyF", "IBSsingleBF",  "IBSmultiBF")

  df$Year <- factor(df$Year, levels = years)
  year_labels <- unlist(paste("Year", years, collapse = "_") %>% stringr::str_split("_"))
  names(year_labels) <- years


  plot <- ggplot(df[df$Rel %in% rel & df$Type %in% type & df$Pops %in% pops, ],
                 aes(x = Value, fill = Pops)) +
          geom_vline(xintercept = c(0, 0.25, 0.5, 0.75), linewidth = 0.25, colour = "grey") +  #remove line if minimum is larger that 0.75 or change the xintercept
          geom_histogram(binwidth = 0.01, position = "identity") +
          facet_grid2(rows = vars(Type), cols = vars(Year), scales = "free_y", independent = "y", labeller = labeller(Type = type_labels, Year = year_labels)) +
          scale_fill_manual("", values = palette, aesthetics = c("colour","fill")) +
          theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "top") +
          scale_x_continuous(breaks=seq(-3, 8, 0.25))
  return(plot)
}

scatterQueens <-  function(df, rel = NULL, type = NULL, pops = NULL, legend.position = NULL, palette = NULL, years = NULL) {
  df$Type <- factor(df$Type, levels = c("IBDe", "IBDr", "IBScolonyF", "IBSsingleBF",  "IBSmultiBF"))
  df$Pops <- factor(df$Pops, levels = c("Mel_MelCross", "Mel_Car", "MelCross_Car", "Car", "Mel", "MelCross"))
  df$Rel  <- factor(df$Rel, levels = c("QQ", "Q"))


  type_labels <- c("IBDe", "IBDr", "IBScolonyF", "IBSsingleBF",  "IBSmultiBF")
  names(type_labels) <- c("IBDe", "IBDr", "IBScolonyF", "IBSsingleBF",  "IBSmultiBF")


  df$Year <- factor(df$Year, levels = years)
  year_labels <- unlist(paste("Year", years, collapse = "_") %>% stringr::str_split("_"))
  names(year_labels) <- years

  a <- pivot_wider(df[df$Rel %in% rel & df$Type %in% type & df$Pops %in% pops, ], names_from = Type, values_from = Value, values_fn = list)
  b <- unnest(a, cols = all_of(type))

  plot <- ggplot(data = b, aes(x = IBDr, y = IBSmultiBF)) +
          geom_point(aes(colour = Pops)) +
          facet_grid2(cols = vars(Year), scales = "free", independent = "y", labeller = labeller(Type = type_labels, Year = year_labels)) +
          scale_fill_manual("", values = palette, aesthetics = c("colour","fill")) +
          theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = legend.position) +
          scale_x_continuous(breaks=seq(-2, 2, 0.25))
  return(plot)
}

plotQueens_heatmapBIND <- function(df, PopIdDF = NULL, years = NULL) {
  df$ID <- as.factor(as.numeric(df$ID))
  df$name <- as.factor(as.numeric(df$name))
  df$Type <- factor(df$Type, levels = c("IBDe", "IBDr", "IBScolonyF", "IBSsingleBF",  "IBSmultiBF"))
  type_labels <- c("IBDe", "IBDr", "IBScolonyF", "IBSsingleBF",  "IBSmultiBF")
  names(type_labels) <- c("IBDe", "IBDr", "IBScolonyF", "IBSsingleBF",  "IBSmultiBF")

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
    facet_grid2(rows = vars(Type), cols = vars(Year), labeller = labeller(Type = type_labels, Year = year_labels), scales = "free", independent = "y")

  return(plot)
}

plotQueens_heatmapSOLO <- function(df, Pop = FALSE, PopIdDF = NULL, legend.position = NULL, years = NULL) {
  df$ID <- as.factor(as.numeric(df$ID))
  df$name <- as.factor(as.numeric(df$name))
  df$Type <- factor(df$Type, levels = c("IBDe", "IBDr", "IBScolonyF", "IBSsingleBF",  "IBSmultiBF"))
  type_labels <- c("IBDe", "IBDr", "IBScolonyF", "IBSsingleBF",  "IBSmultiBF")
  names(type_labels) <- c("IBDe", "IBDr", "IBScolonyF", "IBSsingleBF",  "IBSmultiBF")

  df$Year <- factor(df$Year, levels = years)
  year_labels <- unlist(paste("Year", years, collapse = "_") %>% stringr::str_split("_"))
  names(year_labels) <- years

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
      facet_grid(rows = vars(Type))
  } else {
    plot <- ggplot(data = df, aes(x = ID, y = name, fill = value)) + geom_tile() +
      facet_grid2(rows = vars(Type), labeller = labeller(Type = type_labels), scales = "free")
  }
  return(plot)
}

########################################################
### --- FIGURE 1&2: Pure subspecies (carnica) in years 1 and 10 ---###
########################################################
print("Plot carnica year 1")
load("~/Desktop/GitHub/lstrachan_honeybee_sim/YearCycleSimulation/PlottingData/dataCar.RData")

#Coded out data prep
{
# relCar1 <- prepareDataForPlotting_Colony(ibsMultiDF = colonyCar1$IBSmultiBF,
#                                          ibsSingleDF = colonyCar1$IBSsingleBF,
#                                          ibsColonyDF = colonyCar1$IBScolonyF,
#                                          ibdDF = colonyCar1$IBD,
#                                          Sinv = Sinv,
#                                          idDF = colonyCar1$ID,
#                                          inbreeding = FALSE)
# 
# 
# 
# print("Plot carnica year 10")
# relCar10 <- prepareDataForPlotting_Colony(ibsMultiDF = colonyCar10$IBSmultiBF,
#                                           ibsSingleDF = colonyCar10$IBSsingleBF,
#                                           ibsColonyDF = colonyCar10$IBScolonyF,
#                                           ibdDF = colonyCar10$IBD,
#                                           Sinv = Sinv,
#                                           idDF = colonyCar10$ID,
#                                           inbreeding = FALSE)
# 
# relCar1$Year <- 1
# relCar10$Year <- 10
# dataCar <- rbind(relCar1, relCar10)
# save(... = dataCar, file = "dataCar.RData")
}

CarWWplot<- plotColony(dataCar, type = c("IBDe", "IBDr",  "IBScolonyF", "IBSsingleBF", "IBSmultiBF"), rel = c("WW", "WD", "DD"), years = c(1,10))

CarQWplot <- plotColony(dataCar, type = c("IBDe", "IBDr",  "IBScolonyF", "IBSsingleBF", "IBSmultiBF"),  rel = c("QW", "QD"), years = c(1,10))

#Plot Car inbreeding
load("~/Desktop/GitHub/lstrachan_honeybee_sim/YearCycleSimulation/PlottingData/dataCarF.RData")

# Coded out prep
{
# relCar1F <- prepareDataForPlotting_Colony(ibsMultiDF = colonyCar1$IBSmultiBF,
#                                           ibsSingleDF = colonyCar1$IBSsingleBF,
#                                           ibsColonyDF = colonyCar1$IBScolonyF,
#                                           ibdDF = colonyCar1$IBD,
#                                           Sinv = Sinv,
#                                           idDF = colonyCar1$ID,
#                                           inbreeding = TRUE)
# 
# relCar10F <- prepareDataForPlotting_Colony(ibsMultiDF = colonyCar10$IBSmultiBF,
#                                            ibsSingleDF = colonyCar10$IBSsingleBF,
#                                            ibsColonyDF = colonyCar10$IBScolonyF,
#                                            ibdDF = colonyCar10$IBD,
#                                            Sinv = Sinv,
#                                            idDF = colonyCar10$ID,
#                                            inbreeding = TRUE)
# relCar1F$Year <- 1
# relCar10F$Year <- 10
# 
# dataCarF <- rbind(relCar1F, relCar10F)
# save(dataCarF, file = "dataCarF.Rdata")
}

CarFplot <- plotColony(dataCarF, type =  c("IBDe", "IBDr",  "IBScolonyF", "IBSsingleBF", "IBSmultiBF"), rel = c("WW", "WD", "DD"), years = c(1,10))


#Plot CAR csd Locus
load("~/Desktop/GitHub/lstrachan_honeybee_sim/YearCycleSimulation/PlottingData/dataCarCSDloc.RData")
#Coded out prep
{
# relCar1_csd <- prepareDataForPlotting_Colony(ibsMultiDF = colonyCar1$IBSmultiBFCsdLoci,
#                                              ibsSingleDF = colonyCar1$IBSsingleBFcsdLoci,
#                                              ibdDF = colonyCar1$IBDCsdLoci,
#                                              Sinv = Sinv,
#                                              idDF = colonyCar1$ID)
# 
# 
# relCar10_csd <- prepareDataForPlotting_Colony(ibsMultiDF = colonyCar10$IBSmultiBFCsdLoci,
#                                                  ibsSingleDF = colonyCar10$IBSsingleBFcsdLoci,
#                                                  ibdDF = colonyCar10$IBDCsdLoci,
#                                                  Sinv = Sinv,
#                                                  idDF = colonyCar10$ID)
# relCar1_csd$Year <- 1
# relCar10_csd$Year <- 10
# dataCarCSDloc <- rbind(relCar1_csd, relCar10_csd)
# save(dataCarCSDloc, file = "dataCarCSDloc.Rdata")
}

CarCSDlocplot <- plotColony(dataCarCSDloc, type = c("IBDe", "IBDr", "IBSsingleBF", "IBSmultiBF"), rel = c("WW", "WD", "DD"), years = c(1,10))

#Plot CAR csd Chromosome
load("~/Desktop/GitHub/lstrachan_honeybee_sim/YearCycleSimulation/PlottingData/dataCarCsdChr.RData")
{
# relCar1_csdChr <- prepareDataForPlotting_Colony(ibsMultiDF = colonyCar1$IBSmultiCSDChr,
#                                                    ibsSingleDF = colonyCar1$IBSsingleCSDChr,
#                                                    ibdDF = colonyCar1$IBDcsdChr,
#                                                    Sinv = Sinv,
#                                                    idDF = colonyCar1$ID)
# relCar10_csdChr <- prepareDataForPlotting_Colony(ibsMultiDF = colonyCar10$IBSmultiCSDChr,
#                                                     ibsSingleDF = colonyCar10$IBSsingleCSDChr,
#                                                     ibdDF = colonyCar10$IBDcsdChr,
#                                                     Sinv = Sinv,
#                                                     idDF = colonyCar10$ID)
# relCar1_csdChr$Year <- 1
# relCar10_csdChr$Year <- 10
# 
# dataCarCSDchr <- rbind(relCar1_csdChr, relCar10_csdChr)
# save(... = dataCarCSDchr, file = "dataCarCsdChr.RData")
# 
# #  #filtering
# # relCar10_csdChr <- filter(relCar10_csdChr, Value < 2)
}

CarCSDchr <- plotColony(dataCarCSDchr, rel = c("WW", "WD", "DD"), type = c("IBDe", "IBDr", "IBSsingleBF", "IBSmultiBF"), years = c(1,10))
########################################################
### --- FIGURE 4: Between queens of different subspecies (carnica vs. mellifera) ---###
########################################################
#Plot queens Year 1
print("Plot queens")
load("~/Desktop/GitHub/lstrachan_honeybee_sim/YearCycleSimulation/PlottingData/dataQueens.RData")

#Coded out prep
{
# relQueens1 <- prepareDataForPlotting_Queens(ibsMultiDF =  queens1$IBSmultiBF,
#                                             ibsSingleDF = queens1$IBSsingleBF,
#                                             ibsColonyDF = queens1$IBScolonyF,
#                                             ibdDF = queens1$IBD,
#                                             Sinv = Sinv,
#                                             idPopDF = springerQueensPop1)
# 
# relQueens10 <-  prepareDataForPlotting_Queens(ibsMultiDF =  queens10$IBSmultiBF,
#                                               ibsSingleDF = queens10$IBSsingleBF,
#                                               ibsColonyDF = queens10$IBScolonyF,
#                                               ibdDF = queens10$IBD,
#                                               Sinv = Sinv,
#                                               idPopDF = springerQueensPop10)
# relQueens1$Year <- 1
# relQueens10$Year <- 10
# 
# dataQueens <- rbind(relQueens1, relQueens10)
# save(relQueens10, file = "relQueens10.Rdata")
# save(dataQueens, file = "dataQueens.Rdata")
}

print("Between populations ")
#between populations
QueensQQ <- plotQueens(dataQueens, rel = "QQ", type = c("IBDe", "IBDr", "IBSmultiBF"), pops = c("Mel_MelCross", "Mel_Car", "MelCross_Car"), years = c(1,10), palette = cbPaletteQQ)
QueensQQ

print("Within populations")
#within population - non-diagonal
QueensQ <- plotQueens(dataQueens, rel = "Q", type = c("IBDe", "IBDr","IBSmultiBF"), pops = c("Mel", "Car", "MelCross"), years = c(1,10), palette = cbPaletteQ)
QueensQ
#within population - diagonal
QueensF <- plotQueens(dataQueens, rel = "F", type = c("IBDe", "IBDr","IBSmultiBF"), pops = c("Mel", "Car", "MelCross"), years = c(1,10), palette = cbPaletteQ)
QueensF

#scatter plot looking at Q and QQ
load("~/Desktop/GitHub/lstrachan_honeybee_sim/YearCycleSimulation/PlottingData/relQueens10.RData")

QueensScatterQQ <- scatterQueens(relQueens10, type = c("IBDr", "IBSmultiBF"), rel = "QQ", pops = c("Mel_MelCross", "Mel_Car", "MelCross_Car"), palette = cbPaletteQQ, years = 10 , legend.position = "top")
QueensScatterQ <- scatterQueens(relQueens10, type = c("IBDr", "IBSmultiBF"),  rel = "Q", pops = c("Mel", "Car", "MelCross"), palette = cbPaletteQ, years= 10, legend.position = "top")
QueensScatter <- grid.arrange(QueensScatterQ, QueensScatterQQ, ncol = 2)



print("HeatMaps")
load("~/Desktop/GitHub/lstrachan_honeybee_sim/YearCycleSimulation/PlottingData/dataQheat.RData")

#Coded out prep
{
# relQueens1h <- prepareDataForPlottingHeatMap_Queens(ibsMultiDF =  queens1$IBSmultiBF,
#                                                     ibsSingleDF = queens1$IBSsingleBF,
#                                                     ibsColonyDF = queens1$IBScolonyF,
#                                                     ibdDF = queens1$IBD,
#                                                     Sinv = Sinv,
#                                                     idDF = queens1$ID)
# 
# relQueens10h <- prepareDataForPlottingHeatMap_Queens(ibsMultiDF =  queens10$IBSmultiBF,
#                                                      ibsSingleDF = queens10$IBSsingleBF,
#                                                      ibsColonyDF = queens10$IBScolonyF,
#                                                      ibdDF = queens10$IBD,
#                                                      Sinv = Sinv,
#                                                      idDF = queens10$ID)
# 
# relQueens1h$Year <- 1
# relQueens10h$Year <- 10
# dataQh <- rbind(relQueens1h, relQueens10h)
# save(dataQh, file = "dataQHeat.Rdata")
}

idPopQueens <- rbind(springerQueensPop1, springerQueensPop10)
relQueens1h$Year <- 1
relQueens10h$Year <- 10

plotQueens1h <- plotQueens_heatmapSOLO(relQueens1h, Pop = TRUE, PopIdDF = springerQueensPop1, legend.position = "left")
plotQueens10h <- plotQueens_heatmapSOLO(relQueens10h, Pop = TRUE, PopIdDF = springerQueensPop10, legend.position = "right")
arrangedPlots <- grid.arrange(plotQueens1h, plotQueens10h, ncol = 2)

QueenH <- plotQueens_heatmapBIND(dataQh, PopIdDF = idPopQueens, years = c(1,10))
QueenH

#CSD STUFF
#Plot csd Loc as a histogram (QQ)
load("~/Desktop/GitHub/lstrachan_honeybee_sim/YearCycleSimulation/PlottingData/dataQueensCSDloc.RData")

{
# relQueens1_csdLoc_hist <- prepareDataForPlotting_Queens(ibsMultiDF = springerQueens1$IBSmultiBFCsdLoci,
#                                                         ibsSingleDF = springerQueens1$IBSsingleBFcsdLoci,
#                                                         ibsColonyDF = NULL,
#                                                         ibdDF = springerQueens1$IBDCsdLoci,
#                                                         Sinv = Sinv,
#                                                         idPopDF = springerQueensPop1)
# 
# relQueens10_csdLoc_hist <- prepareDataForPlotting_Queens(ibsMultiDF = springerQueens10$IBSmultiBFCsdLoci,
#                                                             ibsSingleDF = springerQueens10$IBSsingleBFcsdLoci,
#                                                             ibsColonyDF = NULL,
#                                                             ibdDF = queens10$IBDCsdLoci,
#                                                             Sinv = Sinv,
#                                                             idPopDF = springerQueensPop10)
# relQueens1_csdLoc_hist$Year <- 1
# relQueens10_csdLoc_hist$Year <- 10
# dataQueens_csdLoc <- rbind(relQueens1_csdLoc_hist, relQueens10_csdLoc_hist)
# save(dataQueens_csdLoc, file = "dataQueensCSDloc.RData")
}

Queens_csdLoc_hist <- plotQueens(dataQueens_csdLoc, rel = "QQ",  type = c("IBDe", "IBDr", "IBSmultiBF", "IBSsingleBF"), pops = c("Mel_MelCross", "Mel_Car", "MelCross_Car"), years = c(1,10), palette = cbPaletteQQ) 
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
load("~/Desktop/GitHub/lstrachan_honeybee_sim/YearCycleSimulation/PlottingData/dataQueensCSDchr.RData")
{
# relQueen1_csdChr_hist <- prepareDataForPlotting_Queens(ibsMultiDF = springerQueens1$IBSmultiCSDChr,
#                                                         ibsSingleDF = springerQueens1$IBSsingleCSDChr,
#                                                         ibsColonyDF = NULL,
#                                                         ibdDF = springerQueens1$IBDcsdChr,
#                                                         Sinv = Sinv,
#                                                         idPopDF = springerQueensPop1)
# 
# relQueen10_csdChr_hist <- prepareDataForPlotting_Queens(ibsMultiDF = springerQueens10$IBSmultiCSDChr,
#                                                         ibsSingleDF = springerQueens10$IBSsingleCSDChr,
#                                                         ibsColonyDF = NULL,
#                                                         ibdDF = springerQueens10$IBDcsdChr,
#                                                         Sinv = Sinv,
#                                                         idPopDF = springerQueensPop10)
# 
# relQueen1_csdChr_hist$Year <- 1
# relQueen10_csdChr_hist$Year <- 10
# dataQueensCSDchr <- rbind(relQueen1_csdChr_hist, relQueen10_csdChr_hist)
# save(dataQueensCSDchr, file = "dataQueensCSDchr.RData")
}

Queens_csdChr_hist <- plotQueens(dataQueensCSDchr, rel = "QQ",  type = c("IBDe", "IBDr", "IBSmultiBF", "IBSsingleBF"), pops = c("Mel_MelCross", "Mel_Car", "MelCross_Car"), years = c(1,10), palette = cbPaletteQQ)


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


