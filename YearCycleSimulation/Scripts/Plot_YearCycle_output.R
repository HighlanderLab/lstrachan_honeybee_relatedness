rm(list = ls())

#Load libraries
pacman::p_load(tidyverse, Matrix, SIMplyBee, gridExtra, ggh4x, viridisLite)


print("Reading in the data")
data <- load("SpringerSimulation_import.RData") #file created when running YearlyCycle script
Sinv <- IBDe$Sinv

# The data contains two populations - mellifera and carnica
# The carnica stays "pure" throughout the simulation
# Mellifera gets mated with a proportion of carnica drones

#Plotting help
# The colourblind palettes:
cbPaletteQ <- c("#FF6DB6", "#490092", "#6DB6FF")  #Use within the same population
cbPaletteQQ <- c("#993F00", "#FF8E32", "#B6DAFF")  #Use on different populations
paletteViridis <- plasma(n = 3, begin = 0.3, end = 0.9) #Use on Carnica population


print("Assigning objects")
ped <- pedigree
caste <- caste

#Assign Car colony and queens
colonyCar1 <- springerColony1_Car
colonyCar10 <- springerColony10_Car
queens1_MelAF <- springerQueens1_MelAF
queens1_CarAF <- springerQueens1_CarAF
queens10_MelAF <- springerQueens10_MelAF
queens10_CarAF <- springerQueens10_CarAF


getS <- function(Sinv, ids, with = ids, diagOnly = FALSE, vector = FALSE) {
  ids <- as.numeric(ids)
  with <- as.numeric(with)
  x <- sparseMatrix(i = ids, j = 1:length(ids), dims = c(nrow(Sinv), length(ids)))
  M1 <- as(x, "dMatrix")
  Sids <- as.matrix(solve(Sinv, M1)[with,])
  if (dim(Sids)[1] == length(ids)){
    rownames(Sids) <- ids
  } else {
    rownames(Sids) <- with
  }

  if (dim(Sids)[2] == length(ids)){
    colnames(Sids) <- ids
  } else {
    colnames(Sids) <- with
    }
  if (diagOnly) {
    Sids <- diag(Sids)
  }
  if (vector) {
    Sids <- c(as.matrix(Sids))
  }
  return(Sids)
}

determine_sister_type = function(x , df = NULL, WorkersFatherTable = NULL){
  if ((WorkersFatherTable$DPQ[WorkersFatherTable$workers == df$id1[x]] == WorkersFatherTable$DPQ[WorkersFatherTable$workers == df$id2[x]]) &
      (WorkersFatherTable$fathers[WorkersFatherTable$workers == df$id1[x]] != WorkersFatherTable$fathers[WorkersFatherTable$workers == df$id2[x]])){
    ret = "FS"
  } else if (WorkersFatherTable$fathers[WorkersFatherTable$workers == df$id1[x]] == WorkersFatherTable$fathers[WorkersFatherTable$workers == df$id2[x]]){
    ret = "SS"
  } else {
    ret = "HS"
  }
  return(ret)
}

SisterTypeDF <- function(x){ #Doesn't work with merged dataframes, Works with only relCar1/relCar10 then use rbind to merge them
  Ibdr <- x[x$Rel == "WW" & x$Type == "IBDr", ] %>%
          group_by(SisterType) %>%
          reframe(Mean = mean(Value), SD = sd(Value), Type = "IBDr", Year = unique(x$Year))

  Ibde <- x[x$Rel == "WW" & x$Type == "IBDe", ] %>%
          group_by(SisterType) %>%
          reframe(Mean = mean(Value), SD = sd(Value), Type = "IBDe", Year = unique(x$Year))

  IBSmulti <- x[x$Rel == "WW" & x$Type == "IBSmultiAF", ] %>%
              group_by(SisterType) %>%
              reframe(Mean = mean(Value), SD = sd(Value), Type = "IBSmultiAF", Year = unique(x$Year))

  IBSAF0.5 <- x[x$Rel == "WW" & x$Type == "IBSAF0.5", ] %>%
              group_by(SisterType) %>%
              reframe(Mean = mean(Value), SD = sd(Value), Type = "IBSAF0.5", Year = unique(x$Year))

  IBSsingleAF <- x[x$Rel == "WW" & x$Type == "IBSsingleAF", ] %>%
                 group_by(SisterType) %>%
                 reframe(Mean = mean(Value), SD = sd(Value), Type = "IBSsingleAF", Year = unique(x$Year))

  IBScolonyAF <- x[x$Rel == "WW" & x$Type == "IBScolonyAF", ] %>%
                 group_by(SisterType) %>%
                 reframe(Mean = mean(Value), SD = sd(Value), Type = "IBScolonyAF", Year = unique(x$Year))

  df <- list(Ibdr, Ibde, IBSmulti, IBSAF0.5, IBSsingleAF, IBScolonyAF) %>%
        Reduce(function(x, y) merge(x, y, all=TRUE), .)

  return(df)
}

# Plotting functions
prepareDataForPlotting_Colony <- function(ibsMultiDF = NULL, ibsSingleDF = NULL, ibsColonyDF = NULL, ibsAF0.5DF = NULL, ibdDF = NULL, Sinv = NULL, idDF, inbreeding = FALSE, WorkersFatherTable = NULL) {

  #Workers on workers

    print("IBSAF0.5 WW")
    if (inbreeding) {
      ibsAF0.5WW <-  data.frame(Value = diag(ibsAF0.5DF[idDF$workers, idDF$workers]),
                                Rel = "WW",
                                Type = "IBSAF0.5")
    } else {
      x <- ibsAF0.5DF[idDF$workers, idDF$workers]
      ind <- which(lower.tri(x , diag = FALSE) , arr.ind = TRUE )
      ibsAF0.5WW <- data.frame(id1 = dimnames(x)[[2]][ind[,2]] ,
                               id2 = dimnames(x)[[1]][ind[,1]] ,
                               Value = x[ind],
                               Rel = "WW",
                               Type = "IBSAF0.5")
    }
    ibsAF0.5WW$SisterType <- sapply(1:nrow(ibsAF0.5WW), FUN = function(z) determine_sister_type(x = z, df = ibsAF0.5WW, WorkersFatherTable = WorkersFatherTable))



    print("IBSmulti WW")
    if (inbreeding) {
      ibsMultiWW <-  data.frame(Value = diag(ibsMultiDF[idDF$workers, idDF$workers]),
                                Rel = "WW",
                                Type = "IBSmultiAF")
    } else {
      x <- ibsMultiDF[idDF$workers, idDF$workers]
      ind <- which(lower.tri(x , diag = FALSE) , arr.ind = TRUE )
      ibsMultiWW <- data.frame(id1 = dimnames(x)[[2]][ind[,2]] ,
                               id2 = dimnames(x)[[1]][ind[,1]] ,
                               Value = x[ind],
                               Rel = "WW",
                               Type = "IBSmultiAF")
    }
    ibsMultiWW$SisterType <- sapply(1:nrow(ibsMultiWW), FUN = function(x) determine_sister_type(x = x, df = ibsMultiWW, WorkersFatherTable = WorkersFatherTable))


  print("IBSsingle WW")
  if (inbreeding) {
    ibsSingleWW <-  data.frame(Value = diag(ibsSingleDF[idDF$workers, idDF$workers]),
                              Rel = "WW",
                              Type = "IBSsingleAF")
  } else {
    x <- ibsSingleDF[idDF$workers, idDF$workers]
    ind <- which(lower.tri(x , diag = FALSE) , arr.ind = TRUE )
    ibsSingleWW <- data.frame(id1 = dimnames(x)[[2]][ind[,2]] ,
                              id2 = dimnames(x)[[1]][ind[,1]] ,
                              Value = x[ind],
                              Rel = "WW",
                              Type = "IBSsingleAF")
  }
  ibsSingleWW$SisterType <- sapply(1:nrow(ibsSingleWW), FUN =  function(x) determine_sister_type(x = x, df = ibsSingleWW, WorkersFatherTable = WorkersFatherTable))

  if(!is.null(ibsColonyDF)){ #if dataframe is CSD info then ibdOwnDF is NULL
    print("IBScolony WW")
    if (inbreeding) {
      ibsColonyWW <-  data.frame(Value = diag(ibsColonyDF[idDF$workers, idDF$workers]),
                                 Rel = "WW",
                                 Type = "IBScolonyAF")
    } else {
      x <- ibsColonyDF[idDF$workers, idDF$workers]
      ind <- which(lower.tri(x , diag = FALSE) , arr.ind = TRUE )
      ibsColonyWW <- data.frame(id1 = dimnames(x)[[2]][ind[,2]] ,
                                id2 = dimnames(x)[[1]][ind[,1]] ,
                                Value = x[ind],
                                Rel = "WW",
                                Type = "IBScolonyAF")
    }
    ibsColonyWW$SisterType <- sapply(1:nrow(ibsColonyWW), FUN =  function(x) determine_sister_type(x = x, df = ibsColonyWW, WorkersFatherTable = WorkersFatherTable))
  }

  print("IBDr WW")
  if (inbreeding) {
    ibdrWW <-  data.frame(Value = diag(ibdDF[idDF$workers, idDF$workers]),
                               Rel = "WW",
                               Type = "IBDr")
  } else {
    x <- ibdDF[idDF$workers, idDF$workers]
    ind <- which(lower.tri(x , diag = FALSE) , arr.ind = TRUE)
    ibdrWW <- data.frame(id1 = dimnames(x)[[2]][ind[,2]] ,
                              id2 = dimnames(x)[[1]][ind[,1]] ,
                              Value = x[ind],
                              Rel = "WW",
                              Type = "IBDr")
  }
  ibdrWW$SisterType <- sapply(1:nrow(ibdrWW), FUN =  function(x) determine_sister_type(x = x, df = ibdrWW, WorkersFatherTable = WorkersFatherTable))

  if (!is.null(Sinv)) {
    print("IBDe WW")
    if (inbreeding){
      ibdeWW <- data.frame(Value = getS(Sinv, ids = idDF$workers, vector = TRUE, diagOnly = TRUE),
                           Rel = "WW",
                           Type = "IBDe")
    } else {
      x <- getS(Sinv, ids = idDF$workers, vector = FALSE)
      ind <- which(lower.tri(x, diag = FALSE), arr.ind = TRUE)
      ibdeWW <- data.frame(id1 = dimnames(x)[[2]][ind[,2]] ,
                           id2 = dimnames(x)[[1]][ind[,1]] ,
                           Value = x[ind],
                           Rel = "WW",
                           Type = "IBDe")
    }
     ibdeWW$SisterType <- sapply(1:nrow(ibdeWW), FUN = function(x) determine_sister_type(x = x, df = ibdeWW, WorkersFatherTable = WorkersFatherTable))

  }
  # Workers vs drones
  if (isFALSE(inbreeding)) {

      print("IBSAF0.5 WD")
      ibsAF0.5WD <- data.frame(Value = c(ibsAF0.5DF[idDF$workers, idDF$drones]),
                               Rel = "WD",
                               Type = "IBSAF0.5")
    print("IBSmulti WD")
    ibsMultiWD <- data.frame(Value = c(ibsMultiDF[idDF$workers, idDF$drones]),
                             Rel = "WD",
                             Type = "IBSmultiAF")


    print("IBSsingle WD")
    ibsSingleWD <- data.frame(Value =  c(ibsSingleDF[idDF$workers, idDF$drones]),
                              Rel = "WD",
                              Type = "IBSsingleAF")

    if (!is.null(ibsColonyDF)){
      print("IBScolony WD")
      ibsColonyWD <- data.frame(Value = c(ibsColonyDF[idDF$workers, idDF$drones]),
                             Rel = "WD",
                             Type = "IBScolonyAF")
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
      ibsAF0.5WD <- NULL
      ibsSingleWD <- NULL
      ibsColonyWD <- NULL
      ibdrWD <- NULL
      ibdeWD <- NULL
    }

  #Drones vs drones
    print("IBSAF0.5 DD")
    start_time <- Sys.time()
    ibsAF0.5DD <- data.frame(Value = ifelse(inbreeding,
                                            data.frame(diag(ibsAF0.5DF[idDF$drones, idDF$drones])),
                                            ibsAF0.5DF[idDF$drones, idDF$drones] %>% .[lower.tri(., diag = FALSE)] %>% data.frame(.)),
                             Rel = "DD",
                             Type = "IBSAF0.5")
    colnames(ibsAF0.5DD) <- c("Value", "Rel", "Type")

  print("IBSmulti DD")
  start_time <- Sys.time()
  ibsMultiDD <- data.frame(Value = ifelse(inbreeding,
                                          data.frame(diag(ibsMultiDF[idDF$drones, idDF$drones])),
                                          ibsMultiDF[idDF$drones, idDF$drones] %>% .[lower.tri(., diag = FALSE)] %>% data.frame(.)),
                           Rel = "DD",
                           Type = "IBSmultiAF")
  colnames(ibsMultiDD) <- c("Value", "Rel", "Type")

  #Sys.time() - start_time FOR A TIME CHECK

  print("IBSsingle DD")
  start_time <- Sys.time()
  ibsSingleDD <- data.frame(Value = ifelse(inbreeding,
                                           data.frame(diag(ibsSingleDF[idDF$drones, idDF$drones])),
                                           ibsSingleDF[idDF$drones, idDF$drones] %>% .[lower.tri(., diag = FALSE)] %>% data.frame(.)),
                            Rel = "DD",
                            Type = "IBSsingleAF")
  colnames(ibsSingleDD) <- c("Value", "Rel", "Type")
  if(!is.null(ibsColonyDF)){ #if dataframe is CSD info then ibdCdf is NULL
    print("IBScolony DD")
    ibsColonyDD <- data.frame(Value = ifelse(inbreeding,
                                             data.frame(diag(ibsColonyDF[idDF$drones, idDF$drones])),
                                             ibsColonyDF[idDF$drones, idDF$drones] %>% .[lower.tri(., diag = FALSE)] %>% data.frame(.)),
                           Rel = "DD",
                           Type = "IBScolonyAF")
    colnames(ibsColonyDD) <- c("Value", "Rel", "Type")
  }

  print("IBDr DD")
  ibdrDD <- data.frame(Value = ifelse(inbreeding,
                                      data.frame(diag(ibdDF[idDF$drones, idDF$drones])),
                                      ibdDF[idDF$drones, idDF$drones] %>% .[lower.tri(., diag = FALSE)] %>% data.frame(.)),
                       Rel = "DD",
                       Type = "IBDr")
  colnames(ibdrDD) <- c("Value", "Rel", "Type")

  if (!is.null(Sinv)) {
    print("IBDe DD")
    if (inbreeding){
      ibdeDD <- data.frame(Value = getS(Sinv, ids = idDF$drones, vector = TRUE, diagOnly = TRUE),
                           Rel = "DD",
                           Type = "IBDe")
    } else {
      x <- getS(Sinv, ids = idDF$drones, vector = FALSE)
      ind <- which(lower.tri(x, diag = FALSE), arr.ind = TRUE)
      ibdeDD <- data.frame(id1 = dimnames(x)[[2]][ind[,2]] ,
                           id2 = dimnames(x)[[1]][ind[,1]] ,
                           Value = x[ind],
                           Rel = "DD",
                           Type = "IBDe")
  }}

  # Queens vs workers
  if (isFALSE(inbreeding)){

      print("IBSAF0.5 QW")
      ibsAF0.5QW <- data.frame(Value = c(ibsAF0.5DF[idDF$workers, idDF$queen]),
                               Rel = "QW",
                               Type = "IBSAF0.5")
    print("IBSmulti QW")
    ibsMultiQW <- data.frame(Value = c(ibsMultiDF[idDF$workers, idDF$queen]),
                             Rel = "QW",
                             Type = "IBSmultiAF")

    print("IBSsingle QW")
    ibsSingleQW <- data.frame(Value = c(ibsSingleDF[idDF$workers, idDF$queen]),
                              Rel = "QW",
                              Type = "IBSsingleAF")

    if (!is.null(ibsColonyDF)){
      print("IBScolony QW")
      ibsColonyQW <- data.frame(Value = c(ibsColonyDF[idDF$workers, idDF$queen]),
                             Rel = "QW",
                             Type = "IBScolonyAF")
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
      ibsAF0.5QW <- NULL
      ibsSingleQW <- NULL
      ibsColonyQW <- NULL
      ibdrQW <- NULL
      ibdeQW <- NULL
    }
  # Queens vs drones
  if (isFALSE(inbreeding)){
      print("ibsAF0.5 QD")
      ibsAF0.5QD <- data.frame(Value = c(ibsAF0.5DF[idDF$drones, idDF$queen]),
                               Rel = "QD",
                               Type = "IBSAF0.5")

    print("IBSmulti QD")
    ibsMultiQD <- data.frame(Value = c(ibsMultiDF[idDF$drones, idDF$queen]),
                             Rel = "QD",
                             Type = "IBSmultiAF")

    print("IBSsingle QD")
    ibsSingleQD <- data.frame(Value = c(ibsSingleDF[idDF$drones, idDF$queen]),
                              Rel = "QD",
                              Type = "IBSsingleAF")

    if (!is.null(ibsColonyDF)){
      print("IBScolony QD")
      ibsColonyQD <- data.frame(Value = c(ibsColonyDF[idDF$drones, idDF$queen]),
                             Rel = "QD",
                             Type = "IBScolonyAF")
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
      ibsAF0.5QD <- NULL
      ibsSingleQD <- NULL
      ibsColonyQD <- NULL
      ibdrQD <- NULL
      ibdeQD <- NULL
    }
  return(if (is.null(ibsColonyDF)){
    bind_rows(
      ibsMultiWW, ibsMultiWD, ibsMultiDD, ibsMultiQD, ibsMultiQW,
      ibsAF0.5WW, ibsAF0.5WD, ibsAF0.5DD, ibsAF0.5QD, ibsAF0.5QW,
      ibsSingleWW, ibsSingleWD, ibsSingleDD, ibsSingleQD, ibsSingleQW,
      ibsColonyWW = NULL, ibsColonyWD = NULL, ibsColonyDD =  NULL, ibsColonyQD = NULL, ibsColonyQW = NULL,
      ibdrWW, ibdrWD, ibdrDD, ibdrQD, ibdrQW,
      ibdeWW, ibdeWD, ibdeDD, ibdeQD, ibdeQW)
  } else {
    bind_rows(
      ibsMultiWW, ibsMultiWD, ibsMultiDD, ibsMultiQD, ibsMultiQW,
      ibsAF0.5WW, ibsAF0.5WD, ibsAF0.5DD, ibsAF0.5QD, ibsAF0.5QW,
      ibsSingleWW, ibsSingleWD, ibsSingleDD, ibsSingleQD, ibsSingleQW,
      ibsColonyWW, ibsColonyWD, ibsColonyDD, ibsColonyQD, ibsColonyQW,
      ibdrWW, ibdrWD, ibdrDD, ibdrQD, ibdrQW,
      ibdeWW, ibdeWD, ibdeDD, ibdeQD, ibdeQW)
  })
}

prepareDataForPlotting_Queens <- function(ibsMultiDF = NULL, ibsCarDF = NULL, ibsMelDF = NULL, ibsColonyDF = NULL, ibsAF0.5DF = NULL, ibdDF = NULL, Sinv = NULL, idPopDF) {
  #Get populations IDs
  melID <- idPopDF$ID[idPopDF$Pop == "Mel"]
  melCrossID <- idPopDF$ID[idPopDF$Pop == "MelCross"]
  carID <- idPopDF$ID[idPopDF$Pop == "Car"]

  #Compute relationships between queens of different populations (QQ)
  #IBSmulti

    print("IBSAF0.5QQ")
    IBSAF0.5QQ <- rbind(data.frame(Value = as.vector(list(ibsAF0.5DF[melID, melCrossID])[[1]]),
                                   Pops = "Mel_MelCross"),
                        data.frame(Value = as.vector(list(ibsAF0.5DF[melID, carID])[[1]]),
                                   Pops = "Mel_Car"),
                        data.frame(Value = as.vector(list(ibsAF0.5DF[melCrossID, carID])[[1]]),
                                   Pops = "MelCross_Car"))
    IBSAF0.5QQ$Rel = "QQ"
    IBSAF0.5QQ$Type = "IBSAF0.5"

  print("IBSmultiQQ")
  IBSmultiQQ <- rbind(data.frame(Value = as.vector(list(ibsMultiDF[melID, melCrossID])[[1]]),
                                 Pops = "Mel_MelCross"),
                      data.frame(Value = as.vector(list(ibsMultiDF[melID, carID])[[1]]),
                                 Pops = "Mel_Car"),
                      data.frame(Value = as.vector(list(ibsMultiDF[melCrossID, carID])[[1]]),
                                 Pops = "MelCross_Car"))
  IBSmultiQQ$Rel = "QQ"
  IBSmultiQQ$Type = "IBSmultiAF"


  #IBScar - Using the Carnica AF is the single base AF
  print("IBScarQQ")
  IBScarQQ <- rbind(data.frame(Value = as.vector(list(ibsCarDF[melID, melCrossID])[[1]]),
                                  Pops = "Mel_MelCross"),
                       data.frame(Value = as.vector(list(ibsCarDF[melID, carID])[[1]]),
                                  Pops = "Mel_Car"),
                       data.frame(Value = as.vector(list(ibsCarDF[melCrossID, carID])[[1]]),
                                  Pops = "MelCross_Car"))
  IBScarQQ$Rel = "QQ"
  IBScarQQ$Type = "IBScarAF"

  #IBSmel - Using the Mellifera AF as the single base AF
  print("IBSmelQQ")
  IBSmelQQ <- rbind(data.frame(Value = as.vector(list(ibsMelDF[melID, melCrossID])[[1]]),
                               Pops = "Mel_MelCross"),
                    data.frame(Value = as.vector(list(ibsMelDF[melID, carID])[[1]]),
                               Pops = "Mel_Car"),
                    data.frame(Value = as.vector(list(ibsMelDF[melCrossID, carID])[[1]]),
                               Pops = "MelCross_Car"))
  IBSmelQQ$Rel = "QQ"
  IBSmelQQ$Type = "IBSmelAF"

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
    IBScolonyQQ$Type = "IBScolonyAF"
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

    print("IBSAF0.5Q")
    IBSAF0.5Q <- rbind(data.frame(Value = as.vector(list(c(ibsAF0.5DF[melID, melID][lower.tri(ibsAF0.5DF[melID, melID],
                                                                                              diag = FALSE)]))[[1]]),
                                  Pops = "Mel"),
                       data.frame(Value = as.vector(list(c(ibsAF0.5DF[carID, carID][lower.tri(ibsAF0.5DF[carID, carID],
                                                                                              diag = FALSE)]))[[1]]),
                                  Pops = "Car"),
                       data.frame(Value = as.vector(list(c(ibsAF0.5DF[melCrossID, melCrossID][lower.tri(ibsAF0.5DF[melCrossID, melCrossID],
                                                                                                        diag = FALSE)]))[[1]]),
                                  Pops = "MelCross"))
    IBSAF0.5Q$Rel <- "Q"
    IBSAF0.5Q$Type <- "IBSAF0.5"

  #IBSmultiAF
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
  IBSmultiQ$Type <- "IBSmultiAF"


  #IBSCarAF
  print("IBScarQ")
  IBScarQ <- rbind(data.frame(Value = as.vector(list(c(ibsCarDF[melID, melID][lower.tri(ibsCarDF[melID, melID],
                                                                                              diag = FALSE)]))[[1]]),
                                 Pops = "Mel"),
                      data.frame(Value = as.vector(list(c(ibsCarDF[carID, carID][lower.tri(ibsCarDF[carID, carID],
                                                                                              diag = FALSE)]))[[1]]),
                                 Pops = "Car"),
                      data.frame(Value = as.vector(list(c(ibsCarDF[melCrossID, melCrossID][lower.tri(ibsCarDF[melCrossID, melCrossID],
                                                                                                        diag = FALSE)]))[[1]]),
                                 Pops = "MelCross"))
  IBScarQ$Rel <- "Q"
  IBScarQ$Type <- "IBScarAF"

  #IBSMelAF
  print("IBSMelQ")
  IBSmelQ <- rbind(data.frame(Value = as.vector(list(c(ibsMelDF[melID, melID][lower.tri(ibsMelDF[melID, melID],
                                                                                        diag = FALSE)]))[[1]]),
                              Pops = "Mel"),
                   data.frame(Value = as.vector(list(c(ibsMelDF[carID, carID][lower.tri(ibsMelDF[carID, carID],
                                                                                        diag = FALSE)]))[[1]]),
                              Pops = "Car"),
                   data.frame(Value = as.vector(list(c(ibsMelDF[melCrossID, melCrossID][lower.tri(ibsMelDF[melCrossID, melCrossID],
                                                                                                  diag = FALSE)]))[[1]]),
                              Pops = "MelCross"))
  IBSmelQ$Rel <- "Q"
  IBSmelQ$Type <- "IBSmelAF"


  #IBSsingleAF
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
    IBScolonyQ$Type <- "IBScolonyAF"
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

    print("IBSAF0.5F")
    IBSfAF0.5 <- rbind(data.frame(Value = diag(ibsAF0.5DF[melID, melID]),
                                  Pops = "Mel"),
                       data.frame(Value = diag(ibsAF0.5DF[melCrossID, melCrossID]),
                                  Pops = "MelCross"),
                       data.frame(Value = diag(ibsAF0.5DF[carID, carID]),
                                  Pops = "Car"))
    IBSfAF0.5$Type = "IBSAF0.5"
    IBSfAF0.5$Rel = "F"

  #IBSfMulti
  print("IBSfMulti")
  IBSfMulti <- rbind(data.frame(Value = diag(ibsMultiDF[melID, melID]),
                                Pops = "Mel"),
                     data.frame(Value = diag(ibsMultiDF[melCrossID, melCrossID]),
                                Pops = "MelCross"),
                     data.frame(Value = diag(ibsMultiDF[carID, carID]),
                                Pops = "Car"))
  IBSfMulti$Type = "IBSmultiAF"
  IBSfMulti$Rel = "F"


  #IBSfCar
  print("IBSfCar")
  IBSfCar<- rbind(data.frame(Value = diag(ibsCarDF[melID, melID]),
                                 Pops = "Mel"),
                      data.frame(Value = diag(ibsCarDF[melCrossID, melCrossID]),
                                 Pops = "MelCross"),
                      data.frame(Value = diag(ibsCarDF[carID, carID]),
                                 Pops = "Car"))
  IBSfCar$Type = "IBScarAF"
  IBSfCar$Rel = "F"


  #IBSfMel
  print("IBSfMel")
  IBSfMel<- rbind(data.frame(Value = diag(ibsMelDF[melID, melID]),
                             Pops = "Mel"),
                  data.frame(Value = diag(ibsMelDF[melCrossID, melCrossID]),
                             Pops = "MelCross"),
                  data.frame(Value = diag(ibsMelDF[carID, carID]),
                             Pops = "Car"))
  IBSfMel$Type = "IBSmelAF"
  IBSfMel$Rel = "F"

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
    IBSfColony$Type <- "IBScolonyAF"
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
      IBSmultiQQ, IBScarQQ, IBSmelQQ, IBScolonyQQ, IBSAF0.5QQ, IBDrQQ, IBDeQQ,
      IBSmultiQ, IBScarQ, IBSmelQ, IBScolonyQ, IBSAF0.5Q, IBDrQ, IBDeQ,
      IBSfMulti, IBSfCar, IBSfMel, IBSfColony, IBSfAF0.5, IBDrF, IBDeF)
  )
}

plotColony <- function(df, rel = NULL, type = NULL, years = NULL) {
  df$Type <- factor(df$Type, levels = c("IBDe", "IBDr", "IBSsingleAF", "IBScolonyAF", "IBSmultiAF", "IBSAF0.5"))
  df$Rel <- factor(df$Rel, levels= c("WW", "WD", "DD", "QW", "QD"))
  df$SisterType <- factor(df$SisterType, level = c("HS", "FS", "SS"))

  type_labels <- c("IBDe", "IBDr", "IBSsingleAF", "IBScolonyAF", "IBSmultiAF", "IBSAF0.5")
  names(type_labels) <- c("IBDe", "IBDr", "IBSsingleAF", "IBScolonyAF", "IBSmultiAF", "IBSAF0.5")

  df$Year <- factor(df$Year, levels = years)
  year_labels <- unlist(paste("Year", years, collapse = "_") %>% stringr::str_split("_"))
  names(year_labels) <- years

  plot <- ggplot(df[df$Rel %in% rel & df$Type %in% type, ],
              aes(x = Value, fill = Rel)) +
          geom_vline(xintercept = c(0, 0.25, 0.5, 0.75), linewidth = 0.25, colour = "grey") +
          geom_histogram(bins = 200, position = "identity") +
          facet_grid2(rows = vars(Type), cols = vars(Year), scales = "free_y",independent = "y", labeller = labeller(Type = type_labels, Year = year_labels)) +
          scale_fill_manual("", values = paletteViridis, aesthetics = c("colour","fill")) +
          theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "top" , axis.text.y=element_blank(),
                axis.ticks.y=element_blank(), axis.title.y = element_blank()) +
          scale_x_continuous(breaks=seq(-2, 10, 0.25))
  return(plot)
}

plotQueens <- function(df, rel = NULL, type = NULL, pops = NULL, years = NULL, palette = NULL) { #Use palette = cbPaletteQ with same pop and palette = cbPaletteQQ with different pops
  df$Type <- factor(df$Type, levels = c("IBDe", "IBDr", "IBSsingleAF", "IBScolonyAF", "IBSmultiAF", "IBSAF0.5"))
  df$Pops <- factor(df$Pops, levels = c("Mel_MelCross", "Mel_Car", "MelCross_Car", "Mel", "Car", "MelCross"))
  df$Rel  <- factor(df$Rel, levels = c("QQ", "Q", "F"))

  #TODO new type labels required?
  type_labels <- c("IBDe", "IBDr", "IBSsingleAF", "IBScolonyAF", "IBSmultiAF", "IBSAF0.5")
  names(type_labels) <- c("IBDe", "IBDr", "IBSsingleAF", "IBScolonyAF", "IBSmultiAF", "IBSAF0.5")

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

#to change text size axis.text=element_text(size=12), axis.title = element_text(size= 20), strip.text.x = element_text(size = 20))

########################################################
### --- FIGURE 1&2: Pure subspecies (carnica) in years 1 and 10 ---###
########################################################
print("Plot carnica year 1")

#Coded out data prep
{
relCar1 <- prepareDataForPlotting_Colony(ibsMultiDF = colonyCar1$IBSmultiAF,
                                         ibsSingleDF = colonyCar1$IBSsingleAF,
                                         ibsColonyDF = colonyCar1$IBScolonyAF,
                                         ibsAF0.5DF = colonyCar1$IBSAF0.5,
                                         ibdDF = colonyCar1$IBD,
                                         Sinv = Sinv,
                                         idDF = colonyCar1$ID,
                                         inbreeding = FALSE,
                                         WorkersFatherTable = colonyCar1$WorkersFatherTable)



print("Plot carnica year 10")
relCar10 <- prepareDataForPlotting_Colony(ibsMultiDF = colonyCar10$IBSmultiAF,
                                          ibsSingleDF = colonyCar10$IBSsingleAF,
                                          ibsColonyDF = colonyCar10$IBScolonyAF,
                                          ibsAF0.5DF = colonyCar10$IBSAF0.5,
                                          ibdDF = colonyCar10$IBD,
                                          Sinv = Sinv,
                                          idDF = colonyCar10$ID,
                                          inbreeding = FALSE,
                                          WorkersFatherTable = colonyCar10$WorkersFatherTable)

relCar1$Year <- 1
relCar10$Year <- 10
dataCar <- rbind(relCar1, relCar10)
save(... = dataCar, file = "dataCar.RData")
}

CarWWplot<- plotColony(relCar1, type = c("IBDe", "IBDr", "IBSsingleAF"), rel = c("WW", "WD", "DD"), years = c(1))

CarQWplot <- plotColony(dataCar, type = c("IBDe", "IBDr", "IBSsingleAF", "IBScolonyAF", "IBSmultiAF", "IBSAF0.5"),  rel = c("QW", "QD"), years = c(1,10))

dataCar1Sister <- SisterTypeDF(dataCar[dataCar$Year == 1, ])
dataCar10Sister<- SisterTypeDF(dataCar[dataCar$Year == 10, ])

dataCarSisters <- rbind(dataCar1Sister, dataCar10Sister)

Fig1Table <- dataCarSisters[dataCarSisters$Type %in% c("IBDe", "IBDr", "IBSsingleAF"),]
Year1 <- Fig1Table %>% filter(Year == 1 ) %>% pivot_wider(id_cols = c(SisterType), names_from = Type, values_from = Mean) %>% as.data.frame(.) %>% .[, c(1,4,2,3)]
Year10 <- Fig1Table %>% filter(Year == 10 ) %>% pivot_wider(id_cols = c(SisterType), names_from = Type, values_from = Mean) %>% as.data.frame(.)  %>% .[, c(1,4,2,3)]


#Get average mean of SisterType information: e.g using ibdeWW
# tapply(ibdeWW$Value, INDEX= ibdeWW$SisterType, FUN = mean)
#Get mean or sd 
#  tmp<- (relCar1[relCar1$Type == "IBDe",])
# tapply(tmp$Value, INDEX = tmp$Rel, FUN= sd)

#Plot Car inbreeding
# Coded out prep
{
relCar1F <- prepareDataForPlotting_Colony(ibsMultiDF = colonyCar1$IBSmultiAF,
                                          ibsSingleDF = colonyCar1$IBSsingleAF,
                                          ibsColonyDF = colonyCar1$IBScolonyAF,
                                          ibsAF0.5DF = colonyCar1$IBSAF0.5,
                                          ibdDF = colonyCar1$IBD,
                                          Sinv = Sinv,
                                          idDF = colonyCar1$ID,
                                          inbreeding = TRUE,
                                          WorkersFatherTable = colonyCar1$WorkersFatherTable)

relCar10F <- prepareDataForPlotting_Colony(ibsMultiDF = colonyCar10$IBSmultiAF,
                                           ibsSingleDF = colonyCar10$IBSsingleAF,
                                           ibsColonyDF = colonyCar10$IBScolonyAF,
                                           ibsAF0.5DF = colonyCar10$IBSAF0.5,
                                           ibdDF = colonyCar10$IBD,
                                           Sinv = Sinv,
                                           idDF = colonyCar10$ID,
                                           inbreeding = TRUE,
                                           WorkersFatherTable = colonyCar1$WorkersFatherTable)
relCar1F$Year <- 1
relCar10F$Year <- 10

dataCarF <- rbind(relCar1F, relCar10F)
save(dataCarF, file = "dataCarF.Rdata")
}

CarFplot <- plotColony(dataCarF, type =  c("IBDe", "IBDr", "IBSsingleAF", "IBScolonyAF", "IBSmultiAF", "IBSAF0.5"), rel = c("WW", "DD"), years = c(1,10))


#Plot CAR csd Locus
#Coded out prep
{
relCar1_csd <- prepareDataForPlotting_Colony(ibsMultiDF = colonyCar1$IBSmultiAFCsdLocus,
                                             ibsSingleDF = colonyCar1$IBSsingleAFcsdLocus,
                                             ibsColonyDF = colonyCar1$colonyAF_CSDLocus,
                                             ibsAF0.5DF = colonyCar1$IBSAF0.5CSDLocus,
                                             ibdDF = colonyCar1$IBDCsdLocus,
                                             Sinv = Sinv,
                                             idDF = colonyCar1$ID,
                                             WorkersFatherTable = colonyCar1$WorkersFatherTable)


relCar10_csd <- prepareDataForPlotting_Colony(ibsMultiDF = colonyCar10$IBSmultiAFCsdLocus,
                                               ibsSingleDF = colonyCar10$IBSsingleAFcsdLocus,
                                              ibsColonyDF = colonyCar10$colonyAF_CSDLocus,
                                               ibsAF0.5DF = colonyCar10$IBSAF0.5CSDLocus,
                                               ibdDF = colonyCar10$IBDCsdLocus,
                                               Sinv = Sinv,
                                               idDF = colonyCar10$ID,
                                              WorkersFatherTable = colonyCar10$WorkersFatherTable)
relCar1_csd$Year <- 1
relCar10_csd$Year <- 10
dataCarCSDloc <- rbind(relCar1_csd, relCar10_csd)
save(dataCarCSDloc, file = "dataCarCSDloc.Rdata")
}

CarCSDlocplot <- plotColony(dataCarCSDloc, type = c("IBSsingleAF", "IBScolonyAF", "IBSmultiAF", "IBSAF0.5"), rel = c("WW", "WD", "DD"), years = c(1,10))

#Plot CAR csd Chromosome
{
relCar1_csdChr <- prepareDataForPlotting_Colony(ibsMultiDF = colonyCar1$IBSmultiCSDChr,
                                                ibsSingleDF = colonyCar1$IBSsingleCSDChr,
                                                ibsColonyDF = colonyCar1$IBScolonyAF_CSDChr,
                                                ibsAF0.5DF = colonyCar1$IBSAF0.5CSDChr,
                                                ibdDF = colonyCar1$IBDcsdChr,
                                                Sinv = Sinv,
                                                idDF = colonyCar1$ID,
                                                WorkersFatherTable = colonyCar1$WorkersFatherTable)
relCar10_csdChr <- prepareDataForPlotting_Colony(ibsMultiDF = colonyCar10$IBSmultiCSDChr,
                                                  ibsSingleDF = colonyCar10$IBSsingleCSDChr,
                                                 ibsColonyDF = colonyCar10$IBScolonyAF_CSDChr,
                                                  ibsAF0.5DF = colonyCar10$IBSAF0.5CSDChr,
                                                  ibdDF = colonyCar10$IBDcsdChr,
                                                  Sinv = Sinv,
                                                  idDF = colonyCar10$ID,
                                                 WorkersFatherTable = colonyCar10$WorkersFatherTable)
relCar1_csdChr$Year <- 1
relCar10_csdChr$Year <- 10

dataCarCSDchr <- rbind(relCar1_csdChr, relCar10_csdChr)
save(... = dataCarCSDchr, file = "dataCarCsdChr.RData")

}

tmp<- (relCar1_csdChr[relCar1_csdChr$Rel == "WW",])
tapply(tmp$Value, INDEX = tmp$Type, FUN= sd)

CarCSDchr <- plotColony(dataCarCSDchr, rel = c("WW", "WD", "DD"), type = c("IBSsingleAF", "IBScolonyAF", "IBSmultiAF", "IBSAF0.5"), years = c(1, 10))
########################################################
### --- FIGURE 4: Between queens of different subspecies (carnica vs. mellifera) ---###
########################################################
#Plot queens Year 1
print("Plot queens")
#Coded out prep
{
relQueens1_MelAF <- prepareDataForPlotting_Queens(ibsMultiDF =  queens1_MelAF$IBSmultiAF,
                                            ibsSingleDF = queens1_MelAF$IBSsingleAF,
                                            ibsColonyDF = queens1_MelAF$IBScolonyAF,
                                            ibsAF0.5DF = queens1_MelAF$IBSAF0.5,
                                            ibdDF = queens1_MelAF$IBD,
                                            Sinv = Sinv,
                                            idPopDF = springerQueensPop1)

relQueens10_MelAF <-  prepareDataForPlotting_Queens(ibsMultiDF =  queens10_MelAF$IBSmultiAF,
                                              ibsSingleDF = queens10_MelAF$IBSsingleAF,
                                              ibsColonyDF = queens10_MelAF$IBScolonyAF,
                                              ibsAF0.5DF = queens10_MelAF$IBSAF0.5,
                                              ibdDF = queens10_MelAF$IBD,
                                              Sinv = Sinv,
                                              idPopDF = springerQueensPop10)
relQueens1_MelAF$Year <- 1
relQueens10_MelAF$Year <- 10

dataQueens_MelAF <- rbind(relQueens1_MelAF, relQueens10_MelAF)
save(relQueens10_MelAF, file = "relQueens10_MelAF.Rdata")
save(dataQueens_MelAF, file = "dataQueens_MelAF.Rdata")
}
{relQueens1_CarAF <- prepareDataForPlotting_Queens(ibsMultiDF =  queens1_CarAF$IBSmultiAF,
                                                    ibsSingleDF = queens1_CarAF$IBSsingleAF,
                                                    ibsColonyDF = queens1_CarAF$IBScolonyAF,
                                                    ibsAF0.5DF = queens1_CarAF$IBSAF0.5,
                                                    ibdDF = queens1_CarAF$IBD,
                                                    Sinv = Sinv,
                                                    idPopDF = springerQueensPop1)

  relQueens10_CarAF <-  prepareDataForPlotting_Queens(ibsMultiDF =  queens10_CarAF$IBSmultiAF,
                                                      ibsSingleDF = queens10_CarAF$IBSsingleAF,
                                                      ibsColonyDF = queens10_CarAF$IBScolonyAF,
                                                      ibsAF0.5DF = queens10_CarAF$IBSAF0.5,
                                                      ibdDF = queens10_CarAF$IBD,
                                                      Sinv = Sinv,
                                                      idPopDF = springerQueensPop10)
  relQueens1_CarAF$Year <- 1
  relQueens10_CarAF$Year <- 10

  dataQueens_CarAF <- rbind(relQueens1_CarAF, relQueens10_CarAF)
  save(relQueens10_CarAF, file = "relQueens10_CarAF.Rdata")
  save(dataQueens_CarAF, file = "dataQueens_CarAF.Rdata")
}

#Queens with Mel and Car allele frequencies for separate "SingleAF" plots
{
  relQueens1 <- prepareDataForPlotting_Queens(ibsMultiDF = queens1_MelAF$IBSmultiAF, #doesn't matter which since both use same multibase pop
                                            ibsCarDF = queens1_CarAF$IBSsingleAF,
                                            ibsMelDF = queens1_MelAF$IBSsingleAF,
                                            ibsAF0.5DF = queens1_MelAF$IBSAF0.5,
                                            ibdDF = queens1_MelAF$IBD,
                                            Sinv = Sinv,
                                            idPopDF = springerQueensPop1)

relQueens10 <- prepareDataForPlotting_Queens(ibsMultiDF = queens10_MelAF$IBSmultiAF, #doesn't matter which since both use same multibase pop
                                            ibsCarDF = queens10_CarAF$IBSsingleAF,
                                            ibsMelDF = queens10_MelAF$IBSsingleAF,
                                            ibsAF0.5DF = queens10_MelAF$IBSAF0.5,
                                            ibdDF = queens10_MelAF$IBD,
                                            Sinv = Sinv,
                                            idPopDF = springerQueensPop10)

relQueens1$Year <- 1
relQueens10$Year <- 10

dataQueens <- rbind(relQueens1, relQueens10)
save(dataQueens, file = "dataQueens.Rdata")
}
print("Between populations ")
#between populations
QueensQQ <- plotQueens(dataQueens_MelAF, rel = "QQ", type = c("IBDe", "IBDr", "IBSsingleAF", "IBSmultiAF", "IBSAF0.5"), pops = c("Mel_MelCross", "Mel_Car", "MelCross_Car", "IBSAF0.5"), years = c(1,10), palette = cbPaletteQQ)
QueensQQ

print("Within populations")
#within population - non-diagonal
QueensQ <- plotQueens(relQueens10_MelAF, rel = "Q", type = c("IBDe", "IBDr", "IBSmultiAF") , pops = c("Mel", "Car", "MelCross", "IBSAF0.5"), years = c(10), palette = cbPaletteQ)
QueensQ
#within population - diagonal
QueensF <- plotQueens(dataQueens_MelAF, rel = "F", type = c("IBDe", "IBDr", "IBSmultiAF", "IBSsingleAF", "IBSAF0.5"), pops = c("Mel", "Car", "MelCross", "IBSAF0.5"), years = c(1,10), palette = cbPaletteQ)
QueensF



#CSD STUFF
#Plot csd Loc as a histogram (QQ)
{
relQueens1_csdLoc_hist <- prepareDataForPlotting_Queens(ibsMultiDF = springerQueens1_MelAF$IBSmultiAFCsdLocus,
                                                        ibsSingleDF = springerQueens1_MelAF$IBSsingleAFcsdLocus,
                                                        ibsColonyDF = springerQueens1_MelAF$colonyAF_CSDLocus,
                                                        ibsAF0.5DF = springerQueens1_MelAF$IBSAF0.5CSDLocus,
                                                        ibdDF = springerQueens1_MelAF$IBDCsdLocus,
                                                        Sinv = Sinv,
                                                        idPopDF = springerQueensPop1)

relQueens10_csdLoc_hist <- prepareDataForPlotting_Queens(ibsMultiDF = springerQueens10_MelAF$IBSmultiAFCsdLocus,
                                                            ibsSingleDF = springerQueens10_MelAF$IBSsingleAFcsdLocus,
                                                            ibsColonyDF = springerQueens10_MelAF$colonyAF_CSDLocus,
                                                         ibsAF0.5DF = springerQueens10_MelAF$IBSAF0.5CSDLocus,
                                                            ibdDF = springerQueens10_MelAF$IBDCsdLocus,
                                                            Sinv = Sinv,
                                                            idPopDF = springerQueensPop10)
relQueens1_csdLoc_hist$Year <- 1
relQueens10_csdLoc_hist$Year <- 10
dataQueens_csdLoc <- rbind(relQueens1_csdLoc_hist, relQueens10_csdLoc_hist)
save(dataQueens_csdLoc, file = "dataQueensCSDloc.RData")
}

Queens_csdLoc_hist <- plotQueens(dataQueens_csdLoc, rel = "QQ",  type = c("IBDe", "IBDr","IBSsingleAF", "IBSmultiAF", "IBSAF0.5"), pops = c("Mel_MelCross", "Mel_Car", "MelCross_Car"), years = c(1,10), palette = cbPaletteQQ)
Queens_csdLoc_hist


#Plot csd Chr as Histogram
Queens_csdChr_hist <- plotQueens(dataQueensCSDchr, rel = "QQ",  type = c("IBDe", "IBDr", "IBSmultiAF", "IBSsingleAF", "IBSAF0.5"), pops = c("Mel_MelCross", "Mel_Car", "MelCross_Car"), years = c(1,10), palette = cbPaletteQQ)


