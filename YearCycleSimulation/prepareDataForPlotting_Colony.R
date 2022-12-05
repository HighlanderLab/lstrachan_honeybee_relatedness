prepareDataForPlotting_Colony <- function(ibsMultiDF = NULL,
                                          ibsSingleDF = NULL,
                                          ibsOwnDF = NULL,
                                          ibdDF = NULL,
                                          Sinv = NULL,
                                          idDF,
                                          inbreeding = FALSE) {

  #Workers on workers
  print("IBSmulti WW")
  start_time <- Sys.time()
  ibsMultiWW <- data.frame(Value = ifelse(isTRUE(inbreeding),       #NEED TO ADD A ISNULL?
                                     diag(ibsMultiDF[idDF$workers, idDF$workers]),
                                     as.data.frame(ibsMultiDF[idDF$workers, idDF$workers] %>%
                                       c(tmp[lower.tri(tmp, diag = FALSE)]))),
                      Rel = "WW",
                      Type = "IBSmultiBF")
  Sys.time() - start_time #FOR A TIME CHECK

  print("IBSsingle WW")
  start_time <- Sys.time()
  ibsSingleWW <- data.frame(Value = ifelse(isTRUE(inbreeding),
                                     diag(ibsSingleDF[idDF$workers, idDF$workers]),
                                     as.data.frame(ibsSingleDF[idDF$workers, idDF$workers] %>%
                                                     c(tmp[lower.tri(tmp, diag = FALSE)]))),
                      Rel = "WW",
                      Type = "IBSsingleBF")

  if(!is.null(ibsOwnDF)){ #if dataframe is CSD info then ibdOwnDF is NULL
    print("IBSown WW")
    ibsOwnWW <- data.frame(Value = ifelse(isTRUE(inbreeding),
                                        diag(ibsOwnDF[idDF$workers, idDF$workers]),
                                        ibsOwnDF[idDF$workers, idDF$workers] %>%
                                          c(tmp[lower.tri(tmp, diag = FALSE)])),
                         Rel = "WW",
                         Type = "IBSOwnFreq")

  }

  print("IBDr WW")
  ibdrWW <- data.frame(Value = ifelse(isTRUE(inbreeding),
                                      diag(ibdDF[idDF$workers, idDF$workers]),
                                      ibdDF[idDF$workers, idDF$workers] %>%
                                        c(tmp[lower.tri(tmp, diag = FALSE)])),
                       Rel = "WW",
                       Type = "IBDr")

  if (!is.null(Sinv)) {
    print("IBDe WW")
    ibdeWW <- data.frame(Value = ifelse(isTRUE(inbreeding),
                                        getS(Sinv, ids = idDF$workers, vector = TRUE, diagOnly = TRUE),
                                        getS(Sinv, ids = idDF$workers, vector = TRUE) %>%
                                          c(tmp[lower.tri(tmp, diag = FALSE)])),
                         Rel = "WW",
                         Type = "IBDe")
  }
  # Workers vs drones
  print("IBSmulti WD")
  ibsMultiWD <- data.frame(Value = ifelse(isTRUE(inbreeding),
                                     NULL,
                                     c(ibsMultiDF[idDF$workers, idDF$drones])),
                      Rel = "WD",
                      Type = "IBSmultiBF")

  print("IBSsingle WD")
  ibsSingleWD <- data.frame(Value = ifelse(isTRUE(inbreeding),
                                     NULL,
                                     c(ibsSingleDF[idDF$workers, idDF$drones])),
                      Rel = "WD",
                      Type = "IBSsingleBF")

  if (!is.null(ibsOwnDF)){
    print("IBSown WD")
    ibsOwnWD <- data.frame(Value = ifelse(isTRUE(inbreeding),
                                        NULL,
                                        c(ibsOwnDF[idDF$workers, idDF$drones])),
                         Rel = "WD",
                         Type = "IBSOwnFreq")
  }

  print("IBDr WD")
  ibdrWD <- data.frame(Value = ifelse(isTRUE(inbreeding),
                                      NULL,
                                      c(ibdDF[idDF$workers, idDF$drones])),
                       Rel = "WD",
                       Type = "IBDr")

  if (!is.null(Sinv)) {
    print("IBDe WD")
    ibdeWD <- data.frame(Value = ifelse(isTRUE(inbreeding),
                                        NULL,
                                        getS(Sinv, ids = idDF$workers, with = idDF$drones, vector = TRUE)),
                         Rel = "WW",
                         Type = "IBDe")
  }
  #Drones vs drones
  print("IBSmulti DD")
  start_time <- Sys.time()
  ibsMultiDD <- data.frame(Value = ifelse(isTRUE(inbreeding),
                                     diag(ibsMultiDF[idDF$drones, idDF$drones]),
                                     ibsMultiDF[idDF$drones, idDF$drones] %>%
                                       c(tmp[lower.tri(tmp, diag = FALSE)])),
                      Rel = "DD",
                      Type = "IBSmultiBF")
  #Sys.time() - start_time FOR A TIME CHECK

  print("IBSsingle DD")
  start_time <- Sys.time()
  ibsSingleDD <- data.frame(Value = ifelse(isTRUE(inbreeding),
                                     diag(ibsSingleDF[idDF$drones, idDF$drones]),
                                     ibsSingleDF[idDF$drones, idDF$drones] %>%
                                       c(tmp[lower.tri(tmp, diag = FALSE)])),
                      Rel = "DD",
                      Type = "IBSsingleBF")

  if(!is.null(ibsOwnDF)){ #if dataframa is CSD info then ibdCdf is NULL
    print("IBSown DD")
    ibsOwnDD <- data.frame(Value = ifelse(isTRUE(inbreeding),
                                        diag(ibsOwnDF[idDF$drones, idDF$drones]),
                                        ibsOwnDF[idDF$drones, idDF$drones] %>%
                                          c(tmp[lower.tri(tmp, diag = FALSE)])),
                         Rel = "DD",
                         Type = "IBSOwnFreq")

  }

  print("IBDr DD")
  ibdrDD <- data.frame(Value = ifelse(isTRUE(inbreeding),
                                      diag(ibdDF[idDF$drones, idDF$drones]),
                                      ibdDF[idDF$drones, idDF$drones] %>%
                                        c(tmp[lower.tri(tmp, diag = FALSE)])),
                       Rel = "DD",
                       Type = "IBDr")

  if (!is.null(Sinv)) {
    print("IBDe DD")
    ibdeDD <- data.frame(Value = ifelse(isTRUE(inbreeding),
                                        getS(Sinv, ids = idDF$drones, vector = TRUE, diagOnly = TRUE),
                                        getS(Sinv, ids = idDF$drones, vector = TRUE) %>%
                                          c(tmp[lower.tri(tmp, diag = FALSE)])),
                         Rel = "DD",
                         Type = "IBDe")
  }
  # Queens vs workers
  print("IBSmulti QW")
  ibsMultiQW <- data.frame(Value = ifelse(isTRUE(inbreeding),
                                     NULL,
                                     c(ibsMultiDF[idDF$workers, idDF$queen])),
                      Rel = "QW",
                      Type = "IBSmultiBF")

  print("IBSsingle QW")
  ibsSingleQW <- data.frame(Value = ifelse(isTRUE(inbreeding),
                                     NULL,
                                     c(ibsSingleDF[idDF$workers, idDF$queen])),
                      Rel = "QW",
                      Type = "IBSsingleBF")

  if (!is.null(ibsOwnDF)){
    print("IBSown QW")
    ibsOwnQW <- data.frame(Value = ifelse(isTRUE(inbreeding),
                                        NULL,
                                        c(ibsOwnDF[idDF$workers, idDF$queen])),
                         Rel = "QW",
                         Type = "IBSOwnFreq")
  }

  print("IBDr QW")
  ibdrQW <- data.frame(Value = ifelse(isTRUE(inbreeding),
                                      NULL,
                                      c(ibdDF[idDF$workers, idDF$queen])),
                       Rel = "QW",
                       Type = "IBDr")

  if (!is.null(Sinv)) {
    print("IBDe QW")
    ibdeQW <- data.frame(Value = ifelse(isTRUE(inbreeding),
                                        NULL,
                                        getS(Sinv, ids = idDF$workers, with = idDF$queen, vector = TRUE)),
                         Rel = "QW",
                         Type = "IBDe")
  }
  # Queens vs drones
  print("IBSmulti QD")
  ibsMultiQD <- data.frame(Value = ifelse(isTRUE(inbreeding),
                                     NULL,
                                     c(ibsMultiDF[idDF$drones, idDF$queen])),
                      Rel = "QD",
                      Type = "IBS")

  print("IBSsingle QD")
  ibsSingleQD <- data.frame(Value = ifelse(isTRUE(inbreeding),
                                     NULL,
                                     c(ibsSingleDF[idDF$drones, idDF$queen])),
                      Rel = "QD",
                      Type = "IBS")

  if (!is.null(ibsOwnDF)){
    print("IBSown QD")
    ibsOwnQD <- data.frame(Value = ifelse(isTRUE(inbreeding),
                                        NULL,
                                        c(ibsOwnDF[idDF$drones, idDF$queen])),
                         Rel = "QD",
                         Type = "IBSOwnFreq")
  }

  print("IBDr QD")
  ibdrQD <- data.frame(Value = ifelse(isTRUE(inbreeding),
                                      NULL,
                                      c(ibdDF[idDF$drones, idDF$queen])),
                       Rel = "QD",
                       Type = "IBDr")

  if (!is.null(Sinv)) {
    print("IBDe QD")
    ibdeQD <- data.frame(Value = ifelse(isTRUE(inbreeding),
                                        NULL,
                                        getS(Sinv, ids = idDF$drones, with = idDF$queen, vector = TRUE)),
                         Rel = "QD",
                         Type = "IBDe")
  }
  return(if (is.null(ibsOwnDF)){
          rbind(
            ibsMultiWW, ibsMultiWD, ibsMultiDD, ibsMultiQD, ibsMultiQW,
            ibsSingleWW, ibsSingleWD, ibsSingleDD, ibsSingleQD, ibsSingleQW,
            ibsOwnWW = NULL, ibsOwnWD = NULL, ibsOwnDD =  NULL, ibsOwnQD = NULL, ibsOwnQW = NULL,
            ibdrWW, ibdrWD, ibdrDD, ibdrQD, ibdrQW,
            ibdeWW, ibdeWD, ibdeDD, ibdeQD, ibdeQW)
         } else {
           rbind(
             ibsMultiWW, ibsMultiWD, ibsMultiDD, ibsMultiQD, ibsMultiQW,
             ibsSingleWW, ibsSingleWD, ibsSingleDD, ibsSingleQD, ibsSingleQW,
             ibsOwnWW, ibsOwnWD, ibsOwnDD, ibsOwnQD, ibsOwnQW,
             ibdrWW, ibdrWD, ibdrDD, ibdrQD, ibdrQW,
             ibdeWW, ibdeWD, ibdeDD, ibdeQD, ibdeQW)
         })
}

