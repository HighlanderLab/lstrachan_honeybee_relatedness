library(dplyr)
library(tidyr)
library(ggplot2)
library(Matrix)
library(SIMplyBee)

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
cbPalette <- c( "#D55E00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7")
# To use for fills, add
scale_fill_manual(values=cbPalette)

# To use for line and point colors, add
scale_colour_manual(values=cbPalette)

print("Assigning objects")

ped <- pedigree
caste <- caste
# Colony in year 1
# Mellifera
colonyMel1 <- springerColony1_Mel
ibsMel1 <- colonyMel1$IBS
ibsMel1_csdChr <- colonyMel1$IBScsdChr
ibsMel1_csd <- colonyMel1$IBSCsd
ibdMel1 <- colonyMel1$IBD
ibdMel1_csdChr <- colonyMel1$IBDcsdChr
ibdMel1_csd <- colonyMel1$IBDCsd
idMel1 <- colonyMel1$ID
# Mellifera cross
colonyMelCross1 <- springerColony1_MelCross
ibsMelCross1 <- colonyMelCross1$IBS
ibsMelCross1_csdChr <- colonyMelCross1$IBScsdChr
ibsMelCross1_csd <- colonyMelCross1$IBSCsd
ibdMelCross1 <- colonyMelCross1$IBD
ibdMelCross1_csdChr <- colonyMelCross1$IBDcsdChr
ibdMelCross1_csd <- colonyMelCross1$IBDCsd
idMelCross1 <- colonyMelCross1$ID
# Carnica
colonyCar1 <- springerColony1_Car
ibsCar1 <- colonyCar1$IBS
ibsCar1_csdChr <- colonyCar1$IBScsdChr
ibsCar1_csd <- colonyCar1$IBSCsd
ibdCar1 <- colonyCar1$IBD
ibdCar1_csdChr <- colonyCar1$IBDcsdChr
ibdCar1_csd <- colonyCar1$IBDCsd
idCar1 <- colonyCar1$ID
# Queens in year 1
queens1 <- springerQueens1
ibsQueens1 <- queens1$IBS
ibsQueens1_csdChr <- queens1$IBScsdChr
ibsQueens1_csd <- queens1$IBSCsd
ibdQueens1 <- queens1$IBD
ibdQueens1_csdChr <- queens1$IBDcsdChr
ibdQueens1_csd <- queens1$IBDCsd
idQueens1 <- queens1$ID
idPopQueens1 <- springerQueensPop1


# Colony in year 10
# Mellifera
colonyMel10 <- springerColony10_Mel
ibsMel10 <- colonyMel10$IBS
ibsMel10_csdChr <- colonyMel10$IBScsdChr
ibsMel10_csd <- colonyMel10$IBSCsd
ibdMel10 <- colonyMel10$IBD
ibdMel10_csdChr <- colonyMel10$IBDcsdChr
ibdMel10_csd <- colonyMel10$IBDCsd
idMel10 <- colonyMel10$ID
# Mellifera cross
colonyMelCross10 <- springerColony10_MelCross
ibsMelCross10 <- colonyMelCross10$IBS
ibsMelCross10_csdChr <- colonyMelCross10$IBScsdChr
ibsMelCross10_csd <- colonyMelCross10$IBSCsd
ibdMelCross10 <- colonyMelCross10$IBD
ibdMelCross10_csdChr <- colonyMelCross10$IBDcsdChr
ibdMelCross10_csd <- colonyMelCross10$IBDCsd
idMelCross10 <- colonyMelCross10$ID
# Carnica
colonyCar10 <- springerColony10_Car
ibsCar10 <- colonyCar10$IBS
ibsCar10_csdChr <- colonyCar10$IBScsdChr
ibsCar10_csd <- colonyCar10$IBSCsd
ibdCar10 <- colonyCar10$IBD
ibdCar10_csdChr <- colonyCar10$IBDcsdChr
ibdCar10_csd <- colonyCar10$IBDCsd
idCar10 <- colonyCar10$ID
# Queens in year 10
queens10 <- springerQueens10
ibsQueens10 <- queens10$IBS
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
prepareDataForPlotting_Colony <- function(ibsDF = NULL, ibdDF = NULL, Sinv = NULL,  idDF) {
  # workers to workers
  print("IBS WW")
  tmp <- ibsDF[idDF$workers, idDF$workers]
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
  IBS_WD1 <- c(ibsDF[idDF$workers, idDF$drones])
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
  tmp <- ibsDF[idDF$drones, idDF$drones]
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
  print("IBS QD")
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

prepareDataForPlotting_ColonyDiag <- function(ibsDF = NULL, ibdDF = NULL, Sinv = NULL,  idDF) {
  # workers
  print("IBS WW")
  tmp <- diag(ibsDF[idDF$workers, idDF$workers])
  ret <- data.frame(Value = tmp, Rel = "WW", Type = "IBS")

  print("IBDr WW")
  tmp <- diag(ibdDF[idDF$workers, idDF$workers])
  ret <- rbind(ret, data.frame(Value = tmp, Rel = "WW", Type = "IBDr"))

  if (!is.null(Sinv)) {
    print("IBDe WW")
    tmp <- getS(Sinv, ids = idDF$workers, vector = TRUE, diagOnly = T)
    ret <- rbind(ret, data.frame(Value = tmp, Rel = "WW", Type = "IBDe"))
  }
  # drones
  print("IBS DD")
  tmp <- diag(ibsDF[idDF$drones, idDF$drones])
  ret <- rbind(ret, data.frame(Value = tmp, Rel = "DD", Type = "IBS"))

  print("IBDr DD")
  tmp <- diag(ibdDF[idDF$drones, idDF$drones])
  ret <- rbind(ret, data.frame(Value = tmp, Rel = "DD", Type = "IBDr"))

  if (!is.null(Sinv)) {
    print("IBDe DD")
    tmp <- getS(Sinv, ids = idDF$drones, vector = TRUE, diagOnly = T)
    ret <- rbind(ret, data.frame(Value = tmp, Rel = "DD", Type = "IBDe"))
  }
}


prepareDataForPlotting_Queens <- function(ibsDF = NULL, ibdDF = NULL, Sinv = NULL,  idPopDF) {
  #Get populations IDs
  melID <- idPopDF$ID[idPopDF$Pop == "Mel"]
  melCrossID <- idPopDF$ID[idPopDF$Pop == "MelCross"]
  carID <- idPopDF$ID[idPopDF$Pop == "Car"]

  #Compute relationships between queens of different populations (QQ)
  #IBS
  IBSMelMelCross <- data.frame(Value = as.vector(list(ibsDF[melID, melCrossID])[[1]]),
                               Pops = "Mel_MelCross", Rel = "QQ")
  IBSMelCar <- data.frame(Value = as.vector(list(ibsDF[melID, carID])[[1]]),
                          Pops = "Mel_Car", Rel = "QQ")
  IBSMelCrossCar <- data.frame(Value = as.vector(list(ibsDF[melCrossID, carID])[[1]]),
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
  tmp <- ibsDF[melID, melID]
  IBSMel_Mel <- c(tmp[lower.tri(tmp, diag = FALSE)])
  IBSMel_Mel <- data.frame(Value = as.vector(list(IBSMel_Mel)[[1]]), Rel = "Q", Type = "IBS", Pops = "Mel")

  tmp <- ibsDF[carID, carID]
  IBSCar_Car <- c(tmp[lower.tri(tmp, diag = FALSE)])
  IBSCar_Car <-data.frame(Value = as.vector(list(IBSCar_Car)[[1]]), Rel = "Q", Type = "IBS", Pops = "Car")

  tmp <- ibsDF[melCrossID, melCrossID]
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
  inbIBS <- rbind(data.frame(Value = diag(ibsDF[melID, melID]),
                             Pops = "Mel"),
                  data.frame(Value = diag(ibsDF[melCrossID, melCrossID]),
                             Pops = "MelCross"),
                  data.frame(Value = diag(ibsDF[carID, carID]),
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


prepareDataForPlottingHeatMap_Queens <- function(ibsDF = NULL, ibdDF = NULL, Sinv = NULL, idDF = NULL) {
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

plotColony <- function(df, rel = c("WD", "WW", "DD"), type = c("IBDr", "IBDe"), x_axis = NULL) {
  df$Type <- factor(df$Type, levels = c("IBDe", "IBDr", "IBS"))
  df$Rel <- factor(df$Rel, levels= c("WW", "WD", "DD", "QW", "QD"))
  p <- ggplot(df[df$Rel %in% rel & df$Type %in% type, ],
                 aes(x = Value, fill = Rel)) +
    geom_histogram(binwidth = 0.01) + facet_grid(rows = vars(Type), scales = "free") + xlim(x_axis)
  plot <- p + scale_fill_manual("", values=cbPalette, aesthetics = c("colour","fill")) + theme_classic()
  return(plot)
}

plotColonyF <- function(df, rel = c("WW", "DD"), type = c("IBDr", "IBDe"), x_axis = NULL){
  df$Type <- factor(df$Type, levels = c("IBDe", "IBDr", "IBS"))
  df$Rel <- factor(df$Rel, levels= c("WW", "WD", "DD", "QW", "QD"))
  p <- ggplot(df[df$Rel %in% rel & df$Type %in% type, ],
              aes(x = Value, fill = Rel)) +
    geom_histogram(binwidth = 0.01) + facet_wrap(.~Type, scales = "free_y") + xlim(x_axis)
  plot <- p + scale_fill_manual("", values=cbPalette, aesthetics = c("colour","fill")) + theme_classic()
  return(plot)
}

plotQueensQQ <- function(df, rel = c("QQ"), type = c("IBDr", "IBDe", "IBS"), x_axis = c(-1.01, 2.01)) {
  df$Type <- factor(df$Type, levels = c("IBDe", "IBDr", "IBS"))
  df$Pops <- factor(df$Pops, levels = c("Mel_MelCross", "Mel_Car", "MelCross_Car"))
  p <- ggplot(df[df$Rel %in% rel & df$Type %in% type, ],
                   aes(x = Value, fill = Pops)) +
      geom_histogram(binwidth = 0.01) +
    facet_wrap(.~Type, scales = "free_y") + xlim(x_axis) #IBS not to scale so increased scale
  plot <- p + scale_fill_manual("", values=cbPalette, aesthetics = c("colour","fill")) + theme_classic()
  return(plot)
}

scatterQueens <-  function(df, rel = c("QQ"), type = c("IBDr", "IBS"), pops = c("Mel_MelCross", "Mel_Car", "MelCross_Car")) {
  df$Type <- factor(df$Type, levels = c("IBDe", "IBDr", "IBS"))
  df$Pops <- factor(df$Pops, levels = c("Mel_MelCross", "Mel_Car", "MelCross_Car", "Car", "Mel", "MelCross"))
  a <- pivot_wider(df[df$Rel %in% rel & df$Type %in% type & df$Pops %in% pops, ], names_from = Type, values_from = Value, values_fn = list)
  b <- unnest(a, cols = all_of(type))
  c <- ggplot(data = b, aes(x = IBDr, y = IBS)) + geom_point(aes(colour = Pops))
  plot <- c + scale_fill_manual("", values= cbPalette, aesthetics = c("colour","fill")) + theme_classic()
  return(plot)
}

plotQueensQ <- function(df, rel = c("Q", "F"), type = c("IBDr", "IBDe", "IBS"), x_axis = NULL) {
  df$Type <- factor(df$Type, levels = c("IBDe", "IBDr", "IBS"))
  df$Rel <- factor(df$Rel, levels= c("F", "Q"))
  p <- ggplot(df[df$Rel %in% rel & df$Type %in% type, ],
                   aes(x = Value, fill = Pops)) +
      geom_histogram(binwidth = 0.01) +
      facet_wrap(.~Type, scales = "free_y") + xlim(x_axis) #xlim removed
  plot <- p + scale_fill_manual("", values=cbPalette, aesthetics = c("colour","fill")) + theme_classic()
  return(plot)
}

plotQueensF <- function(df, rel = c("Q", "F"), type = c("IBDr", "IBDe", "IBS"), x_axis = c(-0.05, 0.80)) {
  df$Type <- factor(df$Type, levels = c("IBDe", "IBDr", "IBS"))
  df$Rel <- factor(df$Rel, levels= c("F", "Q"))
  p <- ggplot(df[df$Rel %in% rel & df$Type %in% type, ],
              aes(x = Value - 1, fill = Pops)) +
    geom_histogram(binwidth = 0.01) +
    facet_wrap(.~Type, scales = "free_y")+ xlim(x_axis) #xlim removed
  plot <- p + scale_fill_manual("", values=cbPalette, aesthetics = c("colour","fill")) + theme_classic()
  return(plot)
}

plotQueens_heatmap <- function(df, Pop = FALSE, PopIdDF = NULL) {
  df$ID <- as.factor(as.numeric(df$ID))
  df$name <- as.factor(as.numeric(df$name))
  df$Method <- factor(df$Method, levels = c("IBDe", "IBDr", "IBS"))

  if (Pop) {
    df <- merge(df, PopIdDF, by = "ID")
    df <- merge(df, PopIdDF, by.x = "name", by.y = "ID")

    df$PopId1 <- paste0(df$Pop.x, df$ID)
    df$PopId2 <- paste0(df$Pop.y, df$name)

    n <- length(unique(df$name[df$Pop.y == "Mel"]))/2
    breaks = list(df %>% group_by(Pop.y) %>% summarise(mean = unique(name)[n]) %>% unite(PopMean, Pop.y, mean, sep=""))[[1]]$PopMean
    plot <- ggplot(data = df, aes(x = PopId1, y = PopId2, fill = value)) + geom_tile() + scale_fill_gradient(low = "white", high = "blue") + theme_classic() +
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
#Plot carnica year 1
print("Plot carnica year 1")
relCar1 <- prepareDataForPlotting_Colony(ibsDF = ibsCar1, ibdDF = ibdCar1, Sinv = Sinv, idDF = idCar1)

#Carnica WW/WD/DD
plotCar1 <- plotColony(relCar1, type = c("IBDe", "IBDr", "IBS"), rel = c("WW", "WD", "DD"), x_axis = c(-0.01, 1.5))
plotCar1

#Carnica QW/QD
plotCar1Q <- plotColony(relCar1, type = c("IBDe", "IBDr", "IBS"), rel = c("QW", "QD"), x_axis = c(-0.01, 1.5))
plotCar1Q

#Plot Car Year 1 inbreeding
relCar1F <- prepareDataForPlotting_ColonyDiag(ibsDF = ibsCar1, ibdDF = ibdCar1, Sinv = Sinv, idDF = idCar1)
plotCar1F <- plotColonyF(relCar1F, rel = c("WW", "DD"), type = c("IBDe", "IBDr", "IBS"))
plotCar1F

# Plot CAR year 10
print("Plot carnica year 10")
relCar10 <- prepareDataForPlotting_Colony(ibsDF = ibsCar10, ibdDF = ibdCar10, Sinv = Sinv, idDF = idCar10)

#plotCar10 WW/WD/DD
plotCar10 <- plotColony(relCar10, type = c("IBDe", "IBDr","IBS" ), rel = c("WW", "WD", "DD"))
plotCar10

#Carnica Year 10 QW/QD
plotCar10Q <- plotColony(relCar10, type = c("IBDe", "IBDr", "IBS"), rel = c("QW", "QD"))
plotCar10Q

#Plot Car Year 10 inbreeding
relCar10F <- prepareDataForPlotting_ColonyDiag(ibsDF = ibsCar10, ibdDF = ibdCar10, Sinv = Sinv, idDF = idCar10)
plotCar10F <- plotColonyF(relCar10F, rel = c("WW", "DD"), type = c("IBDe", "IBDr", "IBS"))
plotCar10F

#Plot CAR csd Locus
#Year 1
relCar1_csd <- prepareDataForPlotting_Colony(ibsDF = ibsCar1_csd, ibdDF = ibdCar1_csd, Sinv = Sinv, idDF = idCar1)
plotCar1_csdLoc <- plotColony(relCar1_csd, type = c("IBDe", "IBDr", "IBS"))
#plotCar1_csdLoc

#Year 10
relCar10_csd <- prepareDataForPlotting_Colony(ibsDF = ibsCar10_csd, ibdDF = ibdCar10_csd, Sinv = Sinv, idDF = idCar10)
plotCar10_csdLoc <- plotColony(relCar10_csd, type = c("IBDe", "IBDr", "IBS"))
#plotCar10_csdLoc

#Plot CAR csd Chromosome
#Year 1
relCar1_csdChr <- prepareDataForPlotting_Colony(ibsDF = ibsCar1_csdChr, ibdDF = ibdCar1_csdChr, Sinv = Sinv, idDF = idCar1)
plotCar1_csdChr <- plotColony(relCar1_csdChr, type = c("IBDe", "IBDr", "IBS"))
#plotCar1_csdChr

#Year 10
relCar10_csdChr <- prepareDataForPlotting_Colony(ibsDF = ibsCar10_csdChr, ibdDF = ibdCar10_csdChr, Sinv = Sinv, idDF = idCar10)
plotCar10_csdChr <- plotColony(relCar10_csdChr, type = c("IBDe", "IBDr", "IBS"))
#plotCar10_csdChr

#Save pdfs
pdf("Plot_Carnica1.pdf")
plotCar1
dev.off()

pdf("Plot_Carnica1Q.pdf")
plotCar1Q
dev.off()

pdf("Plot_Carnica1F.pdf")
plotCar1F
dev.off()


pdf("Plot_Carnica10.pdf")
plotCar10
dev.off()

pdf("Plot_Carnica10Q.pdf")
plotCar10Q
dev.off()

pdf("Plot_Carnica10F.pdf")
plotCar10F
dev.off()

pdf("Plot_Carnica1_csdLoc.pdf")
plotCar1_csdLoc
dev.off()

pdf("Plot_Carnica10_csdLoc.pdf")
plotCar10_csdLoc
dev.off()

pdf("Plot_Carnica1_csdChr.pdf")
plotCar1_csdChr
dev.off()

pdf("Plot_Carnica10_csdChr.pdf")
plotCar10_csdChr
dev.off()


# #Plot Only workers in year 1 and year 10
# WW1 <- rel1[rel1$Rel %in% c("WW") & rel1$Type %in% c("IBDr", "IBDe"), ]
# WW1$Year <- 1
# WW10 <- rel10[rel10$Rel %in% c("WW") & rel10$Type %in% c("IBDr", "IBDe"), ]
# WW10$Year <- 10
# WW1_10 <- rbind(WW1, WW10)
# plotWW_1_10 <- ggplot(WW1_10,
#                  aes(x = Value, fill = Type)) + geom_histogram(binwidth = 0.01, position = "identity") +
#   facet_grid(cols = vars(Year), scales = "free") + xlim(c(-0.01, 2.01))
# plotWW_1_10

########################################################
### --- FIGURE 4: Between queens of different subspecies (carnica vs. mellifera) ---###
########################################################
#Plot queens Year 1
print("Plot queens year 1")
relQueens1 <- prepareDataForPlotting_Queens(ibsDF = ibsQueens1, ibdDF = ibdQueens1, Sinv = Sinv, idPopDF = idPopQueens1)

plotQueens1QQ <- plotQueensQQ(relQueens1, type = c("IBDr", "IBDe", "IBS"))
#plotQueens1QQ

plotQueensScatter1 <- scatterQueens(relQueens1, type = c("IBDr", "IBS"), pops = c("Mel_MelCross", "Mel_Car", "MelCross_Car"))
#plotQueensScatter1

plotQueens1Q <- plotQueensQ(relQueens1, rel = c("Q"), type = c("IBDr", "IBDe", "IBS"))
#plotQueens1Q

plotQueens1F <- plotQueensF(relQueens1, rel = c("F"), type = c("IBDr", "IBDe", "IBS"))
plotQueens1F

relQueens1h <- prepareDataForPlottingHeatMap_Queens(ibsDF = ibsQueens1, ibdDF = ibdQueens1, Sinv = Sinv, idDF = idQueens1)
plotQueens1h <- plotQueens_heatmap(relQueens1h, Pop = TRUE, PopIdDF = idPopQueens1)
#plotQueens1h

#Plot csd Year 1
relQueens1_csd_hist <- prepareDataForPlotting_Queens(ibsDF = ibsQueens1_csd, ibdDF = ibdQueens1_csd, Sinv = Sinv, idPopDF = idPopQueens1)
plotQueens1_csdLoc<- plotQueensQQ(relQueens1_csd_hist, type = c("IBDr", "IBDe", "IBS"))
#plotQueens1_csdLoc


relQueens1_csd <- prepareDataForPlottingHeatMap_Queens(ibsDF = ibsQueens1_csd, ibdDF = ibdQueens1_csd, Sinv = Sinv, idDF = idQueens1)
plotQueens1_csdLoc_heat <- plotQueens_heatmap(relQueens1_csd, Pop = TRUE, PopIdDF = idPopQueens1)
#plotQueens1_csdLoc_heat

#Plot csd chr Year 1
relQueens1_csdChr <- prepareDataForPlottingHeatMap_Queens(ibsDF = ibsQueens1_csdChr, ibdDF = ibdQueens1_csdChr, Sinv = Sinv, idDF = idQueens1)
plotQueens1_csdChr_heat <- plotQueens_heatmap(relQueens1_csdChr, Pop = TRUE, PopIdDF = idPopQueens1)
#plotQueens1_csdChr_heat



#Save pdfs
pdf("Plot_Queens1QQ.pdf")
plotQueens1QQ
dev.off()

jpeg("PlotQueensScatter1.jpeg")
plotQueensScatter1
dev.off()

pdf("Plot_Queens1Q.pdf")  #Q = non-diagonal
plotQueens1Q
dev.off()

pdf("Plot_Queens1F.pdf")  #F = inbreeding
plotQueens1F
dev.off()

jpeg("Plot_Queens1h.jpeg") #h = heatmaps
plotQueens1h
dev.off()

pdf("PlotQueens1_csdLoc.pdf")
plotQueens1_csdLoc
dev.off()

jpeg("PlotQueens1_csdLoc_heat.jpeg")
plotQueens1_csdLoc_heat
dev.off()

jpeg("Plot_Queens1_csdChr_heat.jpeg")
plotQueens1_csdChr_heat
dev.off()


# --- Year 10 ---#
#Plot queens Year 10
print("Plot queens year 10")
relQueens10 <- prepareDataForPlotting_Queens(ibsDF = ibsQueens10, ibdDF = ibdQueens10, Sinv = Sinv, idPopDF = idPopQueens10)

plotQueens10QQ <- plotQueensQQ(relQueens10, type = c("IBDe", "IBDr", "IBS"))
#plotQueens10QQ

plotQueensScatter10 <- scatterQueens(relQueens10, type = c("IBDr", "IBS"), pops = c("Mel_MelCross", "Mel_Car", "MelCross_Car"))
#plotQueensScatter10

plotQueens10Q <- plotQueensQ(relQueens10, rel = c("Q"), type = c("IBDr", "IBDe", "IBS"))
#plotQueens10Q

plotQueens10F <- plotQueensF(relQueens10, rel = c("F"), type = c("IBDr", "IBDe", "IBS"))
plotQueens10F

relQueens10h <- prepareDataForPlottingHeatMap_Queens(ibsDF = ibsQueens10, ibdDF = ibdQueens10, Sinv = Sinv, idDF = idQueens10)
plotQueens10h <- plotQueens_heatmap(relQueens10h, Pop = T, PopIdDF = idPopQueens10)
#plotQueens10h

#Plot csd Year 10
relQueens10_csd <- prepareDataForPlottingHeatMap_Queens(ibsDF = ibsQueens10_csd, ibdDF = ibdQueens10_csd, Sinv = Sinv, idDF = idQueens10)
plotQueens10_csdLoc_heat <- plotQueens_heatmap(relQueens10_csd, Pop = TRUE, PopIdDF = idPopQueens10)
#plotQueens10_csdLoc_heat

relQueens10_csd_hist <- prepareDataForPlotting_Queens(ibsDF = ibsQueens10_csd, ibdDF = ibdQueens10_csd, Sinv = Sinv, idPopDF = idPopQueens10)
plotQueens10_csdLoc <- plotQueensQQ(relQueens10_csd_hist, type = c("IBDr", "IBDe", "IBS"))
#plotQueens10_csdLoc

#Plot csd chr Year 10
relQueens10_csdChr <- prepareDataForPlottingHeatMap_Queens(ibsDF = ibsQueens10_csdChr, ibdDF = ibdQueens10_csdChr, Sinv = Sinv, idDF = idQueens10)
plotQueens10_csdChr_heat <-plotQueens_heatmap(relQueens10_csdChr, Pop = TRUE, PopIdDF = idPopQueens10)
#plotQueens10_csdChr

#save pdfs
pdf("Plot_Queens10iQQ.pdf")
plotQueens10QQ
dev.off()

jpeg("PlotQueensScatter10.jpeg")
plotQueensScatter10
dev.off()

pdf("Plot_Queens10Q.pdf")  #Q = non-diagonal
plotQueens10Q
dev.off()

pdf("Plot_Queens10F.pdf")  #F = inbreeding
plotQueens10F
dev.off()

jpeg("Plot_Queens10h.jpeg")
plotQueens10h
dev.off()

jpeg("Plot_Queens10_csdLoc_heat.jpeg")
plotQueens10_csdLoc_heat
dev.off()

pdf("PlotQueens10_csdLoc.pdf")
plotQueens10_csdLoc
dev.off()

jpeg("Plot_Queens10_csdChr_heat.jpeg")
plotQueens10_csdChr_heat
dev.off()

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


