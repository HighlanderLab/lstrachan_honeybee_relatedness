library(dplyr)
library(tidyr)
library(ggplot2)
library(Matrix)
library(SIMplyBee)

print("Reading in the data")
#Laura's laptop data 
data <- load("/Users/s2122596/Desktop/GitHub/Data /SpringerSimulation_import_objects.RData")
Sinv <- readMM("/Users/s2122596/Desktop/GitHub/Data /Sinv2.mm")

#Eddie data 
#data <- load("SpringerSimulation_import_objects.RData")
#Sinv <- readMM("Sinv2.mm")
#save.image("~/Documents/")

# The data contains two populations - mellifera and carnica
# The carnica stays "pure" throughout the simulation
# Mellifera gets mated with a proportion of carnica drones

#Plotting help 
# The colourblind palette with grey:
cbPalette <- c("#E69F00", "#0072B2", "#56B4E9", "#009E73", "#F0E442", "#D55E00", "#CC79A7")
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


getS <- function(Sinv, ids, with = ids, diagOnly = F, vector = F) {
  ids <- as.numeric(ids)
  with <- as.numeric(with)
  x <- sparseMatrix(i = ids, j = 1:length(ids), dims = c(nrow(Sinv), length(ids)))
  M1 <- as(x, "dgCMatrix")
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
  print("IBS WW")
  tmp <- ibsDF[idDF$workers, idDF$workers]
  IBS_WW1 <- c(tmp[lower.tri(tmp, diag = TRUE)])
  ret <- data.frame(Value = IBS_WW1, Rel = "WW", Type = "IBS")
  
  print("IBDr WW")
  tmp <- ibdDF[idDF$workers, idDF$workers]
  IBDr_WW1 <- c(tmp[lower.tri(tmp, diag = TRUE)])
  ret <- rbind(ret, data.frame(Value = IBDr_WW1, Rel = "WW", Type = "IBDr"))
  
  if (!is.null(Sinv)) {
    print("IBDe WW")
    ret <- rbind(ret, data.frame(Value = getS(Sinv, ids = idDF$workers, vector = TRUE), 
                                 Rel = "WW", Type = "IBDe"))
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
  
  return(ret)
}

prepareDataForPlotting_Queens <- function(ibsDF = NULL, ibdDF = NULL, Sinv = NULL,  idPopDF) {
  #Get populations IDs
  melID <- idPopDF$ID[idPopDF$Pop == "Mel"]
  melCrossID <- idPopDF$ID[idPopDF$Pop == "MelCross"]
  carID <- idPopDF$ID[idPopDF$Pop == "Car"]
  
  #Compute relationships between queens of different populations (QQ)
  #IBS
  IBSMelMelCross <- data.frame(Value = as.vector(list(ibsDF[melID, melCrossID])[[1]]),
                               Pops = "Mel_MelCross")
  IBSMelCar <- data.frame(Value = as.vector(list(ibsDF[melID, carID])[[1]]),
                          Pops = "Mel_Car")
  IBSMelCrossCar <- data.frame(Value = as.vector(list(ibsDF[melCrossID, carID])[[1]]),
                               Pops = "MelCross_Car")
  IBS <-rbind(IBSMelMelCross, IBSMelCar, IBSMelCrossCar)
  IBS$Type <- "IBS"
  IBS$Rel <- "QQ"
  
  #IBDr
  IBDrMelMelCross <- data.frame(Value = as.vector(list(ibdDF[melID, melCrossID])[[1]]),
                                Pops = "Mel_MelCross")
  IBDrMelCar <- data.frame(Value = as.vector(list(ibdDF[melID, carID])[[1]]),
                           Pops = "Mel_Car")
  IBDrMelCrossCar <- data.frame(Value = as.vector(list(ibdDF[melCrossID, carID])[[1]]),
                                Pops = "MelCross_Car")
  IBDr <- rbind(IBDrMelMelCross, IBDrMelCar, IBDrMelCrossCar)
  IBDr$Type <- "IBDr"
  IBDr$Rel <- "QQ"
  
  #IBDe
  if (!is.null(Sinv)) {
    IBDeMelMelCross <- data.frame(Value = getS(Sinv, ids = melID, with = melCrossID, vector = T),
                                  Pops = "Mel_MelCross")
    IBDeMelCar <- data.frame(Value = getS(Sinv, ids = melID, with = carID, vector = T),
                             Pops = "Mel_Car")
    IBDeMelCrossCar <- data.frame(Value = getS(Sinv, ids = melCrossID, with = carID, vector = T),
                                  Pops = "MelCross_Car")
    IBDe <- rbind(IBDeMelMelCross, IBDeMelCar, IBDeMelCrossCar)
    IBDe$Type <- "IBDe"
    IBDe$Rel <- "QQ"
  }
  
  # Inbreeding (diagonal!!!) if queens (Q)
  inbIBS <- rbind(data.frame(Value = diag(ibsDF[melID, melID]),
                             Pops = "Mel"),
                  data.frame(Value = diag(ibsDF[melCrossID, melCrossID]),
                             Pops = "MelCross"),
                  data.frame(Value = diag(ibsDF[carID, carID]),
                             Pops = "Car"))
  inbIBS$Type = "IBS"
  inbIBS$Rel = "Q"
  IBS <- rbind(IBS, inbIBS)
  
  inbIBDr <- rbind(data.frame(Value = diag(ibdDF[melID, melID]),
                              Pops = "Mel"),
                   data.frame(Value = diag(ibdDF[melCrossID, melCrossID]),
                              Pops = "MelCross"),
                   data.frame(Value = diag(ibdDF[carID, carID]),
                              Pops = "Car"))
  inbIBDr$Type = "IBDr"
  inbIBDr$Rel = "Q"
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
    inbIBDe$Rel = "Q"
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

plotColony <- function(df, rel = c("WD", "WW"), type = c("IBDr", "IBDe")) {
  df$Type <- factor(df$Type, levels = c("IBDe", "IBDr", "IBS"))
  df$Rel <- factor(df$Rel, levels= c("WW", "WD"))
  p <- ggplot(df[df$Rel %in% rel & df$Type %in% type, ],
                 aes(x = Value, fill = Rel)) +
    geom_histogram(binwidth = 0.01, position = "identity") + facet_grid(rows = vars(Type), scales = "free") + xlim(c(-0.01, 2.01))
  plot <- p + scale_fill_manual(values=cbPalette, aesthetics = c("colour","fill")) + theme_classic()
  return(plot)
}

plotQueens <- function(df, rel = c("QQ"), type = c("IBDr", "IBDe"), plot = "histogram") {
  df$Type <- factor(df$Type, levels = c("IBDe", "IBDr", "IBS"))
  if (plot == "histogram") {
    p <- ggplot(df[df$Rel %in% rel & df$Type %in% type, ],
                   aes(x = Value, fill = Rel)) +
      geom_histogram(binwidth = 0.01, position = "identity") +
      facet_grid(cols = vars(Type)) + xlim(c(-0.01, 2.01))
  } else if (plot == "density") {
    p <- ggplot(df[df$Rel %in% rel & df$Type %in% type, ]) + 
      stat_density(aes(x=Value, y=..scaled..), position="dodge", geom="line") +
      facet_grid(cols = vars(Type))
  }
  plot <- p + scale_fill_manual(values=cbPalette, aesthetics = c("colour","fill")) + theme_classic()
  return(plot)
}

plotQueensF <- function(df, rel = "Q", type = c("IBDr", "IBDe"), plot = "histogram") {
  df$Type <- factor(df$Type, levels = c("IBDe", "IBDr", "IBS"))
  if (plot == "histogram") {
    p <- ggplot(df[df$Rel %in% rel & df$Type %in% type, ],
                   aes(x = Value, fill = Rel)) +
      geom_histogram(binwidth = 0.01, position = "identity") +
      facet_grid(cols = vars(Type)) + xlim(c(-0.01, 2.01))
  } else if (plot == "density") {
    p <- ggplot(df[df$Rel %in% rel & df$Type %in% type, ]) + 
      stat_density(aes(x=Value, color=Rel), position="dodge", geom="line") +
      facet_grid(cols = vars(Type)) + theme_by(base_size = 18)
  }
  plot <- p + scale_fill_manual(values=cbPalette, aesthetics = c("colour","fill")) + theme_classic()
  return(plot)
}

plotQueens_heatmap <- function(df, Pop = FALSE, PopIdDF = NULL) {
  df$ID <- as.factor(as.numeric(df$ID))
  df$name <- as.factor(as.numeric(df$name))
  df$Method <- factor(df$Method, levels = c("IBDe", "IBDr", "IBD"))
  
  if (Pop) {
    df <- merge(df, PopIdDF, by = "ID")
    df <- merge(df, PopIdDF, by.x = "name", by.y = "ID")
    
    df$PopId1 <- paste0(df$Pop.x, df$ID)
    df$PopId2 <- paste0(df$Pop.y, df$name)
    
    n <- length(unique(df$name[df$Pop.y == "Mel"]))/2
    breaks = list(df %>% group_by(Pop.y) %>% summarise(mean = unique(name)[n]) %>% unite(PopMean, Pop.y, mean, sep=""))[[1]]$PopMean
    plot <- ggplot(data = df, aes(x = PopId1, y = PopId2, fill = cbPalette)) + geom_tile() +
      scale_x_discrete(breaks = breaks, labels = c("Car", "Mel", "MelCross")) +
      scale_y_discrete(breaks = breaks, labels = c("Car", "Mel", "MelCross")) +
      xlab("") + ylab("") +
      facet_grid(rows = vars(Method))
  } else {
    plot <- ggplot(data = df, aes(x = ID, y = name, fill = cbPalette)) + geom_tile() +
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

plotCar1ibd <- plotColony(relCar1, type = c("IBDe", "IBDr"))
plotCar1ibd

plotCar1ibd_ibs <- plotColony(relCar1, type = c("IBDr", "IBS"))
plotCar1ibd_ibs

# Plot CAR year 10
print("Plot carnica year 10") 
relCar10 <- prepareDataForPlotting_Colony(ibsDF = ibsCar10, ibdDF = ibdCar10, Sinv = Sinv, idDF = idCar10)

plotCar10ibd <- plotColony(relCar1, type = c("IBDe", "IBDr"))
plotCar10ibd

plotCar10ibd_ibs<- plotColony(relCar10, type = c("IBDr", "IBS"))
plotCar10ibd_ibs

#Plot CAR csd Year 10
relCar10_csd <- prepareDataForPlotting_Colony(ibsDF = ibsCar10_csd, ibdDF = ibdCar10_csd, Sinv = Sinv, idDF = idCar10)
plotCar10_csd <- plotColony(relCar10_csd)
plotCar10_csd

#Plot CAR csd chr Year 10
relCar10_csdChr <- prepareDataForPlotting_Colony(ibsDF = ibsCar10_csdChr, ibdDF = ibdCar10_csdChr, Sinv = Sinv, idDF = idCar10)
plotCar10_csdChr <- plotColony(relCar10_csdChr)
plotCar10_csdChr

#Save pdfs
pdf("Plot_Carnica1ibd.pdf")
plotCar1ibd
dev.off()

pdf("Plot_Carnica1ibd_ibs.pdf")
plotCar1ibd_ibs
dev.off()

pdf("Plot_Carnica10ibd.pdf")
plotCar10ibd
dev.off()

pdf("Plot_Carnica10ibd_ibs.pdf")
plotCar10ibd_ibs
dev.off()

pdf("Plot_Carnica10_csd.pdf")
plotCar10_csd
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

plotQueens1ibd <- plotQueens(relQueens1, type = c("IBDr", "IBDe"), plot = "histogram")
plotQueens1ibd

plotQueens1ibs <- plotQueens(relQueens1, type = c("IBS"), plot = "histogram")
plotQueens1ibs

plotQueens1Fibd <- plotQueensF(relQueens1, type = c("IBDr", "IBDe"), plot = "histogram")
plotQueens1Fibd

plotQueens1Fibs <- plotQueensF(relQueens1, type = c("IBS"), plot = "histogram")
plotQueens1Fibs

relQueens1h <- prepareDataForPlottingHeatMap_Queens(ibsDF = ibsQueens1, ibdDF = ibdQueens1, Sinv = Sinv, idDF = idQueens1)
plotQueens1h <- plotQueens_heatmap(relQueens1h, Pop = TRUE, PopIdDF = idPopQueens1)
plotQueens1h

#Plot csd Year 1
relQueens1_csd <- prepareDataForPlottingHeatMap_Queens(ibsDF = ibsQueens1_csd, ibdDF = ibdQueens1_csd, Sinv = Sinv, idDF = idQueens1)
plotQueens1_csd <- plotQueens_heatmap(relQueens1_csd)
plotQueens1_csd

#Plot csd chr Year 1
plotQueens1_csdChr <- plotQueens_heatmap(relQueens1_csdChr)
plotQueens1_csdChr

#save pdfs
#Save pdfs
pdf("Plot_Queens1ibd.pdf")
plotQueens1ibd
dev.off()

pdf("Plot_Queens1ibs.pdf")
plotQueens1ibs
dev.off()

pdf("Plot_Queens1Fibd.pdf")  #F = inbreeding 
plotQueens1Fibd
dev.off()

pdf("Plot_Queens1Fibs.pdf")  
plotQueens1Fibs
dev.off()

pdf("Plot_Queens1h.pdf") #h = heatmaps 
plotQueens1h
dev.off()

pdf("Plot_Queens1_csd.pdf")
plotQueens1_csd
dev.off()

pdf("Plot_Queens1_csdChr.pdf")
plotQueens1_csdChr
dev.off()



# --- Year 10 ---#
#Plot queens Year 10
print("Plot queens year 1") 

relQueens10 <- prepareDataForPlotting_Queens(ibsDF = ibsQueens10, ibdDF = ibdQueens10, Sinv = Sinv, idPopDF = idPopQueens10)
plotQueens10ibd <- plotQueens(relQueens10, type = c("IBDr", "IBDe"), plot = "histogram")
plotQueens10ibd

plotQueens10ibs <- plotQueens(relQueens10, type = c("IBS"), plot = "histogram")
plotQueens10ibs

plotQueens10Fibd <- plotQueensFibd(relQueens10, type = c("IBDr", "IBDe"), plot = "histogram")
plotQueens10Fibd

plotQueens10Fibs <- plotQueensFibs(relQueens10, type = c( "IBS"), plot = "histogram")
plotQueens10Fibs

relQueens10h <- prepareDataForPlottingHeatMap_Queens(ibsDF = ibsQueens10, ibdDF = ibdQueens10, Sinv = Sinv, idDF = idQueens10)
plotQueens10h <- plotQueens_heatmap(relQueens10h, Pop = T, PopIdDF = idPopQueens10)
plotQueens10h

#Plot csd Year 10
relQueens10_csd <- prepareDataForPlottingHeatMap_Queens(ibsDF = ibsQueens10_csd, ibdDF = ibdQueens10_csd, Sinv = Sinv, idDF = idQueens10)
plotQueens10_csd <- plotQueens_heatmap(relQueens10_csd)


#Plot csd chr Year 10
relQueens10_csdChr <- prepareDataForPlottingHeatMap_Queens(ibsDF = ibsQueens10_csdChr, ibdDF = ibdQueens10_csdChr, Sinv = Sinv, idDF = idQueens10)
plotQueens10_csdChr <-plotQueens_heatmap(relQueens10_csdChr)

#save pdfs
pdf("Plot_Queens10ibd.pdf")
plotQueens10ibd
dev.off()

pdf("Plot_Queens10ibs.pdf")
plotQueens10ibs
dev.off()

pdf("Plot_Queens10Fibd.pdf")
plotQueens10Fibd
dev.off()

pdf("Plot_Queens10Fibs.pdf")
plotQueens10Fibs
dev.off()

pdf("Plot_Queens10h.pdf")
plotQueens10h
dev.off()

pdf("Plot_Queens10_csd.pdf")
plotQueens10_csd
dev.off()

pdf("Plot_Queens10_csdChr.pdf")
plotQueens10_csdChr
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


