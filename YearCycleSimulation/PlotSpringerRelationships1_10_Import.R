#Laura's wd :
##setwd("~/Desktop/GitHub/lstrachan_honeybee_sim/YearCycleSimulation")
# Plot the relationships from the springer simulation
#Laura's data file :
#data <- load("~/github//lstrachan_honeybee_sim/YearCycleSimulation/SpringerSimulation.Rdata")
#data <- load("~/EddieDir/YearCycleSimulation/lstrachan_honeybee_sim/YearCycleSimulation/SpringerSimulation_import.Rdata")
print("Reading in the data")
data <- load("SpringerSimulation_import_objects.RData")
#save.image("~/Documents/")

# The data contains two populations - mellifera and carnica
# The carnica stays "pure" throughout the simulation
# Mellifera gets mated with a proportion of carnica drones

print("Assigning objects")
library(ggplot2)
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

# Pedigree
IBDe <- IBDe



# Plotting functions
prepareDataForPlotting_Colony <- function(ibsDF = NULL, ibdDF = NULL, pedDF = NULL,  idDF) {
  #IBS
  tmp <- ibsDF[idDF$workers, idDF$workers]
  IBS_WW1 <- c(tmp[lower.tri(tmp, diag = TRUE)])
  ret <- data.frame(Value = IBS_WW1, Rel = "WW", Type = "IBS")

  tmp <- ibdDF[idDF$workers, idDF$workers]
  IBDr_WW1 <- c(tmp[lower.tri(tmp, diag = TRUE)])
  ret <- rbind(ret, data.frame(Value = IBDr_WW1, Rel = "WW", Type = "IBDr"))

  tmp <- as.matrix(pedDF[idDF$workers, idDF$workers])
  IBDe_WW1 <- c(tmp[lower.tri(tmp, diag = TRUE)])
  ret <- rbind(ret, data.frame(Value = IBDe_WW1, Rel = "WW", Type = "IBDe"))

  # workers vs drones
  IBS_WD1 <- c(ibsDF[idDF$workers, idDF$drones])
  ret <- rbind(ret, data.frame(Value = IBS_WD1, Rel = "WD", Type = "IBS"))

  IBDr_WD1 <- c(ibdDF[idDF$workers, idDF$drones])
  ret <- rbind(ret, data.frame(Value = IBDr_WD1, Rel = "WD", Type = "IBDr"))

  IBDe_WD1 <- c(as.matrix(pedDF[idDF$workers, idDF$drones]))
  ret <- rbind(ret, data.frame(Value = IBDe_WD1, Rel = "WD", Type = "IBDe"))

  return(ret)
}

prepareDataForPlotting_Queens <- function(ibsDF = NULL, ibdDF = NULL, pedDF = NULL,  idDF) {
  #IBS
  IBS <- c(tmp[lower.tri(ibsDF, diag = TRUE)])
  ret <- data.frame(Value = IBS, Rel = "QQ", Type = "IBS")

  IBDr <- c(tmp[lower.tri(ibdDF, diag = TRUE)])
  ret <- rbind(ret, data.frame(Value = IBDr, Rel = "QQ", Type = "IBDr"))

  IBDe <- c(tmp[lower.tri(pedDF, diag = TRUE)])
  ret <- rbind(ret, data.frame(Value = IBDe, Rel = "QQ", Type = "IBDe"))

  return(ret)
}

plotColony <- function(df, rel = c("WD", "WW"), type = c("IBDr", "IBDe")) {
  plot <- ggplot(df[df$Rel %in% rel & df$Type %in% type, ],
         aes(x = Value, fill = Type)) +
    geom_histogram(binwidth = 0.01, position = "identity") +
    facet_grid(cols = vars(Rel)) + xlim(c(-0.01, 2.01))
  return(plot)
}

plotQueens <- function(df, rel = c("QQ"), type = c("IBDr", "IBDe")) {
  plot <- ggplot(df[df$Rel %in% rel & df$Type %in% type, ],
                 aes(x = Value, fill = Type)) +
    geom_histogram(binwidth = 0.01, position = "identity") +
    facet_grid(cols = vars(Rel)) + xlim(c(-0.01, 2.01))
  return(plot)
}

print("Plotting")
########################################################
### --- FIGURE 1&2: Pure subspecies (carnica) in years 1 and 10 ---###
########################################################
# Plot year 1
relCar1 <- prepareDataForPlotting_Colony(ibsDF = ibsCar1, ibdDF = ibdCar1, pedDF = IBDe, idDF = idCar1)
plotCar1 <- plotColony(relCar1, type = "IBDr")
pdf("Plot_Carnica1.pdf")
plotCar1
dev.off()

# Plot year 10
relCar10 <- prepareDataForPlotting_Colony(ibsDF = ibsCar10, ibdDF = ibdCar10, pedDF = IBDe, idDF = idCar10)
plotCar10 <- plotColony(relCar10)
pdf("Plot_Carnica10.pdf")
plotCar10
dev.off()

#
# #Plot Only workers in year 1 and year 10
# WW1 <- rel1[rel1$Rel %in% c("WW") & rel1$Type %in% c("IBDr", "IBDe"), ]
# WW1$Year <- 1
# WW10 <- rel10[rel10$Rel %in% c("WW") & rel10$Type %in% c("IBDr", "IBDe"), ]
# WW10$Year <- 10
# WW1_10 <- rbind(WW1, WW10)
#
# plotWW_1_10 <- ggplot(WW1_10,
#                  aes(x = Value, fill = Type)) + geom_histogram(binwidth = 0.01, position = "identity") +
#   facet_grid(cols = vars(Year), scales = "free") + xlim(c(-0.01, 2.01))
#
# plotWW_1_10

#Plot csd Year 10
relCar10_csd <- prepareDataForPlotting_Colony(ibsDF = ibsCar10_csd, ibdDF = ibdCar10_csd, pedDF = IBDe, idDF = idCar10)
plotCar10_csd <- plotColony(relCar10_csd)
pdf("Plot_Carnica10_csd.pdf")
plotCar10_csd
dev.off()

#Plot csd chr Year 10
relCar10_csdChr <- prepareDataForPlotting_Colony(ibsDF = ibsCar10_csdChr, ibdDF = ibdCar10_csdChr, pedDF = IBDe, idDF = idCar10)
plotCar10_csdChr <- plotColony(relCar10_csdChr)
pdf("Plot_Carnica10_csdChr.pdf")
plotCar10_csdChr
dev.off()

########################################################
### --- FIGURE 4: Between queens of different subspecies (carnica vs. mellifera) ---###
########################################################
#Plot queens Year 10
relQueens1 <- prepareDataForPlotting_Queens(ibsDF = ibsQueens1, ibdDF = ibdQueens1, pedDF = IBDe, idDF = idQueens1)
plotQueens1 <- plotQueens(relQueens1)
pdf("Plot_Queens1.pdf")
plotQueens1
dev.off()

#Plot csd chr Year 10
relQueens10 <- prepareDataForPlotting_Queens(ibsDF = ibsQueens10, ibdDF = ibdQueens10, pedDF = IBDe, idDF = idQueens10)
plotQueens10 <- plotQueens(relQueens10)
plot("Plot_Queens10.pdf")
plotQueens10
dev.off()

