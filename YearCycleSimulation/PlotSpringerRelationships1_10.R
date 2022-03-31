#Laura's wd :
##setwd("~/Desktop/GitHub/lstrachan_honeybee_sim/YearCycleSimulation")
# Plot the relationships from the springer simulation
#Laura's data file :
data <- load("~/github//lstrachan_honeybee_sim/YearCycleSimulation/SpringerSimulation.Rdata")

library(ggplot2)
ped <- SP$pedigree
caste <- SP$caste
# Colony in year 1
ibs1 <- ibs_springerColony1
ibd1 <- ibd_springerColony1
id1 <- springerColony1_id

# Colony in year 10
ibs10 <- ibs_springerColony10
ibd10 <- ibd_springerColony10
id10 <- springerColony10_id

#csd in year 10
ibs10_csd <- ibs_springerColony10csd
ibd10_csd <- ibd_springerColony10csd

#chr3 in year10
ibs10_chr3 <- ibs_springerColony10chr3
ibd10_chr3 <- ibd_springerColony10chr3

# Pedigree
IBDe <- IBDr$S

# Plot year 1
# workers vs workers
tmp <- ibs1[id1$workers, id1$workers]
IBS_WW1 <- c(tmp[lower.tri(tmp, diag = TRUE)])
rel1 <- data.frame(Value = IBS_WW1, Rel = "WW", Type = "IBS")

tmp <- ibd1[id1$workers, id1$workers]
IBDr_WW1 <- c(tmp[lower.tri(tmp, diag = TRUE)])
rel1 <- rbind(rel1, data.frame(Value = IBDr_WW1, Rel = "WW", Type = "IBDr"))

tmp <- as.matrix(IBDe[id1$workers, id1$workers])
IBDe_WW1 <- c(tmp[lower.tri(tmp, diag = TRUE)])
rel1 <- rbind(rel1, data.frame(Value = IBDe_WW1, Rel = "WW", Type = "IBDe"))

# workers vs drones
IBS_WD1 <- c(ibs1[id1$workers, id1$drones])
rel1 <- rbind(rel1, data.frame(Value = IBS_WD1, Rel = "WD", Type = "IBS"))

IBDr_WD1 <- c(ibd1[id1$workers, id1$drones])
rel1 <- rbind(rel1, data.frame(Value = IBDr_WD1, Rel = "WD", Type = "IBDr"))

IBDe_WD1 <- c(as.matrix(IBDe[id1$workers, id1$drones]))
rel1 <- rbind(rel1, data.frame(Value = IBDe_WD1, Rel = "WD", Type = "IBDe"))

#write c("IBDr", "IBDe") if you want both on plot
plot1 <- ggplot(rel1[rel1$Rel %in% c("WD", "WW") & rel1$Type %in% c("IBDe", "IBDr"), ],
                aes(x = Value, fill = Type)) + geom_histogram(binwidth = 0.01) +
  facet_grid(cols = vars(Rel)) + xlim(c(-0.01, 2.01))

plot1

# Plot year 10
# workers vs workers
tmp <- ibs10[id10$workers, id10$workers]
IBS_WW10 <- c(tmp[lower.tri(tmp, diag = TRUE)])
rel10 <- data.frame(Value = IBS_WW10, Rel = "WW", Type = "IBS")

tmp <- ibd10[id10$workers, id10$workers]
IBDr_WW10 <- c(tmp[lower.tri(tmp, diag = TRUE)])
rel10 <- rbind(rel10, data.frame(Value = IBDr_WW10, Rel = "WW", Type = "IBDr"))

tmp <- as.matrix(IBDe[id10$workers, id10$workers])
IBDe_WW10 <- c(tmp[lower.tri(tmp, diag = TRUE)])
rel10 <- rbind(rel10, data.frame(Value = IBDe_WW10, Rel = "WW", Type = "IBDe"))

# workers vs drones
IBS_WD10 <- c(ibs10[id10$workers, id10$drones])
rel10 <- rbind(rel10, data.frame(Value = IBS_WD10, Rel = "WD", Type = "IBS"))

IBDr_WD10 <- c(ibd10[id10$workers, id10$drones])
rel10 <- rbind(rel10, data.frame(Value = IBDr_WD10, Rel = "WD", Type = "IBDr"))

IBDe_WD10 <- c(as.matrix(IBDe[id10$workers, id10$drones]))
rel10 <- rbind(rel10, data.frame(Value = IBDe_WD10, Rel = "WD", Type = "IBDe"))

#write c("IBDr", "IBDe") if you want both on plot
plot10 <- ggplot(rel10[rel10$Rel %in% c("WD", "WW") & rel10$Type %in% c("IBDr", "IBDe"), ],
                aes(x = Value, fill = Type)) + geom_histogram(binwidth = 0.01, position = "identity") +
  facet_grid(cols = vars(Rel)) + xlim(c(-0.01, 2.01))

plot10


#Plot Only workers in year 1 and year 10
WW1 <- rel1[rel1$Rel %in% c("WW") & rel1$Type %in% c("IBDr", "IBDe"), ]
WW1$Year <- 1
WW10 <- rel10[rel10$Rel %in% c("WW") & rel10$Type %in% c("IBDr", "IBDe"), ]
WW10$Year <- 10
WW1_10 <- rbind(WW1, WW10)

plotWW_1_10 <- ggplot(WW1_10,
                 aes(x = Value, fill = Type)) + geom_histogram(binwidth = 0.01, position = "identity") +
  facet_grid(cols = vars(Year), scales = "free") + xlim(c(-0.01, 2.01))

plotWW_1_10

#Plot csd Year 10
# workers vs workers
tmp <- ibs10_csd[id10$workers, id10$workers]
IBS_WW10_csd <- c(tmp[lower.tri(tmp, diag = TRUE)])
rel10_csd <- data.frame(Value = IBS_WW10_csd, Rel = "WW", Type = "IBS")

tmp <- ibd10_csd[id10$workers, id10$workers]
IBDr_WW10_csd <- c(tmp[lower.tri(tmp, diag = TRUE)])
rel10_csd <- rbind(rel10_csd, data.frame(Value = IBDr_WW10_csd, Rel = "WW", Type = "IBDr"))

tmp <- as.matrix(IBDe[id10$workers, id10$workers])
IBDe_WW10_csd <- c(tmp[lower.tri(tmp, diag = TRUE)])
rel10_csd <- rbind(rel10_csd, data.frame(Value = IBDe_WW10_csd, Rel = "WW", Type = "IBDe"))

# workers vs drones
IBS_WD10_csd<- c(ibs10_csd[id10$workers, id10$drones])
rel10_csd <- rbind(rel10_csd, data.frame(Value = IBS_WW10_csd, Rel = "WD", Type = "IBS"))

IBDr_WD10_csd <- c(ibd10_csd[id10$workers, id10$drones])
rel10_csd <- rbind(rel10_csd, data.frame(Value = IBDr_WW10_csd, Rel = "WD", Type = "IBDr"))

IBDe_WD10_csd <- c(as.matrix(IBDe[id10$workers, id10$drones]))
rel10_csd <- rbind(rel10_csd, data.frame(Value = IBDe_WW10_csd, Rel = "WD", Type = "IBDe"))

#write c("A", "IBD") if you want both on plot
plot10_csd <- ggplot(rel10_csd[rel10_csd$Rel %in% c("WD", "WW") & rel10_csd$Type %in% c("IBDr", "IBDe"), ],
                 aes(x = Value, fill = Type)) + geom_histogram(binwidth = 0.01) +
  facet_grid(cols = vars(Rel)) + xlim(c(-0.01, 2.01))

plot10_csd

#Plot chr3 Year 10
# workers vs workers
tmp <- ibs10_chr3[id10$workers, id10$workers]
IBS_WW10_chr3 <- c(tmp[lower.tri(tmp, diag = TRUE)])
rel10_chr3 <- data.frame(Value = IBS_WW10_chr3, Rel = "WW", Type = "IBS")

tmp <- ibd10_chr3[id10$workers, id10$workers]
IBDr_WW10_chr3 <- c(tmp[lower.tri(tmp, diag = TRUE)])
rel10_chr3 <- rbind(rel10_chr3, data.frame(Value = IBDr_WW10_chr3, Rel = "WW", Type = "IBDr"))

tmp <- as.matrix(IBDe[id10$workers, id10$workers])
IBDe_WW10_chr3 <- c(tmp[lower.tri(tmp, diag = TRUE)])
rel10_chr3 <- rbind(rel10_chr3, data.frame(Value = IBDe_WW10_chr3, Rel = "WW", Type = "IBDe"))

# workers vs drones
IBS_WD10_chr3<- c(ibs10_chr3[id10$workers, id10$drones])
rel10_chr3 <- rbind(rel10_chr3, data.frame(Value = IBS_WW10_chr3, Rel = "WD", Type = "IBS"))

IBDr_WD10_chr3 <- c(ibd10_chr3[id10$workers, id10$drones])
rel10_chr3 <- rbind(rel10_chr3, data.frame(Value = IBDr_WW10_chr3, Rel = "WD", Type = "IBDr"))

IBDe_WD10_chr3 <- c(as.matrix(IBDe[id10$workers, id10$drones]))
rel10_chr3 <- rbind(rel10_chr3, data.frame(Value = IBDe_WW10_chr3, Rel = "WD", Type = "IBDe"))

#write c("IBDr", "IBDe") if you want both on plot
plot10_chr3 <- ggplot(rel10_chr3[rel10_chr3$Rel %in% c("WD", "WW") & rel10_chr3$Type %in% c("IBDr", "IBDe"), ],
                     aes(x = Value, fill = Type)) + geom_histogram(binwidth = 0.01) +
  facet_grid(cols = vars(Rel)) + xlim(c(-0.01, 2.01))

plot10_chr3

#Between queens in colonies


