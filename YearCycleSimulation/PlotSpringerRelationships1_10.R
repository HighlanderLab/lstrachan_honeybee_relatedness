#setwd("~/github/lstrachan_honeybee_sim/YearCycleSimulation/")
# Plot the relationships from the springer simulation
#data <- load("SpringerSimulation.Rdata")

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

# Plot year 1
# workers vs workers
tmp <- ibs1[id1$workers, id1$workers]
ibs_WW1 <- c(tmp[lower.tri(tmp, diag = TRUE)])
rel1 <- data.frame(Value = ibs_WW1, Rel = "WW", Type = "IBS")

tmp <- ibd1[id1$workers, id1$workers]
ibd_WW1 <- c(tmp[lower.tri(tmp, diag = TRUE)])
rel1 <- rbind(rel1, data.frame(Value = ibd_WW1, Rel = "WW", Type = "IBD"))

tmp <- as.matrix(S[id1$workers, id1$workers])
A_WW1 <- c(tmp[lower.tri(tmp, diag = TRUE)])
rel1 <- rbind(rel1, data.frame(Value = A_WW1, Rel = "WW", Type = "A"))

# workers vs drones
ibs_WD1 <- c(ibs1[id1$workers, id1$drones])
rel1 <- rbind(rel1, data.frame(Value = ibs_WD1, Rel = "WD", Type = "IBS"))

ibd_WD1 <- c(ibd1[id1$workers, id1$drones])
rel1 <- rbind(rel1, data.frame(Value = ibd_WD1, Rel = "WD", Type = "IBD"))

A_WD1 <- c(as.matrix(S[id1$workers, id1$drones]))
rel1 <- rbind(rel1, data.frame(Value = A_WD1, Rel = "WD", Type = "A"))

#write c("A", "IBD") if you want both on plot 
plot1 <- ggplot(rel1[rel1$Rel %in% c("WD", "WW") & rel1$Type %in% "IBD", ],
                aes(x = Value, fill = Type)) + geom_histogram(binwidth = 0.01) +
  facet_grid(cols = vars(Rel)) + xlim(c(-0.01, 2.01))

plot1

# Plot year 10
# workers vs workers
tmp <- ibs10[id10$workers, id10$workers]
ibs_WW10 <- c(tmp[lower.tri(tmp, diag = TRUE)])
rel10 <- data.frame(Value = ibs_WW10, Rel = "WW", Type = "IBS")

tmp <- ibd10[id10$workers, id10$workers]
ibd_WW10 <- c(tmp[lower.tri(tmp, diag = TRUE)])
rel10 <- rbind(rel10, data.frame(Value = ibd_WW10, Rel = "WW", Type = "IBD"))

tmp <- as.matrix(S[id10$workers, id10$workers])
A_WW10 <- c(tmp[lower.tri(tmp, diag = TRUE)])
rel10 <- rbind(rel10, data.frame(Value = A_WW10, Rel = "WW", Type = "A"))

# workers vs drones
ibs_WD10 <- c(ibs10[id10$workers, id10$drones])
rel10 <- rbind(rel10, data.frame(Value = ibs_WD10, Rel = "WD", Type = "IBS"))

ibd_WD10 <- c(ibd10[id10$workers, id10$drones])
rel10 <- rbind(rel10, data.frame(Value = ibd_WD10, Rel = "WD", Type = "IBD"))

A_WD10 <- c(as.matrix(S[id10$workers, id10$drones]))
rel10 <- rbind(rel10, data.frame(Value = A_WD10, Rel = "WD", Type = "A"))

#write c("A", "IBD") if you want both on plot 
plot10 <- ggplot(rel10[rel10$Rel %in% c("WD", "WW") & rel10$Type %in% c("A", "IBD"), ],
                aes(x = Value, fill = Type)) + geom_histogram(binwidth = 0.01) +
  facet_grid(cols = vars(Rel)) + xlim(c(-0.01, 2.01))

plot10

#Plot csd Year 10
# workers vs workers
tmp <- ibs10_csd[id10$workers, id10$workers]
ibs_WW10_csd <- c(tmp[lower.tri(tmp, diag = TRUE)])
rel10_csd <- data.frame(Value = ibs_WW10_csd, Rel = "WW", Type = "IBS")

tmp <- ibd10_csd[id10$workers, id10$workers]
ibd_WW10_csd <- c(tmp[lower.tri(tmp, diag = TRUE)])
rel10_csd <- rbind(rel10_csd, data.frame(Value = ibd_WW10_csd, Rel = "WW", Type = "IBD"))

tmp <- as.matrix(S[id10$workers, id10$workers])
A_WW10_csd <- c(tmp[lower.tri(tmp, diag = TRUE)])
rel10_csd <- rbind(rel10_csd, data.frame(Value = A_WW10_csd, Rel = "WW", Type = "A"))

# workers vs drones
ibs_WD10_csd<- c(ibs10_csd[id10$workers, id10$drones])
rel10_csd <- rbind(rel10_csd, data.frame(Value = ibs_WD10, Rel = "WD", Type = "IBS"))

ibd_WD10_csd <- c(ibd10_csd[id10$workers, id10$drones])
rel10_csd <- rbind(rel10_csd, data.frame(Value = ibd_WD10_csd, Rel = "WD", Type = "IBD"))

A_WD10_csd <- c(as.matrix(S[id10$workers, id10$drones]))
rel10_csd <- rbind(rel10_csd, data.frame(Value = A_WD10, Rel = "WD", Type = "A"))

#write c("A", "IBD") if you want both on plot 
plot10_csd <- ggplot(rel10_csd[rel10_csd$Rel %in% c("WD", "WW") & rel10_csd$Type %in% c("A", "IBD"), ],
                 aes(x = Value, fill = Type)) + geom_histogram(binwidth = 0.01) +
  facet_grid(cols = vars(Rel)) + xlim(c(-0.01, 2.01))

plot10_csd

#Plot chr3 Year 10
# workers vs workers
tmp <- ibs10_chr3[id10$workers, id10$workers]
ibs_WW10_chr3 <- c(tmp[lower.tri(tmp, diag = TRUE)])
rel10_chr3 <- data.frame(Value = ibs_WW10_chr3, Rel = "WW", Type = "IBS")

tmp <- ibd10_chr3[id10$workers, id10$workers]
ibd_WW10_chr3 <- c(tmp[lower.tri(tmp, diag = TRUE)])
rel10_chr3 <- rbind(rel10_chr3, data.frame(Value = ibd_WW10_chr3, Rel = "WW", Type = "IBD"))

tmp <- as.matrix(S[id10$workers, id10$workers])
A_WW10_chr3 <- c(tmp[lower.tri(tmp, diag = TRUE)])
rel10_chr3 <- rbind(rel10_chr3, data.frame(Value = A_WW10_chr3, Rel = "WW", Type = "A"))

# workers vs drones
ibs_WD10_chr3<- c(ibs10_chr3[id10$workers, id10$drones])
rel10_chr3 <- rbind(rel10_chr3, data.frame(Value = ibs_WD10, Rel = "WD", Type = "IBS"))

ibd_WD10_chr3 <- c(ibd10_chr3[id10$workers, id10$drones])
rel10_chr3 <- rbind(rel10_chr3, data.frame(Value = ibd_WD10_chr3, Rel = "WD", Type = "IBD"))

A_WD10_chr3 <- c(as.matrix(S[id10$workers, id10$drones]))
rel10_chr3 <- rbind(rel10_chr3, data.frame(Value = A_WD10, Rel = "WD", Type = "A"))

#write c("A", "IBD") if you want both on plot 
plot10_chr3 <- ggplot(rel10_chr3[rel10_chr3$Rel %in% c("WD", "WW") & rel10_chr3$Type %in% c("A", "IBD"), ],
                     aes(x = Value, fill = Type)) + geom_histogram(binwidth = 0.01) +
  facet_grid(cols = vars(Rel)) + xlim(c(-0.01, 2.01))

plot10_chr3


