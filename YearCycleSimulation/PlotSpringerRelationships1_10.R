setwd("~/github/lstrachan_honeybee_sim/YearCycleSimulation/")
# Plot the relationships from the springer simulation
data <- load("SpringerSimulation.Rdata")

ped <- SP$pedigree
caste <- SP$caste
# Colony in year 1
ibs1 <- ibs_springerColony1
ibd1 <- ibd_springerColony1
id1 <- springerColony1_id

# Colony in year 10
ibs10 <- ibs_springerColony10
ibd10 <- ibd_springerColony1
id10 <- springerColony10_id
# S matrix
S <- S


# Plot year 1
# Queen vs workers
rel1 <- data.frame(Value = numeric(), Rel = character(), Type = character())

ibs_QW1 <- ibs1[id1$queen, id1$workers]
length(ibs_QW1)
rel1 <- rbind(rel1, data.frame(Value = ibs_QW1, Rel = "QW", Type = "IBS"))

ibd_QW1 <- ibd1[id1$queen, id1$workers]
length(ibd_QW1)
rel1 <- rbind(rel1, data.frame(Value = ibd_QW1, Rel = "QW", Type = "IBD"))

A_QW1 <- S[as.numeric(id1$queen), as.numeric(id1$workers)]
length(A_QW1)
rel1 <- rbind(rel1, data.frame(Value = A_QW1, Rel = "QW", Type = "A"))


# Queen vs drones
ibs_QD1 <- ibs1[id1$queen, id1$drones]
rel1 <- rbind(rel1, data.frame(Value = ibs_QD1, Rel = "QD", Type = "IBS"))

ibd_QD1 <- ibd1[id1$queen, id1$drones]
rel1 <- rbind(rel1, data.frame(Value = ibd_QD1, Rel = "QD", Type = "IBD"))

A_QD1 <- S[as.numeric(id1$queen), as.numeric(id1$drones)]
rel1 <- rbind(rel1, data.frame(Value = A_QD1, Rel = "QD", Type = "A"))

# workers vs workers
ibs_WW1 <- ibs1[id1$queen, id1$drones]
rel1 <- rbind(rel1, data.frame(Value = ibs_QD1, Rel = "WW", Type = "IBS"))

ibd_WW1 <- ibd1[id1$queen, id1$drones]
rel1 <- rbind(rel1, data.frame(Value = ibd_QD1, Rel = "WW", Type = "IBD"))

A_WW <- S[as.numeric(id1$queen), as.numeric(id1$drones)]
rel1 <- rbind(rel1, data.frame(Value = A_QD1, Rel = "WW", Type = "A"))

# workers vs drones
ibs_WD1 <- ibs1[id1$queen, id1$drones]
rel1 <- rbind(rel1, data.frame(Value = ibs_QD1, Rel = "WD", Type = "IBS"))

ibd_WD1 <- ibd1[id1$queen, id1$drones]
rel1 <- rbind(rel1, data.frame(Value = ibd_QD1, Rel = "WD", Type = "IBD"))

A_WD1 <- S[as.numeric(id1$queen), as.numeric(id1$drones)]
rel1 <- rbind(rel1, data.frame(Value = A_QD1, Rel = "WD", Type = "A"))

ggplot(rel1, aes(x = Value, fill = Type)) + geom_histogram(bins = 100) +
  facet_grid(rows = vars(Type),
             cols = )



# Plot year 10
# Queen vs workers
rel10 <- data.frame(Value = numeric(), Rel = character(), Type = character())

ibs_QW10 <- ibs10[id1$queen, id1$workers]
length(ibs_QW1)
rel10 <- rbind(rel10, data.frame(Value = ibs_QW1, Rel = "QW", Type = "IBS"))

ibd_QW10 <- ibd1[id1$queen, id1$workers]
length(ibd_QW10)
rel10 <- rbind(rel10, data.frame(Value = ibd_QW10, Rel = "QW", Type = "IBD"))

A_QW10 <- S[as.numeric(id1$queen), as.numeric(id1$workers)]
length(A_QW10)
rel10 <- rbind(rel1, data.frame(Value = A_QW10, Rel = "QW", Type = "A"))


# Queen vs drones
ibs_QD10 <- ibs10[id1$queen, id1$drones]
rel10 <- rbind(rel10, data.frame(Value = ibs_QD10, Rel = "QD", Type = "IBS"))

ibd_QD10 <- ibd1[id1$queen, id1$drones]
rel10 <- rbind(rel10, data.frame(Value = ibd_QD10, Rel = "QD", Type = "IBD"))

A_QD10 <- S[as.numeric(id1$queen), as.numeric(id1$drones)]
rel10 <- rbind(rel10, data.frame(Value = A_QD10, Rel = "QD", Type = "A"))

# workers vs workers
ibs_WW10 <- ibs10[id1$queen, id1$drones]
rel10 <- rbind(rel10, data.frame(Value = ibs_WW10, Rel = "WW", Type = "IBS"))

ibd_WW10 <- ibd10[id1$queen, id1$drones]
rel10 <- rbind(rel10, data.frame(Value = ibd_WW10, Rel = "WW", Type = "IBD"))

A_WW10<- S[as.numeric(id1$queen), as.numeric(id1$drones)]
rel10 <- rbind(rel1, data.frame(Value = A_WW10, Rel = "WW", Type = "A"))

# workers vs drones
ibs_WD10 <- ibs10[id1$queen, id1$drones]
rel10 <- rbind(rel10, data.frame(Value = ibs_WD10, Rel = "WD", Type = "IBS"))

ibd_WD10 <- ibd10[id1$queen, id1$drones]
rel10 <- rbind(rel10, data.frame(Value = ibd_WD1, Rel = "WD", Type = "IBD"))

A_WD10 <- S[as.numeric(id1$queen), as.numeric(id1$drones)]
rel10 <- rbind(rel10, data.frame(Value = A_WD10, Rel = "WD", Type = "A"))

ggplot(rel1, aes(x = Value, fill = Type)) + geom_histogram(bins = 100) +
  facet_grid(rows = vars(Type),
             cols = )

