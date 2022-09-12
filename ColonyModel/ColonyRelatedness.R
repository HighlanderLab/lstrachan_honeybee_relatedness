setwd("~/Documents/1Projects/SpringerChapter_HoneybeeRelatedness/")

# Colony relatedness
data <- load("SpringerSimulation_import_objects.RData")
pedigree <- as.data.frame(pedigree)
caste <- caste
colonyRecords <- colonyRecords # No duplicated IDs
nrow(colonyRecords)
sum(caste == "queen")
length(unique(pedigree$mother))

# Connect mothers and their YOB
mothersYOB <- colonyRecords[, c("Id", "year"), drop = F]
colnames(mothersYOB) <- c("MId", "yearMother")
colonyRecords <- merge(colonyRecords, mothersYOB, by = "MId", all.x=T)[, union(names(colonyRecords), names(mothersYOB))]
colonyRecords$yearMother[is.na(colonyRecords$yearMother)] <- 0
colonyRecords <- colonyRecords[order(colonyRecords$year), ]
inputPed <- data.frame(yearObs = 0,
                       breeder = 0,
                       testStation = 0,
                       YOBQueen = colonyRecords$year,
                       Queen = colonyRecords$Id,
                       YOBDamOfQueen = colonyRecords$yearMother,
                       DamOfQueen = colonyRecords$MId,
                       YOBSire = 0,
                       Sire = 0,
                       ten = 0,
                       eleven = 0,
                       twelve = 0,
                       thirteen = 0,
                       NS = colonyRecords$nDPQ,
                       ND = colonyRecords$nFathers,
                       YOBDPQ = 0,
                       DPQ = 0,
                       DPQid = 0)

head(inputPed)
write.table(inputPed, "~/Documents/1Projects/SpringerChapter_HoneybeeRelatedness/input-pedigree.txt", quote=F, row.names=F, col.names=F)
