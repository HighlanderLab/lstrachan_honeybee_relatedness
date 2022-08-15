ped <- read.table("~/Documents/1Projects/SpringerChapter_HoneybeeRelatedness/pedigree.txt")
head(ped)
colnames(ped) <- c("ID", "Sex", "Mate", "MatingType", "YOBModified", "DamSeq", "SireSeq", 
                           "ass", "F", "varMS", "testStation", "diagA", "thirteen", "fourteen", 
                           "fifteen", "seq1", "seq2", "YOB", "NS", "ND")
table(ped$Sex)

pedComplete <- read.table("~/Documents/1Projects/SpringerChapter_HoneybeeRelatedness/pedigree_complete.txt")
colnames(pedComplete) <- c("ID", "Sex", "Mate", "MatingType", "YOBModified", "DamSeq", "SireSeq", 
                   "ass", "F", "varMS", "testStation", "diagA", "thirteen", "fourteen", 
                   "fifteen", "seq1", "seq2", "YOB", "NS", "ND")

table(pedComplete$ID)
table(pedComplete$Sex) #1 = dam, 2 = sire, 3 = colony; 11211 matches the number of input queens
table(pedComplete$MatingType)
hist(pedComplete$diagA)
# Which one is the original ID?
length(intersect(inputPed$Queen, pedComplete$ID)) #ID is ID of the queen! The fabricated IDs are for the sire!

ident <- read.table("~/Documents/1Projects/SpringerChapter_HoneybeeRelatedness/ident.txt")

Ainv <- read.table("~/Documents/1Projects/SpringerChapter_HoneybeeRelatedness/AINV.giv")

ggplot(Ainv, aes(x = V1, y = V2, fill = V3)) + geom_tile()
