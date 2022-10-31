CVE <- read.table("/home/audett/scratch/SSD/Analysis/repsMerged/fst/CVE.fst", header = TRUE)
AVE <- read.table("/home/audett/scratch/SSD/Analysis/repsMerged/fst/AVE.fst", header = TRUE)
LVE <- read.table("/home/audett/scratch/SSD/Analysis/repsMerged/fst/LVE.fst", header = TRUE)
LVS <- read.table("/home/audett/scratch/SSD/Analysis/repsMerged/fst/LVS.fst", header = TRUE)


CVE<- CVE[CVE$Means>quantile(CVE$Means, 0.95),]
AVE<- AVE[AVE$Means>quantile(AVE$Means, 0.95),]
LVE<- LVE[LVE$Means>quantile(LVE$Means, 0.95),]
LVS<- LVS[LVS$Means>quantile(LVS$Means, 0.95),]


AVE$chromEnd <- AVE[,2]
CVE$chromEnd <- CVE[,2]
LVE$chromEnd <- LVE[,2]
LVS$chromEnd <- LVS[,2]

AVE <- AVE[,c(1,2,4)]
CVE <- CVE[,c(1,2,4)]
LVE <- LVE[,c(1,2,4)]
LVS <- LVS[,c(1,2,4)]

headers <- c("chrom","chromStart","chromEnd")

colnames(AVE)<-headers
colnames(LVS)<-headers
colnames(CVE)<-headers
colnames(LVE)<-headers



write.table(CVE, file = "/home/audett/scratch/SSD/Analysis/repsMerged/fst/CVE_top5percent.bed", sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(LVS, file = "/home/audett/scratch/SSD/Analysis/repsMerged/fst/LVS_top5percent.bed", sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(AVE, file = "/home/audett/scratch/SSD/Analysis/repsMerged/fst/AVE_top5percent.bed", sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(LVE, file = "/home/audett/scratch/SSD/Analysis/repsMerged/fst/LVE_top5percent.bed", sep = "\t",
            row.names = FALSE, quote = FALSE)

