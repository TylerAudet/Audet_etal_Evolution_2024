
CVE<-read.csv("/home/audett/scratch/SSD/Analysis/sexesMerged/CMH/CVE_pval.csv" ,row.names=FALSE)
LVS<-read.csv("/home/audett/scratch/SSD/Analysis/sexesMerged/CMH/LVS_pval.csv" ,row.names=FALSE) 



data <- na.omit(CVE)
data2 <- na.omit(LVS)

data <- data[data$padj<0.01,]
data2 <- data2[data2$padj<0.01,]

data$chromEnd <- data[,3]
data2$chromEnd <- data2[,3]

data <- data[,c(2,3,7)]
data2 <- data2[,c(2,3,7)]

headers <- c("chrom","chromStart","chromEnd")

colnames(data)<-headers
colnames(data2)<-headers

data$temp <- paste("chr", data$chrom, sep = "")
data2$temp <- paste("chr", data2$chrom, sep = "")

data <- data[,c(4,2,3)]
data2 <- data2[,c(4,2,3)]

colnames(data) <- headers
colnames(data2) <- headers

data$temp <- data$chromStart - 5
data$temp2 <- data$chromEnd + 5
data <- data[c(1,4,5)]

data2$temp <- data2$chromStart - 5
data2$temp2 <- data2$chromEnd + 5
data2 <- data2[c(1,4,5)]

colnames(data) <- headers
colnames(data2) <- headers

write.table(data, file = "/home/audett/scratch/SSD/Analysis/cutoff_10/CVE_signif_CMH.bed", sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(data2, file = "/home/audett/scratch/SSD/Analysis/cutoff_10/LVS_signif_CMH.bed", sep = "\t",
            row.names = FALSE, quote = FALSE)
