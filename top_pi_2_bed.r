
rm(list=ls())

library(dplyr)
library(ggplot2)



pi<-read.table("/home/audett/scratch/SSD/Analysis/sexesMerged/pi/all-pi.txt")


headers<-c("chr","pos","C1","C2","E1","E2","L1","L2","S1","S2")

colnames(pi) <- headers

C1<- pi[,c(1,2,3)]
C2<- pi[,c(1,2,4)]
E1<- pi[,c(1,2,5)]
E2<- pi[,c(1,2,6)]
L1<- pi[,c(1,2,7)]
L2<- pi[,c(1,2,8)]
S1<- pi[,c(1,2,9)]
S2<- pi[,c(1,2,10)]

C1$chromEnd <- C1$pos + 4999
C2$chromEnd <- C2$pos + 4999
E1$chromEnd <- E1$pos + 4999
E2$chromEnd <- E2$pos + 4999
L1$chromEnd <- L1$pos + 4999
L2$chromEnd <- L2$pos + 4999
S1$chromEnd <- S1$pos + 4999
S2$chromEnd <- S2$pos + 4999

C1 <- na.omit(C1)
C2 <- na.omit(C2)
E1 <- na.omit(E1)
E2 <- na.omit(E2)
L1 <- na.omit(L1)
L2 <- na.omit(L2)
S1 <- na.omit(S1)
S2 <- na.omit(S2)

C1 <- C1[C1$C1>quantile(C1$C1, 0.95),]
C2 <- C2[C2$C2>quantile(C2$C2, 0.95),]
E1 <- E1[E1$E1>quantile(E1$E1, 0.95),]
E2 <- E2[E2$E2>quantile(E2$E2, 0.95),]
L1 <- L1[L1$L1>quantile(L1$L1, 0.95),]
L2 <- L2[L2$L2>quantile(L2$L2, 0.95),]
S1 <- S1[S1$S1>quantile(S1$S1, 0.95),]
S2 <- S2[S2$S2>quantile(S2$S2, 0.95),]

C1<- C1[,c(1,2,4)]
C2<- C2[,c(1,2,4)]
E1<- E1[,c(1,2,4)]
E2<- E2[,c(1,2,4)]
L1<- L1[,c(1,2,4)]
L2<- L2[,c(1,2,4)]
S1<- S1[,c(1,2,4)]
S2<- S2[,c(1,2,4)]

headers<- c("chrom","chromStart","chromEnd")

colnames(C1) <- headers
colnames(C2) <- headers
colnames(E1) <- headers
colnames(E2) <- headers
colnames(L1) <- headers
colnames(L2) <- headers
colnames(S1) <- headers
colnames(S2) <- headers

C1$chrom <- sub("^", "chr", C1$chrom )
C2$chrom <- sub("^", "chr", C2$chrom )
E1$chrom <- sub("^", "chr", E1$chrom )
E2$chrom <- sub("^", "chr", E2$chrom )
L1$chrom <- sub("^", "chr", L1$chrom )
L2$chrom <- sub("^", "chr", L2$chrom )
S1$chrom <- sub("^", "chr", S1$chrom )
S2$chrom <- sub("^", "chr", S2$chrom )


write.table(C1, file = "/home/audett/scratch/SSD/Analysis/sexesMerged/pi/C1_top5percent_pi.bed", sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(C2, file = "/home/audett/scratch/SSD/Analysis/sexesMerged/pi/C2_top5percent_pi.bed", sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(E1, file = "/home/audett/scratch/SSD/Analysis/sexesMerged/pi/E1_top5percent_pi.bed", sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(E2, file = "/home/audett/scratch/SSD/Analysis/sexesMerged/pi/E2_top5percent_pi.bed", sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(L1, file = "/home/audett/scratch/SSD/Analysis/sexesMerged/pi/L1_top5percent_pi.bed", sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(L2, file = "/home/audett/scratch/SSD/Analysis/sexesMerged/pi/L2_top5percent_pi.bed", sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(S1, file = "/home/audett/scratch/SSD/Analysis/sexesMerged/pi/S1_top5percent_pi.bed", sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(S2, file = "/home/audett/scratch/SSD/Analysis/sexesMerged/pi/S2_top5percent_pi.bed", sep = "\t",
            row.names = FALSE, quote = FALSE)
