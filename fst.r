rm(list=ls())

library(poolfstat)
library(ggplot2)


psizes <- as.numeric(c('800','800','800','800'))

pnames <- as.character(c('C','E','L','S'))

# Here is where we read the vcf file and call SNPs. 
#1) min.rc =  the minimum # reads that an allele needs to have (across all pools) to be called 
#2) min.cov.per.pool = the minimum allowed read count per pool for SNP to be called
#3) max.cov.per.pool = the maximum read count per pool for SNP to be called 
#4) min.maf = the minimum allele frequency (over all pools) for a SNP to be called (note this is obtained from dividing the read counts for the minor allele over the total read coverage) 
#5) nlines.per.readblock = number of lines in sync file to be read simultaneously 

SG.pooldata <- vcf2pooldata(vcf.file = "/home/audett/scratch/SSD/Analysis/repsMerged/repsMerged_blacklisted.vcf", poolsizes = psizes, poolnames = pnames, min.cov.per.pool = 200, min.maf = 0.05)

#Data consists of 617063 SNPs for 8 Pools


##### And we can compute pairwise FSTs
SG.pair.fst <- compute.pairwiseFST(SG.pooldata, method = "Anova",
                                        output.snp.values = TRUE)

# Extracting fst as a matrix and then making a data.frame with the snp info associated with the fst values
test<-as.matrix(SG.pair.fst@PairwiseSnpFST)
crap <- data.frame(SG.pooldata@snp.info, test)



# Extract the valuess associated with the comparisons we care about
#CVE = controls vs. experimental
#AVE = all samples vs. the experimental
#LVS = large vs. small
#LVE = Large vs. Experimental

CVE <- data.frame(crap[c(1,2,5)])
AVE <- data.frame(crap[c(1,2,5,8,9)])
LVS <- data.frame(crap[c(1,2,10)])
LVE <- data.frame(crap[c(1,2,8)])

AVE$mean <- rowMeans(AVE[c(3:5)], na.rm = TRUE)



AVE <- AVE[c(1,2,6)]


AVE <- na.omit(AVE)
LVS <- na.omit(LVS)
CVE <- na.omit(CVE)
LVE <- na.omit(LVE)


# Changing headers for clarity and to match the rolling average code

headers<-c("ID.Chromosome","ID.Position","Means")
colnames(AVE)<-headers
colnames(LVS)<-headers
colnames(CVE)<-headers
colnames(LVE)<-headers

##AVE<-data.frame(ID=AVE[,c(1:2)], Means=rowMeans(AVE[,-c(1:2)], na.rm=TRUE))

# Make data tables to work with without having to re-run fst code

write.table(CVE, file = "/home/audett/scratch/SSD/Analysis/repsMerged/fst/CVE.fst", sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(LVS, file = "/home/audett/scratch/SSD/Analysis/repsMerged/fst/LVS.fst", sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(AVE, file = "/home/audett/scratch/SSD/Analysis/repsMerged/fst/AVE.fst", sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(LVE, file = "/home/audett/scratch/SSD/Analysis/repsMerged/fst/LVE.fst", sep = "\t",
            row.names = FALSE, quote = FALSE)
