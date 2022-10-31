rm(list=ls())

#install.packages("data.table")
#install.packages("/home/audett/projects/def-idworkin/audett/SSD/scripts/ACER-ACER_1.0.3.tar.gz", repos=NULL, type="source")
#install.packages("matrixStats")
#install.packages("/home/audett/projects/def-idworkin/audett/SSD/scripts/poolSeq-0.3.5.tar.gz", repos=NULL, type="source")

#Loading in the required packages
library(ggplot2)
library(poolSeq)
library(ACER)

reps <- c(1:8)
gen <- rep(0,8)
sync <- read.sync("/home/audett/scratch/SSD/Analysis/sexesMerged/sexesMerged.sync", 
                  gen=gen, repl=reps, 
                  polarization = "minor", 
                  keepOnlyBiallelic = TRUE)


pops <-c('C1','C2','E1','E2','L1','L2','S1','S2')
       

af.mat <- matrix(NA, nrow = nrow(sync@alleles), ncol = 8)
colnames(af.mat) <- pops

for (i in 1:ncol(af.mat)){
  tempdat <- af(sync, repl = i, gen = 0)
  af.mat[,i] <- as.matrix(tempdat)
}

af.mat <- na.omit(af.mat)
head(af.mat)
dim(af.mat)

af.mat2 <- af.mat[,c(1,2,3,4)]
af.mat3 <- af.mat[,c(5,6,7,8)]
dim(af.mat3)

    
#now to make a coverage one. 

cov.mat <- matrix(NA, nrow = nrow(sync@alleles), ncol = 8)
colnames(cov.mat) <- pops

for (i in 1:ncol(cov.mat)){
  tempdat <- coverage(sync, repl = i, gen = 0)
  cov.mat[,i] <- as.matrix(tempdat)
}

crap <- data.frame(cov.mat, sync@alleles[,1:2])
crap[crap==0] <- NA
crap2 <- na.omit(crap)
location <- crap2[,9:10]

cov.mat[cov.mat == 0] <- NA
cov.mat <- na.omit(cov.mat)

dim(cov.mat)

head(cov.mat)

cov.mat2 <- cov.mat[,c(1,2,3,4)]
cov.mat3 <- cov.mat[,c(5,6,7,8)]
dim(cov.mat2)
dim(cov.mat3)

#Now I want to estimate Ne to use below. 
ne <- estimateNe(p0=af.mat2[,"C1"], pt=af.mat2[,"E1"], 
           cov0=cov.mat2[,"C1"], covt=cov.mat2[,"E1"], 
           t=750, method = "P.planII", poolSize=c(200, 200))

# Ne = 184

#Creating the vars for the CMH test 
rep<-c(1,1,2,2) #Number of replicates
Ne<-c(184,184)
tp<-c(0,0,750,750) #Generations of evolution for each sample
ps<-c(400,400,400,400) #Pool size

#cmh test 
#I think I need to include the random pools as a gen 0?

pval <- adapted.cmh.test(freq=af.mat2, coverage=cov.mat2, 
                         Ne=Ne, gen=tp, repl=rep, poolSize=ps)
pval2 <- adapted.cmh.test(freq=af.mat3, coverage=cov.mat3, 
                         Ne=Ne, gen=tp, repl=rep, poolSize=ps)

# Warning messages:
#   1: In adapted.cmh.test(freq = af.mat, coverage = cov.mat, Ne = Ne,  :
#       Ne value(s) which are not integer are converted to integer
#I took care of this (warning 2) by setting left to gen 0 and right to gen 1 
#   2: In adapted.cmh.test(freq = af.mat, coverage = cov.mat, Ne = Ne,  :
#       Value of 'Ne' will be ignored because no random genetic drift is assumed.
#   3: In adapted.cmh.test(freq = af.mat, coverage = cov.mat, Ne = Ne,  :
#       The counts that equal 0 or equal the coverage in all replicates are changed to 1 or to coverage-1 respectively.


padj <- p.adjust(pval, "fdr")
data <- cbind(cov.mat2, pval, padj)
data<-cbind(location,data)
data<-data[,c(1,2,7,8)]
data$neg.log10 <- -log10(data$padj)

padj2 <- p.adjust(pval2, "fdr")
data2 <- cbind(cov.mat3, pval2, padj2)
data2<-cbind(location,data2)
data2<-data2[,c(1,2,7,8)]
data2$neg.log10 <- -log10(data2$padj)






write.csv(data, "/home/audett/scratch/SSD/Analysis/sexesMerged/CMH/CVE_pval.csv")
write.csv(data2, "/home/audett/scratch/SSD/Analysis/sexesMerged/CMH/LVS_pval.csv") 

data <- na.omit(data)
data2 <- na.omit(data2)

data <- data[data$padj<quantile(data$padj,0.01),]
data2 <- data2[data2$padj<quantile(data2$padj,0.01),]

data$chromEnd <- data[,2]
data2$chromEnd <- data2[,2]

data <- data[,c(1,2,6)]
data2 <- data2[,c(1,2,6)]

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

write.table(data, file = "/home/audett/scratch/SSD/Analysis/sexesMerged/CMH/CVE_top1percent.bed", sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(data2, file = "/home/audett/scratch/SSD/Analysis/sexesMerged/CMH/LVS_top1percent.bed", sep = "\t",
            row.names = FALSE, quote = FALSE)
