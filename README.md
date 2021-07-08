# scripts

## run fastqc
QC.sh

## Trim
bbtrim.sh

## QC
QC_after_trim.sh

## Index genome for bwa

bwa index dmel-all-chromosome-r6.23.fasta.gz

## Map genomes

bwa_map.sh

## Sam files converted to bam

sam2bam.sh

## merging first reads with supplimentary reads and then merging sexes
### This requires `mv` to put matching samples in run1 and run2 respectively. The names also must match.
### So to merge E1F and E1M `mv E1F.bam run1/E1.bam` and `mv E1M.bam run2/E1.bam` and then run the merge script.

merge.sh

## Filter for quality over 20

Q20.sh

## Add read-groups because they are necessary for future programs

addrg.sh

## Create an mpileup
### Setting max coverage to 450 to avoid the super high coverage areas causing issues with memory and time
### Maximum coverage of 450 because average coverage should be 400 (males + females at 200 each)
### also swtching back to using info because sharcnet is too time consuming

````
samtools mpileup -Q 20 -q 20 -d 450 \
-f /2/scratch/TylerA/Dmelgenome/gatk/dmel-all-chromosome-r6.23.fa \
./*.bam \
-o Sexes_combined.mpileup
````

## Remove repeat regions

````
perl /home/tylera/bin/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --gtf /2/scratch/TylerA/Dmelgenome/dmel-all-chromosome-r6.23.fasta.out.gff --input ./Sexes_combined.mpileup --output ./Sexes_combined_norepeats.mpileup

/usr/local/RepeatMasker/RepeatMasker -pa 10 -species drosophila -gff dmel-all-chromosome-r6.23.fasta


# Remove suspicious areas

perl /home/tylera/bin/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --gtf /2/scratch/TylerA/SSD/suspicious_coverage.gff --input ./Sexes_combined_norepeats.mpileup --output ./Sexes_combined_norepeats_nosus.mpileup

# ID indels
perl /home/tylera/bin/popoolation2_1201/indel_filtering/identify-indel-regions.pl --input ./Sexes_combined_norepeats_nosus.mpileup --output ./Sexes_combined_norepeat_nosus.gtf --indel-window 10

#Pileup entries processed: 112150681
#Pileup entries containing at least one indel: 1392619
#How many bp of the reference are covered by indel-regions: 24178268


# hard masks indels

perl /home/tylera/bin/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --gtf ./Sexes_combined_norepeat_nosus.gtf --input ./Sexes_combined_norepeats_nosus.mpileup --output ./Sexes_combined_norepeat_nosus_noindel.mpileup

````

## Make a VCF

make_vcf.sh

for sexes combined:
````
java -Xmx32g -jar \
~/bin/VarScan.v2.3.9.jar \
mpileup2indel \
./Sexes_combined_norepeats_nosus.mpileup \
--min-reads2 5 \
--min-coverage 50 \
--p-value 0.1 \
--min-var-freq 0.01 \
--min-freq-for-hom 1 \
--min-avg-qual 20 \
--variants \
--output-vcf 1 \
> ./combined_indels.vcf

java -Xmx32g -jar \
~/bin/VarScan.v2.3.9.jar \
mpileup2snp \
./Sexes_combined_norepeat_nosus_noindel.mpileup \
--min-reads2 5 \
--min-coverage 50 \
--p-value 0.1 \
--min-var-freq 0.01 \
--min-freq-for-hom 1 \
--min-avg-qual 20 \
--variants \
--output-vcf 1 \
> ./combined_variants.vcf

````
for replicates combined min coverage is doubled because coverage should double:

## Make a sync

````
java -ea -jar /usr/local/popoolation/mpileup2sync.jar --threads 16 --input ./Sexes_combined_norepeat_nosus_noindel.mpileup --output ./Sexes_combined_norepeats.sync
````

# Calculate Fst

````
rm(list=ls())

library(poolfstat)
library(WriteXLS)
library(ggplot2)
##### Convert sync file to poolfstat file, and call SNPS

# We first have to give haploid sizes of each pool.
psizes <- as.numeric(c('200','200','200','200','200','200','200','200'))

# Then we give the names of each pool/sample.
pnames <- as.character(c('C1','C2','E1','E2','L1','L2','S1','S2'))

# Here is where we read the sync file and call SNPs. The input file must have the '.sync' extension, and can also be gzipped, like in the example below. The parameters to note are: 

#1) min.rc =  the minimum # reads that an allele needs to have (across all pools) to be called 
#2) min.cov.per.pool = the minimum allowed read count per pool for SNP to be called
#3) max.cov.per.pool = the maximum read count per pool for SNP to be called 
#4) min.maf = the minimum allele frequency (over all pools) for a SNP to be called (note this is obtained from dividing the read counts for the minor allele over the total read coverage) 
#5) nlines.per.readblock = number of lines in sync file to be read simultaneously 

SG.pooldata <- vcf2pooldata(vcf.file = "/2/scratch/TylerA/SSD/bwamap/Sexes_combined_variants.vcf", poolsizes = psizes, poolnames = pnames)

##### And we can compute pairwise FSTs
SG.pair.fst <- computePairwiseFSTmatrix(SG.pooldata, method = "Anova",
                                        output.snp.values = TRUE)

test<-as.matrix(SG.pair.fst$PairwiseSnpFST)
crap <- data.frame(SG.pooldata@snp.info, test[,c(2,3,8,9,24,25,26,27)])


CVE <- data.frame(crap[c(1,2,5,6,7,8)])
LVS <- data.frame(crap[c(1,2,9,10,11,12)])



CVE<-data.frame(ID=CVE[,c(1:2)], Means=rowMeans(CVE[,-c(1:2)], na.rm=TRUE))
LVS<-data.frame(ID=LVS[,c(1:2)], Means=rowMeans(LVS[,-c(1:2)], na.rm=TRUE))

CVE <- CVE[CVE$Means!='NaN',]
LVS <- LVS[LVS$Means!='NaN',]

write.table(CVE, file = "/2/scratch/TylerA/SSD/bwamap/CVE.fst", sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(LVS, file = "/2/scratch/TylerA/SSD/bwamap/LVS.fst", sep = "\t",
            row.names = FALSE, quote = FALSE)
````

## Plot Fst

plot_fst.sh (calls plot_fst.R)



# Checking vcf coverage

````
vcftools --vcf Sexes_combined_variants.vcf --site-mean-depth
````
change name so my plot script can work

````
mv out.ldepth.mean Sexes_combined.coverage
````

````
plot_coverage.sh
````
## Calculate CMH

````
rm(list=ls())
#install.packages("/home/tylera/bin/poolSeq-0.3.5.tar.gz", repos=NULL, type="source")
#install.packages("/home/tylera/bin/ACER-1.0.2.tar.gz")

#Loading in the required packages
library(ggplot2)
library(poolSeq)
library(ACER)

reps <- c(1:8)
gen <- rep(0,8)
sync <- read.sync("/2/scratch/TylerA/SSD/bwamap/Sexes_combined_norepeats.sync", 
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

# Ne = 237

#Creating the vars for the CMH test 
rep<-c(1,1,2,2) #Number of replicates
Ne<-c(237,237)
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


write.csv(data, "/2/scratch/TylerA/SSD/bwamap/CVE_pval.csv")
write.csv(data2, "/2/scratch/TylerA/SSD/bwamap/LVS_pval.csv") 

ddat2<-data


ddat22L <- ddat2[which(ddat2$chr=='2L'),]
ddat22R <- ddat2[which(ddat2$chr=='2R'),]
ddat23L <- ddat2[which(ddat2$chr=='3L'),]
ddat23R <- ddat2[which(ddat2$chr=='3R'),]
ddat24 <- ddat2[which(ddat2$chr=='4'),]
ddat2X <- ddat2[which(ddat2$chr=='X'),]
ddat2 <- rbind(ddat2X, ddat22L, ddat22R, ddat23L, ddat23R, ddat24)


g <- nrow(ddat2[which(ddat2$chr=='2L'),])
h <- nrow(ddat2[which(ddat2$chr=='2R'),])
i <- nrow(ddat2[which(ddat2$chr=='3L'),])
j <- nrow(ddat2[which(ddat2$chr=='3R'),])
k <- nrow(ddat2[which(ddat2$chr=='4'),])
l <- nrow(ddat2[which(ddat2$chr=='X'),])

ddat2$number <-  c((1:l),
                   (l+1):(l+g), 
                   (l+g+1):(l+g+h), 
                   (l+g+h+1):(l+g+h+i),
                   (l+g+h+i+1):(l+g+h+i+j),
                   (l+g+h+i+j+1):(l+g+h+i+j+k))


plot<-ggplot(ddat2, aes(x=number, y=neg.log10, color=chr)) +
  geom_point(size=0.5, show.legend = F, alpha=0.25) +
  theme(panel.background = element_blank()) +
  scale_colour_manual(values=c("seagreen", "darkslateblue", 'darkred', 'darkorchid4', 'darkolivegreen', 'darkblue')) +
  theme(text = element_text(size=20),
        axis.text.x= element_text(size=15), 
        axis.text.y= element_text(size=15))
        
        
png("/2/scratch/TylerA/SSD/bwamap/CVE_pval_plot.png",type="cairo")
plot
dev.off()


ddat2<-data2


ddat22L <- ddat2[which(ddat2$chr=='2L'),]
ddat22R <- ddat2[which(ddat2$chr=='2R'),]
ddat23L <- ddat2[which(ddat2$chr=='3L'),]
ddat23R <- ddat2[which(ddat2$chr=='3R'),]
ddat24 <- ddat2[which(ddat2$chr=='4'),]
ddat2X <- ddat2[which(ddat2$chr=='X'),]
ddat2 <- rbind(ddat2X, ddat22L, ddat22R, ddat23L, ddat23R, ddat24)


g <- nrow(ddat2[which(ddat2$chr=='2L'),])
h <- nrow(ddat2[which(ddat2$chr=='2R'),])
i <- nrow(ddat2[which(ddat2$chr=='3L'),])
j <- nrow(ddat2[which(ddat2$chr=='3R'),])
k <- nrow(ddat2[which(ddat2$chr=='4'),])
l <- nrow(ddat2[which(ddat2$chr=='X'),])

ddat2$number <-  c((1:l),
                   (l+1):(l+g), 
                   (l+g+1):(l+g+h), 
                   (l+g+h+1):(l+g+h+i),
                   (l+g+h+i+1):(l+g+h+i+j),
                   (l+g+h+i+j+1):(l+g+h+i+j+k))


plot<-ggplot(ddat2, aes(x=number, y=neg.log10, color=chr)) +
  geom_point(size=0.5, show.legend = F, alpha=0.25) +
  theme(panel.background = element_blank()) +
  scale_colour_manual(values=c("seagreen", "darkslateblue", 'darkred', 'darkorchid4', 'darkolivegreen', 'darkblue')) +
  theme(text = element_text(size=20),
        axis.text.x= element_text(size=15), 
        axis.text.y= element_text(size=15))
        
        
png("/2/scratch/TylerA/SSD/bwamap/LVS_pval_plot.png",type="cairo")
plot
dev.off()




````
## Plot CMH -log10
# Control vs. SSD-reverse

CVE_pval_plot.png![CVE_pval_plot](https://user-images.githubusercontent.com/77504755/124979733-6d672f00-e001-11eb-9522-9cdc01f1f325.png)

# Large vs. Small

LVS_pval_plot.png![LVS_pval_plot](https://user-images.githubusercontent.com/77504755/124979750-75bf6a00-e001-11eb-8228-10e223540b6a.png)

## Creating a suspicious coverage .gff

The high Fst found in male vs. female comparisons in the conflict_workflow highlights some areas that are most likely not true. So I created a file with those locations with a 10bp before and after window to remove from my genomes.

````
################################# USING POOLFSTAT ###############################
##################################################################################

rm(list=ls())

library(poolfstat)
library(WriteXLS)
library(ggplot2)
library(plyr)
library(dplyr)
##### Convert sync file to poolfstat file, and call SNPS

# We first have to give haploid sizes of each pool.
psizes <- as.numeric(c('100','100','100','100','100','100','100','100','100','100','100','100','100','100','100','100'))

# Then we give the names of each pool/sample.
pnames <- as.character(c('C1F','C1M','C2F','C2M','E1F','E1M','E2F','E2M','L1F','L1M','L2F','L2M','S1F','S1M','S2F','S2M'))

# Here is where we read the sync file and call SNPs. The input file must have the '.sync' extension, and can also be gzipped, like in the example below. The parameters to note are: 

#1) min.rc =  the minimum # reads that an allele needs to have (across all pools) to be called 
#2) min.cov.per.pool = the minimum allowed read count per pool for SNP to be called
#3) max.cov.per.pool = the maximum read count per pool for SNP to be called 
#4) min.maf = the minimum allele frequency (over all pools) for a SNP to be called (note this is obtained from dividing the read counts for the minor allele over the total read coverage) 
#5) nlines.per.readblock = number of lines in sync file to be read simultaneously 

SG.pooldata <- vcf2pooldata(vcf.file = "/2/scratch/TylerA/SSD/Conflict/conflict_norepeat_variants.vcf", poolsizes = psizes, poolnames = pnames,
                                     min.rc = 5, min.cov.per.pool = 50, max.cov.per.pool = 300,
                                     min.maf = 0.001, remove.indels = FALSE, nlines.per.readblock = 1e+06)

#Parsing allele counts
#VarScan like format detected for allele count data: the AD field contains allele depth for the alternate allele and RD field for the reference allele (N.B., positions with more than one alternate allele will be ignored)
#1e+06  lines processed in 0 h  0 m  53 s : 821946 SNPs found
#1429860  lines processed in 0 h  1 m  13 s : 1085681 SNPs found
#Data consists of 1085681 SNPs for 16 Pools



##### From this file we can compute global and per SNP FSTs
SG.snp.fsts <- computeFST(SG.pooldata, method = "Anova", snp.index = NA)

#Assign SNP location for graphing
crap <- data.frame(SG.pooldata@snp.info, SG.snp.fsts$snp.FST)




##### And we can compute pairwise FSTs
SG.pair.fst <- computePairwiseFSTmatrix(SG.pooldata, method = "Anova",
                                        min.cov.per.pool = 50, max.cov.per.pool = 175,
                                        min.maf = 0.001,
                                        output.snp.values = TRUE)

##### If you want to save the pairwise matrix, you can do the following:
#SG.p.fst <- SG.pair.fst$PairwiseFSTmatrix
#SG.p.fst <- as.data.frame(SG.p.fst)
#WriteXLS(SG.p.fst, "SG.p.fst.xls")

test<-as.matrix(SG.pair.fst$PairwiseSnpFST)
crap <- data.frame(SG.pooldata@snp.info, test[,c(1,30,55,76,93,106,115,120)])


C1F_C1M <- data.frame(crap[c(1,2,5)])
C2F_C2M <- data.frame(crap[c(1,2,6)])
E1F_E1M <- data.frame(crap[c(1,2,7)])
E2F_E2M <- data.frame(crap[c(1,2,8)])
L1F_L1M <- data.frame(crap[c(1,2,9)])
L2F_L2M <- data.frame(crap[c(1,2,10)])
S1F_S1M <- data.frame(crap[c(1,2,11)])
S2F_S2M <- data.frame(crap[c(1,2,12)])

C1F_C1M <- na.omit(C1F_C1M)
C2F_C2M <- na.omit(C2F_C2M)
E1F_E1M <- na.omit(E1F_E1M)
E2F_E2M <- na.omit(E2F_E2M)
L1F_L1M <- na.omit(L1F_L1M)
L2F_L2M <- na.omit(L2F_L2M)
S1F_S1M <- na.omit(S1F_S1M)
S2F_S2M <- na.omit(S2F_S2M)

high_C1F_C1M <- C1F_C1M[which(C1F_C1M$C1F_vs_C1M>=0.05),]
high_C2F_C2M <- C2F_C2M[which(C2F_C2M$C2F_vs_C2M>=0.05),]
high_E1F_E1M <- E1F_E1M[which(E1F_E1M$E1F_vs_E1M>=0.05),]
high_E2F_E2M <- E2F_E2M[which(E2F_E2M$E2F_vs_E2M>=0.05),]
high_L1F_L1M <- L1F_L1M[which(L1F_L1M$L1F_vs_L1M>=0.05),]
high_L2F_L2M <- L2F_L2M[which(L2F_L2M$L2F_vs_L2M>=0.05),]
high_S1F_S1M <- S1F_S1M[which(S1F_S1M$S1F_vs_S1M>=0.05),]
high_S2F_S2M <- S2F_S2M[which(S2F_S2M$S2F_vs_S2M>=0.05),]

match<-match_df(high_C1F_C1M, high_C2F_C2M, on = "X2")

test<-match[,1:2]
test$X2<-as.numeric(test$X2)
test$start<-test$X2-10
test$end<-test$X2+10
test$feature<-rep("repeat",length(567))
test$score<-rep("1",length(567))
test$strand<-rep("1",length(567))
test <- data.frame(test[,c(2,1,5,3,4,6,7)])

write.table(test, file = "/2/scratch/TylerA/SSD/suspicious_coverage.gff", sep = "\t",
            row.names = FALSE, quote = FALSE)
````





