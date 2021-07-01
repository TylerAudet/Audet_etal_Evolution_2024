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
### Setting max coverage to 800 to avoid the super high coverage areas causing issues with memory and time
### Maximum coverage of 800 because average coverage should be 400 (males + females at 200 each) and I want to capture possible duplications, but nothing more than that.

mpileup_sexes.sh


## Remove repeat regions

````
perl /home/tylera/bin/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --gtf /2/scratch/TylerA/Dmelgenome/dmel-all-chromosome-r6.23.fasta.out.gff --input ./Sexes_combined.mpileup --output ./Sexes_combined_norepeats.mpileup
````

## Make a VCF

make_vcf.sh

# Calculate Fst

````
rm(list=ls())

library(poolfstat)
library(WriteXLS)
library(ggplot2)
##### Convert sync file to poolfstat file, and call SNPS

# We first have to give haploid sizes of each pool.
psizes <- as.numeric(c('100','100','100','100','100','100','100','100'))

# Then we give the names of each pool/sample.
pnames <- as.character(c('C1','C2','E1','E2','L1','L2','S1','S2'))

# Here is where we read the sync file and call SNPs. The input file must have the '.sync' extension, and can also be gzipped, like in the example below. The parameters to note are: 

#1) min.rc =  the minimum # reads that an allele needs to have (across all pools) to be called 
#2) min.cov.per.pool = the minimum allowed read count per pool for SNP to be called
#3) max.cov.per.pool = the maximum read count per pool for SNP to be called 
#4) min.maf = the minimum allele frequency (over all pools) for a SNP to be called (note this is obtained from dividing the read counts for the minor allele over the total read coverage) 
#5) nlines.per.readblock = number of lines in sync file to be read simultaneously 

SG.pooldata <- vcf2pooldata(vcf.file = "/2/scratch/TylerA/SSD/bwamap/Sexes_combined_variants.vcf", poolsizes = psizes, poolnames = pnames,
                                     min.rc = 5, min.cov.per.pool = 50, max.cov.per.pool = 400,
                                     min.maf = 0.001, remove.indels = FALSE, nlines.per.readblock = 1e+06)


##### From this file we can compute global and per SNP FSTs
SG.snp.fsts <- computeFST(SG.pooldata, method = "Anova", snp.index = NA)

#Assign SNP location for graphing
crap <- data.frame(SG.pooldata@snp.info, SG.snp.fsts$snp.FST)




##### And we can compute pairwise FSTs
SG.pair.fst <- computePairwiseFSTmatrix(SG.pooldata, method = "Anova",
                                        min.cov.per.pool = 50, max.cov.per.pool = 175,
                                        min.maf = 0.001,
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


