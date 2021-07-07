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
````

## Make a VCF

make_vcf.sh

for sexes combined:
````
java -Xmx32g -jar \
~/bin/VarScan.v2.3.9.jar \
mpileup2cns \
./combined_norepeats.mpileup \
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

````
java -Xmx32g -jar \
~/bin/VarScan.v2.3.9.jar \
mpileup2cns \
./combined_norepeats.mpileup \
--min-reads2 5 \
--min-coverage 100 \
--p-value 0.1 \
--min-var-freq 0.01 \
--min-freq-for-hom 1 \
--min-avg-qual 20 \
--variants \
--output-vcf 1 \
> ./combined_variants.vcf
````
## Make a sync

````
java -ea -jar /usr/local/popoolation/mpileup2sync.jar --threads 16 --input ./Sexes_combined_norepeats.mpileup --output ./Sexes_combined_norepeats.sync
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

# Control vs. SSD-reverse
CVE.png![CVE](https://user-images.githubusercontent.com/77504755/124128284-a116e780-da4a-11eb-8c2f-906695b1f090.png)

# Large vs. Small

LVS.png![LVS](https://user-images.githubusercontent.com/77504755/124128316-aaa04f80-da4a-11eb-9923-1e9e5c91e144.png)

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
#ne <- estimateNe(p0=af.mat2[,"CMO.L"], pt=af.mat2[,"CMO.R"], 
#           cov0=cov.mat2[,"CMO.L"], covt=cov.mat2[,"CMO.R"], 
#           t=0, method = "P.planII", poolSize=c(150, 150))

#Creating the vars for the CMH test 
rep<-c(1,1,2,2) #Number of replicates
Ne<-c(160,160)
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
````
## Plot CMH -log10








