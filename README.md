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



