
# 1) Trimming was done with bbduk (bbmap v. 38.86)
```
bbduk.sh \
in1=R1.fastq \
in2=R2.fastq \
out1=R1_trimmed.fastq \
out2=R2_trimmed.fastq \
ref=AllAdapters.fa \
threads=32 ftr=149 ktrim=r k=23 mink=7 hdist=1 tpe tbo \
qtrim=rl trimq=20 minlength=36
```
# 2) Genome was mapped using bwa-mem v. 0.7.17 to Drosophila reference genome 6.23
```
bwa mem -t 32 \
-M ref.fa \
R1_trimmed.fastq \
R2_trimmed.fastq \
> mapped.sam
```

# 3) Sam files were converted to bam files using Samtools v. 1.12 at the same time the core genome was extracted to filter out other contigs including reads that mapped to commensal genomes
```
samtools view -h -@ 32 in.sam | \
awk '{if ($3 == "2L" || $3 == "2R" || $3 == "3L" || $3 == "3R" || $3 == "4" || $3 == "X" || \
$2 == "SN:2L" || $2 == "SN:2R" || $2 == "SN:3L" || $2 == "SN:3R" || $2 == "SN:4" || $2 == "SN:X") {print $0}}' | \
samtools view -b -@ 32 -o out.bam
```


# 4) Supplimentary reads from an additional run of sexuencing were merged together with samtoold v. 1.12

This was done by first seperating each run of sequencing in to two directories called run1 and run2, and finally merging these directories in to a final directory called merged. This method was used to merge suppletary runs of sequencing together as well as to merge sexes for analyses where sexes are pooled together.
```
samtools merge merged.bam \
run1/*.bam \
run2/*.bam
```
# 5) Mark and then remove read groups using Samtools v. 1.12

Done in four stages. Files are first sorted by name, then fixmate is used to add quality tags to reads. Then files are sorted by coordinate and mardup is used to mark the duplicates with the -r flags to remove those duplicate reads.

```
samtools sort -n -@ 32 -o out.bam in.bam

samtools fixmate -m -u -@ 32 in.bam out.bam

samtools sort -@ 32 -o out.bam in.bam

samtools markdup -l 150 -r -s -f stats.txt -d 2500 -@ 32 in.bam out.bam
```

# 6) Read-groups were added using picard v 2.26 and then indels were marked and realigned around using GATK v. 3.8

```
java -jar -Xmx10g /picard.jar AddOrReplaceReadGroups \
INPUT=in.bam \
OUTPUT=out_RG.bam \
SORT_ORDER=coordinate \
RGID=library \
RGLB=library \
RGPL=illumina \
RGSM=Stewart \
RGPU=library \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=SILENT
```
```
java -Xmx32g -jar /GenomeAnalysisTK.jar -I in_RG.bam \
-R /ref/dmel-all-chromosome-r6.23.fasta \
-T RealignerTargetCreator \
-o out.intervals
```
```
java -Xmx10g -jar /GenomeAnalysisTK.jar -I in_RG.bam \
-R /ref/dmel-all-chromosome-r6.23.fasta \
-T IndelRealigner -targetIntervals in.intervals \
-o out.bam
```

# 7) Samtools v. 1.12 was used to create an mpileup with sequence and SNP quality thresholds set to 20 and maximum depth set to 1.5x expected depth to remove suspiciously high coverage areas.
```
samtools mpileup -Q 20 -q 20 -d 300 \
-f /ref/dmel-all-chromosome-r6.23.fasta \
in.bam \
-o out.mpileup
```
mpileup files were created for 1) treatments where sexes were combined 2) treatments where sexes were kept seperately 3) 8 individual bam files created by merging all sequences together and randomly sampling to creat 8 'null' populations.

# 8) Repeteive regions were removed using popoolation v. 1.2.2.
These regions were the known transposable elements in the reference genome version 6.23, other "blacklisted" regions that have been shown to cause issues in SNP calling in the drosophila genome (Amemiya et al. 2019; https://github.com/Boyle-Lab/Blacklist/blob/master/lists/dm6-blacklist.v2.bed.gz), and a created bedfile to isolate regions that show suspiciously high inter-sex Fst and were verified to be transposable or repetative elements.
```
curl -O http://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.23_FB2018_04/fasta/dmel-all-transposon-r6.23.fasta.gz
```

## First Identify repeats using RepeatMasker v. 4
```
/path/to/RepeatMasker \
-pa 20 \
--lib /path/to/transposons/dmel-all-transposon-r6.23.fasta \
--gff /path/to/reference/genome/dmel-all-chromosome-r6.23.fasta	
```

## Next remove repetetive regions using popoolation v. 1.2.2.
```
perl /popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl \
--gtf /ReapeatMasker/output/Dmelgenome/dmel-all-chromosome-#r6.23.fasta.out.gff \
--input in.mpileup \
--output out.mpileup
```
# 9) Identify indels using popoolation 2 v. 1.2 and remove them using popoolation v 1.2 
```
perl /path/to/popoolation2_1201/indel_filtering/identify-indel-regions.pl \
--input in.mpileup \
--output out.gtf \
--indel-window 15
```
```
perl /path/to/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl \
--gtf in.gtf \
--input in.mpileup \
--output out_noindel.mpileup
```

# 10) SNP calling was performed using poolSNP v. 1

This was done chromosome by chromosome to save memory. So `awk` was first used to break mpileup in to chromosomes.

```
bash /PoolSNP-master/PoolSNP.sh \
mpileup=in.mpileup \
reference=/ref/all_ref.fa \
names=C1,C2,E1,E2,L1,L2,S1,S2 \
max-cov=0.98 \
min-cov=30 \
min-count=10 \
min-freq=0.01 \
miss-frac=0.2 \
jobs=32 \
BS=1 \
output=out.vcf
```

# 11) This VCF was filtered for the ENCODE blacklist using bedtools v. 2.3
```
bedtools intersect -v -a in.vcf \
-b /blacklist/dm6-blacklist.v2.bed \
> out.vcf
```















## Make a VCF with poolSNP

for sexes combined:
````

# Need to pull out chromosome by chromosome because it is very memory intensive to ru poolSNP on everything at once
awk '{if ($1 == "2L") {print $0}}' /scratch/audett/SSD/Stewart/Stewart_repeatmasked_indelmasked.mpileup > /scratch/audett/SSD/Stewart/2L.mpileup
awk '{if ($1 == "2R") {print $0}}' /scratch/audett/SSD/Stewart/Stewart_repeatmasked_indelmasked.mpileup > /scratch/audett/SSD/Stewart/2R.mpileup
awk '{if ($1 == "3L") {print $0}}' /scratch/audett/SSD/Stewart/Stewart_repeatmasked_indelmasked.mpileup > /scratch/audett/SSD/Stewart/3L.mpileup
awk '{if ($1 == "3R") {print $0}}' /scratch/audett/SSD/Stewart/Stewart_repeatmasked_indelmasked.mpileup > /scratch/audett/SSD/Stewart/3R.mpileup
awk '{if ($1 == "X") {print $0}}' /scratch/audett/SSD/Stewart.Stewart_repeatmasked_indelmasked.mpileup > /scratch/audett/SSD/Stewart/X.mpileup
awk '{if ($1 == "4") {print $0}}' /scratch/audett/SSD/Stewart/Stewart_repeatmasked_indelmasked.mpileup > /scratch/audett/SSD/Stewart/4.mpileup

# Run poolSNP one chromosome at a time
bash /home/tylera/bin/PoolSNP-master/PoolSNP.sh \
mpileup=/2/scratch/TylerA/SSD/merged/2L.mpileup \
reference=/2/scratch/TylerA/Dmelgenome/dmel-all-chromosome-r6.23.fasta \
names=C,E,L,S \
max-cov=0.98 \
min-cov=40 \
min-count=20 \
min-freq=0.01 \
miss-frac=0.2 \
jobs=16 \
BS=1 \
output=/2/scratch/TylerA/SSD/poolSNP/poolSNP_variants_2L

bash /home/tylera/bin/PoolSNP-master/PoolSNP.sh \
mpileup=/2/scratch/TylerA/SSD/merged/2R.mpileup \
reference=/2/scratch/TylerA/Dmelgenome/dmel-all-chromosome-r6.23.fasta \
names=C,E,L,S \
max-cov=0.98 \
min-cov=40 \
min-count=20 \
min-freq=0.01 \
miss-frac=0.2 \
jobs=16 \
BS=1 \
output=/2/scratch/TylerA/SSD/poolSNP/poolSNP_variants_2R

bash /home/tylera/bin/PoolSNP-master/PoolSNP.sh \
mpileup=/2/scratch/TylerA/SSD/merged/3L.mpileup \
reference=/2/scratch/TylerA/Dmelgenome/dmel-all-chromosome-r6.23.fasta \
names=C,E,L,S \
max-cov=0.98 \
min-cov=40 \
min-count=20 \
min-freq=0.01 \
miss-frac=0.2 \
jobs=16 \
BS=1 \
output=/2/scratch/TylerA/SSD/poolSNP/poolSNP_variants_3L

bash /home/tylera/bin/PoolSNP-master/PoolSNP.sh \
mpileup=/2/scratch/TylerA/SSD/merged/3R.mpileup \
reference=/2/scratch/TylerA/Dmelgenome/dmel-all-chromosome-r6.23.fasta \
names=C,E,L,S \
max-cov=0.98 \
min-cov=40 \
min-count=20 \
min-freq=0.01 \
miss-frac=0.2 \
jobs=16 \
BS=1 \
output=/2/scratch/TylerA/SSD/poolSNP/poolSNP_variants_3R

bash /home/tylera/bin/PoolSNP-master/PoolSNP.sh \
mpileup=/2/scratch/TylerA/SSD/merged/4.mpileup \
reference=/2/scratch/TylerA/Dmelgenome/dmel-all-chromosome-r6.23.fasta \
names=C,E,L,S \
max-cov=0.98 \
min-cov=40 \
min-count=20 \
min-freq=0.01 \
miss-frac=0.2 \
jobs=16 \
BS=1 \
output=/2/scratch/TylerA/SSD/poolSNP/poolSNP_variants_4

bash /home/tylera/bin/PoolSNP-master/PoolSNP.sh \
mpileup=/2/scratch/TylerA/SSD/merged/X.mpileup \
reference=/2/scratch/TylerA/Dmelgenome/dmel-all-chromosome-r6.23.fasta \
names=C,E,L,S \
max-cov=0.98 \
min-cov=40 \
min-count=20 \
min-freq=0.01 \
miss-frac=0.2 \
jobs=16 \
BS=1 \
output=/2/scratch/TylerA/SSD/poolSNP/poolSNP_variants_X

# Defaults to gzip files, so all files need to be gunzipped to sed and cat them
gunzip poolSNP_variants_2L.vcf.gz
gunzip poolSNP_variants_2R.vcf.gz
gunzip poolSNP_variants_3L.vcf.gz
gunzip poolSNP_variants_3R.vcf.gz
gunzip poolSNP_variants_4.vcf.gz
gunzip poolSNP_variants_X.vcf.gz

# Remove the first 18 lines which are summaries from all files other than 2L so when they are cat together they don't have summaries stuck in the middle
sed '1,18d' poolSNP_variants_2R.vcf > 2R.vcf
sed '1,18d' poolSNP_variants_3L.vcf > 3L.vcf
sed '1,18d' poolSNP_variants_3R.vcf > 3R.vcf
sed '1,18d' poolSNP_variants_X.vcf > X.vcf
sed '1,18d' poolSNP_variants_4.vcf > 4.vcf

# Cat all back together to get a completed VCF
cat poolSNP_variants_2L.vcf 2R.vcf 3L.vcf 3R.vcf 4.vcf X.vcf > merged_poolSNP.vcf

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

# We first have to give haploid sizes of each pool.
psizes <- as.numeric(c('400','400','400','400'))

# Then we give the names of each pool/sample.
pnames <- as.character(c('C','E','L','S'))

# Here is where we read the vcf file and call SNPs. 
#1) min.rc =  the minimum # reads that an allele needs to have (across all pools) to be called 
#2) min.cov.per.pool = the minimum allowed read count per pool for SNP to be called
#3) max.cov.per.pool = the maximum read count per pool for SNP to be called 
#4) min.maf = the minimum allele frequency (over all pools) for a SNP to be called (note this is obtained from dividing the read counts for the minor allele over the total read coverage) 
#5) nlines.per.readblock = number of lines in sync file to be read simultaneously 

SG.pooldata <- vcf2pooldata(vcf.file = "/2/scratch/TylerA/SSD/poolSNP/merged_poolSNP.vcf", poolsizes = psizes, poolnames = pnames, min.cov.per.pool = 50, min.maf = 0.05)

#Data consists of 872491 SNPs for 4 Pools


##### And we can compute pairwise FSTs
SG.pair.fst <- compute.pairwiseFST(SG.pooldata, method = "Anova",
                                        output.snp.values = TRUE)

# Extracting fst as a matrix and then making a data.frame with the snp info associated with the fst values
test<-as.matrix(SG.pair.fst@PairwiseSnpFST)
crap <- data.frame(SG.pooldata@snp.info, test[,c(1,4,5,6)])

# Extract the valuess associated with the comparisons we care about
#CVE = controls vs. experimental
#AVE = all samples vs. the experimental
#LVS = large vs. small

CVE <- data.frame(crap[c(1,2,5)])
AVE <- data.frame(crap[c(1,2,5,6,7)])
LVS <- data.frame(crap[c(1,2,8)])

#AVE has 872491 snps

AVE <- AVE[AVE$X1!='NaN',]
#818705  SNPs after removing NaNs from column X1

AVE <- AVE[AVE$X2!='NaN',]
#745111

AVE <- AVE[AVE$X3!='NaN',]
#706251

# LVS has 872491

LVS <- LVS[LVS$X4!='NaN',]
#802721
        
CVE<-na.omit(CVE)
LVS<-na.omit(LVS)
AVE<-na.omit(AVE)

# Changing headers for clarity and to match the rolling average code

headers<-c("ID.Chromosome","ID.Position","Means")
colnames(LVS)<-headers


AVE<-data.frame(ID=AVE[,c(1:2)], Means=rowMeans(AVE[,-c(1:2)], na.rm=TRUE))

# Make data tables to work with without having to re-run fst code

write.table(CVE, file = "/2/scratch/TylerA/SSD/poolSNP/CVE.fst", sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(LVS, file = "/2/scratch/TylerA/SSD/poolSNP/LVS.fst", sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(AVE, file = "/2/scratch/TylerA/SSD/poolSNP/AVE.fst", sep = "\t",
            row.names = FALSE, quote = FALSE)
 
            
# Rolling average code
# This requires changing the name of the file to data, in the future I'm going to make this a shell script. Also requires doing it chromosome by chromosome, this also needs to be fixed in the future.

#Reading table on local
data<-read.table("/2/scratch/TylerA/SSD/poolSNP/AVE.fst",header=TRUE)

# OR #

#data<-read.table("/2/scratch/TylerA/SSD/poolSNP/LVS.fst",header=TRUE)

# Break everything in to chromosomes again so that two different chromosomes are not included in the same window

data2L<-data[data$ID.Chromosome=="2L",]
data2R<-data[data$ID.Chromosome=="2R",]
data3L<-data[data$ID.Chromosome=="3L",]
data3R<-data[data$ID.Chromosome=="3R",]
data4<-data[data$ID.Chromosome=="4",]
dataX<-data[data$ID.Chromosome=="X",]

# Sliding window function

sliding_window<-function(data,window){
temp.vec=vector("list",0)
results.vec=vector("list",0)
x=0
y=1000
i=1
while (i <= length(data$Means)){
 if (data$ID.Position[i] >= x & data$ID.Position[i] < y) {
    temp.vec<-append(temp.vec,data$Means[i])
    i<-i+1
  } else {
    if (data$ID.Position[i] <y) {
      temp.vec<-append(temp.vec,NA)
      i<-i+1
      } else {
    temp.vec<-as.numeric(temp.vec)
    results.vec <- append(results.vec,mean(temp.vec, na.omit=TRUE))
    temp.vec=vector("list",0)
    x <- x + window
    y <- y + window
      }
  }
}
print(results.vec)
}

# Running sliding window function on all chromosomes

fst2L<-sliding_window(data2L,1000)
fst2R<-sliding_window(data2R,1000)
fst3L<-sliding_window(data3L,1000)
fst3R<-sliding_window(data3R,1000)
fst4<-sliding_window(data4,1000)
fstX<-sliding_window(dataX,1000)

# Making matrices which are easier to work with

fst2L <- as.matrix(fst2L)
fst2R <- as.matrix(fst2R)
fst3L <- as.matrix(fst3L)
fst3R <- as.matrix(fst3R)
fst4 <- as.matrix(fst4)
fstX <- as.matrix(fstX)

fst2L <- fst2L[which(fst2L!="NaN"),]
fst2R <- fst2R[which(fst2R!="NaN"),]
fst3L <- fst3L[which(fst3L!="NaN"),]
fst3R <- fst3R[which(fst3R!="NaN"),]
fst4 <- fst4[which(fst4!="NaN"),]
fstX <- fstX[which(fstX!="NaN"),]

# Merging all in to one file

fst2L<-cbind(rep("2L",length(fst2L)),fst2L)
fst2R<-cbind(rep("2R",length(fst2R)),fst2R)
fst3L<-cbind(rep("3L",length(fst3L)),fst3L)
fst3R<-cbind(rep("3R",length(fst3R)),fst3R)
fst4<-cbind(rep("4",length(fst4)),fst4)
fstX<-cbind(rep("X",length(fstX)),fstX)
fst_data<-rbind(fstX,fst2L,fst2R,fst3L,fst3R,fst4)

# Changing column names

crap<-c("chr","fst")
colnames(fst_data)<-crap
fst_data<-as.data.frame(fst_data)

# This code adds numeric indicators for each row to make sure the order in the plots are correct

g <- nrow(fst_data[which(fst_data$chr=='2L'),])
h <- nrow(fst_data[which(fst_data$chr=='2R'),])
i <- nrow(fst_data[which(fst_data$chr=='3L'),])
j <- nrow(fst_data[which(fst_data$chr=='3R'),])
k <- nrow(fst_data[which(fst_data$chr=='4'),])
l <- nrow(fst_data[which(fst_data$chr=='X'),])

fst_data$number <-  c((1:l),
                   (l+1):(l+g), 
                   (l+g+1):(l+g+h), 
                   (l+g+h+1):(l+g+h+i),
                   (l+g+h+i+1):(l+g+h+i+j),
                   (l+g+h+i+j+1):(l+g+h+i+j+k))

fst_data <- as.data.frame(lapply(fst_data, unlist))
fst_data$number<-as.numeric(fst_data$number)
fst_data$chr<-as.factor(fst_data$chr)

# Plotting

plot<-ggplot(fst_data, aes(x=as.numeric(number), y=as.numeric(fst), color=as.factor(chr))) +
  geom_point(size=0.5, show.legend = F, alpha = 0.2) +
  theme(panel.background = element_blank()) +
  scale_colour_manual(values=c("seagreen", "darkslateblue", 'darkred', 'darkorchid4', 'darkolivegreen', 'darkblue')) +
  theme(text = element_text(size=20),
        axis.text.x= element_text(size=15), 
        axis.text.y= element_text(size=15))
        
        
# Need to change file name !!!!!!!!!

png("/2/scratch/TylerA/SSD/poolSNP/AVE_Fst.png",type="cairo")
plot
dev.off()

````

## Plot Fst

# Experimental vs. All other treatments

AVE_Fst.png![AVE_Fst](https://user-images.githubusercontent.com/77504755/127500637-8264e8d5-9e5d-4122-9a0f-a07affd8b8a7.png)

# Checking vcf coverage

````
vcftools --vcf /2/scratch/TylerA/SSD/poolSNP/merged_poolSNP.vcf --site-mean-depth

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
#cov.mat3 <- cov.mat[,c(5,6,7,8)]
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

# Same plotting code as done in the fst scripts

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

CVE_pval_plot.png![CVE_pval_plot](https://user-images.githubusercontent.com/77504755/125775510-a4b78588-172b-404f-8d80-94c8a6275473.png)


# Large vs. Small

LVS_pval_plot.png![LVS_pval_plot](https://user-images.githubusercontent.com/77504755/125775533-a6b0fbc0-0f2a-4f32-bbec-5f08b17446cc.png)

## Creating a suspicious coverage .gff

A high Fst found in male vs. female control comparisons in the conflict_workflow highlights some areas that are most likely not true. So I created a file with those locations with an fst >0.1 in male vs. female comparisons in control populations with a 10bp before and after window to remove from my genomes.

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

I then use this suspicious_coverage.gff and repeat masker to mask these suspicious areas. This was done at the same stage of the original repeat masking and all steps afterwards were repeated with these suspicious areas removed.


# Comparing areas with high Fst and statistically significant CMH values

I want to pull out all spots in merged_AVE with an fst > 0.75 and then compare them to spots with a CMH padj < 0.01

````

awk '$3 > 0.75  {print $1, $2, $3}' /2/scratch/TylerA/SSD/poolSNP/AVE.fst > /2/scratch/TylerA/SSD/results/highfst_AVE.fst

awk -F "," '$5 < 0.01 {print $2, $3, $5}' /2/scratch/TylerA/SSD/bwamap/CVE_pval.csv > /2/scratch/TylerA/SSD/results/highpval_CVE.csv

sed "s/\"//g;s/,/\t/g" /2/scratch/TylerA/SSD/results/highpval_CVE.csv > /2/scratch/TylerA/SSD/results/highpval_CVE.cmh

````

There are no matching SNPs between these two files. So I'm going to make a bed file with a +/- 10bp buffer for each to look for close together SNPs

````
rm(list=ls())

data<-read.table("/2/scratch/TylerA/SSD/results/highfst_AVE.fst")

fst<-data

fst$V2<-as.numeric(fst$V2)

fst<-fst[-1,]

fst$start<-fst$V2-10
fst$end<-fst$V2+10

fst<-data.frame(fst[c(1,4,5)])

fst$V1 <- sub("^", "chr", fst$V1 )

headers<-c("chrom","chromStart","chromEnd")
colnames(fst)<-headers



write.table(fst,file="/2/scratch/TylerA/SSD/results/highfst_AVE.bed",sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)


cmh<-read.table("/2/scratch/TylerA/SSD/results/highpval_CVE.cmh")


cmh$V2<-as.numeric(cmh$V2)

cmh<-cmh[-1,]
cmh<-cmh[-1,]
cmh<-cmh[-1,]
cmh<-cmh[-1,]

cmh$start<-cmh$V2-10
cmh$end<-cmh$V2+10

cmh<-data.frame(cmh[c(1,4,5)])

cmh$V1 <- sub("^", "chr", cmh$V1 )

colnames(cmh)<-headers


write.table(cmh,file="/2/scratch/TylerA/SSD/results/highpval_CVE.bed",sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)


````

Now I want to look for overlap in these bed files

````
bedtools intersect -a /2/scratch/TylerA/SSD/results/highpval_CVE.bed -b /2/scratch/TylerA/SSD/results/highfst_AVE.bed > /2/scratch/TylerA/SSD/results/matches.bed

````

There are 462 overlapping regions

````
bedtools intersect -a /2/scratch/TylerA/Dmelgenome/dmel-all-r6.23.gtf -b /2/scratch/TylerA/SSD/results/matches.bed

sed 's/chr//g' /2/scratch/TylerA/SSD/results/matches.bed > /2/scratch/TylerA/SSD/results/matches-vcf.bed

vcftools --vcf /2/scratch/TylerA/SSD/poolSNP/merged_poolSNP.vcf --bed /2/scratch/TylerA/SSD/results/matches-vcf.bed --out /2/scratch/TylerA/SSD/results/genes-of-interest --recode --keep-INFO-all

java -Xmx32g -jar /usr/local/gatk/GenomeAnalysisTK.jar -R /2/scratch/TylerA/Dmelgenome/gatk/dmel-all-chromosome-r6.23.fasta -V genes-of-interest.recode.vcf  -T VariantsToTable -F CHROM -F POS -F TYPE -F REF -F ALT -o interesting_loci.table

java -Xmx8g -jar ~/bin/snpEff/snpEff.jar -ud 0 Drosophila_melanogaster interesting_loci.table > loci.ann.table

grep -v 'NO_VARIATION'  loci.ann.table > loci.ann.filtered.table
grep -v 'intron_variant'  loci.ann.filtered.table > loci.ann.filtered.twice.table

awk '{print $8}'  > loci.ann.filtered.twice.table test.table

awk -F '[\|]' '{ print $4 }' test.table > gene_ids.table


awk '{print $1 " " $2}' loci.ann.filtered.twice.table test.table > locations.table

#First 6 lines are uninformative or white space

sed -i '1,6d' locations.table
sed -i '1,6d' gene_ids.table
````
This gives me 130 genes that have both an fst > 0.75 and a CMH padj < 0.01 with SNPs in a genic region.

# Looking for selective sweeps with Tajima's D

````

awk '{print $1,$2,$3,$4,$5,$6}' Sexes_combined_norepeats_nosus.mpileup > C1.pileup
awk '{print $1,$2,$3,$7,$8,$9}' Sexes_combined_norepeats_nosus.mpileup > C2.pileup
awk '{print $1,$2,$3,$10,$11,$12}' Sexes_combined_norepeats_nosus.mpileup > E1.pileup
awk '{print $1,$2,$3,$13,$14,$15}' Sexes_combined_norepeats_nosus.mpileup > E2.pileup

awk '{print $1,$2,$3,$16,$17,$18}' Sexes_combined_norepeats_nosus.mpileup > L1.pileup
awk '{print $1,$2,$3,$19,$20,$21}' Sexes_combined_norepeats_nosus.mpileup > L2.pileup
awk '{print $1,$2,$3,$22,$23,$24}' Sexes_combined_norepeats_nosus.mpileup > S1.pileup
awk '{print $1,$2,$3,$25,$26,$27}' Sexes_combined_norepeats_nosus.mpileup > S2.pileup


python ~/bin/1.4.4/pool-hmm.py -f /2/scratch/TylerA/SSD/bwamap/E2 -n 400 -R 2L -a unknown -P 8 -p -k 0.001 --theta 0.005
python ~/bin/1.4.4/pool-hmm.py -f /2/scratch/TylerA/SSD/bwamap/E2 -n 400 -R 2R -a unknown -P 8 -p -k 0.001 --theta 0.005
python ~/bin/1.4.4/pool-hmm.py -f /2/scratch/TylerA/SSD/bwamap/E2 -n 400 -R 3L -a unknown -P 8 -p -k 0.001 --theta 0.005
python ~/bin/1.4.4/pool-hmm.py -f /2/scratch/TylerA/SSD/bwamap/E2 -n 400 -R 3R -a unknown -P 8 -p -k 0.001 --theta 0.005
python ~/bin/1.4.4/pool-hmm.py -f /2/scratch/TylerA/SSD/bwamap/E2 -n 300 -R X -a unknown -P 8 -p -k 0.001 --theta 0.005
python ~/bin/1.4.4/pool-hmm.py -f /2/scratch/TylerA/SSD/bwamap/E2 -n 400 -R 4 -a unknown -P 8 -p -k 0.001 --theta 0.005

python ~/bin/1.4.4/pool-hmm.py -f /2/scratch/TylerA/SSD/bwamap/E1 -n 400 -R 2L -a unknown -P 8 -p -k 0.001 --theta 0.005
python ~/bin/1.4.4/pool-hmm.py -f /2/scratch/TylerA/SSD/bwamap/E1 -n 400 -R 2R -a unknown -P 8 -p -k 0.001 --theta 0.005
python ~/bin/1.4.4/pool-hmm.py -f /2/scratch/TylerA/SSD/bwamap/E1 -n 400 -R 3L -a unknown -P 8 -p -k 0.001 --theta 0.005
python ~/bin/1.4.4/pool-hmm.py -f /2/scratch/TylerA/SSD/bwamap/E1 -n 400 -R 3R -a unknown -P 8 -p -k 0.001 --theta 0.005
python ~/bin/1.4.4/pool-hmm.py -f /2/scratch/TylerA/SSD/bwamap/E1 -n 300 -R X -a unknown -P 8 -p -k 0.001 --theta 0.005
python ~/bin/1.4.4/pool-hmm.py -f /2/scratch/TylerA/SSD/bwamap/E1 -n 400 -R 4 -a unknown -P 8 -p -k 0.001 --theta 0.005



~~~~~

sed 1d E1_2L.stat > E1_2L.txt
sed 1d E1_2R.stat > E1_2R.txt
sed 1d E1_3L.stat > E1_3L.txt
sed 1d E1_3R.stat > E1_3R.txt
sed 1d E1_4.stat > E1_4.txt
sed 1d E1_X.stat > E1_X.txt

sed 1d E2_2L.stat > E2_2L.txt
sed 1d E2_2R.stat > E2_2R.txt
sed 1d E2_3L.stat > E2_3L.txt
sed 1d E2_3R.stat > E2_3R.txt
sed 1d E2_4.stat > E2_4.txt
sed 1d E2_X.stat > E2_X.txt



####### To R #######

E1_2L<-read.table("/2/scratch/TylerA/SSD/results/E1_2L.txt")
E1_2R<-read.table("/2/scratch/TylerA/SSD/results/E1_2R.txt")
E1_3L<-read.table("/2/scratch/TylerA/SSD/results/E1_3L.txt")
E1_3R<-read.table("/2/scratch/TylerA/SSD/results/E1_3R.txt")
E1_4<-read.table("/2/scratch/TylerA/SSD/results/E1_4.txt")
E1_X<-read.table("/2/scratch/TylerA/SSD/results/E1_X.txt")



headers<-c("chrom","chromStart","chromEnd")

E1_2L$chrom<-c("chr2L")
E1_2R$chrom<-c("chr2R")
E1_3L$chrom<-c("chr3L")
E1_3R$chrom<-c("chr3R")
E1_4$chrom<-c("chr4")
E1_X$chrom<-c("chrX")

E1_2L <- subset(E1_2L, select=c(chrom,V1,V2))
E1_2R <- subset(E1_2R, select=c(chrom,V1,V2))
E1_3L <- subset(E1_3L, select=c(chrom,V1,V2))
E1_3R <- subset(E1_3R, select=c(chrom,V1,V2))
E1_4 <- subset(E1_4, select=c(chrom,V1,V2))
E1_X <- subset(E1_X, select=c(chrom,V1,V2))


colnames(E1_2L)<-headers
colnames(E1_2R)<-headers
colnames(E1_3L)<-headers
colnames(E1_3R)<-headers
colnames(E1_4)<-headers
colnames(E1_X)<-headers

E1<-rbind(E1_2L,E1_2R,E1_3L,E1_3R,E1_4,E1_X)

write.table(E1,file="/2/scratch/TylerA/SSD/results/E1_sweep.bed",sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)


E2_2L<-read.table("/2/scratch/TylerA/SSD/results/E2_2L.txt")
E2_2R<-read.table("/2/scratch/TylerA/SSD/results/E2_2R.txt")
E2_3L<-read.table("/2/scratch/TylerA/SSD/results/E2_3L.txt")
E2_3R<-read.table("/2/scratch/TylerA/SSD/results/E2_3R.txt")
E2_4<-read.table("/2/scratch/TylerA/SSD/results/E2_4.txt")
E2_X<-read.table("/2/scratch/TylerA/SSD/results/E2_X.txt")

E2_2L$chrom<-c("chr2L")
E2_2R$chrom<-c("chr2R")
E2_3L$chrom<-c("chr3L")
E2_3R$chrom<-c("chr3R")
E2_4$chrom<-c("chr4")
E2_X$chrom<-c("chrX")

E2_2L <- subset(E2_2L, select=c(chrom,V1,V2))
E2_2R <- subset(E2_2R, select=c(chrom,V1,V2))
E2_3L <- subset(E2_3L, select=c(chrom,V1,V2))
E2_3R <- subset(E2_3R, select=c(chrom,V1,V2))
E2_4 <- subset(E2_4, select=c(chrom,V1,V2))
E2_X <- subset(E2_X, select=c(chrom,V1,V2))


colnames(E2_2L)<-headers
colnames(E2_2R)<-headers
colnames(E2_3L)<-headers
colnames(E2_3R)<-headers
colnames(E2_4)<-headers
colnames(E2_X)<-headers

E2<-rbind(E2_2L,E2_2R,E2_3L,E2_3R,E2_4,E2_X)

write.table(E2,file="/2/scratch/TylerA/SSD/results/E2_sweep.bed",sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)



# Back to terminal

bedtools intersect -a matches.bed -b E1_sweep.bed > matched_sweeps.bed
bedtools intersect -a matched_sweeps.bed -b E2_sweep.bed > sweeps.bed






pool-hmm code: python ~/bin/1.4.4/pool-hmm.py -f /2/scratch/TylerA/SSD/bwamap/E1 -n 6 -R 2L -a unknown -P 8 -p -k 0.001 --theta 0.005
-n: number of chromosomes
-R: region or chromosome to look at (need to do 1 at a time because of memory intensiveness)
-a: site frequency spectrum, can be provided or calculated first if unknown
-P: threads
-p: tells it to actually give you results (predict selective sweeps)
-k: per site transition probability between hidden states. Used the number listed in the example in the README.txt because they used Drosophila as the example data
-theta: scaled mutation rate. used the example number for the same reason as -k




````




# Calculataing Tajima's D

````


~~~

python /2/scratch/TylerA/SSD/scripts/VCF2sync.py \
--vcf poolSNP_reps_2L.vcf \
> 2L.sync


python /2/scratch/TylerA/SSD/scripts/VCF2sync.py \
--vcf poolSNP_reps_2R.vcf \
> 2R.sync

python /2/scratch/TylerA/SSD/scripts/VCF2sync.py \
--vcf poolSNP_reps_3L.vcf \
> 3L.sync

python /2/scratch/TylerA/SSD/scripts/VCF2sync.py \
--vcf poolSNP_reps_3R.vcf \
> 3R.sync

python /2/scratch/TylerA/SSD/scripts/VCF2sync.py \
--vcf poolSNP_reps_4.vcf \
> 4.sync

python /2/scratch/TylerA/SSD/scripts/VCF2sync.py \
--vcf poolSNP_reps_X.vcf \
> X.sync

#resample

python /2/scratch/TylerA/SSD/scripts/SubsampleSync.py \
--sync 2L.sync \
--target-cov 300 \
--min-cov 50 \
> subsample_2L.sync

python /2/scratch/TylerA/SSD/scripts/SubsampleSync.py \
--sync 2R.sync \
--target-cov 300 \
--min-cov 50 \
> subsample_2R.sync

python /2/scratch/TylerA/SSD/scripts/SubsampleSync.py \
--sync 3L.sync \
--target-cov 300 \
--min-cov 50 \
> subsample_3L.sync

python /2/scratch/TylerA/SSD/scripts/SubsampleSync.py \
--sync 3R.sync \
--target-cov 300 \
--min-cov 50 \
> subsample_3R.sync

python /2/scratch/TylerA/SSD/scripts/SubsampleSync.py \
--sync 4.sync \
--target-cov 300 \
--min-cov 50 \
> subsample_4.sync

python /2/scratch/TylerA/SSD/scripts/SubsampleSync.py \
--sync X.sync \
--target-cov 300 \
--min-cov 50 \
> subsample_X.sync

#calculate true window size

#awk '{if ($1 == "2L") {print $0}}' /2/scratch/TylerA/SSD/bwamap/Sexes_combined_norepeat_nosus.gtf > ./2L_indels.gtf

#awk '{if ($1 == "2L") {print $0}}' /2/scratch/TylerA/Dmelgenome/dmel-all-chromosome-r6.23.fasta.out.gff > ./2L_transposons.gtf




python /2/scratch/TylerA/SSD/scripts/TrueWindows.py \
--badcov poolSNP_reps_2L_BS.txt \
--window 10000 \
--step 10000 \
--chromosomes 2L:23513712 \
--output truewindows_2L

python /2/scratch/TylerA/SSD/scripts/TrueWindows.py \
--badcov poolSNP_reps_2R_BS.txt \
--window 10000 \
--step 10000 \
--chromosomes 2R:25286936 \
--output truewindows_2R


python /2/scratch/TylerA/SSD/scripts/TrueWindows.py \
--badcov poolSNP_reps_3L_BS.txt \
--window 10000 \
--step 10000 \
--chromosomes 3L:28110227 \
--output truewindows_3L


python /2/scratch/TylerA/SSD/scripts/TrueWindows.py \
--badcov poolSNP_reps_3R_BS.txt \
--window 10000 \
--step 10000 \
--chromosomes 3R:32079331 \
--output truewindows_3R

python /2/scratch/TylerA/SSD/scripts/TrueWindows.py \
--badcov poolSNP_reps_4_BS.txt \
--window 10000 \
--step 10000 \
--chromosomes 4:1348131 \
--output truewindows_4


python /2/scratch/TylerA/SSD/scripts/TrueWindows.py \
--badcov poolSNP_reps_X_BS.txt \
--window 10000 \
--step 10000 \
--chromosomes X:23542271 \
--output truewindows_X


# Calculate pop parameters

python /2/scratch/TylerA/SSD/scripts/PoolGen_var.py \
--input 200_subsample_2L.sync \
--pool-size 200,200,200,200,200,200,200,200 \
--min-count 2 \
--window 10000 \
--step 10000 \
--sitecount truewindows_2L-10000-10000.txt \
--min-sites-frac 0.75 \
--output 200_2L_D

python /2/scratch/TylerA/SSD/scripts/PoolGen_var.py \
--input subsample_2R.sync \
--pool-size 200,200,200,200,200,200,200,200 \
--min-count 2 \
--window 10000 \
--step 10000 \
--sitecount truewindows_2R-10000-10000.txt \
--min-sites-frac 0.75 \
--output 2R_D

python /2/scratch/TylerA/SSD/scripts/PoolGen_var.py \
--input subsample_3L.sync \
--pool-size 200,200,200,200,200,200,200,200 \
--min-count 2 \
--window 10000 \
--step 10000 \
--sitecount truewindows_3L-10000-10000.txt \
--min-sites-frac 0.75 \
--output 3L_D

python /2/scratch/TylerA/SSD/scripts/PoolGen_var.py \
--input subsample_3R.sync \
--pool-size 200,200,200,200,200,200,200,200 \
--min-count 2 \
--window 10000 \
--step 10000 \
--sitecount truewindows_3R-10000-10000.txt \
--min-sites-frac 0.75 \
--output 3R_D

python /2/scratch/TylerA/SSD/scripts/PoolGen_var.py \
--input subsample_4.sync \
--pool-size 200,200,200,200,200,200,200,200 \
--min-count 2 \
--window 10000 \
--step 10000 \
--sitecount truewindows_4-10000-10000.txt \
--min-sites-frac 0.75 \
--output 4_D

python /2/scratch/TylerA/SSD/scripts/PoolGen_var.py \
--input subsample_X.sync \
--pool-size 200,200,200,200,200,200,200,200 \
--min-count 2 \
--window 10000 \
--step 10000 \
--sitecount truewindows_X-10000-10000.txt \
--min-sites-frac 0.75 \
--output X_D

~~~~


cat 2L_D_10000_10000.D 2R_D_10000_10000.D 3L_D_10000_10000.D 3R_D_10000_10000.D 4_D_10000_10000.D X_D_10000_10000.D > all.D

cat 2L_D_10000_10000.pi 2R_D_10000_10000.pi 3L_D_10000_10000.pi 3R_D_10000_10000.pi 4_D_10000_10000.pi X_D_10000_10000.pi > all.pi

~~~
````

Now graphing in R and extracting interesting outliers (top and bottom 1%)

````
rm(list=ls())

library(dplyr)
library(ggplot2)


D<-read.table("~/Desktop/all.D")
pi<-read.table("~/Desktop/all.pi")


ddat2<-D

headers<-c("chr","pos","C1","C2","E1","E2","L1","L2","S1","S2")

colnames(ddat2)<-headers

# Convert na in column D to 0

ddat2 <- data.frame(lapply(ddat2, function(x) {
  gsub("na", "0", x)
}))

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

ddat2$E1<-as.numeric(ddat2$E1)
ddat2$E2<-as.numeric(ddat2$E2)
ddat2$C1<-as.numeric(ddat2$C1)
ddat2$C2<-as.numeric(ddat2$C2)
ddat2$L1<-as.numeric(ddat2$L1)
ddat2$L2<-as.numeric(ddat2$L2)
ddat2$S1<-as.numeric(ddat2$S1)
ddat2$S2<-as.numeric(ddat2$S2)


#filter(ddat2, chr == "2L") %>%

  ggplot(ddat2, aes(x=number,y=E1)) +
  geom_point(size=0.5, show.legend = F, alpha=0.25) +
  facet_grid(vars(chr)) +
  theme(panel.background = element_blank()) +
  geom_smooth(aes(y=E1), method="auto", col="blue")

filter(ddat2, chr == "2R") %>%
  ggplot(aes(x=number,y=E1)) +
  geom_point(size=0.5, show.legend = F, alpha=0.25) +
  theme(panel.background = element_blank()) +
  geom_smooth(aes(y=E1), method="auto", col="blue")

filter(ddat2, chr == "3L") %>%
  ggplot(aes(x=number,y=E1)) +
  geom_point(size=0.5, show.legend = F, alpha=0.25) +
  theme(panel.background = element_blank()) +
  geom_smooth(aes(y=E1), method="auto", col="blue")

filter(ddat2, chr == "3R") %>%
  ggplot(aes(x=number,y=E1)) +
  geom_point(size=0.5, show.legend = F, alpha=0.25) +
  theme(panel.background = element_blank()) +
  geom_smooth(aes(y=E1), method="auto", col="blue")





plot_E2<-ggplot(ddat2,aes(x=number,y=E2,color=chr)) +
  geom_point(size=0.5, show.legend = F, alpha=0.25) +
  scale_colour_manual(values=c("seagreen", "darkslateblue", 'darkred', 'darkorchid4', 'darkolivegreen', 'darkblue')) +
  theme(panel.background = element_blank()) +
  geom_smooth(aes(y=E2), method="auto", col="blue")

plot_C1<-ggplot(ddat2,aes(x=number,y=C1,color=chr)) +
  geom_point(size=0.5, show.legend = F, alpha=0.25) +
  theme(panel.background = element_blank()) +
  geom_smooth(aes(y=C1), method="auto", col="blue")

plot_C2<-ggplot(ddat2,aes(x=number,y=C2,color=chr)) +
  geom_point(size=0.5, show.legend = F, alpha=0.25) +
  theme(panel.background = element_blank()) +
  geom_smooth(aes(y=C2), method="auto", col="blue")

plot_L1<-ggplot(ddat2,aes(x=number,y=L1,color=chr)) +
  geom_point(size=0.5, show.legend = F, alpha=0.25) +
  theme(panel.background = element_blank()) +
  geom_smooth(aes(y=L1), method="auto", col="blue")

plot_L2<-ggplot(ddat2,aes(x=number,y=L2,color=chr)) +
  geom_point(size=0.5, show.legend = F, alpha=0.25) +
  theme(panel.background = element_blank()) +
  geom_smooth(aes(y=L2), method="auto", col="blue")

plot_S1<-ggplot(ddat2,aes(x=number,y=S1,color=chr)) +
  geom_point(size=0.5, show.legend = F, alpha=0.25) +
  theme(panel.background = element_blank()) +
  geom_smooth(aes(y=S1), method="auto", col="blue")

plot_S2<-ggplot(ddat2,aes(x=number,y=S2,color=chr)) +
  geom_point(size=0.5, show.legend = F, alpha=0.25) +
  theme(panel.background = element_blank()) +
  geom_smooth(aes(y=S2), method="auto", col="blue")

##


with(ddat2, tapply(E1, chr,
                   quantile, probs = c(0, 0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999, 1),
                   na.rm = TRUE))

# 1% for E1
#2l: -2.589, 5.147
#2R: -2.419, 5.021
#3L: -2.724, 4.723
#3R: -2.266, 4.96
#4: -1.93, 1.307
#X: -2.135, 4.016


with(ddat2, tapply(E2, chr,
                   quantile, probs = c(0, 0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999, 1),
                   na.rm = TRUE))

# 1 for E2
#2L: -2.68, 5.12
#2R: -2.525, 4.91 
#3L: -2.45, 5.03
#3R: -2.46, 4.98
#4: -1.86, 0.76
#X: -2.46, 4.31


### Make bed files for D

library(dplyr)


E1_D <- data.frame(D[c(1,2,5)])
E2_D <- data.frame(D[c(1,2,6)])

E1_D$V1 <- sub("^", "chr", E1_D$V1)
E2_D$V1 <- sub("^", "chr", E2_D$V1)

E1_D$V2<-as.numeric(E1_D$V2)
E2_D$V2<-as.numeric(E2_D$V2)

E1_D$start<-E1_D$V2-10000
E1_D$end<-E1_D$V2+10000

E2_D$start<-E2_D$V2-10000
E2_D$end<-E2_D$V2+10000

E1_D_2L <- filter(E1_D, V1 == "chr2L") %>%
  filter(V5 <= -2.589 | V5 >=  5.147)

E1_D_2R <- filter(E1_D, V1 == "chr2R") %>%
  filter(V5 <= -2.419 | V5 >=  5.021)

E1_D_3L <- filter(E1_D, V1 == "chr3L") %>%
  filter(V5 <= -2.724 | V5 >=  4.723)

E1_D_3R <- filter(E1_D, V1 == "chr3R") %>%
  filter(V5 <= -2.266 | V5 >=  4.96)

E1_D_4 <- filter(E1_D, V1 == "chr4") %>%
  filter(V5 <= -1.93 | V5 >=  1.307)

E1_D_X <- filter(E1_D, V1 == "chrX") %>%
  filter(V5 <= -2.135 | V5 >=  4.016)

# 1% for E1
#2l: -2.589, 5.147
#2R: -2.419, 5.021
#3L: -2.724, 4.723
#3R: -2.266, 4.96
#4: -1.93, 1.307
#X: -2.135, 4.016

E2_D_2L <- filter(E2_D, V1 == "chr2L") %>%
  filter(V6 <= -2.68 | V6 >= 5.12)

E2_D_2R <- filter(E2_D, V1 == "chr2R") %>%
  filter(V6 <= -2.53 | V6 >= 4.91)

E2_D_3L <- filter(E2_D, V1 == "chr3L") %>%
  filter(V6 <= -2.45 | V6 >= 5.03)

E2_D_3R <- filter(E2_D, V1 == "chr3R") %>%
  filter(V6 <= -2.46 | V6 >= 4.98)

E2_D_4 <- filter(E2_D, V1 == "chr4") %>%
  filter(V6 <= -1.86 | V6 >= 0.76)

E2_D_X <- filter(E2_D, V1 == "chrX") %>%
  filter(V6 <= -2.46 | V6 >= 4.31)

# 1 for E2
#2L: -2.68, 5.12
#2R: -2.525, 4.91 
#3L: -2.45, 5.03
#3R: -2.46, 4.98
#4: -1.86, 0.76
#X: -2.46, 4.31

E1_D_2L <- E1_D_2L[,c(1,4,5)]
E1_D_2R <- E1_D_2R[,c(1,4,5)]
E1_D_3L <- E1_D_3L[,c(1,4,5)]
E1_D_3R <- E1_D_3R[,c(1,4,5)]
E1_D_4 <- E1_D_4[,c(1,4,5)]
E1_D_X <- E1_D_X[,c(1,4,5)]

E2_D_2L <- E2_D_2L[,c(1,4,5)]
E2_D_2R <- E2_D_2R[,c(1,4,5)]
E2_D_3L <- E2_D_3L[,c(1,4,5)]
E2_D_3R <- E2_D_3R[,c(1,4,5)]
E2_D_4 <- E2_D_4[,c(1,4,5)]
E2_D_X <- E2_D_X[,c(1,4,5)]

E1_D_final <- rbind(E1_D_2L, E1_D_2R, E1_D_3L, E1_D_3R, E1_D_4, E1_D_X)
E2_D_final <- rbind(E2_D_2L, E2_D_2R, E2_D_3L, E2_D_3R, E2_D_4, E2_D_X)


headers<-c("chrom","chromStart","chromEnd")
colnames(E1_D_final)<-headers
colnames(E2_D_final)<-headers

write.table(E1_D_final,file="~/Desktop/E1_D.bed",sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)
write.table(E2_D_final,file="~/Desktop/E2_D.bed",sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)

````


# Running pool-hmm to find homologous areas




awk '{print $1,$2,$3,$4,$5,$6}' Sexes_combined_norepeats_nosus.mpileup > C1.pileup
awk '{print $1,$2,$3,$7,$8,$9}' Sexes_combined_norepeats_nosus.mpileup > C2.pileup
awk '{print $1,$2,$3,$10,$11,$12}' Sexes_combined_norepeats_nosus.mpileup > E1.pileup
awk '{print $1,$2,$3,$13,$14,$15}' Sexes_combined_norepeats_nosus.mpileup > E2.pileup

awk '{print $1,$2,$3,$16,$17,$18}' Sexes_combined_norepeats_nosus.mpileup > L1.pileup
awk '{print $1,$2,$3,$19,$20,$21}' Sexes_combined_norepeats_nosus.mpileup > L2.pileup
awk '{print $1,$2,$3,$22,$23,$24}' Sexes_combined_norepeats_nosus.mpileup > S1.pileup
awk '{print $1,$2,$3,$25,$26,$27}' Sexes_combined_norepeats_nosus.mpileup > S2.pileup






python ~/bin/1.4.4/pool-hmm.py -f /2/scratch/TylerA/SSD/bwamap/E2 -n 400 -R 2L -a unknown -P 8 -p -k 0.001 --theta 0.005
python ~/bin/1.4.4/pool-hmm.py -f /2/scratch/TylerA/SSD/bwamap/E2 -n 400 -R 2R -a unknown -P 8 -p -k 0.001 --theta 0.005
python ~/bin/1.4.4/pool-hmm.py -f /2/scratch/TylerA/SSD/bwamap/E2 -n 400 -R 3L -a unknown -P 8 -p -k 0.001 --theta 0.005
python ~/bin/1.4.4/pool-hmm.py -f /2/scratch/TylerA/SSD/bwamap/E2 -n 400 -R 3R -a unknown -P 8 -p -k 0.001 --theta 0.005
python ~/bin/1.4.4/pool-hmm.py -f /2/scratch/TylerA/SSD/bwamap/E2 -n 300 -R X -a unknown -P 8 -p -k 0.001 --theta 0.005
python ~/bin/1.4.4/pool-hmm.py -f /2/scratch/TylerA/SSD/bwamap/E2 -n 400 -R 4 -a unknown -P 8 -p -k 0.001 --theta 0.005

python ~/bin/1.4.4/pool-hmm.py -f /2/scratch/TylerA/SSD/bwamap/E1 -n 400 -R 2L -a unknown -P 8 -p -k 0.001 --theta 0.005
python ~/bin/1.4.4/pool-hmm.py -f /2/scratch/TylerA/SSD/bwamap/E1 -n 400 -R 2R -a unknown -P 8 -p -k 0.001 --theta 0.005
python ~/bin/1.4.4/pool-hmm.py -f /2/scratch/TylerA/SSD/bwamap/E1 -n 400 -R 3L -a unknown -P 8 -p -k 0.001 --theta 0.005
python ~/bin/1.4.4/pool-hmm.py -f /2/scratch/TylerA/SSD/bwamap/E1 -n 400 -R 3R -a unknown -P 8 -p -k 0.001 --theta 0.005
python ~/bin/1.4.4/pool-hmm.py -f /2/scratch/TylerA/SSD/bwamap/E1 -n 300 -R X -a unknown -P 8 -p -k 0.001 --theta 0.005
python ~/bin/1.4.4/pool-hmm.py -f /2/scratch/TylerA/SSD/bwamap/E1 -n 400 -R 4 -a unknown -P 8 -p -k 0.001 --theta 0.005



python ~/bin/1.4.4/pool-hmm.py -f /2/scratch/TylerA/SSD/bwamap/C2 -n 400 -R 2L -a unknown -P 8 -p -k 0.001 --theta 0.005
python ~/bin/1.4.4/pool-hmm.py -f /2/scratch/TylerA/SSD/bwamap/C2 -n 400 -R 2R -a unknown -P 8 -p -k 0.001 --theta 0.005
python ~/bin/1.4.4/pool-hmm.py -f /2/scratch/TylerA/SSD/bwamap/C2 -n 400 -R 3L -a unknown -P 8 -p -k 0.001 --theta 0.005
python ~/bin/1.4.4/pool-hmm.py -f /2/scratch/TylerA/SSD/bwamap/C2 -n 400 -R 3R -a unknown -P 8 -p -k 0.001 --theta 0.005
python ~/bin/1.4.4/pool-hmm.py -f /2/scratch/TylerA/SSD/bwamap/C2 -n 300 -R X -a unknown -P 8 -p -k 0.001 --theta 0.005
python ~/bin/1.4.4/pool-hmm.py -f /2/scratch/TylerA/SSD/bwamap/C2 -n 400 -R 4 -a unknown -P 8 -p -k 0.001 --theta 0.005

python ~/bin/1.4.4/pool-hmm.py -f /2/scratch/TylerA/SSD/bwamap/C1 -n 400 -R 2L -a unknown -P 8 -p -k 0.001 --theta 0.005
python ~/bin/1.4.4/pool-hmm.py -f /2/scratch/TylerA/SSD/bwamap/C1 -n 400 -R 2R -a unknown -P 8 -p -k 0.001 --theta 0.005

~~~~~~
python ~/bin/1.4.4/pool-hmm.py -f /2/scratch/TylerA/SSD/bwamap/C1 -n 400 -R 3L -a unknown -P 8 -p -k 0.001 --theta 0.005
python ~/bin/1.4.4/pool-hmm.py -f /2/scratch/TylerA/SSD/bwamap/C1 -n 400 -R 3R -a unknown -P 8 -p -k 0.001 --theta 0.005
python ~/bin/1.4.4/pool-hmm.py -f /2/scratch/TylerA/SSD/bwamap/C1 -n 300 -R X -a unknown -P 8 -p -k 0.001 --theta 0.005
python ~/bin/1.4.4/pool-hmm.py -f /2/scratch/TylerA/SSD/bwamap/C1 -n 400 -R 4 -a unknown -P 8 -p -k 0.001 --theta 0.005


~~~~~

sed 1d E1_2L.stat > E1_2L.txt
sed 1d E1_2R.stat > E1_2R.txt
sed 1d E1_3L.stat > E1_3L.txt
sed 1d E1_3R.stat > E1_3R.txt
sed 1d E1_4.stat > E1_4.txt
sed 1d E1_X.stat > E1_X.txt

sed 1d E2_2L.stat > E2_2L.txt
sed 1d E2_2R.stat > E2_2R.txt
sed 1d E2_3L.stat > E2_3L.txt
sed 1d E2_3R.stat > E2_3R.txt
sed 1d E2_4.stat > E2_4.txt
sed 1d E2_X.stat > E2_X.txt



####### To R #######

E1_2L<-read.table("/2/scratch/TylerA/SSD/results/E1_2L.txt")
E1_2R<-read.table("/2/scratch/TylerA/SSD/results/E1_2R.txt")
E1_3L<-read.table("/2/scratch/TylerA/SSD/results/E1_3L.txt")
E1_3R<-read.table("/2/scratch/TylerA/SSD/results/E1_3R.txt")
E1_4<-read.table("/2/scratch/TylerA/SSD/results/E1_4.txt")
E1_X<-read.table("/2/scratch/TylerA/SSD/results/E1_X.txt")



headers<-c("chrom","chromStart","chromEnd")

E1_2L$chrom<-c("chr2L")
E1_2R$chrom<-c("chr2R")
E1_3L$chrom<-c("chr3L")
E1_3R$chrom<-c("chr3R")
E1_4$chrom<-c("chr4")
E1_X$chrom<-c("chrX")

E1_2L <- subset(E1_2L, select=c(chrom,V1,V2))
E1_2R <- subset(E1_2R, select=c(chrom,V1,V2))
E1_3L <- subset(E1_3L, select=c(chrom,V1,V2))
E1_3R <- subset(E1_3R, select=c(chrom,V1,V2))
E1_4 <- subset(E1_4, select=c(chrom,V1,V2))
E1_X <- subset(E1_X, select=c(chrom,V1,V2))


colnames(E1_2L)<-headers
colnames(E1_2R)<-headers
colnames(E1_3L)<-headers
colnames(E1_3R)<-headers
colnames(E1_4)<-headers
colnames(E1_X)<-headers

E1<-rbind(E1_2L,E1_2R,E1_3L,E1_3R,E1_4,E1_X)

write.table(E1,file="/2/scratch/TylerA/SSD/results/E1_sweep.bed",sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)


E2_2L<-read.table("/2/scratch/TylerA/SSD/results/E2_2L.txt")
E2_2R<-read.table("/2/scratch/TylerA/SSD/results/E2_2R.txt")
E2_3L<-read.table("/2/scratch/TylerA/SSD/results/E2_3L.txt")
E2_3R<-read.table("/2/scratch/TylerA/SSD/results/E2_3R.txt")
E2_4<-read.table("/2/scratch/TylerA/SSD/results/E2_4.txt")
E2_X<-read.table("/2/scratch/TylerA/SSD/results/E2_X.txt")

E2_2L$chrom<-c("chr2L")
E2_2R$chrom<-c("chr2R")
E2_3L$chrom<-c("chr3L")
E2_3R$chrom<-c("chr3R")
E2_4$chrom<-c("chr4")
E2_X$chrom<-c("chrX")

E2_2L <- subset(E2_2L, select=c(chrom,V1,V2))
E2_2R <- subset(E2_2R, select=c(chrom,V1,V2))
E2_3L <- subset(E2_3L, select=c(chrom,V1,V2))
E2_3R <- subset(E2_3R, select=c(chrom,V1,V2))
E2_4 <- subset(E2_4, select=c(chrom,V1,V2))
E2_X <- subset(E2_X, select=c(chrom,V1,V2))


colnames(E2_2L)<-headers
colnames(E2_2R)<-headers
colnames(E2_3L)<-headers
colnames(E2_3R)<-headers
colnames(E2_4)<-headers
colnames(E2_X)<-headers

E2<-rbind(E2_2L,E2_2R,E2_3L,E2_3R,E2_4,E2_X)

write.table(E2,file="/2/scratch/TylerA/SSD/results/E2_sweep.bed",sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)





bedtools intersect -a matches.bed -b E1_sweep.bed > matched_sweeps.bed
bedtools intersect -a matched_sweeps.bed -b E2_sweep.bed > sweeps.bed






pool-hmm code: python ~/bin/1.4.4/pool-hmm.py -f /2/scratch/TylerA/SSD/bwamap/E1 -n 6 -R 2L -a unknown -P 8 -p -k 0.001 --theta 0.005
-n: number of chromosomes
-R: region or chromosome to look at (need to do 1 at a time because of memory intensiveness)
-a: site frequency spectrum, can be provided or calculated first if unknown
-P: threads
-p: tells it to actually give you results (predict selective sweeps)
-k: per site transition probability between hidden states. Used the number listed in the example in the README.txt because they used Drosophila as the example data
-theta: scaled mutation rate. used the example number for the same reason as -k










