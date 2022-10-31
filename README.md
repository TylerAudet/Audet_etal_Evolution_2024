
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
--indel-window 10
```
```
perl /path/to/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl \
--gtf in.gtf \
--input in.mpileup \
--output out_noindel.mpileup
```

# 10) SNP calling was performed using poolSNP v. 1

This was done chromosome by chromosome to save memory. So `awk` was first used to break mpileup in to chromosomes.

This was done with: 1) all samples seperate (to look for sex specific SNPs), 2) sexes merged (For CMH testing), 3) replicates and sexes merged (for Fst tests)

```
bash /PoolSNP-master/PoolSNP.sh \
mpileup=in.mpileup \
reference=/ref/all_ref.fa \
names=C,E,L,S \
max-cov=0.98 \
min-cov=60 \
min-count=20 \
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



## Make sync files from the SNP calls
Script borrowed from Martin Kapun


````

python /scripts/VCF2sync.py \
--vcf /chrom.vcf \
> /chrom.sync

````
## Calculate Fst

## Calculate CMH


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










