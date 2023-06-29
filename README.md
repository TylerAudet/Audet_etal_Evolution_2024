
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

using grenedalf from the mpileups

### Seperate Fst out by sample

order: 1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.3,2.4,2.5,2.6,2.7,2.8,3.4,3.5,3.6,3.7,3.8,4.5,4.6,4.7,4.8,5.6,5.7,5.8,6.7,6.8,7.8
        C1C2,C1E1,C1E2,C1L1,C1L2,C1S1,C1S2,C2E1,C2E2,C2L1,C2L2,C2S1,C2S2,E1E2,E1L1,E1L2,E1S1,E1S2,L1L2,L1S1,L1S2,S1S2

awk -F "," '{print $1, $2, $3, $6}' sexesMerged_fst.csv > C1E1.fst

## Calculate CMH

using ACER in R

# Comparing areas with high Fst and statistically significant CMH values

Now I want to look for overlap in these bed files

````
bedtools intersect -a /home/audett/scratch/SSD/Analysis/sexesMerged/CMH/CVE_top1percent.bed -b /home/audett/scratch/SSD/Analysis/repsMerged/fst/CVE_top5percent.bed > /highfst_lowCMH.bed

````
There are 405 overlapping regions
````
# pull the overlapping areas out of vcf
vcftools --vcf /home/audett/scratch/SSD/Analysis/sexesMerged/sexesMerged_clean.vcf --bed /home/audett/scratch/SSD/Analysis/highfst_lowCMH.bed --out /home/audett/scratch/SSD/Analysis/top5fst_bottom1CMH --recode --keep-INFO-all

grep -v 'NO_VARIATION\|intron_variant'  /home/audett/scratch/SSD/Analysis/top5fst_bottom1CMH_annotated.txt > /home/audett/scratch/SSD/Analysis/top5fst_bottom1CMH_interesting.txt

awk '{print $8}' /home/audett/scratch/SSD/Analysis/top5fst_bottom1CMH_interesting.txt | awk -F '[\|]' '{if (NR>20) print $4}' | 
awk '{ a[$1]++ } END { for (b in a) { print b } }' > genes_highfst_lowcmh.txt

````




````


# Calculataing Tajima's D


# resample chromosomes to the same depth

python /home/audett/projects/def-idworkin/audett/SSD/scripts/SubsampleSync.py \
--sync /home/audett/scratch/SSD/Analysis/sexesMerged/syncs/2L.sync \
--target-cov 300 \
--min-cov 50 \
> /home/audett/scratch/SSD/Analysis/sexesMerged/syncs/subsample_2L.sync

python /home/audett/projects/def-idworkin/audett/SSD/scripts/SubsampleSync.py \
--sync /home/audett/scratch/SSD/Analysis/sexesMerged/syncs/2R.sync \
--target-cov 300 \
--min-cov 50 \
> /home/audett/scratch/SSD/Analysis/sexesMerged/syncs/subsample_2R.sync

python /home/audett/projects/def-idworkin/audett/SSD/scripts/SubsampleSync.py \
--sync /home/audett/scratch/SSD/Analysis/sexesMerged/syncs/3L.sync \
--target-cov 300 \
--min-cov 50 \
> /home/audett/scratch/SSD/Analysis/sexesMerged/syncs/subsample_3L.sync

python /home/audett/projects/def-idworkin/audett/SSD/scripts/SubsampleSync.py \
--sync /home/audett/scratch/SSD/Analysis/sexesMerged/syncs/3R.sync \
--target-cov 300 \
--min-cov 50 \
> /home/audett/scratch/SSD/Analysis/sexesMerged/syncs/subsample_3R.sync

python /home/audett/projects/def-idworkin/audett/SSD/scripts/SubsampleSync.py \
--sync /home/audett/scratch/SSD/Analysis/sexesMerged/syncs/4.sync \
--target-cov 300 \
--min-cov 50 \
> /home/audett/scratch/SSD/Analysis/sexesMerged/syncs/subsample_4.sync

python /home/audett/projects/def-idworkin/audett/SSD/scripts/SubsampleSync.py \
--sync /home/audett/scratch/SSD/Analysis/sexesMerged/syncs/X.sync \
--target-cov 300 \
--min-cov 50 \
> /home/audett/scratch/SSD/Analysis/sexesMerged/syncs/subsample_X.sync




# True windows


python /2/scratch/TylerA/SSD/scripts/TrueWindows.py \
--badcov /home/audett/scratch/SSD/Analysis/sexesMerged/SNPcalls/2L_BS.txt \
--window 5000 \
--step 5000 \
--chromosomes 2L:23513712 \
--output /home/audett/scratch/SSD/Analysis/sexesMerged/pi/truewindows_2L

python /2/scratch/TylerA/SSD/scripts/TrueWindows.py \
--badcov /home/audett/scratch/SSD/Analysis/sexesMerged/SNPcalls/2R_BS.txt \
--window 5000 \
--step 5000 \
--chromosomes 2R:25286936 \
--output /home/audett/scratch/SSD/Analysis/sexesMerged/pi/truewindows_2R


python /2/scratch/TylerA/SSD/scripts/TrueWindows.py \
--badcov /home/audett/scratch/SSD/Analysis/sexesMerged/SNPcalls/3L_BS.txt \
--window 5000 \
--step 5000 \
--chromosomes 3L:28110227 \
--output /home/audett/scratch/SSD/Analysis/sexesMerged/pi/truewindows_3L


python /2/scratch/TylerA/SSD/scripts/TrueWindows.py \
--badcov /home/audett/scratch/SSD/Analysis/sexesMerged/SNPcalls/3R_BS.txt \
--window 5000 \
--step 5000 \
--chromosomes 3R:32079331 \
--output /home/audett/scratch/SSD/Analysis/sexesMerged/pi/truewindows_3R

python /2/scratch/TylerA/SSD/scripts/TrueWindows.py \
--badcov /home/audett/scratch/SSD/Analysis/sexesMerged/SNPcalls/4_BS.txt \
--window 5000 \
--step 5000 \
--chromosomes 4:1348131 \
--output /home/audett/scratch/SSD/Analysis/sexesMerged/pi/truewindows_4


python /2/scratch/TylerA/SSD/scripts/TrueWindows.py \
--badcov /home/audett/scratch/SSD/Analysis/sexesMerged/SNPcalls/X_BS.txt \
--window 5000 \
--step 5000 \
--chromosomes X:23542271 \
--output /home/audett/scratch/SSD/Analysis/sexesMerged/pi/truewindows_X


# Calculate pop parameters

python /home/audett/projects/def-idworkin/audett/SSD/scripts/PoolGen_var.py \
--input /home/audett/scratch/SSD/Analysis/sexesMerged/syncs/subsample_2L.sync \
--pool-size 200,200,200,200,200,200,200,200 \
--min-count 2 \
--window 5000 \
--step 5000 \
--sitecount /home/audett/scratch/SSD/Analysis/sexesMerged/pi/truewindows_2L-5000-5000.txt \
--min-sites-frac 0.75 \
--output window_5000_subsample_300_2L

python /home/audett/projects/def-idworkin/audett/SSD/scripts/PoolGen_var.py \
--input /home/audett/scratch/SSD/Analysis/sexesMerged/syncs/subsample_2R.sync \
--pool-size 200,200,200,200,200,200,200,200 \
--min-count 2 \
--window 5000 \
--step 5000 \
--sitecount /home/audett/scratch/SSD/Analysis/sexesMerged/pi/truewindows_2R-5000-5000.txt \
--min-sites-frac 0.75 \
--output window_5000_subsample_300_2R

python /home/audett/projects/def-idworkin/audett/SSD/scripts/PoolGen_var.py \
--input /home/audett/scratch/SSD/Analysis/sexesMerged/syncs/subsample_3L.sync \
--pool-size 200,200,200,200,200,200,200,200 \
--min-count 2 \
--window 5000 \
--step 5000 \
--sitecount /home/audett/scratch/SSD/Analysis/sexesMerged/pi/truewindows_3L-5000-5000.txt \
--min-sites-frac 0.75 \
--output window_5000_subsample_300_3L

python /home/audett/projects/def-idworkin/audett/SSD/scripts/PoolGen_var.py \
--input /home/audett/scratch/SSD/Analysis/sexesMerged/syncs/subsample_3R.sync \
--pool-size 200,200,200,200,200,200,200,200 \
--min-count 2 \
--window 5000 \
--step 5000 \
--sitecount /home/audett/scratch/SSD/Analysis/sexesMerged/pi/truewindows_3R-5000-5000.txt \
--min-sites-frac 0.75 \
--output window_5000_subsample_300_3R

python /home/audett/projects/def-idworkin/audett/SSD/scripts/PoolGen_var.py \
--input /home/audett/scratch/SSD/Analysis/sexesMerged/syncs/subsample_4.sync \
--pool-size 200,200,200,200,200,200,200,200 \
--min-count 2 \
--window 5000 \
--step 5000 \
--sitecount /home/audett/scratch/SSD/Analysis/sexesMerged/pi/truewindows_4-5000-5000.txt \
--min-sites-frac 0.75 \
--output window_5000_subsample_300_4

python /home/audett/projects/def-idworkin/audett/SSD/scripts/PoolGen_var.py \
--input /home/audett/scratch/SSD/Analysis/sexesMerged/syncs/subsample_X.sync \
--pool-size 200,200,200,200,200,200,200,200 \
--min-count 2 \
--window 5000 \
--step 5000 \
--sitecount /home/audett/scratch/SSD/Analysis/sexesMerged/pi/truewindows_X-5000-5000.txt \
--min-sites-frac 0.75 \
--output window_5000_subsample_300_X

~~~~


cat *D > all_D.txt

cat *pi > all_pi.txt

~~~
````

Now graphing in R and extracting interesting outliers (top and bottom 1%)

````
vcftools --vcf /home/audett/scratch/SSD/Analysis/sexesMerged/SNPcalls/sexesMerged_clean.vcf --bed /home/audett/scratch/SSD/Analysis/regions_of_interest_C1.bed --out /home/audett/scratch/SSD/Analysis/allcriteria_C1 --recode --keep-INFO-all

vcftools --vcf /home/audett/scratch/SSD/Analysis/sexesMerged/SNPcalls/sexesMerged_clean.vcf --bed /home/audett/scratch/SSD/Analysis/regions_of_interest_C2.bed --out /home/audett/scratch/SSD/Analysis/allcriteria_C2 --recode --keep-INFO-all

vcftools --vcf /home/audett/scratch/SSD/Analysis/sexesMerged/SNPcalls/sexesMerged_clean.vcf --bed /home/audett/scratch/SSD/Analysis/regions_of_interest_E1.bed --out /home/audett/scratch/SSD/Analysis/allcriteria_E1 --recode --keep-INFO-all

vcftools --vcf /home/audett/scratch/SSD/Analysis/sexesMerged/SNPcalls/sexesMerged_clean.vcf --bed /home/audett/scratch/SSD/Analysis/regions_of_interest_E2.bed --out /home/audett/scratch/SSD/Analysis/allcriteria_E2 --recode --keep-INFO-all

vcftools --vcf /home/audett/scratch/SSD/Analysis/sexesMerged/SNPcalls/sexesMerged_clean.vcf --bed /home/audett/scratch/SSD/Analysis/regions_of_interest_L1.bed --out /home/audett/scratch/SSD/Analysis/allcriteria_L1 --recode --keep-INFO-all

vcftools --vcf /home/audett/scratch/SSD/Analysis/sexesMerged/SNPcalls/sexesMerged_clean.vcf --bed /home/audett/scratch/SSD/Analysis/regions_of_interest_L2.bed --out /home/audett/scratch/SSD/Analysis/allcriteria_L2 --recode --keep-INFO-all

vcftools --vcf /home/audett/scratch/SSD/Analysis/sexesMerged/SNPcalls/sexesMerged_clean.vcf --bed /home/audett/scratch/SSD/Analysis/regions_of_interest_S1.bed --out /home/audett/scratch/SSD/Analysis/allcriteria_S1 --recode --keep-INFO-all

vcftools --vcf /home/audett/scratch/SSD/Analysis/sexesMerged/SNPcalls/sexesMerged_clean.vcf --bed /home/audett/scratch/SSD/Analysis/regions_of_interest_S2.bed --out /home/audett/scratch/SSD/Analysis/allcriteria_S2 --recode --keep-INFO-all



grep -v 'NO_VARIATION\|intron_variant'  /home/audett/scratch/SSD/Analysis/C1_annotated.txt > /home/audett/scratch/SSD/Analysis/C1_interesting.txt
grep -v 'NO_VARIATION\|intron_variant'  /home/audett/scratch/SSD/Analysis/C2_annotated.txt > /home/audett/scratch/SSD/Analysis/C2_interesting.txt
grep -v 'NO_VARIATION\|intron_variant'  /home/audett/scratch/SSD/Analysis/E1_annotated.txt > /home/audett/scratch/SSD/Analysis/E1_interesting.txt
grep -v 'NO_VARIATION\|intron_variant'  /home/audett/scratch/SSD/Analysis/E2_annotated.txt > /home/audett/scratch/SSD/Analysis/E2_interesting.txt
grep -v 'NO_VARIATION\|intron_variant'  /home/audett/scratch/SSD/Analysis/L1_annotated.txt > /home/audett/scratch/SSD/Analysis/L1_interesting.txt
grep -v 'NO_VARIATION\|intron_variant'  /home/audett/scratch/SSD/Analysis/L2_annotated.txt > /home/audett/scratch/SSD/Analysis/L2_interesting.txt
grep -v 'NO_VARIATION\|intron_variant'  /home/audett/scratch/SSD/Analysis/S1_annotated.txt > /home/audett/scratch/SSD/Analysis/S1_interesting.txt
grep -v 'NO_VARIATION\|intron_variant'  /home/audett/scratch/SSD/Analysis/S2_annotated.txt > /home/audett/scratch/SSD/Analysis/S2_interesting.txt




awk '{print $8}' /home/audett/scratch/SSD/Analysis/top5fst_bottom1CMH_interesting.txt | awk -F '[\|]' '{if (NR>20) print $4}' | 
awk '{ a[$1]++ } END { for (b in a) { print b } }' > genes_highfst_lowcmh.txt





````



