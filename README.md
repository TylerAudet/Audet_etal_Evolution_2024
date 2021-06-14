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

merge.sh

## Filter for quality over 20

Q20.sh

## Add read-groups because they are necessary for future programs

addrg.sh

## Create dictionary file for gatk

dict.sh

## mark indels using gatk and then realign around them

gatk_indel_target.sh

gatk_realign.sh 
