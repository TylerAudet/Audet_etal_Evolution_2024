# scripts

# run fastqc
QC.sh

# Trim
bbtrim.sh

# QC
QC_after_trim.sh

# Index genome for bwa

bwa index dmel-all-chromosome-r6.23.fasta.gz

# Map genomes

bwa_map.sh

# Sam files converted to bam

sam2bam.sh


