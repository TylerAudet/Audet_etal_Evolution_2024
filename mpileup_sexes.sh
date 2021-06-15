#!/bin/bash
#SBATCH -t 8:00:00
#SBATCH -A def-idworkin
#SBATCH --mem=32G
#SBATCH --mail-user=audett@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module load samtools

samtools mpileup -Q 20 -q 20 -d 600 \
-f /home/audett/projects/def-idworkin/audett/SSD/ref/dmel-all-chromosome-r6.23.fasta \
/home/audett/projects/def-idworkin/audett/SSD/dedup/*_RG_dedup.bam \
-o Sexes_combined.mpileup

#-Q minimum base quality
#-q minimum mapping quality
#-d max coverage
