#!/bin/bash
#SBATCH -t 00-06:00:00
#SBATCH -A def-idworkin
#SBATCH --cpus-per-task=32
#SBATCH --mem=0 # reserve the whole node
#SBATCH --mail-user=audett@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module load samtools
parallel --willcite --results {/.}.bam samtools view -b -q 20 -@ 16 \
::: /home/audett/projects/def-idworkin/audett/SSD/Sexes_merged/*.bam
