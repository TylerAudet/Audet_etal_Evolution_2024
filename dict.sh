#!/bin/bash
#SBATCH -t 0:30:00
#SBATCH --mem=32G
#SBATCH -A def-idworkin
#SBATCH --mail-user=audett@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

picard/2.23.2

java -Xmx32g -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary R=./dmel-all-chromosome-r6.23.fasta O=dmel-all-chromosome-r6.23.dict
