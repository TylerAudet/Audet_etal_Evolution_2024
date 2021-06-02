#!/bin/bash
#SBATCH -t 4:00:00 # I am unsure if that is how long it will take
#SBATCH -A def-idworkin
#SBATCH --array=0-32
#SBATCH --cpus-per-task=32
#SBATCH --mem=6G
#SBATCH --mail-user=audett@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module load bwa/0.7.17

# make the common path available in variable (this is optional, and just for the looks
BASE=/home/audett/projects/def-idworkin/audett/SSD

# make it arrays so you can index them
declare -a forward=( ${BASE}/trimmed/*_R1_trimmed.fastq )

# get the individual reads.
R1=${forward[${SLURM_ARRAY_TASK_ID}]}
 
# get the prefix of the reads
prefix=${R1%%_R1_trimmed.fastq}
 
bwa mem -t ${SLURM_CPUS_PER_TASK} -M ${BASE}/ref/dmel-all-chromosome-r6.23.fasta.gz ${prefix}_R1_trimmed.fastq \
   ${prefix}_R2_trimmed.fastq > ${BASE}/mapped${prefix}.sam

# -M Mark shorter split hits as secondary (for Picard compatibility).
# -t used to tell BWA to use 32 threads to speed things up

