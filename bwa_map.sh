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
dir=/home/audett/projects/def-idworkin/audett/SSD/trimmed/
out=/home/audett/projects/def-idworkin/audett/SSD/mapped/

# make it arrays so you can index them
declare -a forward=( /home/audett/projects/def-idworkin/audett/SSD/trimmed/*_R1_trimmed.fastq )

# get the individual reads.
R1=${forward[${SLURM_ARRAY_TASK_ID}]}
 
# get the prefix of the reads
#prefix=${R1%%_R1_trimmed.fastq}
 
name=${R1}
base=`basename ${name} _R1_trimmed.fastq`
 
bwa mem -t ${SLURM_CPUS_PER_TASK} -M /home/audett/projects/def-idworkin/audett/SSD/ref/dmel-all-chromosome-r6.23.fasta.gz ${dir}${base}_R1_trimmed.fastq \
   ${dir}${base}_R2_trimmed.fastq > ${out}${base}.sam

# -M Mark shorter split hits as secondary (for Picard compatibility).
# -t used to tell BWA to use 32 threads to speed things up

