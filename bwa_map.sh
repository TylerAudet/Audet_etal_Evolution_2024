#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH -A def-idworkin
#SBATCH -n 32
#SBATCH --mem-per-cpu=10G
#SBATCH --mail-user=audett@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


# directory of processed sequences with BBduk
trim_dir=/home/audett/projects/def-idworkin/audett/SSD/trimmed

# variable for the reference genome
refGenome=/home/audett/projects/def-idworkin/audett/SSD/genomes/dmel-all-chromosome-r6.23.fasta.gz

# make output directory from mapping outputs
output=/home/audett/projects/def-idworkin/audett/SSD/mapped

#list all files to be read (this selects the left end from each PE pair)
files=(${trim_dir}/*_R1_trimmed.fastq)

#For loop over every file
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _R1_trimmed.fastq`

bwa mem -t 32 -M ${refGenome} \
${trim_dir}/${base}_R1_trimmed.fastq \
${trim_dir}/${base}_R2_trimmed.fastq \
> ${output}/${base}_bwa_PE.SAM
done

# -M Mark shorter split hits as secondary (for Picard compatibility).
# -t used to tell BWA to use 32 threads to speed things up

