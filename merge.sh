#!/bin/bash
#SBATCH -t 3:00:00
#SBATCH -A def-idworkin
#SBATCH --array=0-28
#SBATCH --cpus-per-task=32
#SBATCH --mem=10G
#SBATCH --mail-user=audett@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module load samtools/1.12

run1=/home/audett/projects/def-idworkin/audett/SSD/run1
run2=/home/audett/projects/def-idworkin/audett/SSD/run2

merged_dir=/home/audett/projects/def-idworkin/audett/SSD/Sexes_merged

declare -a forward=( ${run1}/*.bam )

R1=${forward[${SLURM_ARRAY_TASK_ID}]}

name=${R1}
base=`basename ${name} .bam`


samtools merge ${merged_dir}/${base}.bam \
${run1}/${base}.bam \
${run2}/${base}.bam
