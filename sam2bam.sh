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

in=/home/audett/projects/def-idworkin/audett/SSD/mapped

declare -a forward=( ${in}/*.sam )

R1=${forward[${SLURM_ARRAY_TASK_ID}]}

name=${R1}
base=`basename ${name} .sam`

samtools view -b -@32 ${in}/${base}.sam | samtools sort -o ${in}/${base}.bam
