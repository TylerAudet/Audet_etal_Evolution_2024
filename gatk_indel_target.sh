#!/bin/bash
#SBATCH -t 3:00:00
#SBATCH -A def-idworkin
#SBATCH --array=0-7
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --mail-user=audett@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module load nixpkgs/16.09
module load gatk/3.8

#Path to input directory
in_dir=/home/audett/projects/def-idworkin/audett/SSD/dedup

#Path to output directory
out_dir=/home/audett/projects/def-idworkin/audett/SSD/gatk


declare -a forward=( /home/audett/projects/def-idworkin/audett/SSD/dedup/*_RG_dedup.bam )

R1=${forward[${SLURM_ARRAY_TASK_ID}]}

name=${R1}
base=`basename ${name} _RG_dedup.bam`

java -Xmx32g -jar "${EBROOTGATK}"/GenomeAnalysisTK.jar -I ${in_dir}/${base}_RG_dedup.bam \
-R /home/audett/projects/def-idworkin/audett/SSD/ref/dmel-all-chromosome-r6.23.fasta \
-T RealignerTargetCreator \
-o ${out_dir}/${base}.intervals
