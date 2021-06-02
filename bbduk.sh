#!/bin/bash
#SBATCH -t 4:00:00
#SBATCH -A def-idworkin
#SBATCH --array=0-57
#SBATCH --cpus-per-task=32
#SBATCH --mem=10G
#SBATCH --mail-user=audett@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module load bbmap/38.86

in=/home/audett/projects/def-idworkin/audett/SSD/genomes
out=/home/audett/projects/def-idworkin/audett/SSD/trimmed

declare -a forward=( ${in}/*_R1.fastq )

R1=${forward[${SLURM_ARRAY_TASK_ID}]}

name=${R1}
base=`basename ${name} _R1.fastq`

bbduk.sh \
in1=${in}/${base}_R1.fastq \
in2=${in}/${base}_R2.fastq \
out1=${out}/${base}_R1_trimmed.fastq \
out2=${out}/${base}_R2_trimmed.fastq \
ref=/home/audett/projects/def-idworkin/audett/SSD/trimmed
threads=32 ftr=149 ktrim=r k=23 mink=10 hdist=1 tpe tbo \
qtrim=rl trimq=20 minlength=36 2> /home/audett/projects/def-idworkin/audett/SSD/trimmed/log

done
