#!/bin/bash
#SBATCH -t 6:00:00
#SBATCH -A def-idworkin
#SBATCH -n 32
#SBATCH --mem-per-cpu=4G
#SBATCH --mail-user=audett@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

raw_dir=/home/audett/projects/def-idworkin/audett/SSD/genomes

trim_dir=/home/audett/projects/def-idworkin/audett/SSD/trimmed

files=(${raw_dir}/*_R1.fastq)

for file in ${files[@]}
do
name=${file}
base=`basename ${name} _R1.fastq`

bbduk.sh \
in1=${raw_dir}/${base}_R1.fastq \
in2=${raw_dir}/${base}_R2.fastq \
out1=${trim_dir}/${base}_R1_trimmed.fastq \
out2=${trim_dir}/${base}_R2_trimmed.fastq \
ref=/home/audett/projects/def-idworkin/audett/SSD/trimmed/AllAdapters.fa \
threads=32 ftr=149 ktrim=r k=23 mink=6 hdist=1 tpe tbo \
qtrim=rl trimq=20 minlength=36 2> /home/audett/projects/def-idworkin/audett/SSD/trimmed/log/${base}.log

done
