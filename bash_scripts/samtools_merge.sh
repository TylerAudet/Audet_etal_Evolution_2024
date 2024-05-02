

module load StdEnv/2020
module load samtools/1.15.1

run1=/home/audett/scratch/SSD/F100/run1
run2=/home/audett/scratch/SSD/F100/run2

merged_dir=/home/audett/scratch/SSD/F100/repsMerged

declare -a forward=( ${run1}/*.bam )

R1=${forward[${SLURM_ARRAY_TASK_ID}]}

name=${R1}
base=`basename ${name} .bam`


samtools merge ${merged_dir}/${base}.bam \
${run1}/${base}.bam \
${run2}/${base}.bam
