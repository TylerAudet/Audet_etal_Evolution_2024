

module load samtools/1.15

in=/scratch/audett/SSD/Stewart/trash_can/realigned
out=/scratch/audett/SSD/Stewart/indel_test

declare -a forward=( ${in}/*.bam )

R1=${forward[${SLURM_ARRAY_TASK_ID}]}

name=${R1}
base=`basename ${name} .bam`

samtools view -@ 32 -b ${in}/${base}.bam 2L > ${out}/${base}_2L.bam
