

module load samtools/1.15

in=/home/audett/scratch/SSD/F100/fixmated
out=/home/audett/scratch/SSD/F100/coord_sorted

declare -a forward=( ${in}/*.bam )

R1=${forward[${SLURM_ARRAY_TASK_ID}]}

name=${R1}
base=`basename ${name} .bam`

samtools sort -@ 32 -o ${out}/${base}_coordsorted.bam ${in}/${base}.bam

