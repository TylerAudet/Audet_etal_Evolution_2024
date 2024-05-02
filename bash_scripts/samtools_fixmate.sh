

module load samtools/1.15

in=/home/audett/scratch/SSD/F100/name_sorted_array
out=/home/audett/scratch/SSD/F100/fixmated

declare -a forward=( ${in}/*.bam )

R1=${forward[${SLURM_ARRAY_TASK_ID}]}

name=${R1}
base=`basename ${name} .bam`

samtools fixmate -m -u -@ 32 -O bam ${in}/${base}.bam ${out}/${base}_fixmated.bam
