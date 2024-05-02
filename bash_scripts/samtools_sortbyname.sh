

module load samtools/1.15

in=/home/audett/scratch/SSD/F100/raw_bams
out=/home/audett/scratch/SSD/F100/name_sorted_array
declare -a forward=( ${in}/*.bam )

R1=${forward[${SLURM_ARRAY_TASK_ID}]}

name=${R1}
base=`basename ${name} .bam`

samtools sort -n -@ 32 -o ${out}/${base}_namesorted.bam ${in}/${base}.bam
