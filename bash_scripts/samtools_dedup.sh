

module load samtools/1.15

in=/home/audett/scratch/SSD/F100/coord_sorted
out=/home/audett/scratch/SSD/F100/dedup

declare -a forward=( ${in}/*_namesorted_fixmated_coordsorted.bam )

R1=${forward[${SLURM_ARRAY_TASK_ID}]}

name=${R1}
base=`basename ${name} _namesorted_fixmated_coordsorted.bam`

samtools markdup -l 75 -r -s -f stats.txt -d 100 -@ 32 ${in}/${base}_namesorted_fixmated_coordsorted.bam ${out}/${base}.bam

