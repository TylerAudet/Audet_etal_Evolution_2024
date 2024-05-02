

module load samtools/1.15

in=/scratch/audett/SSD/dedup2/

declare -a forward=( ${in}*.bam )

R1=${forward[${SLURM_ARRAY_TASK_ID}]}

name=${R1}
base=`basename ${name} .bam`

samtools quickcheck -v ${in}${base}.bam > /scratch/audett/SSD/dedup2/check.txt
