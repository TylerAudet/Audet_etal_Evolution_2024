

module load samtools/1.15

in=/home/audett/scratch/SSD/Stewart

declare -a forward=( ${in}/*.sam )

R1=${forward[${SLURM_ARRAY_TASK_ID}]}

name=${R1}
base=`basename ${name} .sam`

samtools view -b -@32 ${in}/${base}.sam | samtools sort ${in}/${base}.bam | samtools view -b -q 20 -@ 32 -o ${in}/${base}.bam
