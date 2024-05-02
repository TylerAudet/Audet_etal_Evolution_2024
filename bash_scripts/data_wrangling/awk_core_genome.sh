

module load samtools/1.15

in=/home/audett/scratch/SSD/F100/mapped
out=/home/audett/scratch/SSD/F100/raw_bams

declare -a forward=( ${in}/*.sam )

R1=${forward[${SLURM_ARRAY_TASK_ID}]}

name=${R1}
base=`basename ${name} .sam`

samtools view -h -@ 32 ${in}/${base}.sam | \
awk '{if ($3 == "2L" || $3 == "2R" || $3 == "3L" || $3 == "3R" || $3 == "4" || $3 == "X" || \
$2 == "SN:2L" || $2 == "SN:2R" || $2 == "SN:3L" || $2 == "SN:3R" || $2 == "SN:4" || $2 == "SN:X") {print $0}}' | \
samtools view -b -@ 32 -o ${out}/${base}.bam -

