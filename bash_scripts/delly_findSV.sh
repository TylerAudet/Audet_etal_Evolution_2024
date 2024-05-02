

module load StdEnv/2023
module load gcc/12.3
module load delly/1.1.8

declare -a forward=( /home/audett/scratch/SSD/Stewart/realigned/*.bam )

R1=${forward[${SLURM_ARRAY_TASK_ID}]}

name=${R1}
base=`basename ${name} .bam`

delly call -g /home/audett/projects/def-idworkin/audett/SSD/ref/dmel-all-chromosome-r6.23.fasta \
/home/audett/scratch/SSD/Stewart/realigned/${base}.bam \
> /home/audett/scratch/SSD/Stewart/delly_output/${base}_delly.vcf


