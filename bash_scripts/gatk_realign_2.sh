

module load nixpkgs/16.09
module load gatk/3.8

#Path to input directory
in_dir=/home/audett/scratch/SSD/Analysis/addRG

#Path to output directory
out_dir=/home/audett/scratch/SSD/Analysis/realigned

declare -a forward=( /home/audett/scratch/SSD/Analysis/addRG/*_RG.bam )

R1=${forward[${SLURM_ARRAY_TASK_ID}]}

name=${R1}
base=`basename ${name} _RG.bam`

java -Xmx10g -jar ${EBROOTGATK}/GenomeAnalysisTK.jar -I ${in_dir}/${base}_RG.bam \
-R /home/audett/projects/def-idworkin/audett/SSD/ref/dmel-all-chromosome-r6.23.fasta \
-T IndelRealigner -targetIntervals /home/audett/scratch/SSD/Analysis/indel_targets/${base}.intervals \
-o ${out_dir}/${base}.bam
