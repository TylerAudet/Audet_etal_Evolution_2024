

module load nixpkgs/16.09
module load gatk/3.8

#Path to input directory
in_dir=/home/audett/scratch/SSD/F100/addRG

#Path to output directory
out_dir=/home/audett/scratch/SSD/F100/realigned

declare -a forward=( ${in_dir}/*_RG.bam )

R1=${forward[${SLURM_ARRAY_TASK_ID}]}

name=${R1}
base=`basename ${name} _RG.bam`

java -Xmx10g -jar ${EBROOTGATK}/GenomeAnalysisTK.jar --fix_misencoded_quality_scores \
-I ${in_dir}/${base}_RG.bam \
-R /home/audett/projects/def-idworkin/audett/SSD/ref/dmel-all-chromosome-r6.23.fasta \
-T IndelRealigner -targetIntervals /home/audett/scratch/SSD/F100/indel_targets/${base}.intervals \
-o ${out_dir}/${base}.bam
