#! /bin/bash
#!/bin/bash
#SBATCH -t 3:00:00
#SBATCH -A def-idworkin
#SBATCH --array=0-7
#SBATCH --mem=10G
#SBATCH --mail-user=audett@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


module load picard/2.23.3

#Variable for project:
project_dir=/home/audett/projects/def-idworkin/audett/SSD/Q20

declare -a forward=( /home/audett/projects/def-idworkin/audett/SSD/Q20/*.bam )

# get the individual reads.
R1=${forward[${SLURM_ARRAY_TASK_ID}]}


name=${R1}
base=`basename ${name} .bam`


java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=${project_dir}/${base}.bam \
  O=${project_dir}/${base}_RG.bam \
  RGID=1_2 \
  RGLB=library1 \
  RGPL=illumina \
  RGPU=None \
  RGSM=${base}
