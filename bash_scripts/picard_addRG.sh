

module load picard/2.26.3

in=/home/audett/scratch/SSD/F100/dedup
out=/home/audett/scratch/SSD/F100/addRG

declare -a forward=( ${in}/*.bam )

R1=${forward[${SLURM_ARRAY_TASK_ID}]}

name=${R1}
base=`basename ${name} .bam`

java -jar -Xmx10g $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
INPUT=${in}/${base}.bam \
OUTPUT=${out}/${base}_RG.bam \
SORT_ORDER=coordinate \
RGID=library \
RGLB=library \
RGPL=illumina \
RGSM=Stewart \
RGPU=library \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=SILENT
