

# make the common path available in variable (this is optional, and just for the looks
dir=/home/audett/scratch/SSD/Analysis/repsMerged/
out=/home/audett/scratch/SSD/Analysis/repsMerged/

# make it arrays so you can index them
declare -a forward=( /home/audett/scratch/SSD/Analysis/repsMerged/*.mpileup )

# get the individual reads.
R1=${forward[${SLURM_ARRAY_TASK_ID}]}

name=${R1}
base=`basename ${name} .mpileup`

bash /home/audett/projects/def-idworkin/audett/SSD/scripts/PoolSNP-master/PoolSNP.sh \
mpileup=${dir}${base}.mpileup \
reference=/home/audett/projects/def-idworkin/audett/SSD/ref/dmel-all-chromosome-r6.23.fasta \
names=C,E,L,S \
max-cov=0.98 \
min-cov=120 \
min-count=40 \
min-freq=0.01 \
miss-frac=0.2 \
jobs=32 \
BS=1 \
output=${dir}${base}
