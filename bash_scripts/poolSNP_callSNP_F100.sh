


# make the common path available in variable (this is optional, and just for the looks
dir=/home/audett/scratch/SSD/F100/SNP_calling/

# make it arrays so you can index them
declare -a forward=( ${dir}*.mpileup )

# get the individual reads.
R1=${forward[${SLURM_ARRAY_TASK_ID}]}

name=${R1}
base=`basename ${name} .mpileup`

bash /home/audett/projects/def-idworkin/audett/SSD/scripts/PoolSNP-master/PoolSNP.sh \
mpileup=${dir}${base}.mpileup \
reference=/home/audett/projects/def-idworkin/audett/SSD/ref/dmel-all-chromosome-r6.23.fasta \
names=C1,C2,L1,L2,S1,S2 \
max-cov=0.98 \
min-cov=10 \
min-count=2 \
min-freq=0.01 \
miss-frac=0.2 \
jobs=32 \
BS=1 \
output=/home/audett/scratch/SSD/F100/SNP_calling/${base}
