

module load vcftools


bed_dir=/home/audett/scratch/SSD/Analysis/cutoff_10

# make it arrays so you can index them
declare -a forward=(${bed_dir}/*.bed)

# get the individual reads.
R1=${forward[${SLURM_ARRAY_TASK_ID}]}

name=${R1}
base=`basename ${name} .bed`

vcftools --vcf /home/audett/scratch/SSD/Analysis/repsMerged/repsMerged_blacklisted.vcf --bed ${bed_dir}/${base}.bed --out ${bed_dir}/${base} --recode --keep-INFO-all

