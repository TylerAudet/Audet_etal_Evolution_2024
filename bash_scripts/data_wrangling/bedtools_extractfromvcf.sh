

module load bedtools


bed_dir=/home/audett/scratch/SSD/repsMerged/find_top_sites

#make it arrays so you can index them
declare -a forward=(${bed_dir}/*.bed)

# get the individual reads.
R1=${forward[${SLURM_ARRAY_TASK_ID}]}

name=${R1}
base=`basename ${name} .bed`

bedtools intersect -header -a /home/audett/scratch/SSD/repsMerged/repsMerged_blacklisted.vcf \
-b ${bed_dir}/${base}.bed > ${bed_dir}/${base}_sites.vcf



