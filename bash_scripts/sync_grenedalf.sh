


dir=/home/audett/scratch/SSD/repsMerged/find_top_sites

#make it arrays so you can index them
declare -a forward=(${dir}/*_sites.vcf)

# get the individual reads.
R1=${forward[${SLURM_ARRAY_TASK_ID}]}

name=${R1}
base=`basename ${name} _sites.vcf`

/home/audett/projects/def-idworkin/audett/SSD/scripts/grenedalf/bin/grenedalf sync \
--vcf-path ${dir}/${base}_sites.vcf \
--threads 16 \
--out-dir ${dir}/ \
--file-prefix ${base}
#--sample-name-list C1,C2,L1,L2,S1,S2
