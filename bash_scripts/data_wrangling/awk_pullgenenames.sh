

dir=/home/audett/scratch/SSD/Analysis/EVC

# make it arrays so you can index them
#declare -a forward=(${dir}/*.txt)

# get the individual reads.
#R1=${forward[${SLURM_ARRAY_TASK_ID}]}

#name=${R1}
#base=`basename ${name} .txt`

awk '{print $8}' ${dir}/final_EVCgenes_updatedCMH.txt | \
awk -F '|' '{print $4, $5}' | \
sort | uniq > ${dir}/final_EVCgenes_updatedCMH_names.txt

