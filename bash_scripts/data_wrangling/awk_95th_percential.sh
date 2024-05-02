

dir=/home/audett/scratch/SSD/Stewart/fst

declare -a forward=( ${dir}/*.fst )

R1=${forward[${SLURM_ARRAY_TASK_ID}]}

name=${R1}
base=`basename ${name} .fst`

sort -n ${dir}/${base}.fst  | awk '{all[NR] = $3} END{print all[int(NR*0.99)]}'\
 > ${dir}/${base}_top1_fst.txt
