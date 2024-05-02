

dir=/home/audett/scratch/SSD/Stewart/fst

declare -a forward=( ${dir}/*.fst )

R1=${forward[${SLURM_ARRAY_TASK_ID}]}

name=${R1}
base=`basename ${name} .fst`

awk '{if ($3 >= 0.05) {print $1,$2-1,$2+1}}' ${dir}/${base}.fst | sed "1d" \
| sed 1i"chrom\tchromStart\tchromEnd" | \
sed 's/ /\t/g' \
> ${dir}/high_${base}.bed
