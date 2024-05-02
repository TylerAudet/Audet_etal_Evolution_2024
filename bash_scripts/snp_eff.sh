

module load snpeff/5.0

dir=/home/audett/scratch/SSD/Analysis/E1

# make it arrays so you can index them
#declare -a forward=(${dir}/*.vcf)

# get the individual reads.
#R1=${forward[${SLURM_ARRAY_TASK_ID}]}

#name=${R1}
#base=`basename ${name} .vcf`

java -Xmx8g -jar /home/audett/projects/def-idworkin/audett/SSD/scripts/snpEff/snpEff.jar \
-ud -download Drosophila_melanogaster ${dir}/high_E1_3L.vcf > ${dir}/high_E1_3L.txt
