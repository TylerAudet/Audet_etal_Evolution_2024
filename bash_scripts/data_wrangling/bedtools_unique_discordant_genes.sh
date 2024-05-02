

module load bedtools/2.30.0

dir=/home/audett/scratch/SSD/Stewart/fst

bedtools intersect -wa -a ${dir}/sex_E1.bed -b ${dir}/sex_E2.bed | \
bedtools intersect -v -a stdin \
-b ${dir}/sex_C1.bed \
${dir}/sex_C2.bed \
${dir}/sex_L1.bed \
${dir}/sex_L2.bed \
${dir}/sex_S1.bed \
${dir}/sex_S2.bed \
| sed 1i"chrom\tchromStart\tchromEnd" > ${dir}/sexy_genes.bed
