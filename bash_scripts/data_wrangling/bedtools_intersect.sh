

module load bedtools/2.30.0

bed_dir=/home/audett/scratch/SSD/Analysis/EVC/

#make it arrays so you can index them
#files=(${bed_dir}/*_2L.bam)

# get the individual reads.


#for file in ${files[@]}
#do
#name=${file}
#base=`basename ${name} _2L.bam`


#echo "directory"
#echo ${bed_dir}
#echo "file"
#echo ${file}
#echo "name"
#echo ${name}


bedtools intersect -sorted -header -a /home/audett/scratch/SSD/Analysis/EVC/EvC_regions_of_interest_withHeader.vcf \
-b /home/audett/scratch/SSD/Analysis/E1/high_E1_3L.bed \
> /home/audett/scratch/SSD/Analysis/E1/overlap_with_EvC.vcf

#done


