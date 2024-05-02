

module load StdEnv/2020
module load bedtools/2.30.0

bedtools intersect -header -v -a /home/audett/scratch/SSD/Analysis/repsMerged/repsMerged.vcf \
-b /home/audett/projects/def-idworkin/audett/SSD/ref/blacklist/dm6-blacklist.v2.bed \
> /home/audett/scratch/SSD/Analysis/repsMerged/repsMerged_blacklisted.vcf
