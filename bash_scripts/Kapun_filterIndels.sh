

module load python/2.7.18

python2.7 /home/audett/projects/def-idworkin/audett/SSD/scripts/FilterPosFromVCF.py \
--indel /home/audett/scratch/SSD/Stewart/Stewart_repeatANDindelsmasked.indels \
--te /home/audett/projects/def-idworkin/audett/SSD/ref/dmel-all-chromosome-r6.23.fasta.out.gff \
--vcf /home/audett/scratch/SSD/Stewart/Stewart_blacklisted.vcf \
> Stewart_SNPs_clean.vcf
