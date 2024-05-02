

module load vcftools


vcftools --vcf /home/audett/scratch/SSD/Analysis/sexesMerged/SNPcalls/sexesMerged_clean.vcf --bed /home/audett/scratch/SSD/Analysis/regions_of_interest_C1.bed --out /home/audett/scratch/SSD/Analysis/allcriteria_C1 --recode --keep-INFO-all

vcftools --vcf /home/audett/scratch/SSD/Analysis/sexesMerged/SNPcalls/sexesMerged_clean.vcf --bed /home/audett/scratch/SSD/Analysis/regions_of_interest_C2.bed --out /home/audett/scratch/SSD/Analysis/allcriteria_C2 --recode --keep-INFO-all

vcftools --vcf /home/audett/scratch/SSD/Analysis/sexesMerged/SNPcalls/sexesMerged_clean.vcf --bed /home/audett/scratch/SSD/Analysis/regions_of_interest_E1.bed --out /home/audett/scratch/SSD/Analysis/allcriteria_E1 --recode --keep-INFO-all

vcftools --vcf /home/audett/scratch/SSD/Analysis/sexesMerged/SNPcalls/sexesMerged_clean.vcf --bed /home/audett/scratch/SSD/Analysis/regions_of_interest_E2.bed --out /home/audett/scratch/SSD/Analysis/allcriteria_E2 --recode --keep-INFO-all

vcftools --vcf /home/audett/scratch/SSD/Analysis/sexesMerged/SNPcalls/sexesMerged_clean.vcf --bed /home/audett/scratch/SSD/Analysis/regions_of_interest_L1.bed --out /home/audett/scratch/SSD/Analysis/allcriteria_L1 --recode --keep-INFO-all

vcftools --vcf /home/audett/scratch/SSD/Analysis/sexesMerged/SNPcalls/sexesMerged_clean.vcf --bed /home/audett/scratch/SSD/Analysis/regions_of_interest_L2.bed --out /home/audett/scratch/SSD/Analysis/allcriteria_L2 --recode --keep-INFO-all

vcftools --vcf /home/audett/scratch/SSD/Analysis/sexesMerged/SNPcalls/sexesMerged_clean.vcf --bed /home/audett/scratch/SSD/Analysis/regions_of_interest_S1.bed --out /home/audett/scratch/SSD/Analysis/allcriteria_S1 --recode --keep-INFO-all

vcftools --vcf /home/audett/scratch/SSD/Analysis/sexesMerged/SNPcalls/sexesMerged_clean.vcf --bed /home/audett/scratch/SSD/Analysis/regions_of_interest_S2.bed --out /home/audett/scratch/SSD/Analysis/allcriteria_S2 --recode --keep-INFO-all



grep -v 'NO_VARIATION\|intron_variant'  /home/audett/scratch/SSD/Analysis/C1_annotated.txt > /home/audett/scratch/SSD/Analysis/C1_interesting.txt
grep -v 'NO_VARIATION\|intron_variant'  /home/audett/scratch/SSD/Analysis/C2_annotated.txt > /home/audett/scratch/SSD/Analysis/C2_interesting.txt
grep -v 'NO_VARIATION\|intron_variant'  /home/audett/scratch/SSD/Analysis/E1_annotated.txt > /home/audett/scratch/SSD/Analysis/E1_interesting.txt
grep -v 'NO_VARIATION\|intron_variant'  /home/audett/scratch/SSD/Analysis/E2_annotated.txt > /home/audett/scratch/SSD/Analysis/E2_interesting.txt
grep -v 'NO_VARIATION\|intron_variant'  /home/audett/scratch/SSD/Analysis/L1_annotated.txt > /home/audett/scratch/SSD/Analysis/L1_interesting.txt
grep -v 'NO_VARIATION\|intron_variant'  /home/audett/scratch/SSD/Analysis/L2_annotated.txt > /home/audett/scratch/SSD/Analysis/L2_interesting.txt
grep -v 'NO_VARIATION\|intron_variant'  /home/audett/scratch/SSD/Analysis/S1_annotated.txt > /home/audett/scratch/SSD/Analysis/S1_interesting.txt
grep -v 'NO_VARIATION\|intron_variant'  /home/audett/scratch/SSD/Analysis/S2_annotated.txt > /home/audett/scratch/SSD/Analysis/S2_interesting.txt


