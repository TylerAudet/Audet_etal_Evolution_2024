

java -Xmx8g -jar /home/audett/projects/def-idworkin/audett/SSD/scripts/Gowinda-1.12.jar \
--annotation-file /home/audett/projects/def-idworkin/audett/SSD/ref/chr_dmel-all-r6.23.gtf \
--gene-set-file /home/audett/projects/def-idworkin/audett/SSD/ref/FBgn_funcassociate_go_associations.txt \
--snp-file /home/audett/scratch/SSD/sexesMerged/sexesMerged_clean.vcf \
--candidate-snp-file /home/audett/scratch/SSD/Analysis/LVS/Small_sites_CMH_CvSfst_LvSfst_pi_modelled.vcf \
--output-file /home/audett/scratch/SSD/Analysis/LVS/GO_output_Small.txt \
--mode gene \
--gene-definition updownstream5000 \
--simulations 1000000 \
--threads 16 \
--detailed-log
