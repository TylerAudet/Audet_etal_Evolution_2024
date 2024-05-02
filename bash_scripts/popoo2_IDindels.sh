

module load nixpkgs/16.09
module load perl/5.22.4

perl /home/audett/projects/def-idworkin/audett/SSD/scripts/popoolation2_1201/indel_filtering/identify-indel-regions.pl \
--input /scratch/audett/SSD/Analysis/repsMerged/repsMerged_repeatmasked.mpileup \
--output /home/audett/projects/def-idworkin/audett/SSD/ref/repsMerged.gtf \
--indel-window 10

