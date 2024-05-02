

module load nixpkgs/16.09
module load perl/5.22.4

perl /home/audett/projects/def-idworkin/audett/SSD/scripts/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl \
--gtf /home/audett/projects/def-idworkin/audett/SSD/ref/dmel-all-chromosome-r6.23.fasta.out.gff \
--input /scratch/audett/SSD/Analysis/repsMerged/repsMerged.mpileup \
--output /scratch/audett/SSD/Analysis/repsMerged/repsMerged_repeatmasked.mpileup

