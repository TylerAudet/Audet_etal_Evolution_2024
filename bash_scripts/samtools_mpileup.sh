

module load samtools/1.15

samtools mpileup -6 -q 20 \
-f /home/audett/projects/def-idworkin/audett/SSD/ref/dmel-all-chromosome-r6.23.fasta \
/home/audett/scratch/SSD/F100/repsMerged/*.bam \
-o /home/audett/scratch/SSD/F100/repsMerged/repsMerged_F100.mpileup
