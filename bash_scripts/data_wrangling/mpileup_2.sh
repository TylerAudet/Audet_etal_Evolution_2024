

module load samtools/1.15

samtools mpileup -Q 20 -q 20 -d 800 \
-f /home/audett/projects/def-idworkin/audett/SSD/ref/dmel-all-chromosome-r6.23.fasta \
/scratch/audett/SSD/null/realigned/*.bam \
-o /scratch/audett/SSD/null/null.mpileup
