

module load samtools
parallel --willcite --results {/.}.bam samtools view -b -q 20 -@ 16 \
::: /home/audett/projects/def-idworkin/audett/SSD/Sexes_merged/*.bam
