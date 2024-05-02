

module load samtools
parallel --willcite samtools index -@ 8 {} \
::: /home/audett/scratch/SSD/F100/raw_bams/*.bam
