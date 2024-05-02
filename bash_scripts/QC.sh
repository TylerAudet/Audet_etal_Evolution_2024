

module load fastqc/0.11.9
parallel --will-cite --jobs $SLURM_CPUS_PER_TASK fastqc -t 1  \
-o /home/audett/projects/def-idworkin/audett/SSD/genomes/F100/QC {} \
::: /home/audett/projects/def-idworkin/audett/SSD/genomes/F100/*.txt
