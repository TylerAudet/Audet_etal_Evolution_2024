
#!/bin/bash
#SBATCH -t 00-01:00:00
#SBATCH -A def-idworkin
#SBATCH --cpus-per-task=32
#SBATCH --mem=0 # reserve the whole node
#SBATCH --mail-user=audett@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module load fastqc/0.11.9
parallel --will-cite --jobs $SLURM_CPUS_PER_TASK fastqc -t 1  \
-o /home/audett/projects/def-idworkin/audett/SSD/QC_after_trim {} \
::: /home/audett/projects/def-idworkin/audett/SSD/trimmed/*fastq
