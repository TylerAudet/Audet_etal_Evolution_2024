#!/bin/bash
#SBATCH -t 3-00:00:00
#SBATCH -A def-idworkin
#SBATCH --mem=32G
#SBATCH --mail-user=audett@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module load StdEnv/2020
module load varscan/2.4.2

java -Xmx32g -jar \
$EBROOTVARSCAN/VarScan.v2.4.2.jar \
mpileup2cns \
/home/audett/projects/def-idworkin/audett/SSD/Sexes_combined_norepeats.mpileup \
--min-reads2 5 \
--min-coverage 50
--p-value 0.1 \
--min-var-freq 0.01 \
--min-freq-for-hom 1 \
--min-avg-qual 20 \
--variants \
--output-vcf 1 \
| bgzip > /home/audett/projects/def-idworkin/audett/SSD/Sexes_combined_variants.vcf.bgz
