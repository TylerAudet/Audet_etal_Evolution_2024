

module load StdEnv/2020
module load java/17.0.2

#dir=/home/audett/scratch/SSD/Analysis/EVC

java -Xmx8g -jar /home/audett/projects/def-idworkin/audett/SSD/scripts/snpEff/snpEff.jar \
-ud -download Drosophila_melanogaster /home/audett/scratch/SSD/Analysis/EVC/EvC_sig.vcf \
> /home/audett/scratch/SSD/Analysis/EVC/EvC_sig_annotated.txt
