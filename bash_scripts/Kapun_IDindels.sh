

module load python/2.7.18

python /home/audett/projects/def-idworkin/audett/SSD/scripts/DetectIndels.py \
--mpileup /home/audett/scratch/SSD/Stewart/Stewart_repeatmasked/Stewart_repeatmasked.mpileup \
--minimum-count 1 \
--mask 10 \
> Stewart_repeatANDindelsmasked.mpileup
