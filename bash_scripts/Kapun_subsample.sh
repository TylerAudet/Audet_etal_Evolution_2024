

module load python/2.7.18

python2.7 /home/audett/projects/def-idworkin/audett/SSD/scripts/SubsampleSync.py \
--sync /home/audett/scratch/SSD/Stewart/Stewart_good.sync \
--target-cov 150 \
--min-cov 100 \
> /home/audett/scratch/SSD/Stewart/subsample150_Stewart_good.sync
