

module load python/2.7.18


python /home/audett/projects/def-idworkin/audett/SSD/scripts/SubsampleSync.py \
--sync /home/audett/scratch/SSD/Analysis/sexesMerged/syncs/2L.sync \
--target-cov 300 \
--min-cov 50 \
> /home/audett/scratch/SSD/Analysis/sexesMerged/syncs/subsample_2L.sync

python /home/audett/projects/def-idworkin/audett/SSD/scripts/SubsampleSync.py \
--sync /home/audett/scratch/SSD/Analysis/sexesMerged/syncs/2R.sync \
--target-cov 300 \
--min-cov 50 \
> /home/audett/scratch/SSD/Analysis/sexesMerged/syncs/subsample_2R.sync

python /home/audett/projects/def-idworkin/audett/SSD/scripts/SubsampleSync.py \
--sync /home/audett/scratch/SSD/Analysis/sexesMerged/syncs/3L.sync \
--target-cov 300 \
--min-cov 50 \
> /home/audett/scratch/SSD/Analysis/sexesMerged/syncs/subsample_3L.sync

python /home/audett/projects/def-idworkin/audett/SSD/scripts/SubsampleSync.py \
--sync /home/audett/scratch/SSD/Analysis/sexesMerged/syncs/3R.sync \
--target-cov 300 \
--min-cov 50 \
> /home/audett/scratch/SSD/Analysis/sexesMerged/syncs/subsample_3R.sync

python /home/audett/projects/def-idworkin/audett/SSD/scripts/SubsampleSync.py \
--sync /home/audett/scratch/SSD/Analysis/sexesMerged/syncs/4.sync \
--target-cov 300 \
--min-cov 50 \
> /home/audett/scratch/SSD/Analysis/sexesMerged/syncs/subsample_4.sync

python /home/audett/projects/def-idworkin/audett/SSD/scripts/SubsampleSync.py \
--sync /home/audett/scratch/SSD/Analysis/sexesMerged/syncs/X.sync \
--target-cov 300 \
--min-cov 50 \
> /home/audett/scratch/SSD/Analysis/sexesMerged/syncs/subsample_X.sync
