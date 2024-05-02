

module load python/2.7.18

python /home/audett/projects/def-idworkin/audett/SSD/scripts/PoolGen_var.py \
--input /home/audett/scratch/SSD/Analysis/sexesMerged/syncs/subsample_2L.sync \
--pool-size 200,200,200,200,200,200,200,200 \
--min-count 2 \
--window 5000 \
--step 5000 \
--sitecount /home/audett/scratch/SSD/Analysis/sexesMerged/pi/truewindows_2L-5000-5000.txt \
--min-sites-frac 0.75 \
--output window_5000_subsample_300_2L

python /home/audett/projects/def-idworkin/audett/SSD/scripts/PoolGen_var.py \
--input /home/audett/scratch/SSD/Analysis/sexesMerged/syncs/subsample_2R.sync \
--pool-size 200,200,200,200,200,200,200,200 \
--min-count 2 \
--window 5000 \
--step 5000 \
--sitecount /home/audett/scratch/SSD/Analysis/sexesMerged/pi/truewindows_2R-5000-5000.txt \
--min-sites-frac 0.75 \
--output window_5000_subsample_300_2R

python /home/audett/projects/def-idworkin/audett/SSD/scripts/PoolGen_var.py \
--input /home/audett/scratch/SSD/Analysis/sexesMerged/syncs/subsample_3L.sync \
--pool-size 200,200,200,200,200,200,200,200 \
--min-count 2 \
--window 5000 \
--step 5000 \
--sitecount /home/audett/scratch/SSD/Analysis/sexesMerged/pi/truewindows_3L-5000-5000.txt \
--min-sites-frac 0.75 \
--output window_5000_subsample_300_3L

python /home/audett/projects/def-idworkin/audett/SSD/scripts/PoolGen_var.py \
--input /home/audett/scratch/SSD/Analysis/sexesMerged/syncs/subsample_3R.sync \
--pool-size 200,200,200,200,200,200,200,200 \
--min-count 2 \
--window 5000 \
--step 5000 \
--sitecount /home/audett/scratch/SSD/Analysis/sexesMerged/pi/truewindows_3R-5000-5000.txt \
--min-sites-frac 0.75 \
--output window_5000_subsample_300_3R

python /home/audett/projects/def-idworkin/audett/SSD/scripts/PoolGen_var.py \
--input /home/audett/scratch/SSD/Analysis/sexesMerged/syncs/subsample_4.sync \
--pool-size 200,200,200,200,200,200,200,200 \
--min-count 2 \
--window 5000 \
--step 5000 \
--sitecount /home/audett/scratch/SSD/Analysis/sexesMerged/pi/truewindows_4-5000-5000.txt \
--min-sites-frac 0.75 \
--output window_5000_subsample_300_4

python /home/audett/projects/def-idworkin/audett/SSD/scripts/PoolGen_var.py \
--input /home/audett/scratch/SSD/Analysis/sexesMerged/syncs/subsample_X.sync \
--pool-size 200,200,200,200,200,200,200,200 \
--min-count 2 \
--window 5000 \
--step 5000 \
--sitecount /home/audett/scratch/SSD/Analysis/sexesMerged/pi/truewindows_X-5000-5000.txt \
--min-sites-frac 0.75 \
--output window_5000_subsample_300_X
