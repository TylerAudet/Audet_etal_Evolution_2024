

module load python/2.7.18

python /home/audett/projects/def-idworkin/audett/SSD/scripts/TrueWindows.py \
--badcov /home/audett/scratch/SSD/Analysis/sexesMerged/SNPcalls/2L_BS.txt \
--window 5000 \
--step 5000 \
--chromosomes 2L:23513712 \
--output /home/audett/scratch/SSD/Analysis/sexesMerged/pi/truewindows_2L

python /home/audett/projects/def-idworkin/audett/SSD/scripts/TrueWindows.py \
--badcov /home/audett/scratch/SSD/Analysis/sexesMerged/SNPcalls/2R_BS.txt \
--window 5000 \
--step 5000 \
--chromosomes 2R:25286936 \
--output /home/audett/scratch/SSD/Analysis/sexesMerged/pi/truewindows_2R


python /home/audett/projects/def-idworkin/audett/SSD/scripts/TrueWindows.py \
--badcov /home/audett/scratch/SSD/Analysis/sexesMerged/SNPcalls/3L_BS.txt \
--window 5000 \
--step 5000 \
--chromosomes 3L:28110227 \
--output /home/audett/scratch/SSD/Analysis/sexesMerged/pi/truewindows_3L


python /home/audett/projects/def-idworkin/audett/SSD/scripts/TrueWindows.py \
--badcov /home/audett/scratch/SSD/Analysis/sexesMerged/SNPcalls/3R_BS.txt \
--window 5000 \
--step 5000 \
--chromosomes 3R:32079331 \
--output /home/audett/scratch/SSD/Analysis/sexesMerged/pi/truewindows_3R

python /home/audett/projects/def-idworkin/audett/SSD/scripts/TrueWindows.py \
--badcov /home/audett/scratch/SSD/Analysis/sexesMerged/SNPcalls/4_BS.txt \
--window 5000 \
--step 5000 \
--chromosomes 4:1348131 \
--output /home/audett/scratch/SSD/Analysis/sexesMerged/pi/truewindows_4


python /home/audett/projects/def-idworkin/audett/SSD/scripts/TrueWindows.py \
--badcov /home/audett/scratch/SSD/Analysis/sexesMerged/SNPcalls/X_BS.txt \
--window 5000 \
--step 5000 \
--chromosomes X:23542271 \
--output /home/audett/scratch/SSD/Analysis/sexesMerged/pi/truewindows_X
