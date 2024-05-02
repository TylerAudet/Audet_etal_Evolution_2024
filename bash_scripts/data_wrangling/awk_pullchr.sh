


awk '{if ($1 == "2L") {print $0}}' /scratch/audett/SSD/Analysis/repsMerged/repsMerged_indelsmasked.mpileup > \
/scratch/audett/SSD/Analysis/repsMerged/2L.mpileup
awk '{if ($1 == "2R") {print $0}}' /scratch/audett/SSD/Analysis/repsMerged/repsMerged_indelsmasked.mpileup > \
/scratch/audett/SSD/Analysis/repsMerged/2R.mpileup
awk '{if ($1 == "3L") {print $0}}' /scratch/audett/SSD/Analysis/repsMerged/repsMerged_indelsmasked.mpileup > \
/scratch/audett/SSD/Analysis/repsMerged/3L.mpileup
awk '{if ($1 == "3R") {print $0}}' /scratch/audett/SSD/Analysis/repsMerged/repsMerged_indelsmasked.mpileup > \
/scratch/audett/SSD/Analysis/repsMerged/3R.mpileup
awk '{if ($1 == "X") {print $0}}' /scratch/audett/SSD/Analysis/repsMerged/repsMerged_indelsmasked.mpileup > \
/scratch/audett/SSD/Analysis/repsMerged/X.mpileup
awk '{if ($1 == "4") {print $0}}' /scratch/audett/SSD/Analysis/repsMerged/repsMerged_indelsmasked.mpileup > \
/scratch/audett/SSD/Analysis/repsMerged/4.mpileup
