

/home/audett/projects/def-idworkin/audett/SSD/scripts/grenedalf/bin/grenedalf fst \
--window-type sliding \
--window-sliding-width 1 \
--method unbiased-nei \
--pool-sizes 400 \
--threads 16 \
--sync-path /home/audett/scratch/SSD/repsMerged/repsMerged_subsetted_sync.sync \
--sample-name-list C,E,L,S \
--omit-na-windows \
--out-dir /home/audett/scratch/SSD/Stewart/repsMerged \
--file-prefix Stewart_from_filteredSync_no_windows
