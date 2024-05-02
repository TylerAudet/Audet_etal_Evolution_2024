

/home/audett/projects/def-idworkin/audett/SSD/scripts/grenedalf/bin/grenedalf diversity \
--filter-sample-min-count 5 \
--measure theta-pi \
--window-type sliding \
--window-sliding-width 10000 \
--pool-sizes 200 \
--threads 16 \
--sync-path /home/audett/scratch/SSD/sexesMerged/sexesMerged_subsetted_sync.sync \
--out-dir /home/audett/scratch/SSD/sexesMerged \
--file-prefix sexesMerged_10000_pi \
--sample-name-list C1,C2,E1,E2,L1,L2,S1,S2
