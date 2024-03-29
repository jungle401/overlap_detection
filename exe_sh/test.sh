./app \
--reads_fasta /mnt/es/ness/johnson/thesis/data/Ecoli/reads/aligned_and_sampled/minimap2/DevNet-P6C4/shrinked/head_5000.fasta \
--fmi_src build \
--fmi_dir /mnt/es/ness/johnson/thesis/overlaps/overlap_ver4/batchReads/midFiles/toy \
--kmer_size 15 \
--dist_hori_filter 5 \
--bin_width 500 \
--thres_one_bin_least_score 30 \
--thres_bin_score 80 \
--numBin_sliding 3 \
--thrsDcnt_toGraphBin 10 \
--min_antiDiag_space 40 \
--output_dir ../output/toy
