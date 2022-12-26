./app \
--reads_fasta /mnt/es/ness/johnson/thesis/data/Ecoli/reads/aligned_and_sampled/minimap2/DevNet-P6C4/rngAll/aln.lsmq10.fasta \
--fmi_src build \
--fmi_dir ../mid_files/devnet \
--kmer_size 16 \
--dist_hori_filter 800 \
--bin_width 400 \
--thres_one_bin_least_score 80 \
--maxNumSeedAnchors 100 \
--thres_bin_score 40 \
--numBin_sliding 3 \
--thrsDcnt_toGraphBin 10 \
--min_antiDiag_space 80 \
--batchSize_readsNum 70000 \
--output_dir ../output/devnet
# ~/thesis/evaluation/src/exe_shs/minOlen_2000/ecoli/rngAll/ecoli-devNet_mime.sh
