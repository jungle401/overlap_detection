./app \
--reads_fasta /mnt/es/ness/johnson/thesis/data/Scerevisiae/reads/processed/chrAll/aln.lsmq10.fasta \
--fmi_src build \
--batchSize_readsNum 230000 \
--fmi_dir ../midFiles/scere \
--kmer_size 16 \
--dist_hori_filter 16 \
--bin_width 600 \
--thres_one_bin_least_score 40 \
--thres_bin_score 80 \
--maxNumSeedAnchors 200 \
--numBin_sliding 3 \
--thrsDcnt_toGraphBin 10 \
--min_antiDiag_space 40 \
--output_dir ../output/scere
~/thesis/evaluation/src/exe_shs/minOlen_2000/scere/chrAll/mime.sh
