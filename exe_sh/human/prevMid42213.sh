./app \
--reads_fasta /mnt/es/ness/johnson/thesis/data/Human/reads/processed/prev_mid_42213/chr1/aln.lsmq10.fasta \
--fmi_src load \
--batchSize_readsNum 230000 \
--fmi_dir ../midFiles/human/prev_mid_42213/ \
--kmer_size 16 \
--dist_hori_filter 800 \
--bin_width 400 \
--thres_one_bin_least_score 40 \
--thres_bin_score 80 \
--maxNumSeedAnchors 1000 \
--numBin_sliding 3 \
--thrsDcnt_toGraphBin 10 \
--min_antiDiag_space 40 \
--output_dir ../output/human/prev_mid_42213/
# ~/thesis/evaluation/src/exe_shs/minOlen_2000/human/prev_mid_42213/mime.sh
