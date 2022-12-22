./app \
--reads_fasta /mnt/es/ness/johnson/thesis/data/Human/reads/processed/m13101_67/chr1/aln.lsmq10.fasta\
--fmi_src build \
--batchSize_readsNum 230000 \
--fmi_dir ../midFiles/human/m13101_67/ \
--kmer_size 16 \
--dist_hori_filter 0 \
--bin_width 400 \
--thres_one_bin_least_score 40 \
--thres_bin_score 80 \
--maxNumSeedAnchors 50 \
--numBin_sliding 3 \
--thrsDcnt_toGraphBin 10 \
--min_antiDiag_space 40 \
--output_dir ../output/human/m13101_067/
~/thesis/evaluation/src/exe_shs/minOlen_2000/human/m13101_67/chr1/mime.sh
