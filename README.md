# overlap_detection
An overlap detection specific to detecting overlaps among PacBio long reads, which have high error rate but long sequencing length. The longer sequencing length is beneficial for detecting false alignment or interference caused by repeats. Here, this overlap detection tool is able to distinguish false alignments, which may contribute to particular downstream analyses that value the precision of alignments, such as variant calling, error detection, etc.

Install:
```
git clone git@github.com:jungle401/overlap_detection.git
cd REAL/build
make real -j
```

Usage:
```
/path/to/real\
  --reads_fasta input.fasta\
  --sequencing_depth 20\
  --output_dir ./output
```

The parameter `sequencing_depth` is set to manditory for evaluating false overlaps caused by interspersed repeats.

The output file would be in `.ovl` format: (specific to this project)
```
read1_id  read2_id  read1_start  read1_end  read2_start  read2_end
```