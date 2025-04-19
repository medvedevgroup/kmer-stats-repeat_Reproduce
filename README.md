# kmer-stats-repeat_Reproduce

Reproduce the figures in the paper.

## Dependencies

	1. Snakemake
	1. Matplotlib

## kmer_info

1. Run following command, to generate the statistics information of sequences in the paper

   ```
   snakemake --core [N]
   ```

   

2. The results are in corresponding folders, like `./chr21_centromere`

3. For occurrence kmer count, we list the number of occurrence $i$ and the counts of kmers with occurence $i$ in `occ_counts.txt`.  The contents are shown in the form "#occurrence : #kmer counts". We also draw a plot `kmer_distribution_histogram.png` to show this information. 

4. For overlapped kmers, we list the value of $(occ-sep)$ and the counts of kmers with this difference in `diff_counts.txt`. The contents are shown in the form "oct-sep : #kmer counts". (The definitions of $occ$ and $sep$ for a kmer can be found in the paper).

