# kmer-stats-repeat_Reproduce

Reproduce the figures in the paper.

## Dependencies

	1. Snakemake
	1. Matplotlib

## kmer_info

1. Use the following command to compile

   ```
   g++ -std=c++17 kmer_info.cpp -o kmer_info
   ```

   

2. Run following command, to generate the statistics information of sequences in the paper

   ```
   snakemake --core [N]
   ```

3. The results are in corresponding folders, like `./chr21_centromere`

4. For occurrence kmer count, we list the number of occurrence $i$ and the counts of kmers with occurence $i$ in `occ_counts.txt`.  The contents are shown in the form "#occurrence : #kmer counts". We also draw a plot `kmer_distribution_histogram.png` to show this information. 

5. For overlapped kmers, we list the value of $(occ-sep)$ and the counts of kmers with this difference in `diff_counts.txt`. The contents are shown in the form "oct-sep : #kmer counts". (The definitions of $occ$ and $sep$ for a kmer can be found in the paper).

## Comparasion_of_three_estimators

1. `cd ./Comparasion_of_three_estimators` and run `make` in current folder to compile.

2. Users can modify settings in config.yaml. Here is an example

   The range of $k,r$ a hard code in snakefile.

   ```config.yaml
   analyses:
     - name: "chr6_RepeatMasker" # Experiment name
       output_dir: "chr6_RepeatMasker" # Experiment output folder
       fasta: "/research/hvw5426/splitted_fasta/chr6.fasta" # input fasta sequence
       simulation_params:
         length: 100000 # number of kmers selected 
         start_pos: 559707 # start position 
         precision: 1e-10 # precision for Newton Method
         replicates: 100 # number of replicates for mutation process
       plot_params:
         k_variation:
           thresholds: [32, 630] # threshold to split the boxplot panels.
           ylimits: [0.8, 0.02, 1.0] # max scale for each panel
         r_variation:
           threshold: 0.25 # threshold to split the boxplot panels.
   ```

3. un following command, to generate the statistics information of sequences in the paper

   ```
   snakemake --core [N]
   ```

4. The results are in your specified folders, like `./chr6_RepeatMasker`. The simulation outputs for each $r,k$ settings can be found in folders `./specified/varing_k` and `./specified/varing_r`.  The output files are in csv form, three columns represent r_strong, r_weak and r_mash. The boxplots are also stored in the specified folder. 

## Comparasion_of_rstrong_and_sketch

1. `cd ./ Comparasion_of_rstrong_and_sketch/` and run `make` in current folder to compile.

2. Users can modify settings in config.yaml. Here is an example

   The range of $k,r$ a hard code in snakefile.

   ```config.yaml
   analyses:
     - name: "chr6_RepeatMasker_sketch_100" # Experiment name
       output_dir: "chr6_RepeatMasker_sketch_100" # Experiment output folder
       fasta: "/research/hvw5426/splitted_fasta/chr6.fasta" # input fasta sequence
       simulation_params:
         length: 100000 # number of kmers selected 
         start_pos: 559707 # start position 
         precision: 1e-10 # precision for Newton Method
         replicates: 100 # number of replicates for mutation process
         sketch_repeats: 100 # number of replicates for sketch process
       plot_params:
         k_variation:
           thresholds: [32, 630] # threshold to split the boxplot panels.
           ylimits: [0.8, 0.02, 1.0] # max scale for each panel
         r_variation:
           threshold: 0.25 # threshold to split the boxplot panels.
   ```

3. un following command, to generate the statistics information of sequences in the paper

   ```
   snakemake --core [N]
   ```

4. The results are in your specified folders, like `./chr6_RepeatMasker_sketch_100`. The simulation outputs for each $r,k$ settings can be found in folders `./specified/varing_k` and `./specified/varing_r`.  The output files are in csv form, three columns represent r_strong, r_sketch ($\theta=0.1$) and r_sketch ($\theta=0.01$​) . The boxplots are also stored in the specified folder. 

## Sketch_bias_plot/

1. `cd ./Sketch_bias_plot//` and run `make` in current folder to compile.

2. Users can modify settings in config.yaml. Here is an example

   The range of $k,r$ a hard code in snakefile. Because we hope to check the bias of only sketch process, we mutate the string only once (then we remove the `replicates` parameter from config.yaml). 

   ```config.yaml
   analyses:
     - name: "chr6_RepeatMasker" # Experiment name
       output_dir: "chr6_RepeatMasker" # Experiment output folder
       fasta: "/research/hvw5426/splitted_fasta/chr6.fasta" # input fasta sequence
       simulation_params:
         length: 100000 # number of kmers selected 
         start_pos: 559707 # start position 
         precision: 1e-10 # precision for Newton Method
         sketch_repeats: 100 # number of replicates for sketch process
       plot_params:
         k_variation:
           thresholds: [32, 630] # threshold to split the boxplot panels.
           ylimits: [0.8, 0.02, 1.0] # max scale for each panel
         r_variation:
           threshold: 0.151 # threshold to split the boxplot panels.
           ylimits: [0.2, 1.1] # max scale for each panel
   ```

3. un following command, to generate the statistics information of sequences in the paper

   ```
   snakemake --core [N]
   ```

4. The results are in your specified folders, like `./chr6_RepeatMasker`. The simulation outputs for each $r,k$ settings can be found in folders `./specified/varing_k` and `./specified/varing_r`.  The output files are in csv form, three columns represent r_strong (does not change in a output file), r_sketch ($\theta=0.1$) and r_sketch ($\theta=0.01$) . The boxplots are also stored in the specified folder. 