# kmer-stats-repeat_Reproduce

Reproduce the empirical results  in the paper.

## Dependencies

1. Snakemake
2. Matplotlib

## Datasets

We use following four datasets in our paper. The four complete chromosomes are in folder `./sequence_data/`

| Dataset name               | starting position | Chromosome |
| -------------------------- | ----------------- | ---------- |
| D-easy/chr6_RepeatMasker   | 559,707           | Chr6       |
| D-med/chrY_RBMY1A1         | 22,410,155        | ChrY       |
| D-hard/chrY_SimpleRepeat   | 22,420,785        | ChrY       |
| D-hardest/chr21_centromere | 10,962,854        | Chr21      |

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

5. For overlapped kmers, we list the value of $(occ-sep)$ and the counts of kmers with this difference in `diff_counts.txt`. The contents are shown in the form "occ-sep : #kmer counts". (The definitions of $occ$ and $sep$​ for a kmer can be found in the paper).

6. For the hamming distance distribution among kmer pairs, we list the hamming distance and the counts of kmer pairs in `HD_counts.txt`. We also draw a plot `hamming_distance_distribution_histogram.png` to show this information. 

## Comparison of two estimators

1. This folder is used to compare the performance of $\hat{r}$ and $r_{mash}$. 

2. `cd ./Comparasion_of_two_estimators` and run `make` in current folder to compile.

3. Users can modify settings in config.yaml. Here is an example

   The range of $k,r$ are hard code in snakefile.

   ```config.yaml
   analyses:
     - name: "chr6_RepeatMasker" # Experiment name
       output_dir: "chr6_RepeatMasker" # Experiment output folder
       fasta: "../sequence_data/chr6.fasta" # input fasta sequence
       simulation_params:
         length: 100000 # number of kmers selected 
         start_pos: 559707 # start position 
         precision: 1e-10 # precision for Newton Method
         replicates: 100 # number of replicates for mutation process
         fixed_r: 0.01 # set fixed r
       plot_params:
         k_variation:
           thresholds: [32, 630] # threshold to split the boxplot panels.
           ylimits: [0.8, 0.02, 1.0] # max scale for each panel
         r_variation:
           threshold: 0.25 # threshold to split the boxplot panels.
           fixed_k: 20 # set fixed k
   ```

4. Run following command, to generate the comparisons of three estimators and draw the corresponding box plots in the paper

   ```
   snakemake --core [N]
   ```

5. The results are in your specified folders, like `./chr6_RepeatMasker`. The simulation outputs for each $r,k$​ settings can be found in folders `./specified/varing_k` and `./specified/varing_r`.  The output files are in csv form, three columns represent r_strong, r_weak and r_mash. The boxplots are also stored in the specified folder. 

## Comparison of two estimators (relative error)

1. This folder is used to compare the performance of $\hat{r}$ and $r_{mash}$ based on relative error. 

2. `cd ./Relative_Comparasion_of_two_estimators` and run `make` in current folder to compile.

3. Users can modify settings in config.yaml. Here is an example

   The range of $k,r$ a hard code in snakefile.

   ```config.yaml
   analyses:
     - name: "chr6_RepeatMasker" # Experiment name
       output_dir: "chr6_RepeatMasker" # Experiment output folder
       fasta: "../sequence_data/chr6.fasta" # input fasta sequence
       simulation_params:
         length: 100000 # number of kmers selected 
         start_pos: 559707 # start position 
         precision: 1e-10 # precision for Newton Method
         replicates: 100 # number of replicates for mutation process
         fixed_r: 0.01 # set fixed r
       plot_params:
         k_variation:
           thresholds: [32, 630] # threshold to split the boxplot panels.
           ylimits: [0.8, 0.02, 1.0] # max scale for each panel
         r_variation:
           threshold: 0.25 # threshold to split the boxplot panels.
           fixed_k: 20 # set fixed k
   ```

4. Run following command, to generate the comparisons of three estimators and draw the corresponding box plots in the paper

   ```
   snakemake --core [N]
   ```

5. The results are in your specified folders, like `./chr6_RepeatMasker`. The simulation outputs for each $r,k$​ settings can be found in folders `./specified/varing_k` and `./specified/varing_r`.  The output files are in csv form, three columns represent r_strong, r_weak and r_mash. The boxplots are also stored in the specified folder. 

## Heatmap

1. This folder is used to generate the deviation of $\hat{r}$ on each $r,k$ setting and visualize them in a heatmap. 

2. `cd ./Heatmap/` and run `make` in current folder to compile.

3. Users can modify settings in config.yaml. Here is an example

   ```config.yaml
   analyses:
     - name: "chr21_centromere" # Experiment name
       output_dir: "chr21_centromere" # Experiment output folder
       fasta: "../sequence_data/chr21.fasta" # input fasta sequence
       simulation_params:
         length: 100000 # number of kmers selected
         start_pos: 10962854 # start position 
         precision: 1e-10 # precision for Newton Method
         replicates: 100 # number of replicates for substitution process
         # settings for r values
         r_values: [0.001, 0.011, 0.021, 0.031, 0.041, 0.051, 0.061, 0.071, 0.081, 0.091, 0.101, 0.151, 0.201, 0.251, 0.301]
         # setting for k values
         k_values: [4, 8, 16, 32, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360, 390, 420, 450, 480, 510, 540, 570, 600]
       plot_params:
          heatmap:
           title: "Estimate Performance"
           cmap: "coolwarm"
           annotate: True
           figsize: "14,10"
   ```

4. Run following command, to show the bias of estimator through heat map.

   ```
   snakemake --core [N]
   ```

5. The results are in your specified folders, like `./chr21_centromere/`. The simulation outputs for each $r,k$ settings can be found in folders `./specified_folder/results`. The boxplots are also stored in the specified folder. 

## Bias of estimator based on sketching

1. We use following experiments to show the bias of $r_{sketch}$

2. `cd ./Comparasion_of_rstrong_and_sketch/` and run `make` in current folder to compile.

3. Users can modify settings in config.yaml. Here is an example

   The range of $k,r$ a hard code in snakefile.

   ```config.yaml
   analyses:
     - name: "chr6_RepeatMasker_sketch_100" # Experiment name
       output_dir: "chr6_RepeatMasker_sketch_100" # Experiment output folder
       fasta: "../sequence_data/chr6.fasta" # input fasta sequence
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
           fixed_r: 0.01 # set fixed r
         r_variation:
           threshold: 0.25 # threshold to split the boxplot panels.
           fixed_k: 20 # set fixed k
   ```

4. Run following command, to generate the comparisons of $r_{strong}$ and $r_{sketch}$ for different settings and draw the corresponding box plots in the paper

   ```
   snakemake --core [N]
   ```

5. The results are in your specified folders, like `./chr6_RepeatMasker_sketch_100`. The simulation outputs for each $r,k$ settings can be found in folders `./specified/varing_k` and `./specified/varing_r`.  The output files are in csv form, three columns represent r_strong, r_sketch ($\theta=0.1$) and r_sketch ($\theta=0.01$​) . The boxplots are also stored in the specified folder. 

## Sketch_bias_plot/

1. We check the bias of only sketch process using following experiments, we mutate the string only once and replicate sketching process multiple times.  

2. `cd ./Sketch_bias_plot/` and run `make` in current folder to compile.

3. Users can modify settings in config.yaml. Here is an example

   The range of $k,r$ a hard code in snakefile. 

   ```config.yaml
   analyses:
     - name: "chr6_RepeatMasker" # Experiment name
       output_dir: "chr6_RepeatMasker" # Experiment output folder
       fasta: "../sequence_data/chr6.fasta" # input fasta sequence
       simulation_params:
         length: 100000 # number of kmers selected 
         start_pos: 559707 # start position 
         precision: 1e-10 # precision for Newton Method
         sketch_repeats: 100 # number of replicates for sketch process
       plot_params:
         k_variation:
           thresholds: [32, 630] # threshold to split the boxplot panels.
           ylimits: [0.8, 0.02, 1.0] # max scale for each panel
           fixed_r: 0.01 # set fixed r
         r_variation:
           threshold: 0.151 # threshold to split the boxplot panels.
           ylimits: [0.2, 1.1] # max scale for each panel
           fixed_k: 20 # set fixed k
   ```

4. Run following command, to show the bias of sketching process and draw the corresponding box plots in the paper

   ```
   snakemake --core [N]
   ```

5. The results are in your specified folders, like `./chr6_RepeatMasker`. The simulation outputs for each $r,k$ settings can be found in folders `./specified/varing_k` and `./specified/varing_r`.  The output files are in csv form, three columns represent r_strong (does not change in a output file), r_sketch ($\theta=0.1$) and r_sketch ($\theta=0.01$) . The boxplots are also stored in the specified folder. 

## Error bounds of q

1. We use this experiment to show the theoretical bounds of $q$. 

2. `cd ./Error_bounds_of_q/` and run `make` in current folder to compile.

3. Users can modify settings in config.yaml. Here is an example


```
analyses:
  - name: "chr21_centromere" # Experiment name
    output_dir: "chr21_centromere" # Experiment output folder
    fasta: "../sequence_data/chr21.fasta" # input fasta sequence
    simulation_params:
      length: 100000 # number of kmers selected 
      start_pos: 10962854 # start position 
      precision: 1e-8 # precision for Newton Method
      replicates: 100 # number of replicates
    plot_params:
      r_variation:
        fixed_k: 30 # set fixed k
```

4. The results are in your specified folders, like `./chr21_centromere/`. The simulation outputs for each $r$ settings can be found in folders `./specified/varing_r`.  The first row of output file is in the form of `lower bound, upper bound`, the remaining raws are the value of $\hat{r}$ in each replicates.
