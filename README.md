# kmer-stats-repeat: Reproduction Framework

This repository contains code to reproduce the empirical results from our paper.

## Dependencies

- Snakemake
- Matplotlib
- C++17

## Datasets

Four genomic datasets used in our paper are in the `./sequence_data/` directory. Their detailed information is shown below. 

| Dataset Name               | Chromosome | Starting Position | Complexity Level |
|----------------------------|------------|-------------------|------------------|
| D-easy/chr6_RepeatMasker   | Chr6       | 559,707           | Easy             |
| D-med/chrY_RBMY1A1         | ChrY       | 22,410,155        | Medium           |
| D-hard/chrY_SimpleRepeat   | ChrY       | 22,420,785        | Hard             |
| D-hardest/chr21_centromere | Chr21      | 10,962,854        | Hardest          |

## $K$-mer statistics of datasets

The `./kmer_info` analyzes k-mer statistics from sequence data.

### Compilation

```bash
g++ -std=c++17 kmer_info.cpp -o kmer_info
```

### Running the Analysis

```bash
snakemake --core [N]
```

Where `[N]` is the number of CPU cores to utilize.

### Results

Results are stored in dataset-specific folders (e.g., `./chr21_centromere/`):

1. **Occurrence k-mer count**: `occ_counts.txt` lists the number of occurrences `i` and the count of k-mers with occurrence `i` in the format "occurrence : count". A visual representation is provided in `kmer_distribution_histogram.png`.

2. **Overlapped k-mers**: `diff_counts.txt` contains the value of `(occ-sep)` and the count of k-mers with this difference in the format "occ-sep : count". (See paper for definitions of `occ` and `sep`).

3. **Hamming distance distribution**: `HD_counts.txt` shows the Hamming distance and the count of k-mer pairs. This information is also visualized in `hamming_distance_distribution_histogram.png`.

## Comparison of Estimators

This module compares the performance of estimators $\hat{r}$ and $r_{mash}$.

### Setup and Configuration

1. Navigate to the comparison directory:
   ```bash
   cd ./Comparison_of_two_estimators
   make
   ```

2. Modify `config.yaml` to specify your experiments:
   ```yaml
   analyses:
     - name: "chr6_RepeatMasker"
       output_dir: "chr6_RepeatMasker"
       fasta: "../sequence_data/chr6.fasta"
       simulation_params:
         length: 100000
         start_pos: 559707
         precision: 1e-10
         replicates: 100
         fixed_r: 0.01
       plot_params:
         k_variation:
           thresholds: [32, 630]
           ylimits: [0.8, 0.02, 1.0]
         r_variation:
           threshold: 0.25
           fixed_k: 20
   ```

3. Run the analysis:
   ```bash
   snakemake --core [N]
   ```

4. Results are stored in the specified output directory (e.g., `./chr6_RepeatMasker/`) with simulation data in `./specified/varying_k` and `./specified/varying_r` subdirectories. The output files are in CSV format with three columns representing $\hat{r}$ and $r_{mash}$ for each simulation replicate. Boxplots visualizing the comparison between these estimators are generated in the specified output directory, with separate panels based on the threshold parameters defined in the configuration file.

## Relative Error Comparison

This module is similar to the previous one but evaluates estimator performance using relative error metrics.

### Setup and Configuration

1. Navigate to the relative comparison directory:
   ```bash
   cd ./Relative_Comparison_of_two_estimators
   make
   ```

2. Modify `config.yaml` to specify your experiments (similar format to the previous module):
   ```yaml
   analyses:
     - name: "chr6_RepeatMasker"
       output_dir: "chr6_RepeatMasker"
       fasta: "../sequence_data/chr6.fasta"
       simulation_params:
         length: 100000
         start_pos: 559707
         precision: 1e-10
         replicates: 100
         fixed_r: 0.01
       plot_params:
         k_variation:
           thresholds: [32, 630]
           ylimits: [0.8, 0.02, 1.0]
         r_variation:
           threshold: 0.25
           fixed_k: 20
   ```

3. Run the analysis:
   ```bash
   snakemake --core [N]
   ```

4. Results are stored in the specified output directory (e.g., `./chr6_RepeatMasker/`) with simulation data in `./specified/varying_k` and `./specified/varying_r` subdirectories. The output files are in CSV format with three columns representing the relative errors of  $\hat{r}$ and $r_{mash}$ for each simulation replicate. Boxplots visualizing these relative errors are generated in the specified output directory, with separate panels based on the threshold parameters defined in the configuration.

## Heatmap Analysis

This module generates a heatmap showing the deviation of $\hat{r}$ across different $(r,k)$ parameter combinations.

### Setup and Configuration

1. Navigate to the heatmap directory:
   ```bash
   cd ./Heatmap/
   make
   ```

2. Modify `config.yaml`:
   ```yaml
   analyses:
     - name: "chr21_centromere"
       output_dir: "chr21_centromere"
       fasta: "../sequence_data/chr21.fasta"
       simulation_params:
         length: 100000
         start_pos: 10962854
         precision: 1e-10
         replicates: 100
         r_values: [0.001, 0.011, 0.021, 0.031, 0.041, 0.051, 0.061, 0.071, 0.081, 0.091, 0.101, 0.151, 0.201, 0.251, 0.301]
         k_values: [4, 8, 16, 32, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360, 390, 420, 450, 480, 510, 540, 570, 600]
       plot_params:
          heatmap:
           title: "Estimate Performance"
           cmap: "coolwarm"
           annotate: True
           figsize: "14,10"
   ```

3. Run the analysis:
   ```bash
   snakemake --core [N]
   ```

4. Results are stored in the specified output directory (e.g., `./chr21_centromere/`), with raw simulation data in `./specified_folder/results`. The primary output is a comprehensive heatmap visualization (`heatmap_deviation.png`) showing the bias of $\hat{r}$ estimator across the specified r and k value combinations. The heatmap is colored according to the selected colormap, with annotations displaying the actual values when enabled in the configuration. The raw data for each $(r,k)$ combination are stored in separate CSV files in the results directory, containing the estimator values for each replicate.

## Sketching Bias Analysis

This module examines the bias of $r_{sketch}$ compared to $r_{strong}$.

### Setup and Configuration

1. Navigate to the comparison directory:
   ```bash
   cd ./Comparison_of_rstrong_and_sketch/
   make
   ```

2. Modify `config.yaml` to specify your experiments:
   ```yaml
   analyses:
     - name: "chr6_RepeatMasker_sketch_100"
       output_dir: "chr6_RepeatMasker_sketch_100"
       fasta: "../sequence_data/chr6.fasta"
       simulation_params:
         length: 100000
         start_pos: 559707
         precision: 1e-10
         replicates: 100
         sketch_repeats: 100
       plot_params:
         k_variation:
           thresholds: [32, 630]
           ylimits: [0.8, 0.02, 1.0]
           fixed_r: 0.01
         r_variation:
           threshold: 0.25
           fixed_k: 20
   ```

3. Run the analysis:
   ```bash
   snakemake --core [N]
   ```

4. Results are stored in the specified output directory (e.g., `./chr6_RepeatMasker_sketch_100/`) with simulation data in `./specified/varying_k` and `./specified/varying_r` subdirectories. Output files are in CSV format with three columns representing r_strong, r_sketch ($\theta=0.1$), and r_sketch ($\theta=0.01$). Visualization is provided through boxplots in the same directory.

## Single-Mutation Sketching Bias

This module focuses on isolating the bias from the sketching process by applying only one mutation and replicating the sketching process multiple times.

1. Navigate to the sketching bias directory:
   ```bash
   cd ./Sketch_bias_plot/
   make
   ```

2. Modify `config.yaml` to specify your experiments:
   ```yaml
   analyses:
     - name: "chr6_RepeatMasker"
       output_dir: "chr6_RepeatMasker"
       fasta: "../sequence_data/chr6.fasta"
       simulation_params:
         length: 100000
         start_pos: 559707
         precision: 1e-10
         sketch_repeats: 100
       plot_params:
         k_variation:
           thresholds: [32, 630]
           ylimits: [0.8, 0.02, 1.0]
           fixed_r: 0.01
         r_variation:
           threshold: 0.151
           ylimits: [0.2, 1.1]
           fixed_k: 20
   ```

3. Run the analysis:
   ```bash
   snakemake --core [N]
   ```

4. Results are stored in the specified output directory (e.g., `./chr6_RepeatMasker/`) with simulation data in `./specified/varying_k` and `./specified/varying_r` subdirectories. The output files are in CSV format with three columns representing r_strong (constant within each output file), r_sketch ($\theta=0.1$), and r_sketch ($\theta=0.01$). Boxplots visualizing these results are also provided in the same directory.

## Error Bounds of q

This experiment demonstrates the theoretical bounds of parameter $q$.

### Setup and Configuration

1. Navigate to the error bounds directory:
   ```bash
   cd ./Error_bounds_of_q/
   make
   ```

2. Configure with `config.yaml`:
   ```yaml
   analyses:
     - name: "chr21_centromere"
       output_dir: "chr21_centromere"
       fasta: "../sequence_data/chr21.fasta"
       simulation_params:
         length: 100000
         start_pos: 10962854
         precision: 1e-8
         replicates: 100
       plot_params:
         r_variation:
           fixed_k: 30
   ```

3. Run the analysis:
   ```bash
   snakemake --core [N]
   ```

4. Results are stored in the specified output directory (e.g., `./chr21_centromere/`) with detailed data in `./specified/varying_r`. Each output file corresponds to a specific r value and contains:
   - The first row showing the lower and upper theoretical bounds in the format "lower_bound, upper_bound"
   - Subsequent rows showing the $\hat{r}$ values for each simulation replicate
   
   A visualization comparing the theoretical bounds with the empirical results of $\hat{r}$ values is also generated, allowing for visual confirmation of the theoretical predictions across different $r$ values while keeping $k$ fixed as specified in the configuration.

## Relation of $P_{empty}$ and Unstableness

This experiment demonstrates $P_{empty}$ can be a diagnostic criterion for determining the unstableness for large $r$ and $k$​.

### Setup and Configuration

1. Navigate to the error bounds directory:
   ```bash
   cd ./P_empty_and_unstableness/
   make
   ```

2. Configure with `config.yaml`:
   ```yaml
   analyses:
     - name: "chr21_centromere"
       output_dir: "chr21_centromere"
       fasta: "../sequence_data/D-hardest.fasta"
       simulation_params:
         length: 100000
         start_pos: 0
         precision: 1e-8
         replicates: 100
       plot_params:
         r_variation:
           threshold: 0.201
           fixed_k: 30
   ```
   
3. Run the analysis:
   ```bash
   snakemake --core [N]
   ```

4. Results are stored in the specified output directory (e.g., `./chr21_centromere/`) with detailed data in `./specified/varying_r`. Each output file corresponds to a specific r value and contains:
   - The first row showing the value of $P_{empty}$
   - Subsequent rows showing the $\hat{r}$ values for each simulation replicate
   
   A visualization empirical results of $\hat{r}$ values and the curve of $P_{empty}$ is also generated.

## $P_{empty}$ Heatmap

This folder is used to generate the heat map of $P_{empty}$ for the settings of different $k$ and $r$.

### Setup and Configuration

1. Navigate to the error bounds directory:
   ```bash
   cd ./P_empty_and_unstableness/
   make
   ```

2. Configure with `config.yaml`:
   ```yaml
   analyses:  
     - name: "chrY_SimpleRepeat" # chrY:22420785-22423059,
       output_dir: "chrY_SimpleRepeat"
       fasta: "../sequence_data/D-hard.fasta"
       simulation_params:
         length: 2264
         start_pos: 0
         precision: 1e-10
         replicates: 100
         r_values: [0.001, 0.031, 0.061, 0.091, 0.121, 0.151, 0.181, 0.211, 0.241, 0.271, 0.301, 0.331, 0.361]
         k_values: [8,10,12,14,16,18,20,22,24,26,28,30,32]
       plot_params:
          heatmap:
           title: "HARD"
           cmap: "Greys"
           annotate: True
           figsize: "10,10"
   ```
   
3. Run the analysis:
   ```bash
   snakemake --core [N]
   ```

4. Results are stored in the specified output directory (e.g., `./chrY_SimpleRepeat/`), with raw simulation data in `./specified_folder/results`. The primary output is a comprehensive heatmap visualization (`heatmap.png`) showing the bias of $\hat{r}$ estimator across the specified r and k value combinations. The heatmap is colored according to the selected colormap, with annotations displaying the actual values when enabled in the configuration. The $P_{empty}$ for each $(r,k)$ combination are stored in separate CSV files in the results directory.

