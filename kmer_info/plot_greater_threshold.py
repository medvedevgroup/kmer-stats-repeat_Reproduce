import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import sys
import os

def read_data_from_file(file_path):
    """
    Read data from a file
    
    Parameters:
    file_path (str): Path to the data file
    
    Returns:
    str: File contents
    """
    try:
        with open(file_path, 'r') as file:
            return file.read()
    except FileNotFoundError:
        print(f"File {file_path} not found")
        return ""
    except IOError:
        print(f"Error reading file {file_path}")
        return ""

def process_kmer_data(data):
    """
    Process k-mer data
    
    Parameters:
    data (str): Raw data string
    
    Returns:
    DataFrame: Processed data
    """
    # Process data lines
    lines = [line.strip().split(' : ') for line in data.strip().split('\n')]
    df = pd.DataFrame(lines, columns=['copy_number', 'kmer_count'])
    df['copy_number'] = df['copy_number'].astype(int)
    df['kmer_count'] = df['kmer_count'].astype(int)
    
    return df

def analyze_kmer_distribution(df, thresholds):
    """
    Analyze k-mer distribution
    
    Parameters:
    df (DataFrame): K-mer data
    thresholds (list): List of thresholds to analyze
    
    Returns:
    DataFrame: Analysis results
    """
    # Calculate total k-mer count
    total_kmers = df['kmer_count'].sum()
    
    results = {}
    
    for threshold in thresholds:
        # Calculate percentage of k-mers above each threshold
        kmers_above_threshold = df[df['copy_number'] >= threshold]['kmer_count'].sum()
        percentage = (kmers_above_threshold / total_kmers) * 100
        results[threshold] = percentage
    
    # Create results DataFrame
    results_df = pd.DataFrame({
        'Threshold (≥)': list(results.keys()),
        'Percentage of kmers': list(results.values())
    })
    
    # Print detailed results
    print("Total number of kmers:", total_kmers)
    print("\nPercentage of kmers with copy numbers above thresholds:")
    for index, row in results_df.iterrows():
        print(f"≥ {row['Threshold (≥)']: >3}: {row['Percentage of kmers']:.2f}%")
    
    return results_df

def plot_kmer_distribution(results_df, output_dir):
    """
    Create k-mer distribution plot
    
    Parameters:
    results_df (DataFrame): Analysis results DataFrame
    output_dir (str): Directory to save the plot
    """
    # Set plot style
    plt.style.use('default')
    mpl.rcParams['axes.linewidth'] = 1.5
    mpl.rcParams['axes.edgecolor'] = '#333333'
    
    # Create figure
    plt.figure(figsize=(3.5, 2), dpi=100)
    
    # Sort the dataframe by threshold in descending order for proper y-axis ordering
    results_df = results_df.sort_values('Threshold (≥)', ascending=False)
    
    # Create custom y positions with reduced spacing
    num_bars = len(results_df)
    # Use smaller spacing value (0.7 instead of the default 1.0)
    positions = np.arange(num_bars) * 0.3
    
    # Draw bar plot with custom positions and narrower bars
    bar_color = '#555555'
    bars = plt.barh(positions, 
                   results_df['Percentage of kmers'], 
                   color=bar_color, 
                   height=0.2,  # Narrower bars
                   edgecolor='#333333',
                   linewidth=1)
    
    # Set custom y-tick positions and labels
    plt.yticks(positions, results_df['Threshold (≥)'].astype(str))
    
    # Set background
    plt.gca().set_facecolor('white')
    plt.gcf().set_facecolor('white')
    
    # Add grid lines
    plt.grid(axis='x', linestyle='--', alpha=0.3, color='gray')
    
    # Adjust layout and limits
    plt.subplots_adjust(left=0.12)
    plt.xlim(0, max(results_df['Percentage of kmers'])*1.15)
    
    # Set y-axis limits to fit the custom positions
    plt.ylim(-0.2, max(positions) + 0.2)
    
    # Move x-axis to the top of the plot
    ax = plt.gca()
    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()
    
    # Set the x-axis label at the top
    plt.xlabel('Percentage of kmers (%)', fontsize=16, fontweight='bold', color='#333333')
    
    # Set tick parameters
    plt.tick_params(axis='both', which='major', labelsize=10, colors='#333333', width=1.5)
    plt.xticks(fontweight='bold', fontsize=16)
    plt.yticks(fontweight='bold', fontsize=14)
    
    # Final layout adjustment
    plt.tight_layout(pad=1.5)
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Save figures
    output_path = os.path.join(output_dir, 'kmer_distribution_histogram.png')
    plt.savefig(output_path, dpi=100, bbox_inches='tight')
    plt.close()
    
    print(f"Plot saved to {output_path}")

def main():
    # Check if correct number of arguments is provided
    if len(sys.argv) < 3:
        print("Usage: python plot_greater_threshold.py <input_file> <output_dir> <thresholds>")
        print("Example: python plot_greater_threshold.py input.txt ./output 2 3 4 5")
        sys.exit(1)
    
    # Parse arguments
    input_file = sys.argv[1]
    output_dir = sys.argv[2]
    thresholds = [int(t) for t in sys.argv[3:]]
    
    # Read data
    data = read_data_from_file(input_file)
    
    if data:
        # Process data
        df = process_kmer_data(data)
        
        # Analyze data
        results_df = analyze_kmer_distribution(df, thresholds)
        
        # Create plot
        plot_kmer_distribution(results_df, output_dir)
    else:
        print("Unable to read data, program terminated")

if __name__ == '__main__':
    main()