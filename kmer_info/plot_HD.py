import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import pandas as pd
import sys
import os
from scipy.special import comb

def read_data_from_file(file_path):
    """
    Read data from a file containing Hamming distance data
    
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

def process_hamming_data(data):
    """
    Process Hamming distance data
    
    Parameters:
    data (str): Raw data string
    
    Returns:
    DataFrame: Processed data with Hamming distances and counts
    """
    # Process data lines
    lines = [line.strip().split(' : ') for line in data.strip().split('\n')]
    df = pd.DataFrame(lines, columns=['hamming_distance', 'pair_count'])
    df['hamming_distance'] = df['hamming_distance'].astype(int)
    df['pair_count'] = df['pair_count'].astype(int)
    
    return df

def calculate_theoretical_distribution(k, total_pairs):
    """
    Calculate theoretical Hamming distance distribution for random k-mers
    
    Parameters:
    k (int): Length of k-mers
    total_pairs (int): Total number of k-mer pairs
    
    Returns:
    DataFrame: Theoretical distribution data
    """
    # For random k-mers over {A,C,G,T}, probability of mismatch at each position is 3/4
    # Probability of having exactly h mismatches follows binomial distribution
    
    distances = range(k + 1)
    probabilities = [comb(k, h) * (0.75**h) * (0.25**(k-h)) for h in distances]
    
    # Calculate expected counts based on total pairs
    expected_counts = [p * total_pairs for p in probabilities]
    
    # Create DataFrame
    df_theory = pd.DataFrame({
        'hamming_distance': distances,
        'expected_count': expected_counts
    })
    
    return df_theory

def plot_hamming_histogram(df, output_dir, use_log_scale=True):
    """
    Create Hamming distance histogram plot with precise control
    
    Parameters:
    df (DataFrame): Hamming distance data
    output_dir (str): Directory to save the plot
    use_log_scale (bool): Whether to use logarithmic scale for x-axis
    """
    # Infer k value from data
    k_length = df['hamming_distance'].max()
    
    # Calculate theoretical distribution
    total_pairs = df['pair_count'].sum()
    df_theory = calculate_theoretical_distribution(k_length, total_pairs)
    
    # Create figure with object-oriented approach for better control
    # NOTE: For horizontal bars with x-axis at top, we need to invert the usual width/height
    fig, ax = plt.subplots(figsize=(8, 12), dpi=100)
    
    # Set style parameters
    mpl.rcParams['axes.linewidth'] = 1.5
    mpl.rcParams['axes.edgecolor'] = '#333333'
    ax.set_facecolor('white')
    fig.set_facecolor('white')
    
    # Determine max value for x-axis
    max_count = max(df['pair_count'].max(), df_theory['expected_count'].max())
    
    # Set up x-axis scale and limits
    if use_log_scale:
        # Convert to logarithmic scale
        ax.set_xscale('log', base=10)
        
        # Calculate max power for setting range
        max_power = int(np.ceil(np.log10(max_count)))
        
        # Set x-axis limits FIRST (this is critical)
        # Start exactly at 1.0 (10^0) and end slightly beyond max value
        ax.set_xlim(1.0, 10**(max_power + 0.1))
        
        # Use fixed step of 3 to create nice evenly spaced ticks like 10^0, 10^3, 10^6
        step = 3
        
        # Create powers list starting from 0 with step of 3
        powers = list(range(0, max_power + step, step))
        
        # Trim the list if the last power exceeds max_power significantly
        if powers[-1] > max_power + 1:
            powers = powers[:-1]
            
        # Create tick positions and labels
        tick_positions = [10**p for p in powers]
        tick_labels = [f"$10^{{{p}}}$" for p in powers]
    else:
        # Linear scale
        ax.set_xlim(0, max_count * 1.1)
        
        # Create 4 evenly spaced ticks
        num_ticks = 4
        tick_positions = np.linspace(0, max_count, num_ticks)
        tick_labels = [f"{int(x):,}" for x in tick_positions]
    
    # Set tick positions and labels - with extra bold weight
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, fontweight='extra bold', fontsize=50)
    
    # Plot bars AFTER setting scales and limits
    # Use the 'align' parameter and left edge exactly at 1.0 for log scale
    bars = ax.barh(
        df['hamming_distance'],
        df['pair_count'],
        height=0.8,
        color='#555555',
        edgecolor='#333333',
        linewidth=1,
        left=1.0 if use_log_scale else 0,  # Critical: start exactly at 10^0
        align='center'
    )
    
    # Plot theoretical distribution
    ax.plot(
        df_theory['expected_count'],
        df_theory['hamming_distance'],
        'r--',
        linewidth=2.5,
        alpha=0.8
    )
    
    # Set y-axis ticks with extra bold weight
    min_dist = min(df['hamming_distance'])
    max_dist = max(df['hamming_distance'])
    num_y_ticks = 4
    ytick_values = np.linspace(min_dist, max_dist, num_y_ticks, dtype=int)
    ax.set_yticks(ytick_values)
    ax.set_yticklabels(ytick_values, fontweight='extra bold', fontsize=50)
    
    # Add grid lines
    ax.grid(axis='x', linestyle='--', alpha=0.3, color='gray')
    
    # Add labels
    ax.set_ylabel('Hamming Distance', fontsize=50, fontweight='bold', color='#333333')
    
    # Move x-axis to top and make sure it's visible
    # This ensures the x-axis ticks and labels appear at the top
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    
    # Ensure the top spine is visible
    ax.spines['top'].set_visible(True)
    ax.spines['bottom'].set_visible(False)  # Hide bottom spine
    
    # Set tick parameters
    ax.tick_params(axis='both', which='major', labelsize=30, colors='#333333', width=1.5)
    ax.tick_params(axis='x', which='major', top=True, bottom=False)  # Ensure x ticks are at top
    
    # Ensure bars align with y-axis by moving left spine to x=1
    if use_log_scale:
        # Create a secondary axis with the exact same limits for handling the spine
        # This is a workaround for precise spine positioning in log scale
        ax2 = ax.twiny()
        ax2.set_xlim(ax.get_xlim())
        ax2.set_xscale('log', base=10)
        ax2.set_xticks([])
        ax2.set_xticklabels([])
        ax2.spines['bottom'].set_visible(False)
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        
        # Position the left spine at exactly x=1
        ax2.spines['left'].set_position(('data', 1.0))
        ax2.spines['left'].set_linewidth(1.5)
        ax2.spines['left'].set_edgecolor('#333333')
        
        # Hide the original left spine
        ax.spines['left'].set_visible(False)
    
    # Final layout adjustment
    plt.tight_layout(pad=2.0)
    
    # Manually adjust the position of the x-axis tick labels to be ABOVE the plot
    for tick in ax.get_xticklabels():
        tick.set_y(1.05)  # Move labels up above the plot
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Save figure
    output_path = os.path.join(output_dir, 'hamming_distance_histogram.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    print(f"Plot saved to {output_path}")

def analyze_hamming_data(df):
    """
    Analyze Hamming distance data and print statistics
    
    Parameters:
    df (DataFrame): Hamming distance data
    """
    # Calculate total pairs
    total_pairs = df['pair_count'].sum()
    
    # Find most common Hamming distance
    most_common_idx = df['pair_count'].idxmax()
    most_common_distance = df.loc[most_common_idx, 'hamming_distance']
    most_common_count = df.loc[most_common_idx, 'pair_count']
    
    # Calculate average Hamming distance
    weighted_sum = sum(df['hamming_distance'] * df['pair_count'])
    avg_distance = weighted_sum / total_pairs
    
    # Print statistics
    print("\nHamming Distance Analysis:")
    print(f"Total number of k-mer pairs: {total_pairs:,}")
    print(f"Most common Hamming distance: {most_common_distance} ({most_common_count:,} pairs)")
    print(f"Average Hamming distance: {avg_distance:.2f}")

def main():
    # Check if correct number of arguments is provided
    if len(sys.argv) < 3:
        print("Usage: python plot_HD.py <input_file> <output_dir> [--no-log]")
        print("Example: python plot_HD.py HD_counts.txt ./output")
        sys.exit(1)
    
    # Parse arguments
    input_file = sys.argv[1]
    output_dir = sys.argv[2]
    
    # Check for log scale option
    use_log_scale = True
    if len(sys.argv) > 3 and sys.argv[3] == "--no-log":
        use_log_scale = False
    
    # Read data
    data = read_data_from_file(input_file)
    
    if data:
        # Process data
        df = process_hamming_data(data)
        
        # Analyze data
        analyze_hamming_data(df)
        
        # Create plot
        plot_hamming_histogram(df, output_dir, use_log_scale)
    else:
        print("Unable to read data, program terminated")

if __name__ == '__main__':
    main()