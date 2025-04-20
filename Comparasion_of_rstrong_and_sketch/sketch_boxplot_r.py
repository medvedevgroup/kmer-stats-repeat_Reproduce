import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from glob import glob
import argparse

# Default k value
k = 30

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Sketch R-value Variation Boxplot')
    parser.add_argument('-d', '--directory', type=str, default='./Estimate_sketch_r',
                        help='Input directory with sketch output files')
    parser.add_argument('-t', '--threshold', type=float, default=None, 
                        help='Threshold for splitting r values (default: median of data)')
    parser.add_argument('-k', '--kmer', type=int, default=30,
                        help='k-mer size to analyze (default: 30)')
    args = parser.parse_args()
    
    # Update global k value
    global k
    k = args.kmer
    
    # Process data files
    all_data = process_all_files(args.directory)
    
    if all_data is not None:
        print(f"Processed data for {len(all_data['true_r'].unique())} different r values")
        
        # Create the dual subplot boxplot with specified threshold
        create_dual_subplot_boxplots(all_data, r_threshold=args.threshold)
    else:
        print("No data to visualize")

# Function to extract r value from filename
def extract_r_value(filename):
    try:
        # Extract the r value from a filename like "r0.231_k30.output"
        parts = filename.split('_')
        r_part = parts[0]
        r_value = float(r_part[1:])
        
        return r_value
    except (ValueError, IndexError) as e:
        print(f"Error extracting r value from filename {filename}: {e}")
        return None

# Function to process a single file
def process_file(file_path):
    try:
        # Extract true r value from filename
        true_r = extract_r_value(os.path.basename(file_path))
        if true_r is None:
            return pd.DataFrame()
        
        # Print debugging info
        print(f"Processing file: {file_path}, extracted r value: {true_r}")
        
        # Try to read with specific "comma + space" delimiter
        try:
            # Use pandas with specific delimiter ", " (comma followed by space)
            df = pd.read_csv(file_path, header=None, names=['r_strong', 'r_sketch'],
                            sep=", ", engine='python')
        except Exception as e:
            print(f"Reading with 'comma + space' delimiter failed: {e}")
            
            # Fallback: Try manually parsing
            print(f"Attempt 2: Manual parsing of {file_path}")
            rows = []
            with open(file_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    # Split by comma+space
                    parts = line.split(", ")
                    if len(parts) >= 2:
                        try:
                            rows.append([float(parts[0]), float(parts[1])])
                        except ValueError:
                            print(f"  Warning: Could not convert line to floats: {line}")
                            continue
            
            if rows:
                df = pd.DataFrame(rows, columns=['r_strong', 'r_sketch'])
            else:
                raise ValueError("No valid rows found after manual parsing")
        
        # Add the true r value to the DataFrame
        df['true_r'] = true_r
        # Hard-code two theta values
        df['theta'] = [0.1, 0.01] * (len(df) // 2)
        if len(df) % 2 == 1:  # If odd number of rows, add one more theta value
            df.loc[len(df)-1, 'theta'] = 0.1
        
        # Print a sample of the processed data
        print(f"Sample of processed data: {df.head(2)}")
        
        return df
        
    except Exception as e:
        print(f"Error processing file {file_path}: {e}")
        # Return an empty DataFrame with the correct structure
        return pd.DataFrame(columns=['r_strong', 'r_sketch', 'true_r', 'theta'])

# Function to process all files
def process_all_files(directory='./Estimate_sketch_r'):
    # Try multiple directories
    directories = [directory, '.', './output', '../output']
    
    file_list = []
    for dir_path in directories:
        pattern = os.path.join(dir_path, f"r*_k{k}.output")
        current_files = glob(pattern)
        if current_files:
            print(f"Found {len(current_files)} files in {dir_path}")
            file_list = current_files
            break
    
    if not file_list:
        print(f"No files found matching the pattern r*_k{k}.output in any directory")
        return None
    
    print(f"Found {len(file_list)} files matching the pattern")
    
    # Process each file and combine the results
    dataframes = []
    for file in file_list:
        df = process_file(file)
        if not df.empty:
            dataframes.append(df)
    
    if not dataframes:
        print("No valid data found in any of the files")
        return None
        
    all_data = pd.concat(dataframes, ignore_index=True)
    
    return all_data

# Function to create a melted DataFrame for plotting
def prepare_data_for_plotting(data):
    # Sort the dataframe by true_r
    data = data.sort_values('true_r')
    
    # Get unique r values (sorted)
    true_r_values = sorted(data['true_r'].unique())
    
    # Get unique theta values (sorted)
    theta_values = sorted(data['theta'].unique(), reverse=True)
    
    print(f"True r values: {true_r_values}")
    print(f"Theta values: {theta_values}")
    
    # Create a new column for method, combining r_strong and r_sketch with theta
    melted_data = []
    
    # First, add the r_strong data (same for all theta values)
    for r in true_r_values:
        r_data = data[data['true_r'] == r]
        for _, row in r_data.iterrows():
            melted_data.append({
                'true_r': row['true_r'],
                'estimated_r': row['r_strong'],
                'method': 'r_strong'
            })
    
    # Then, add r_sketch data for each theta
    for theta in theta_values:
        theta_data = data[data['theta'] == theta]
        for r in true_r_values:
            r_theta_data = theta_data[theta_data['true_r'] == r]
            for _, row in r_theta_data.iterrows():
                melted_data.append({
                    'true_r': row['true_r'],
                    'estimated_r': row['r_sketch'],
                    'method': f'r_sketch (θ={theta})'
                })
    
    return pd.DataFrame(melted_data), true_r_values, theta_values

# Function to create boxplots with two subplots based on r value
def create_dual_subplot_boxplots(data, r_threshold=None):
    # Prepare data for plotting
    melted_data, true_r_values, theta_values = prepare_data_for_plotting(data)
    
    # Custom color palette
    methods = ['r_strong'] + [f'r_sketch (θ={theta})' for theta in theta_values]
    colors = ['#FF0000', '#0000FF', '#FFD700', '#00FF00', '#800080']  # Add more colors if needed
    custom_palette = {method: color for method, color in zip(methods, colors)}
    
    # Calculate threshold for splitting the data (if not provided)
    if r_threshold is None:
        r_threshold = np.percentile(true_r_values, 50)
        print(f"Using median as threshold: {r_threshold:.3f}")
    
    # Split data into two groups based on r value
    small_r_data = melted_data[melted_data['true_r'] < r_threshold]
    large_r_data = melted_data[melted_data['true_r'] >= r_threshold]
    
    # Get unique r values for each group (sorted)
    small_r_values = sorted([r for r in true_r_values if r < r_threshold])
    large_r_values = sorted([r for r in true_r_values if r >= r_threshold])
    
    print(f"Small r values (< {r_threshold}): {small_r_values}")
    print(f"Large r values (>= {r_threshold}): {large_r_values}")
    
    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(36, 15), gridspec_kw={'width_ratios': [1.1, 1]})
    
    # Function to plot on a specific axis with specific data
    def plot_on_axis(ax, data, r_values, title):
        # Create the boxplot
        sns.set(style="whitegrid", font_scale=1.2)
        sns_plot = sns.boxplot(x='true_r', y='estimated_r', hue='method', data=data,
                         palette=custom_palette,
                         width=0.9,
                         linewidth=1.0,
                         fliersize=5,
                         ax=ax)
        
        # Add a line showing the true r values (baseline)
        x_positions = range(len(r_values))
        ax.plot(x_positions, r_values, 'o-', color='gray', alpha=0.8, markersize=8, 
                linestyle='--', linewidth=2.0, label='True r')
        
        # Set x-axis ticks to show the true r values
        ax.set_xticks(x_positions)
        ax.set_xticklabels([f"{r:.3f}" for r in r_values], rotation=45)
        
        # Set title and labels
        ax.set_title(title, fontsize=24)
        ax.set_xlabel('True r Value', fontsize=20)
        ax.set_ylabel('Estimated r Value', fontsize=20)
        
        # Add grid for easier reading
        ax.grid(True, axis='y', linestyle='--', alpha=0.7)
        
        return sns_plot
    
    # Plot on first subplot (r < threshold)
    if len(small_r_values) > 0:
        plot1 = plot_on_axis(ax1, small_r_data, small_r_values, f'r < {r_threshold:.3f}')
    else:
        ax1.text(0.5, 0.5, 'No data in this range', 
                horizontalalignment='center', verticalalignment='center',
                transform=ax1.transAxes, fontsize=20)
        ax1.set_title(f'r < {r_threshold:.3f}', fontsize=24)
    
    # Plot on second subplot (r >= threshold)
    if len(large_r_values) > 0:
        plot2 = plot_on_axis(ax2, large_r_data, large_r_values, f'r ≥ {r_threshold:.3f}')
    else:
        ax2.text(0.5, 0.5, 'No data in this range', 
                horizontalalignment='center', verticalalignment='center',
                transform=ax2.transAxes, fontsize=20)
        ax2.set_title(f'r ≥ {r_threshold:.3f}', fontsize=24)
    
    # Handle legends
    if ax2.get_legend():
        ax2.get_legend().remove()
    
    if ax1.get_legend():
        ax1.legend(fontsize=16, frameon=True, framealpha=0.9, loc='upper left')
    
    # Add overall title
    plt.suptitle(f'Comparison of r_strong and r_sketch with different θ values (k={k})', fontsize=28)
    
    # Adjust layout
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    
    # Save the plot
    output_filename = f"dual_sketch_boxplot_k{k}.png"
    plt.savefig(output_filename, dpi=300)
    print(f"Plot saved as '{output_filename}'")

# Main execution
if __name__ == "__main__":
    main()