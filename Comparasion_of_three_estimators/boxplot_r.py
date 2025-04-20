import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from glob import glob

def main():
    parser = argparse.ArgumentParser(description='R-value Variation Boxplot')
    parser.add_argument('-d', '--directory', type=str, help='Input directory', required=True)
    parser.add_argument('-t', '--threshold', type=float, default=0.25, 
                        help='Threshold for splitting r values')
    parser.add_argument('-k', '--kmer', type=int, default=30, 
                        help='k-mer size to analyze')
    args = parser.parse_args()

    global k
    k = args.kmer

    all_data = process_all_files(args.directory)
    
    if all_data is not None:
        create_dual_subplot_boxplots(
            all_data, 
            r_threshold=args.threshold
        )

# Function to extract r value from filename
def extract_r_value(filename):
    try:
        # Extract the r value from a filename like "r0.001_k30.output"
        r_part = filename.split('_')[0]
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
            df = pd.read_csv(file_path, header=None, names=['r_strong', 'r_weak', 'r_Mash'],
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
                    if len(parts) >= 3:
                        try:
                            rows.append([float(parts[0]), float(parts[1]), float(parts[2])])
                        except ValueError:
                            print(f"  Warning: Could not convert line to floats: {line}")
                            continue
            
            if rows:
                df = pd.DataFrame(rows, columns=['r_strong', 'r_weak', 'r_Mash'])
            else:
                raise ValueError("No valid rows found after manual parsing")
        
        # Add the true r value to the DataFrame
        df['true_r'] = true_r
        
        # Print a sample of the processed data
        print(f"Sample of processed data: {df.head(2)}")
        
        return df
        
    except Exception as e:
        print(f"Error processing file {file_path}: {e}")
        # Return an empty DataFrame with the correct structure
        return pd.DataFrame(columns=['r_strong', 'r_weak', 'r_Mash', 'true_r'])

# Function to process all files
def process_all_files(directory='.'):
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

# Function to create boxplots with two subplots based on r value
def create_dual_subplot_boxplots(data, r_threshold=0.25):
    # Sort the dataframe by true_r
    data = data.sort_values('true_r')
    
    # Split data into two groups based on r value
    small_r_data = data[data['true_r'] < r_threshold]
    large_r_data = data[data['true_r'] >= r_threshold]
    
    # Get unique r values for each group (sorted)
    small_r_values = sorted(small_r_data['true_r'].unique())
    large_r_values = sorted(large_r_data['true_r'].unique())
    
    print(f"Small r values (< {r_threshold}): {small_r_values}")
    print(f"Large r values (>= {r_threshold}): {large_r_values}")
    
    # Custom color palette with high contrast colors
    custom_palette = {"r_strong": "#FF0000", "r_weak": "#0000FF", "r_Mash": "#FFD700"}
    
    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(24, 12), gridspec_kw={'width_ratios': [3, 1]})
    
    # Function to plot on a specific axis with specific data
    def plot_on_axis(ax, data, r_values, title):
        # Melt the dataframe to prepare for boxplot
        melted_data = pd.melt(data, 
                            id_vars=['true_r'], 
                            value_vars=['r_strong', 'r_weak', 'r_Mash'],
                            var_name='method', 
                            value_name='estimated_r')
        
        # Create the boxplot
        sns.set(style="whitegrid", font_scale=1.2)
        sns_plot = sns.boxplot(x='true_r', y='estimated_r', hue='method', data=melted_data,
                         palette=custom_palette,
                         width=0.9,
                         linewidth=1.0,
                         fliersize=5,
                         ax=ax)
        
        # Add a line showing the true r values (baseline)
        for i, r in enumerate(r_values):
            ax.plot(i, r, 'o', color='gray', alpha=0.8, markersize=8)
        
        # Connect the true r values with a line
        ax.plot(range(len(r_values)), r_values, color='gray', alpha=0.6, linestyle='--', 
                linewidth=2.0, label='True r')
        
        # Set x-axis ticks to show the true r values
        ax.set_xticks(range(len(r_values)))
        ax.set_xticklabels([f"{r:.3f}" for r in r_values], rotation=45)
        
        # Set title and labels
        ax.set_title(title, fontsize=24)
        ax.set_xlabel('True r Value', fontsize=20)
        ax.set_ylabel('Estimated r Value', fontsize=20)
        
        # Add grid for easier reading
        ax.grid(True, axis='y', linestyle='--', alpha=0.7)
        
        # Handle the case where we need to adjust y-axis limits
        if len(data) > 0:
            ymin = data[['r_strong', 'r_weak', 'r_Mash']].min().min() * 0.9
            ymax = data[['r_strong', 'r_weak', 'r_Mash']].max().max() * 1.1
            
            # Ensure the true r values are visible
            ymin = min(ymin, min(r_values) * 0.9)
            ymax = max(ymax, max(r_values) * 1.1)
            
            # Set y-axis limits
            ax.set_ylim(max(0, ymin), ymax)
        
        return sns_plot
    
    # Plot on first subplot (r < threshold)
    if len(small_r_values) > 0:
        plot1 = plot_on_axis(ax1, small_r_data, small_r_values, f'r < {r_threshold}')
    
    # Plot on second subplot (r >= threshold)
    if len(large_r_values) > 0:
        plot2 = plot_on_axis(ax2, large_r_data, large_r_values, f'r ≥ {r_threshold}')
    
    # Create common legend for both subplots
    handles, labels = ax1.get_legend_handles_labels()
    
    # Keep the first subplot's legend and remove the second subplot's legend
    if ax2.get_legend():
        ax2.get_legend().remove()
    
    # Keep or customize the legend in the first subplot
    if ax1.get_legend():
        # Improve the first subplot's legend
        ax1.legend(fontsize=16, frameon=True, framealpha=0.9, loc='upper left')
    
    # Add overall title
    plt.suptitle(f'Comparison of three estimators (k={k})', fontsize=28)
    
    # Adjust layout to make room for the titles
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    
    # Save the plot
    output_filename = f"dual_boxplot_k{k}.png"
    plt.savefig(output_filename, dpi=300)
    print(f"Plot saved as '{output_filename}'")

# Main execution
if __name__ == "__main__":
    main()