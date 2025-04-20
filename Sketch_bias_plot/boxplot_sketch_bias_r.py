import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from glob import glob
import argparse

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Sketch R-value Variation Boxplot')
    parser.add_argument('-d', '--directory', type=str, default='./Estimate_sketch_r_mut1',
                        help='Input directory with sketch output files')
    parser.add_argument('-k', '--k_value', type=int, default=30,
                        help='k value used in the filenames')
    parser.add_argument('--r_threshold', type=float, default=None,
                        help='Threshold for splitting r values (default: median of r values)')
    parser.add_argument('--ylim', nargs='+', type=float, default=[0.2, 1.1],
                        help='Maximum y-axis limits for the two panels [small_r, large_r]')
    args = parser.parse_args()
    
    # Process all the data files
    all_data = process_all_files(args.directory, args.k_value)
    
    if all_data is not None:
        print(f"Processed data for {len(all_data['true_r'].unique())} different r values")
        
        # Create the dual-panel boxplot
        create_dual_subplot_boxplots(all_data, args.r_threshold, args.ylim, args.k_value)
    else:
        print("No data to visualize")

# Function to extract r value from filename
def extract_r_from_filename(filename):
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
        true_r = extract_r_from_filename(os.path.basename(file_path))
        if true_r is None:
            return pd.DataFrame()
        
        # Print debugging info
        print(f"Processing file: {file_path}, extracted r value: {true_r}")
        
        # Try to read the CSV file
        try:
            # The file has three columns: r_strong, r_sketch (theta=0.1), r_sketch (theta=0.01)
            df = pd.read_csv(file_path, header=None, 
                           names=['r_strong', 'r_sketch_theta_0.1', 'r_sketch_theta_0.01'],
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
                    if len(parts) >= 3:  # We expect 3 columns now
                        try:
                            rows.append([float(parts[0]), float(parts[1]), float(parts[2])])
                        except ValueError:
                            print(f"  Warning: Could not convert line to floats: {line}")
                            continue
            
            if rows:
                df = pd.DataFrame(rows, columns=['r_strong', 'r_sketch_theta_0.1', 'r_sketch_theta_0.01'])
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
        return pd.DataFrame(columns=['r_strong', 'r_sketch_theta_0.1', 'r_sketch_theta_0.01', 'true_r'])

# Function to process all files
def process_all_files(directory='./Estimate_sketch_r_mut1', k_value=30):
    # Define the pattern for files
    pattern = os.path.join(directory, f"r*_k{k_value}.output")
    file_list = glob(pattern)
    
    if not file_list:
        print(f"No files found matching the pattern r*_k{k_value}.output in {directory}")
        # Try alternative directories
        alt_directories = ['.', './output', '../output']
        for alt_dir in alt_directories:
            pattern = os.path.join(alt_dir, f"r*_k{k_value}.output")
            file_list = glob(pattern)
            if file_list:
                print(f"Found files in alternative directory: {alt_dir}")
                break
        
        if not file_list:
            print(f"No files found in any directory")
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

# Function to prepare data for plotting
def prepare_data_for_plotting(data):
    # Sort the dataframe by true_r
    data = data.sort_values('true_r')
    
    # Get unique r values (sorted)
    true_r_values = sorted(data['true_r'].unique())
    
    print(f"True r values: {true_r_values}")
    
    # Create a new column for method, combining r_sketch with theta
    melted_data = []
    
    # Add r_sketch data for each theta (0.1 and 0.01)
    for r in true_r_values:
        r_data = data[data['true_r'] == r]
        
        # Get r_strong values (same for all rows with the same true_r)
        r_strong = r_data['r_strong'].iloc[0]
        
        for _, row in r_data.iterrows():
            # Add r_sketch (theta = 0.1)
            melted_data.append({
                'true_r': row['true_r'],
                'estimated_r': row['r_sketch_theta_0.1'],
                'method': 'r_sketch (θ=0.1)',
                'r_strong': r_strong
            })
            
            # Add r_sketch (theta = 0.01)
            melted_data.append({
                'true_r': row['true_r'],
                'estimated_r': row['r_sketch_theta_0.01'],
                'method': 'r_sketch (θ=0.01)',
                'r_strong': r_strong
            })
    
    return pd.DataFrame(melted_data), true_r_values

# Function to create boxplots with two subplots based on r value
def create_dual_subplot_boxplots(data, r_threshold=None, ylimits=[0.2, 1.1], k_value=30):
    # Prepare data for plotting
    melted_data, true_r_values = prepare_data_for_plotting(data)
    
    # Custom color palette for r_sketch (theta = 0.1) and r_sketch (theta = 0.01)
    custom_palette = {
        'r_sketch (θ=0.1)': '#0000FF',  # Blue
        'r_sketch (θ=0.01)': '#FFD700'   # Yellow
    }
    
    # Calculate threshold for splitting the data
    if r_threshold is None:
        r_threshold = np.percentile(true_r_values, 50)  # Default to median
        print(f"Using median as r threshold: {r_threshold}")
    else:
        print(f"Using provided r threshold: {r_threshold}")
    
    # Ensure we have two y-limits
    if len(ylimits) < 2:
        ylimits = ylimits + [0.2] * (2 - len(ylimits))
    
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
    def plot_on_axis(ax, data, r_values, title, ylim=None):
        if len(data) == 0 or len(r_values) == 0:
            ax.text(0.5, 0.5, 'No data available', 
                    horizontalalignment='center', verticalalignment='center')
            return None, None, None
        
        # Create the boxplot for r_sketch values
        sns.set(style="whitegrid", font_scale=1.2)
        sns_plot = sns.boxplot(x='true_r', y='estimated_r', hue='method', data=data,
                         palette=custom_palette,
                         width=0.9,
                         linewidth=1.0,
                         fliersize=5,
                         ax=ax)
        
        # Add r_strong values as horizontal lines across each true_r group
        x_positions = range(len(r_values))
        for i, r in enumerate(r_values):
            # Get the r_strong value for this true_r
            r_data = data[data['true_r'] == r]
            if not r_data.empty:
                r_strong = r_data['r_strong'].iloc[0]
                
                # Calculate the width of each group (for the horizontal line)
                width = 0.4  # Half the default width of boxplot groups
                
                # Draw horizontal line for r_strong (in red)
                ax.plot([i-width, i+width], [r_strong, r_strong], '-', color='#FF0000', 
                        linewidth=2.0, solid_capstyle='round')
        
        # Add a line showing the true r values (baseline)
        # Remove the label parameter to avoid duplicate "True r" in legend
        ax.plot(x_positions, r_values, 'o-', color='gray', alpha=0.8, markersize=8, 
                linestyle='--', linewidth=2.0)
        
        # Set x-axis ticks to show the true r values
        ax.set_xticks(x_positions)
        ax.set_xticklabels([f"{r:.3f}" for r in r_values], rotation=45)
        
        # Set title and labels
        ax.set_title(title, fontsize=24)
        ax.set_xlabel('True r Value', fontsize=20)
        ax.set_ylabel('Estimated r Value', fontsize=20)
        
        # Set y-axis limits if provided
        # if ylim is not None:
        #     max_y = max(data['estimated_r'].max(), max(r_values) * 1.1)
        #     ax.set_ylim(0, min(max_y, ylim))
        # Set y-axis limits if provided
        if ylim is not None:
            data_max = max(data['estimated_r'].max(), max(r_values))
            upper_limit = max(ylim, data_max * 1.1) 
            ax.set_ylim(0, upper_limit)
        
        # Add grid for easier reading
        ax.grid(True, axis='y', linestyle='--', alpha=0.7)
        
        # Add a manually created legend entry for r_strong and true r
        from matplotlib.lines import Line2D
        custom_lines = [
            Line2D([0], [0], color='#FF0000', lw=2),
            Line2D([0], [0], color='gray', linestyle='--', marker='o', markersize=6, lw=2)
        ]
        custom_labels = ['r_strong', 'True r']
        
        # Get the existing legend handles and labels
        handles, labels = ax.get_legend_handles_labels()
        
        # Combine with the custom legend entries
        all_handles = handles + custom_lines
        all_labels = labels + custom_labels
        
        return sns_plot, all_handles, all_labels
    
    # Plot on first subplot (r < threshold)
    plot1, handles1, labels1 = plot_on_axis(ax1, small_r_data, small_r_values, 
                                          f'r < {r_threshold:.3f}', ylim=ylimits[0])
    
    # Plot on second subplot (r >= threshold)
    plot2, handles2, labels2 = plot_on_axis(ax2, large_r_data, large_r_values, 
                                          f'r ≥ {r_threshold:.3f}', ylim=ylimits[1])
    
    # Remove both legends first
    for ax in [ax1, ax2]:
        if ax.get_legend():
            ax.get_legend().remove()
    
    # Create a single legend on the first subplot with all entries
    if plot1 is not None and handles1 is not None:
        ax1.legend(handles=handles1, labels=labels1, fontsize=16, frameon=True, framealpha=0.9, loc='upper left')
    elif plot2 is not None and handles2 is not None:
        ax2.legend(handles=handles2, labels=labels2, fontsize=16, frameon=True, framealpha=0.9, loc='upper left')
    
    # Add overall title
    plt.suptitle(f'Comparison of r_strong and r_sketch with different θ values (k={k_value})', fontsize=28)
    
    # Adjust layout
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    
    # Save the plot
    plt.savefig(f"dual_sketch_boxplot_k{k_value}.png", dpi=300)
    print(f"Plot saved as 'dual_sketch_boxplot_k{k_value}.png'")

# Main execution
if __name__ == "__main__":
    main()