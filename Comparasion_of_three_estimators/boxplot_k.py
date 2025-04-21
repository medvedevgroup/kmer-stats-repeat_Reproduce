import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from glob import glob

def main():
    parser = argparse.ArgumentParser(description='K-mer Variation Boxplot')
    parser.add_argument('-d', '--directory', type=str, help='Input directory', required=True)
    parser.add_argument('-t', '--thresholds', nargs='+', type=float, default=[32, 630], 
                        help='Thresholds for splitting k values')
    parser.add_argument('--ylim', nargs='+', type=float, default=[0.8, 0.05, 1.0],
                        help='Maximum y-axis limits for the three panels [small_k, medium_k, large_k]')
    parser.add_argument('-r', '--fixed_r', type=float, default=0.01,
                        help='Fixed r value used in the filenames')
    args = parser.parse_args()

    all_data = process_all_files(args.directory, args.fixed_r)
    
    if all_data is not None:
        create_triple_subplot_boxplot(
            all_data, 
            thresholds=args.thresholds,
            ylimits=args.ylim,
            fixed_r=args.fixed_r
        )

def extract_k_value(filename):
    try:
        # Extract the k part from filename (e.g., k4_r0.01.output -> k4)
        k_part = filename.split('_')[0]
        # Convert the numeric part to integer (e.g., k4 -> 4)
        k_value = int(k_part[1:])
        return k_value
    except (ValueError, IndexError) as e:
        print(f"Error extracting k value from filename {filename}: {e}")
        return None

def process_file(file_path):
    try:
        # Get the k value from the filename
        true_k = extract_k_value(os.path.basename(file_path))
        if true_k is None:
            return pd.DataFrame()

        # Manual parsing - most reliable method for parsing the file content
        rows = []
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                parts = line.split(", ")
                if len(parts) >= 3:
                    try:
                        rows.append([float(parts[0]), float(parts[1]), float(parts[2])])
                    except ValueError:
                        continue
        
        if rows:
            df = pd.DataFrame(rows, columns=['r_strong', 'r_weak', 'r_Mash'])
            df['k_value'] = true_k
            return df
        else:
            return pd.DataFrame()
        
    except Exception as e:
        print(f"Error processing file {file_path}: {e}")
        return pd.DataFrame()

def process_all_files(directory='.', fixed_r=0.01):
    # Print the directory we're searching in for debugging
    print(f"Searching for files in directory: {directory}")
    
    # List all files in the directory to verify contents
    try:
        all_files = os.listdir(directory)
        print(f"Files in directory: {all_files}")
    except Exception as e:
        print(f"Error listing directory: {e}")
        all_files = []
    
    # Create the pattern for finding files
    pattern = os.path.join(directory, f"k*_r{fixed_r}.output")
    print(f"Looking for files matching pattern: {pattern}")
    
    # Use glob to find matching files
    file_list = glob(pattern)
    
    if file_list:
        print(f"Found {len(file_list)} matching files:")
        for f in file_list:
            print(f"  - {f}")
    else:
        print("No files found with the specific pattern.")
        
        # Try with a more general pattern
        general_pattern = os.path.join(directory, "k*.output")
        general_files = glob(general_pattern)
        if general_files:
            print(f"Found {len(general_files)} files with general pattern '{general_pattern}':")
            for f in general_files:
                print(f"  - {f}")
            
            # Use the general files if we need to
            file_list = general_files
    
    if not file_list:
        print("Could not find any matching files")
        return None
    
    # Process each file and combine results
    dataframes = []
    for file in file_list:
        df = process_file(file)
        if not df.empty:
            dataframes.append(df)
    
    if not dataframes:
        print("No valid data found in any files")
        return None
        
    all_data = pd.concat(dataframes, ignore_index=True)
    return all_data

def create_triple_subplot_boxplot(data, thresholds=[32, 630], ylimits=[0.8, 0.05, 1.0], fixed_r=0.01):
    # Convert thresholds to integers for cleaner display
    thresholds = [int(t) for t in thresholds]
    
    # Make sure we have 3 y-limits
    if len(ylimits) < 3:
        ylimits = ylimits + [1.0] * (3 - len(ylimits))
    
    # Split data into three groups based on the provided thresholds
    small_k_data = data[data['k_value'] <= thresholds[0]]
    medium_k_data = data[(data['k_value'] > thresholds[0]) & (data['k_value'] <= thresholds[1])]
    large_k_data = data[data['k_value'] > thresholds[1]]
    
    # Get unique k values for each group (sorted)
    small_k_values = sorted(small_k_data['k_value'].unique())
    medium_k_values = sorted(medium_k_data['k_value'].unique())
    large_k_values = sorted(large_k_data['k_value'].unique())
    
    # Print summary of k values in each range
    print(f"Small k values (≤ {thresholds[0]}): {small_k_values}")
    print(f"Medium k values ({thresholds[0]} < k ≤ {thresholds[1]}): {medium_k_values}")
    print(f"Large k values (> {thresholds[1]}): {large_k_values}")
    
    # Create figure with three subplots
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(24, 10), 
                                        gridspec_kw={'width_ratios': [1, 2, 1]})
    
    # Style the axes
    for ax in [ax1, ax2, ax3]:
        ax.spines['top'].set_alpha(0.3)
        ax.spines['right'].set_alpha(0.3)
        ax.spines['bottom'].set_alpha(0.3)
        ax.spines['left'].set_alpha(0.3)
        ax.spines['top'].set_linewidth(0.5)
        ax.spines['right'].set_linewidth(0.5)
        ax.spines['bottom'].set_linewidth(0.5)
        ax.spines['left'].set_linewidth(0.5)
    
    # Define colors for different methods
    colors = {"r_strong": "#D00000", "r_weak": "#0000D0", "r_Mash": "#D0A000"}
    
    # Function to plot boxplots on a given axis
    def plot_on_axis(ax, data, k_values, panel_index):
        # If there's no data for this panel, add a message
        if len(k_values) == 0:
            ax.text(0.5, 0.5, 'No data in this range', 
                    horizontalalignment='center',
                    verticalalignment='center',
                    transform=ax.transAxes,
                    fontsize=16)
            # Set specific title based on panel
            if panel_index == 0:
                ax.set_title(f'k ≤ {thresholds[0]}', fontsize=20)
            elif panel_index == 1:
                ax.set_title(f'{thresholds[0]} < k ≤ {thresholds[1]}', fontsize=20)
            else:
                ax.set_title(f'k > {thresholds[1]}', fontsize=20)
            return []
        
        positions = {k: i for i, k in enumerate(k_values)}
        
        legend_handles = []
        
        # Offset for each method's boxplot
        offset = [-0.25, 0, 0.25]
        
        # Create boxplots for each method
        for i, method in enumerate(['r_strong', 'r_weak', 'r_Mash']):
            boxplot_data = []
            pos = []
            
            for k in k_values:
                values = data[data['k_value'] == k][method].dropna()
                if len(values) > 0:
                    boxplot_data.append(values)
                    pos.append(positions[k] + offset[i])
            
            if boxplot_data:
                bp = ax.boxplot(
                    boxplot_data, 
                    positions=pos,
                    widths=0.15,
                    patch_artist=True,
                    boxprops={'facecolor': colors[method], 'color': colors[method], 'alpha': 0.7, 'linewidth': 1.5},
                    medianprops={'color': colors[method], 'linewidth': 1.5},
                    whiskerprops={'color': colors[method], 'linewidth': 1.2, 'alpha': 0.9},
                    capprops={'color': colors[method], 'linewidth': 1.2, 'alpha': 0.9},
                    showfliers=True  # Always show outliers
                )
                
                legend_handles.append(plt.Rectangle((0,0), 1, 1, fc=colors[method], alpha=0.7))
        
        # Add horizontal line for true r value
        line = ax.axhline(y=fixed_r, color='grey', linestyle='--', linewidth=2)
        legend_handles.append(line)
        
        # Set labels and title
        ax.set_xlabel('k Value', fontsize=18)
        
        # Set specific title based on panel
        if panel_index == 0:
            ax.set_title(f'k ≤ {thresholds[0]}', fontsize=20)
        elif panel_index == 1:
            ax.set_title(f'{thresholds[0]} < k ≤ {thresholds[1]}', fontsize=20)
        else:
            ax.set_title(f'k > {thresholds[1]}', fontsize=20)
        
        # Limit number of x-ticks to avoid overcrowding
        if len(k_values) > 15:
            step = max(1, len(k_values) // 15)
            shown_k = k_values[::step]
        else:
            shown_k = k_values
        
        ax.set_xticks([positions[k] for k in shown_k])
        ax.set_xticklabels([str(k) for k in shown_k], rotation=45)
        
        # Add grid for better readability
        ax.grid(axis='y', linestyle='--', alpha=0.7)
        
        return legend_handles
    
    # Create plots for each group of k values
    legend_handles1 = plot_on_axis(ax1, small_k_data, small_k_values, 0)
    legend_handles2 = plot_on_axis(ax2, medium_k_data, medium_k_values, 1)
    legend_handles3 = plot_on_axis(ax3, large_k_data, large_k_values, 2)
    
    # Set appropriate y-axis limits for each panel
    def set_panel_y_limits(ax, data_subset, max_limit, panel_name):
        if not data_subset.empty:
            min_val = data_subset[['r_strong', 'r_weak', 'r_Mash']].min().min()
            max_val = data_subset[['r_strong', 'r_weak', 'r_Mash']].max().max()
            
            # Add padding but cap at the specified maximum
            y_min = max(0, min_val * 0.9)
            y_max = min(max_limit, max_val * 1.1)
            
            # Make sure fixed_r is visible
            if fixed_r < y_max and fixed_r > y_min:
                # If fixed_r is within range, don't change limits
                pass
            elif fixed_r <= y_min:
                # If fixed_r is below range, extend lower limit
                y_min = fixed_r * 0.9
            elif fixed_r >= y_max:
                # If fixed_r is above range, extend upper limit
                y_max = fixed_r * 1.1
            
            # Set the limits
            ax.set_ylim(y_min, y_max)
            print(f"{panel_name} panel y-axis range: {y_min:.5f} to {y_max:.5f}")
        else:
            # For empty panels, set default range
            ax.set_ylim(0, max_limit)
            print(f"{panel_name} panel: No data, using default y-axis")
    
    # Set limits for each panel with user-specified maximum values
    set_panel_y_limits(ax1, small_k_data, ylimits[0], "Small k")    # Small k values
    set_panel_y_limits(ax2, medium_k_data, ylimits[1], "Medium k")  # Medium k values 
    set_panel_y_limits(ax3, large_k_data, ylimits[2], "Large k")    # Large k values
    
    # Add y-axis label to the first subplot only
    ax1.set_ylabel('Estimated r Value', fontsize=18)
    
    # Add legend to the first subplot if it has data
    if legend_handles1:
        ax1.legend(
            legend_handles1, 
            ['r_strong', 'r_weak', 'r_Mash', f'True r={fixed_r}'],
            fontsize=14,
            loc='upper right'
        )
    # If first panel has no data but second does, add legend there
    elif legend_handles2:
        ax2.legend(
            legend_handles2, 
            ['r_strong', 'r_weak', 'r_Mash', f'True r={fixed_r}'],
            fontsize=14,
            loc='upper right'
        )
    
    # Add title to the entire figure
    plt.suptitle(
        f'Comparison of Three Estimators (r={fixed_r})', 
        fontsize=24
    )
    
    # Adjust layout
    plt.tight_layout()
    plt.subplots_adjust(top=0.9)
    
    # Save the figure with both names
    output_filename = f"triple_boxplot_r{fixed_r}.png"
    plt.savefig(output_filename, dpi=300)
    print(f"Plot saved as '{output_filename}'")

# Main execution
if __name__ == "__main__":
    main()