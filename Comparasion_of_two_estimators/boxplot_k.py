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
        create_dynamic_subplot_boxplot(
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
                if len(parts) >= 2:  # Changed from 3 to 2 since we're removing r_weak
                    try:
                        # Only keeping r_hat (formerly r_strong) and r_Mash
                        rows.append([float(parts[0]), float(parts[1])])
                    except ValueError:
                        continue
        
        if rows:
            df = pd.DataFrame(rows, columns=['r_hat', 'r_Mash'])
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

def plot_on_axis(ax, data, k_values, panel_index, fixed_r):
    # Define colors for different methods
    colors = {"r_hat": "#D00000", "r_Mash": "#D0A000"}  # D0A000
    
    # If there's no data for this panel, add a message
    if len(k_values) == 0:
        ax.text(0.5, 0.5, 'No data in this range', 
                horizontalalignment='center',
                verticalalignment='center',
                transform=ax.transAxes,
                fontsize=16)
        # Keep panel title without text
        return []
    
    positions = {k: i for i, k in enumerate(k_values)}
    
    legend_handles = []
    
    # Offset for each method's boxplot
    offset = [-0.15, 0.15]
    
    # Create boxplots for each method
    for i, method in enumerate(['r_hat', 'r_Mash']):
        boxplot_data = []
        pos = []
        
        for k in k_values:
            values = data[data['k_value'] == k][method].dropna()
            if len(values) > 0:
                boxplot_data.append(values)
                pos.append(positions[k] + offset[i])
        
        if boxplot_data:
            # Increase linewidth for all boxplot elements - more significant increase for crowded areas
            linewidth_multiplier = 2.0
            if len(boxplot_data) > 10:
                # For more crowded panels, use thicker lines
                linewidth_multiplier = 0.5
            
            bp = ax.boxplot(
                boxplot_data, 
                positions=pos,
                widths=0.15,
                patch_artist=True,
                boxprops={'facecolor': colors[method], 
                          'color': colors[method], 
                          'alpha': 0.7, 
                          'linewidth': 2.5 * linewidth_multiplier},  # Increased from 1.5 to 2.5
                medianprops={'color': colors[method], 
                             'linewidth': 3.0 * linewidth_multiplier},  # Increased from 1.5 to 3.0
                whiskerprops={'color': colors[method], 
                              'linewidth': 2.0 * linewidth_multiplier,  # Increased from 1.2 to 2.0
                              'alpha': 0.9},
                capprops={'color': colors[method], 
                          'linewidth': 2.0 * linewidth_multiplier,  # Increased from 1.2 to 2.0
                          'alpha': 0.9},
                # Make sure we're including all data in the box plots
                showfliers=False  # Remove outliers
            )
            
            legend_handles.append(plt.Rectangle((0,0), 1, 1, fc=colors[method], alpha=0.7))
    
    # Add horizontal line for true r value - increased thickness
    line = ax.axhline(y=fixed_r, color='grey', linestyle='--', linewidth=3.0)  # Increased from 2.0 to 3.0
    legend_handles.append(line)
    
    # Set labels with simplified x-axis label
    ax.set_xlabel('k', fontsize=30, fontweight='bold')  # Added bold
    
    # Determine which k values to show based on panel index and number of data points
    if len(k_values) > 10:
        # For panels with many data points, show fewer ticks
        step = max(2, len(k_values) // 8)  # Reduced from 10 to 8 ticks max for more labels
        shown_k = k_values[::step]
    else:
        # For panels with fewer data points, show all or every other tick
        step = 2 if len(k_values) > 5 else 1
        shown_k = k_values[::step]
        
    ax.set_xticks([positions[k] for k in shown_k])
    ax.set_xticklabels([str(k) for k in shown_k], rotation=45, fontsize=30, fontweight='bold')  # Added bold
    
    # Increase y-axis tick font size and make bold
    ax.tick_params(axis='y', labelsize=30)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontweight('bold')
    
    # Add grid for better readability - made grid lines bolder
    ax.grid(axis='y', linestyle='--', alpha=0.7, linewidth=1.0)  # Increased grid linewidth
    
    # Style the axis - increased linewidth for all spines
    for spine in ['top', 'right', 'bottom', 'left']:
        ax.spines[spine].set_alpha(0.5)  # Increased from 0.3 to 0.5
        ax.spines[spine].set_linewidth(1.0)  # Increased from 0.5 to 1.0
    
    return legend_handles
    
    # Add horizontal line for true r value - increased thickness
    line = ax.axhline(y=fixed_r, color='grey', linestyle='--', linewidth=3.0)  # Increased from 2.0 to 3.0
    legend_handles.append(line)
    
    # Set labels with simplified x-axis label
    ax.set_xlabel('k', fontsize=30, fontweight='bold')  # Added bold
    
    # Determine which k values to show based on panel index and number of data points
    if len(k_values) > 10:
        # For panels with many data points, show fewer ticks
        step = max(2, len(k_values) // 8)  # Reduced from 10 to 8 ticks max for more labels
        shown_k = k_values[::step]
    else:
        # For panels with fewer data points, show all or every other tick
        step = 2 if len(k_values) > 5 else 1
        shown_k = k_values[::step]
        
    ax.set_xticks([positions[k] for k in shown_k])
    ax.set_xticklabels([str(k) for k in shown_k], rotation=45, fontsize=30, fontweight='bold')  # Added bold
    
    # Increase y-axis tick font size and make bold
    ax.tick_params(axis='y', labelsize=30)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontweight('bold')
    
    # Add grid for better readability - made grid lines bolder
    ax.grid(axis='y', linestyle='--', alpha=0.7, linewidth=1.0)  # Increased grid linewidth
    
    # Style the axis - increased linewidth for all spines
    for spine in ['top', 'right', 'bottom', 'left']:
        ax.spines[spine].set_alpha(0.5)  # Increased from 0.3 to 0.5
        ax.spines[spine].set_linewidth(1.0)  # Increased from 0.5 to 1.0
    
    return legend_handles

def set_panel_y_limits(ax, data_subset, max_limit, fixed_r, panel_name):
    if not data_subset.empty:
        # Calculate boxplot boundaries for each method and each k value
        all_lower_bounds = []
        all_upper_bounds = []
        
        for method in ['r_hat', 'r_Mash']:
            for k in data_subset['k_value'].unique():
                values = data_subset[data_subset['k_value'] == k][method].dropna()
                if len(values) >= 5:  # Need enough data for meaningful quartiles
                    q1 = np.percentile(values, 25)
                    q3 = np.percentile(values, 75)
                    iqr = q3 - q1
                    # Calculate whisker boundaries (which are what we see in boxplot)
                    lower_bound = max(0, q1 - 1.5 * iqr)
                    upper_bound = q3 + 1.5 * iqr
                    all_lower_bounds.append(lower_bound)
                    all_upper_bounds.append(upper_bound)
        
        # If we have calculated boundaries, use them
        if all_lower_bounds and all_upper_bounds:
            min_val = min(all_lower_bounds)
            max_val = max(all_upper_bounds)
        else:
            # Fallback to actual data range
            min_val = data_subset[['r_hat', 'r_Mash']].min().min()
            max_val = data_subset[['r_hat', 'r_Mash']].max().max()
        
        # Make sure fixed_r is within range if needed
        min_val = min(min_val, fixed_r)
        max_val = max(max_val, fixed_r)
        
        # Add minimal padding
        y_min = max(0, min_val - (max_val - min_val) * 0.05)  # 5% padding below
        y_max = max_val + (max_val - min_val) * 0.05  # 5% padding above
        
        # Make sure y_max doesn't exceed max_limit
        y_max = min(y_max, max_limit) 
        
        # Set the limits
        ax.set_ylim(y_min, y_max)
        
        # Set nice round y-tick values - approximately 6 ticks
        range_size = y_max - y_min
        
        # Determine appropriate tick increment based on the range
        if range_size <= 0.05:
            increment = 0.01
        elif range_size <= 0.1:
            increment = 0.02
        elif range_size <= 0.2:
            increment = 0.05
        elif range_size <= 0.5:
            increment = 0.1
        elif range_size <= 1:
            increment = 0.2
        else:
            # For larger ranges, use a power of 10
            power = np.floor(np.log10(range_size))
            increment = 10 ** power / 5
        
        # Calculate nice round starting point that includes y_min
        start = np.floor(y_min / increment) * increment
        
        # Create array of tick positions
        yticks = np.arange(start, y_max + 0.1 * increment, increment)
        
        # If we have too many ticks, increase the step size
        if len(yticks) > 6:
            step = max(2, len(yticks) // 5)  # Aim for 5-6 ticks
            yticks = yticks[::step]
        
        # Apply the ticks to the axis
        ax.set_yticks(yticks)
        
        # Format tick labels to avoid too many decimal places
        if increment < 0.01:
            ax.set_yticklabels([f"{x:.3f}" for x in yticks])
        elif increment < 0.1:
            ax.set_yticklabels([f"{x:.2f}" for x in yticks])
        elif increment < 1:
            ax.set_yticklabels([f"{x:.1f}" for x in yticks])
        else:
            ax.set_yticklabels([f"{int(x)}" if x == int(x) else f"{x:.1f}" for x in yticks])
        
        print(f"{panel_name} panel y-axis range: {y_min:.5f} to {y_max:.5f} with {len(yticks)} ticks")
        print(f"  Based on boxplot boundaries: min={min_val:.5f}, max={max_val:.5f}")
    else:
        # For empty panels, set default range
        ax.set_ylim(0, max_limit)
        print(f"{panel_name} panel: No data, using default y-axis")

def create_dynamic_subplot_boxplot(data, thresholds=[32, 630], ylimits=[0.8, 0.05, 1.0], fixed_r=0.01):
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
    
    # Count how many non-empty panels we have and calculate width ratios
    panel_data = []
    width_ratios = []
    panel_types = []  # Keep track of panel types (small, medium, large)
    
    # Calculate total number of data points for proportional sizing
    total_k_values = len(small_k_values) + len(medium_k_values) + len(large_k_values)
    
    if len(small_k_values) > 0:
        # Calculate width ratio based on data proportion
        ratio = max(1, round(len(small_k_values) / total_k_values * 10)) * 2
        panel_data.append((small_k_data, small_k_values, ylimits[0], "Small k"))
        width_ratios.append(ratio)
        panel_types.append("small")
    
    if len(medium_k_values) > 0:
        # Calculate width ratio based on data proportion
        ratio = max(1, round(len(medium_k_values) / total_k_values * 10))
        panel_data.append((medium_k_data, medium_k_values, ylimits[1], "Medium k"))
        width_ratios.append(ratio)
        panel_types.append("medium")
    
    if len(large_k_values) > 0:
        # Calculate width ratio based on data proportion
        ratio = max(1, round(len(large_k_values) / total_k_values * 10))
        panel_data.append((large_k_data, large_k_values, ylimits[2], "Large k"))
        width_ratios.append(ratio)
        panel_types.append("large")
    
    # Number of non-empty panels
    non_empty_panels = len(panel_data)
    
    print(f"Creating {non_empty_panels} panel(s) with width ratios {width_ratios}")
    
    if non_empty_panels == 0:
        print("No data to plot in any range")
        return
    
    # Create figure with appropriate number of subplots and width ratios
    fig, axes = plt.subplots(1, non_empty_panels, figsize=(24, 10), 
                            gridspec_kw={'width_ratios': width_ratios}) #
    
    # Handle the case of single subplot (axes is not an array in that case)
    if non_empty_panels == 1:
        axes = [axes]  # Convert to list for consistent indexing
    
    # Plot each non-empty panel
    legend_handles_list = []
    for i, (panel_subset, k_vals, y_limit, panel_name) in enumerate(panel_data):
        # Plot on current axis
        legend_handles = plot_on_axis(axes[i], panel_subset, k_vals, i, fixed_r)
        legend_handles_list.append(legend_handles)
        
        # Set y-axis limits for current panel
        set_panel_y_limits(axes[i], panel_subset, y_limit, fixed_r, panel_name)
    
    # Add y-axis label to the first subplot only - made bold
    axes[0].set_ylabel('Estimated r', fontsize=35, fontweight='bold')  # Added bold 
    
    # Determine which panel should have the legend
    legend_panel_index = -1  # Default to last panel
    
    # If we have at least 2 panels and medium panel exists, put legend on medium panel
    if non_empty_panels >= 2:
        try:
            # Find the index of the medium panel if it exists
            medium_panel_index = panel_types.index("medium")
            legend_panel_index = medium_panel_index
        except ValueError:
            # Medium panel doesn't exist, use last panel
            legend_panel_index = non_empty_panels - 1
    
    # Add legend to the selected panel - increased font weight
    if legend_handles_list and legend_panel_index >= 0:
        legend = axes[legend_panel_index].legend(
            legend_handles_list[legend_panel_index], 
            [r"$\hat{r}$", "r_Mash", f'True r={fixed_r}'],
            fontsize=24,
            loc='upper right',
            frameon=True,  # Add frame around legend
            framealpha=0.8,  # Make frame slightly transparent
            edgecolor='black'  # Add border to the legend
        )
        # Bold the legend text
        for text in legend.get_texts():
            text.set_fontweight('bold')
    
    # Adjust layout
    plt.tight_layout()
    
    # Save the figure
    output_filename = f"triple_boxplot_r{fixed_r}.png"
    plt.savefig(output_filename, dpi=100)  # Increased DPI from 100 to 150
    print(f"Plot saved as '{output_filename}'")

# Main execution
if __name__ == "__main__":
    main()