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
        
        # Create the single panel boxplot
        create_single_panel_boxplot(all_data)
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
        
        print(f"Processing file: {file_path}, extracted r value: {true_r}")
        
        # Read the file and parse the data
        rows = []
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                # Split by comma
                parts = line.split(",")
                if len(parts) >= 3:  # We expect at least 3 columns
                    try:
                        # Extract r_strong and both theta values
                        r_strong = float(parts[0])
                        r_sketch_theta_01 = float(parts[1])
                        r_sketch_theta_001 = float(parts[2])
                        
                        # Add two rows to the dataframe - one for each theta value
                        rows.append([r_strong, r_sketch_theta_01, true_r, 0.1])
                        rows.append([r_strong, r_sketch_theta_001, true_r, 0.01])
                    except ValueError:
                        print(f"  Warning: Could not convert line to floats: {line}")
                        continue
        
        if rows:
            df = pd.DataFrame(rows, columns=['r_strong', 'r_sketch', 'true_r', 'theta'])
            # Print a sample of the processed data
            print(f"Sample of processed data: {df.head(2)}")
            return df
        else:
            print(f"No valid data found in {file_path}")
            return pd.DataFrame(columns=['r_strong', 'r_sketch', 'true_r', 'theta'])
        
    except Exception as e:
        print(f"Error processing file {file_path}: {e}")
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
                'method': r'$\hat{r}$'  # Changed from 'r_strong' to '$\hat{r}$'
            })
    
    # Then, add r_sketch data for each theta
    for theta in theta_values:
        theta_data = data[data['theta'] == theta]
        for r in true_r_values:
            r_theta_data = theta_data[theta_data['true_r'] == r]
            for _, row in r_theta_data.iterrows():
                # MODIFIED HERE: Changed label format for r_sketch with LaTeX notation
                if theta == 0.1:
                    method_label = r'$\hat{r}$ using $\theta = 0.1$'
                elif theta == 0.01:
                    method_label = r'$\hat{r}$ using $\theta = 0.01$'
                else:
                    method_label = f'$\\hat{{r}}$ using $\\theta = {theta}$'
                    
                melted_data.append({
                    'true_r': row['true_r'],
                    'estimated_r': row['r_sketch'],
                    'method': method_label
                })
    
    return pd.DataFrame(melted_data), true_r_values, theta_values

# Function to calculate non-outlier min and max for each group
def calculate_non_outlier_range(data, true_r_values, methods):
    # Store min and max for each group
    min_vals = []
    max_vals = []
    
    # For each r value and method combination
    for r in true_r_values:
        for method in methods:
            # Get the data for this group
            group_data = data[(data['true_r'] == r) & (data['method'] == method)]['estimated_r']
            
            if not group_data.empty:
                # Calculate IQR (Interquartile Range)
                q1 = group_data.quantile(0.25)
                q3 = group_data.quantile(0.75)
                iqr = q3 - q1
                
                # Calculate lower and upper bounds (1.5 * IQR)
                lower_bound = q1 - 1.5 * iqr
                upper_bound = q3 + 1.5 * iqr
                
                # Filter out outliers
                non_outliers = group_data[(group_data >= lower_bound) & (group_data <= upper_bound)]
                
                if not non_outliers.empty:
                    min_vals.append(non_outliers.min())
                    max_vals.append(non_outliers.max())
    
    # Also include true r values in the range
    min_vals.extend(true_r_values)
    max_vals.extend(true_r_values)
    
    # Return min and max of all non-outlier values
    if min_vals and max_vals:
        return min(min_vals), max(max_vals)
    else:
        return None, None

# Function to create a single panel boxplot
def create_single_panel_boxplot(data):
    # Prepare data for plotting
    melted_data, true_r_values, theta_values = prepare_data_for_plotting(data)
    
    # MODIFIED HERE: Updated method labels for the legend with LaTeX notation
    methods = [r'$\hat{r}$']
    for theta in theta_values:
        if theta == 0.1:
            methods.append(r'$\hat{r}$ using $\theta = 0.1$')
        elif theta == 0.01:
            methods.append(r'$\hat{r}$ using $\theta = 0.01$')
        else:
            methods.append(f'$\\hat{{r}}$ using $\\theta = {theta}$')
    
    # Custom color palette
    colors = ['#FF0000', '#0000FF', '#FFD700', '#00FF00', '#800080']  # Add more colors if needed
    custom_palette = {method: color for method, color in zip(methods, colors)}
    
    # Calculate non-outlier range
    y_min_non_outlier, y_max_non_outlier = calculate_non_outlier_range(melted_data, true_r_values, methods)
    
    # Create figure with a single panel
    fig, ax = plt.subplots(figsize=(24, 15))
    
    # Create the boxplot without outliers (showfliers=False)
    sns.set(style="whitegrid", font_scale=1.5)  # Increased font scale
    sns_plot = sns.boxplot(x='true_r', y='estimated_r', hue='method', data=melted_data,
                     palette=custom_palette,
                     width=0.9,
                     linewidth=1.0,
                     showfliers=False,  # Remove outliers
                     ax=ax)
    
    # Add a line showing the true r values (baseline)
    x_positions = range(len(true_r_values))
    ax.plot(x_positions, true_r_values, 'o-', color='gray', alpha=0.8, markersize=8, 
            linestyle='--', linewidth=2.0, label='True r')
    
    # Set x-axis ticks - show only every third tick
    visible_ticks = []
    tick_labels = []
    
    for i, r in enumerate(true_r_values):
        if i % 3 == 0:  # Show every 3rd tick (skip 2)
            visible_ticks.append(i)
            tick_labels.append(f"{r:.3f}")
    
    ax.set_xticks(visible_ticks)
    ax.set_xticklabels(tick_labels, rotation=45, fontsize=55)  # Increased tick font size
    
    # Set simplified axis labels with larger font
    ax.set_xlabel('r', fontsize=60)  # Changed from 'True r Value' to 'r'
    ax.set_ylabel('Estimated r', fontsize=60)  # Changed from 'Estimated r Value' to 'estimated r'
    
    # Increase y-axis tick font size
    ax.tick_params(axis='y', labelsize=55)
    
    # Make y-axis ticks more sparse (around 5 ticks)
    # This code was moved to the y-axis limits section
    
    # Add grid for easier reading
    ax.grid(True, axis='y', linestyle='--', alpha=0.7)
    
    # Handle legend with larger font
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, fontsize=40, frameon=True, framealpha=0.9, loc='upper left')
    
    # Set y-axis limits based on non-outlier range
    if y_min_non_outlier is not None and y_max_non_outlier is not None:
        # Add a small padding
        padding = (y_max_non_outlier - y_min_non_outlier) * 0.05
        y_min = max(0, y_min_non_outlier - padding)  # Ensure y_min is not negative
        y_max = y_max_non_outlier + padding
        
        # Set the limits
        ax.set_ylim(y_min, y_max)
        print(f"Setting y-axis range to: {y_min:.6f} - {y_max:.6f} (based on non-outlier data)")
        
        # Create nice, round ticks (0.05, 0.10, etc.)
        # Find appropriate step size to get about 5-6 ticks
        range_size = y_max - y_min
        
        # Determine the appropriate step size for nice round numbers
        if range_size <= 0.25:
            step = 0.05  # For smaller ranges: 0.05, 0.10, 0.15...
        elif range_size <= 0.5:
            step = 0.1   # 0.1, 0.2, 0.3...
        elif range_size <= 1.0:
            step = 0.2   # 0.2, 0.4, 0.6...
        elif range_size <= 2.0:
            step = 0.5   # 0.5, 1.0, 1.5...
        else:
            step = 1.0   # 1.0, 2.0, 3.0...
        
        # Create round-numbered ticks
        # Start from a nice round number
        start = np.floor(y_min / step) * step
        y_ticks = np.arange(start, y_max + step, step)
        
        # Filter out any ticks outside our actual range
        y_ticks = y_ticks[(y_ticks >= y_min) & (y_ticks <= y_max)]
        
        ax.set_yticks(y_ticks)
        
        # Format y-tick labels
        if step >= 1.0:
            ax.set_yticklabels([f"{y:.0f}" for y in y_ticks])
        elif step >= 0.1:
            ax.set_yticklabels([f"{y:.1f}" for y in y_ticks])
        else:
            ax.set_yticklabels([f"{y:.2f}" for y in y_ticks])
    else:
        # Fallback to a reasonable range if calculation failed
        y_min = 0
        y_max = max(true_r_values) * 2
        ax.set_ylim(y_min, y_max)
        print(f"Using fallback y-axis range: {y_min} - {y_max}")
        
        # Create nice round ticks even in fallback mode
        range_size = y_max - y_min
        
        # Determine step size for nice round numbers
        if range_size <= 0.25:
            step = 0.05
        elif range_size <= 0.5:
            step = 0.1
        elif range_size <= 1.0:
            step = 0.2
        elif range_size <= 2.0:
            step = 0.5
        else:
            step = 1.0
            
        # Create round-numbered ticks
        start = np.floor(y_min / step) * step
        y_ticks = np.arange(start, y_max + step, step)
        y_ticks = y_ticks[(y_ticks >= y_min) & (y_ticks <= y_max)]
        
        ax.set_yticks(y_ticks)
        
        # Format y-tick labels based on step size
        if step >= 1.0:
            ax.set_yticklabels([f"{y:.0f}" for y in y_ticks])
        elif step >= 0.1:
            ax.set_yticklabels([f"{y:.1f}" for y in y_ticks])
        else:
            ax.set_yticklabels([f"{y:.2f}" for y in y_ticks])
    
    # Adjust layout
    plt.tight_layout()
    
    # Save the plot
    output_filename = f"single_panel_sketch_r_k{k}.png"
    plt.savefig(output_filename, dpi=100)
    print(f"Plot saved as '{output_filename}'")

# Main execution
if __name__ == "__main__":
    main()