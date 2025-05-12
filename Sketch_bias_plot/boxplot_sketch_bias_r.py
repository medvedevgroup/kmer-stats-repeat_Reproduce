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
    args = parser.parse_args()
    
    # Process all the data files
    all_data = process_all_files(args.directory, args.k_value)
    
    if all_data is not None:
        print(f"Processed data for {len(all_data['true_r'].unique())} different r values")
        
        # Create the single panel boxplot without outliers
        create_single_panel_boxplot(all_data, args.k_value)
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
    
    # Get all r_strong values as well
    r_strong_vals = data['r_strong'].unique()
    max_vals.extend(r_strong_vals)
    
    # Return min and max of all non-outlier values
    if min_vals and max_vals:
        return min(min_vals), max(max_vals)
    else:
        return None, None

# Function to create a single panel boxplot
def create_single_panel_boxplot(data, k_value=30):
    # Prepare data for plotting
    melted_data, true_r_values = prepare_data_for_plotting(data)
    
    # Custom color palette for r_sketch (theta = 0.1) and r_sketch (theta = 0.01)
    custom_palette = {
        'r_sketch (θ=0.1)': '#0000FF',  # Blue
        'r_sketch (θ=0.01)': '#FFD700'   # Yellow
    }
    
    # Calculate non-outlier range
    methods = ['r_sketch (θ=0.1)', 'r_sketch (θ=0.01)']
    y_min_non_outlier, y_max_non_outlier = calculate_non_outlier_range(melted_data, true_r_values, methods)
    
    # Create figure with a single panel
    fig, ax = plt.subplots(figsize=(24, 15))
    
    # Create the boxplot for r_sketch values without outliers
    sns.set(style="whitegrid", font_scale=1.5)  # Increased font scale
    sns_plot = sns.boxplot(x='true_r', y='estimated_r', hue='method', data=melted_data,
                     palette=custom_palette,
                     width=0.9,
                     linewidth=1.0,
                     showfliers=False,  # Remove outliers
                     ax=ax)
    
    # Add r_strong values as horizontal lines across each true_r group
    x_positions = range(len(true_r_values))
    for i, r in enumerate(true_r_values):
        # Get the r_strong value for this true_r
        r_data = melted_data[melted_data['true_r'] == r]
        if not r_data.empty:
            r_strong = r_data['r_strong'].iloc[0]
            
            # Calculate the width of each group (for the horizontal line)
            width = 0.4  # Half the default width of boxplot groups
            
            # Draw horizontal line for r_strong (in red)
            ax.plot([i-width, i+width], [r_strong, r_strong], '-', color='#FF0000', 
                    linewidth=2.0, solid_capstyle='round')
    
    # Set x-axis ticks - show only every 3rd tick
    visible_ticks = []
    tick_labels = []
    
    for i, r in enumerate(true_r_values):
        if i % 3 == 0:  # Every 3rd tick (skip 2)
            visible_ticks.append(i)
            tick_labels.append(f"{r:.3f}")
    
    ax.set_xticks(visible_ticks)
    ax.set_xticklabels(tick_labels, rotation=45, fontsize=60)  # Increased tick label font size
    
    # Set simplified axis labels with larger font
    ax.set_xlabel('r', fontsize=60)  # Changed from 'True r Value' to 'r'
    ax.set_ylabel('Estimated r', fontsize=60)  # Changed from 'Estimated r Value' to 'estimated r'
    
    # Set y-axis limits based on non-outlier range
    if y_min_non_outlier is not None and y_max_non_outlier is not None:
        # Add a small padding
        padding = (y_max_non_outlier - y_min_non_outlier) * 0.05
        y_min = max(0, y_min_non_outlier - padding)  # Ensure y_min is not negative
        y_max = y_max_non_outlier + padding
        
        # Set the limits
        ax.set_ylim(y_min, y_max)
        print(f"Setting y-axis range to: {y_min:.6f} - {y_max:.6f} (based on non-outlier data)")
    else:
        # Fallback to a reasonable range
        y_min = 0
        y_max = max(true_r_values) * 2
        ax.set_ylim(y_min, y_max)
        print(f"Using fallback y-axis range: {y_min} - {y_max}")
    
    # Reduce Y-axis ticks density - show fewer ticks
    # First, get the current ticks
    current_yticks = ax.get_yticks()
    
    # Only keep every 3rd tick to make them more sparse
    reduced_yticks = current_yticks[::3]
    
    # Apply the reduced ticks
    ax.set_yticks(reduced_yticks)
    
    # Increase y-tick font size
    ax.tick_params(axis='y', labelsize=60)
    
    # Add grid for easier reading
    ax.grid(True, axis='y', linestyle='--', alpha=0.7)
    
    # Create a manually crafted legend with just r_strong (removed True r)
    from matplotlib.lines import Line2D
    custom_lines = [
        Line2D([0], [0], color='#FF0000', lw=2)
    ]
    custom_labels = ['$\hat{r}$']
    
    # Get the existing legend handles and labels
    handles, labels = ax.get_legend_handles_labels()
    
    # Combine with the custom legend entries
    all_handles = handles + custom_lines
    all_labels = labels + custom_labels
    
    # Add the combined legend with larger font
    ax.legend(handles=all_handles, labels=all_labels, fontsize=45, frameon=True, framealpha=0.9, loc='upper left')
    
    # Adjust layout
    plt.tight_layout()
    
    # Save the plot
    output_filename = f"single_panel_sketch_r_k{k_value}.png"
    plt.savefig(output_filename, dpi=100)
    print(f"Plot saved as '{output_filename}'")

# Main execution
if __name__ == "__main__":
    main()