import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from glob import glob
import argparse

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='K-mer Variation Boxplot')
    parser.add_argument('-d', '--directory', type=str, default='./Estimate_sketch_r_mut1',
                        help='Input directory with sketch output files')
    parser.add_argument('-t', '--thresholds', nargs='+', type=int, default=[32, 570], 
                        help='Thresholds for splitting k values [small_max, medium_max]')
    parser.add_argument('--ylim', nargs='+', type=float, default=[0.8, 0.035, 1.0],
                        help='Y-axis limits for the three panels [small_k, medium_k, large_k]')
    # Add the parameters that are being passed in the command
    parser.add_argument('-k', nargs='+', type=int, 
                        help='Thresholds for splitting k values (alternative to -t)')
    parser.add_argument('--theta', nargs='+', type=float, 
                        help='Theta values for r_sketch calculation')
    parser.add_argument('-r', '--fixed_r', type=float,
                        help='Fixed r value to use for filename pattern and visualization')
    args = parser.parse_args()
    
    # Use -k values if provided, otherwise use -t thresholds
    if args.k and len(args.k) >= 2:
        args.thresholds = args.k[:2]  # Take the first two values if more are provided
        print(f"Using k thresholds from -k argument: {args.thresholds}")
    
    # Ensure we have exactly 2 thresholds
    if len(args.thresholds) != 2:
        print(f"Warning: Expected 2 threshold values, got {len(args.thresholds)}. Using default values [32, 570].")
        args.thresholds = [32, 570]
    
    # Ensure we have exactly 3 y-limits
    if len(args.ylim) != 3:
        print(f"Warning: Expected 3 y-limit values, got {len(args.ylim)}. Using default values [0.8, 0.035, 1.0].")
        args.ylim = [0.8, 0.035, 1.0]
    
    # Determine which r values to process
    if args.fixed_r is not None:
        # If -r argument was provided, use that specific r value
        r_values = [args.fixed_r]
        print(f"Using fixed r value from command line: {args.fixed_r}")
    else:
        # Otherwise, find all r values from filenames
        r_values = find_r_values(args.directory)
        if not r_values:
            print("No valid r values found in files.")
            return
        print(f"Found r values from filenames: {r_values}")
    
    # Get theta values if provided
    theta_values = args.theta if args.theta else None
    if theta_values:
        print(f"Using custom theta values: {theta_values}")
    
    # Process each r value
    for r_value in r_values:
        print(f"\nProcessing data files for r = {r_value}")
        all_data = process_all_files(args.directory, r_value, theta_values)
        
        if all_data is not None and not all_data.empty:
            print(f"Processed data for {len(all_data['k_value'].unique())} different k values with r = {r_value}")
            
            # Create the triple subplot boxplot with the specified parameters
            create_triple_subplot_boxplot(all_data, args.thresholds, args.ylim, r_value, theta_values)
        else:
            print(f"No data to visualize for r = {r_value}")

# Function to find all r values from filenames in the directory
def find_r_values(directory):
    # Get all output files
    file_list = glob(os.path.join(directory, "k*_r*.output"))
    
    if not file_list:
        print(f"No files found matching the pattern k*_r*.output in {directory}")
        return []
    
    # Extract r values from filenames
    r_values = set()
    for file in file_list:
        r_value = extract_r_from_filename(os.path.basename(file))
        if r_value is not None:
            r_values.add(r_value)
    
    return sorted(r_values)

# Function to extract k value from filename
def extract_k_from_filename(filename):
    try:
        # Extract the k value from a filename like "k30_r0.01.output"
        parts = filename.split('_')
        k_part = parts[0]
        k_value = int(k_part[1:])
        
        return k_value
    except (ValueError, IndexError) as e:
        print(f"Error extracting k value from filename {filename}: {e}")
        return None

# Function to extract r value from filename
def extract_r_from_filename(filename):
    try:
        # Extract the r value from a filename like "k30_r0.01.output"
        parts = filename.split('_')
        for part in parts:
            if part.startswith('r') and part.endswith('.output'):
                # Extract just the numerical part without '.output'
                r_value = float(part[1:-7])
                return r_value
            elif part.startswith('r'):
                # If the format is different
                r_value = float(part[1:])
                return r_value
        return None
    except (ValueError, IndexError) as e:
        print(f"Error extracting r value from filename {filename}: {e}")
        return None

# Function to process a single file
def process_file(file_path, theta_values=None):
    try:
        # Extract k and r values from filename
        filename = os.path.basename(file_path)
        k_value = extract_k_from_filename(filename)
        r_value = extract_r_from_filename(filename)
        
        if k_value is None or r_value is None:
            return pd.DataFrame()
        
        # Print debugging info
        print(f"Processing file: {file_path}, extracted k value: {k_value}, r value: {r_value}")
        
        # Try to read the CSV file
        try:
            # The file has three columns: r_strong, r_sketch (theta=0.1), r_sketch (theta=0.01)
            # If custom theta values are provided, adjust column names
            if theta_values is not None and len(theta_values) >= 2:
                column_names = ['r_strong', 
                               f'r_sketch_theta_{theta_values[0]}', 
                               f'r_sketch_theta_{theta_values[1]}']
            else:
                column_names = ['r_strong', 'r_sketch_theta_0.1', 'r_sketch_theta_0.01']
                
            df = pd.read_csv(file_path, header=None, 
                           names=column_names,
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
                if theta_values is not None and len(theta_values) >= 2:
                    column_names = ['r_strong', 
                                   f'r_sketch_theta_{theta_values[0]}', 
                                   f'r_sketch_theta_{theta_values[1]}']
                else:
                    column_names = ['r_strong', 'r_sketch_theta_0.1', 'r_sketch_theta_0.01']
                    
                df = pd.DataFrame(rows, columns=column_names)
            else:
                raise ValueError("No valid rows found after manual parsing")
        
        # Add the k value and r value to the DataFrame
        df['k_value'] = k_value
        df['r_value'] = r_value
        
        # Print a sample of the processed data
        print(f"Sample of processed data: {df.head(2)}")
        
        return df
        
    except Exception as e:
        print(f"Error processing file {file_path}: {e}")
        # Return an empty DataFrame with the correct structure
        return pd.DataFrame(columns=['r_strong', 'r_sketch_theta_0.1', 'r_sketch_theta_0.01', 'k_value', 'r_value'])

# Function to process all files with a specific r value
def process_all_files(directory='./Estimate_sketch_r_mut1', r_value=0.01, theta_values=None):
    # Define the pattern for files
    pattern = os.path.join(directory, f"k*_r{r_value}.output")
    file_list = glob(pattern)
    
    if not file_list:
        print(f"No files found matching the pattern k*_r{r_value}.output in {directory}")
        return None
    
    print(f"Found {len(file_list)} files matching the pattern for r = {r_value}")
    
    # Process each file and combine the results
    dataframes = []
    for file in file_list:
        df = process_file(file, theta_values)
        if not df.empty:
            dataframes.append(df)
    
    if not dataframes:
        print(f"No valid data found in any of the files for r = {r_value}")
        return None
        
    all_data = pd.concat(dataframes, ignore_index=True)
    
    return all_data

# Function to create boxplots with three subplots (small, medium, large k values)
def create_triple_subplot_boxplot(data, thresholds=[32, 570], ylimits=[0.8, 0.035, 1.0], r_value=0.01, theta_values=None):
    # Convert k_value to numeric and sort
    data['k_value'] = pd.to_numeric(data['k_value'])
    
    # Get unique k values
    k_values = sorted(data['k_value'].unique())
    print(f"Unique k values: {k_values}")
    
    if len(k_values) == 0:
        print("No data available")
        return
    
    # Unpack thresholds
    small_threshold, medium_threshold = thresholds
    
    # Split data into three groups based on the thresholds
    small_k_data = data[data['k_value'] <= small_threshold]
    medium_k_data = data[(data['k_value'] > small_threshold) & (data['k_value'] <= medium_threshold)]
    large_k_data = data[data['k_value'] > medium_threshold]
    
    # Get unique k values for each group
    small_k_values = sorted(small_k_data['k_value'].unique())
    medium_k_values = sorted(medium_k_data['k_value'].unique())
    large_k_values = sorted(large_k_data['k_value'].unique())
    
    print(f"Small k values (≤ {small_threshold}): {small_k_values}")
    print(f"Medium k values ({small_threshold} < k ≤ {medium_threshold}): {medium_k_values}")
    print(f"Large k values (> {medium_threshold}): {large_k_values}")
    
    # Create melted data for plotting
    melted_data = []
    
    # Use the provided theta values or default to 0.1 and 0.01
    if theta_values and len(theta_values) >= 2:
        theta1 = theta_values[0]
        theta2 = theta_values[1]
        col1 = f'r_sketch_theta_{theta1}'
        col2 = f'r_sketch_theta_{theta2}'
        method1 = f'r_sketch (θ={theta1})'
        method2 = f'r_sketch (θ={theta2})'
    else:
        theta1 = 0.1
        theta2 = 0.01
        col1 = 'r_sketch_theta_0.1'
        col2 = 'r_sketch_theta_0.01'
        method1 = 'r_sketch (θ=0.1)'
        method2 = 'r_sketch (θ=0.01)'
    
    # Add r_sketch data for each theta value
    for k in k_values:
        k_data = data[data['k_value'] == k]
        
        # Get r_strong value (same for all rows with the same k_value)
        r_strong = k_data['r_strong'].iloc[0]
        
        for _, row in k_data.iterrows():
            # Add r_sketch for first theta value
            melted_data.append({
                'k_value': k,
                'estimated_r': row[col1],
                'method': method1,
                'r_strong': r_strong
            })
            
            # Add r_sketch for second theta value
            melted_data.append({
                'k_value': k,
                'estimated_r': row[col2],
                'method': method2,
                'r_strong': r_strong
            })
    
    melted_df = pd.DataFrame(melted_data)
    
    # Custom color palette for r_sketch (first and second theta values)
    custom_palette = {
        method1: '#0000FF',  # Blue
        method2: '#FFD700'   # Yellow
    }
    
    # Create figure with three subplots
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(24, 10), 
                                       gridspec_kw={'width_ratios': [1, 2, 1]})
    
    # Unpack y-limits
    small_ylim, medium_ylim, large_ylim = ylimits
    
    # Function to plot on a specific axis
    def plot_on_axis(ax, data, k_values, title, ylim=None):
        if len(data) == 0 or len(k_values) == 0:
            ax.text(0.5, 0.5, 'No data available', 
                    horizontalalignment='center', verticalalignment='center')
            return None, None, None
        
        # Create the boxplot for r_sketch values
        sns.set(style="whitegrid", font_scale=1.2)
        sns_plot = sns.boxplot(x='k_value', y='estimated_r', hue='method', data=data,
                         palette=custom_palette,
                         width=0.8,
                         linewidth=1.0,
                         fliersize=5,
                         ax=ax)
        
        # Add r_strong values as horizontal lines across each k_value group
        x_positions = range(len(k_values))
        for i, k in enumerate(k_values):
            # Get the r_strong value for this k_value
            k_data = data[data['k_value'] == k]
            if not k_data.empty:
                r_strong = k_data['r_strong'].iloc[0]
                
                # Calculate the width of each group (for the horizontal line)
                width = 0.4  # Half the default width of boxplot groups
                
                # Draw horizontal line for r_strong (in red)
                ax.plot([i-width, i+width], [r_strong, r_strong], '-', color='#FF0000', 
                        linewidth=2.0, solid_capstyle='round')
        
        # Add a line showing the true r value (r_value)
        ax.axhline(y=r_value, color='gray', alpha=0.8, linestyle='--', 
                  linewidth=2.0)
        
        # Set x-axis ticks to show the k values
        if len(k_values) > 10:
            # Show subset of k values to avoid crowding
            step = max(1, len(k_values) // 10)
            shown_k = k_values[::step]
            ax.set_xticks([k_values.index(k) for k in shown_k])
            ax.set_xticklabels([str(k) for k in shown_k], rotation=45)
        else:
            ax.set_xticks(range(len(k_values)))
            ax.set_xticklabels([str(k) for k in k_values], rotation=45)
        
        # Set title and labels
        ax.set_title(title, fontsize=20)
        ax.set_xlabel('k Value', fontsize=18)
        
        # Set y-axis limits if provided
        if ylim is not None:
            ax.set_ylim(0, ylim)
        
        # Add grid for easier reading
        ax.grid(True, axis='y', linestyle='--', alpha=0.7)
        
        # Add manually created legend entries
        from matplotlib.lines import Line2D
        custom_lines = [
            Line2D([0], [0], color='#FF0000', lw=2),
            Line2D([0], [0], color='gray', linestyle='--', lw=2)
        ]
        custom_labels = ['r_strong', f'True r={r_value}']
        
        # Get the existing legend handles and labels
        handles, labels = ax.get_legend_handles_labels()
        
        # Combine with the custom legend entries
        all_handles = handles + custom_lines
        all_labels = labels + custom_labels
        
        return sns_plot, all_handles, all_labels
    
    # Plot on each subplot
    plot1, handles1, labels1 = plot_on_axis(ax1, 
                                          melted_df[melted_df['k_value'] <= small_threshold], 
                                          small_k_values, 
                                          f'k ≤ {small_threshold}',
                                          ylim=small_ylim)
    
    plot2, handles2, labels2 = plot_on_axis(ax2, 
                                          melted_df[(melted_df['k_value'] > small_threshold) & 
                                                  (melted_df['k_value'] <= medium_threshold)], 
                                          medium_k_values, 
                                          f'{small_threshold} < k ≤ {medium_threshold}',
                                          ylim=medium_ylim)
    
    plot3, handles3, labels3 = plot_on_axis(ax3, 
                                          melted_df[melted_df['k_value'] > medium_threshold], 
                                          large_k_values, 
                                          f'k > {medium_threshold}',
                                          ylim=large_ylim)
    
    # Remove all legends first
    for ax in [ax1, ax2, ax3]:
        if ax.get_legend():
            ax.get_legend().remove()
    
    # Create a legend on the middle panel (upper right corner)
    # Use handles/labels from whichever panel has data
    legend_handles = None
    legend_labels = None
    
    if handles2 is not None:  # Prefer middle panel if it has data
        legend_handles = handles2
        legend_labels = labels2
    elif handles1 is not None:  # Use left panel if middle has no data
        legend_handles = handles1
        legend_labels = labels1
    elif handles3 is not None:  # Use right panel if others have no data
        legend_handles = handles3
        legend_labels = labels3
        
    if legend_handles is not None:
        ax2.legend(handles=legend_handles, labels=legend_labels, 
                  fontsize=16, frameon=True, framealpha=0.9, 
                  loc='upper right')
    
    # Add y-axis label only to the first subplot
    ax1.set_ylabel('Estimated r Value', fontsize=18)
    
    # Add overall title
    plt.suptitle(f'Comparison of r_strong and r_sketch with different θ values (r={r_value})', fontsize=24)
    
    # Adjust layout
    plt.tight_layout()
    plt.subplots_adjust(top=0.9)  # Make room for the suptitle
    
    # Save the plot
    plt.savefig(f"boxplot_triple_k_comparison_r{r_value}.png", dpi=300)
    print(f"Plot saved as 'boxplot_triple_k_comparison_r{r_value}.png'")

# Main execution
if __name__ == "__main__":
    main()