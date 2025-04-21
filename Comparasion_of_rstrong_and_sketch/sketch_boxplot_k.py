import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from glob import glob
import argparse

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Sketch K-mer Variation Boxplot')
    parser.add_argument('-d', '--directory', type=str, default='./Estimate_sketch_r',
                        help='Input directory with sketch output files')
    parser.add_argument('-t', '--thresholds', nargs='+', type=float, default=[32, 570], 
                        help='Thresholds for splitting k values [small_max, medium_max]')
    parser.add_argument('--ylim', nargs='+', type=float, default=[0.8, 0.035, 1.0],
                        help='Maximum y-axis limits for the three panels [small_k, medium_k, large_k]')
    parser.add_argument('--theta', nargs='+', type=float, default=[0.1, 0.01],
                        help='Theta values to display in the plot, in the order they should appear')
    parser.add_argument('-r', '--fixed_r', type=float, default=0.01,
                        help='Fixed r value used in the filenames')
    args = parser.parse_args()
    
    # Process data files
    all_data = process_all_files(args.directory, args.fixed_r)
    
    if all_data is not None:
        print(f"Processed data for {len(all_data['k_value'].unique())} different k values")
        
        # Create the triple subplot boxplot with specified parameters
        create_triple_subplot_boxplot(
            all_data, 
            selected_theta_values=args.theta,
            thresholds=args.thresholds,
            ylimits=args.ylim,
            fixed_r=args.fixed_r
        )
    else:
        print("No data to visualize")

# Function to extract k value from filename
def extract_k_value(filename):
    try:
        # Extract k value from a filename like "k30_r0.01.output"
        parts = filename.split('_')
        k_part = parts[0]
        k_value = int(k_part[1:])
        
        return k_value
    except (ValueError, IndexError) as e:
        print(f"Error extracting k value from filename {filename}: {e}")
        return None

# Function to process a single file
def process_file(file_path, fixed_r):
    try:
        # Extract k value from filename
        k_value = extract_k_value(os.path.basename(file_path))
        if k_value is None:
            return pd.DataFrame()
        
        # Print debugging info
        print(f"Processing file: {file_path}, extracted k: {k_value}")
        
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
        
        # Add metadata to the DataFrame
        df['k_value'] = k_value
        df['r_value'] = fixed_r
        
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
        return pd.DataFrame(columns=['r_strong', 'r_sketch', 'k_value', 'r_value', 'theta'])

# Function to process all files
def process_all_files(directory='./Estimate_sketch_r', fixed_r=0.01):
   # Try multiple directories
    directories = [directory, '.', './output', '../output']
    
    file_list = []
    for dir_path in directories:
        pattern = os.path.join(dir_path, f"k*_r{fixed_r}.output")
        current_files = glob(pattern)
        if current_files:
            print(f"Found {len(current_files)} files in {dir_path}")
            file_list = current_files
            break
    
    if not file_list:
        print(f"No files found matching the pattern k*_r{fixed_r}.output in any directory")
        return None
    
    print(f"Found {len(file_list)} files matching the pattern")
    
    # Process each file and combine the results
    dataframes = []
    for file in file_list:
        df = process_file(file, fixed_r)
        if not df.empty:
            dataframes.append(df)
    
    if not dataframes:
        print("No valid data found in any of the files")
        return None
        
    all_data = pd.concat(dataframes, ignore_index=True)
    
    return all_data

# Function to create boxplots with three subplots (small, medium, large k values)
def create_triple_subplot_boxplot(data, selected_theta_values=[0.1, 0.01], thresholds=[32, 570], ylimits=[0.8, 0.035, 1.0], fixed_r=0.01):
    # Ensure we have two thresholds
    if len(thresholds) < 2:
        thresholds = [32, 570]  # Default values
    
    # Ensure we have three y-limits
    if len(ylimits) < 3:
        ylimits = ylimits + [1.0] * (3 - len(ylimits))
    
    # Filter data for selected theta values
    filtered_data = data[data['theta'].isin(selected_theta_values)]
    
    # Convert k_value to numeric and sort
    filtered_data['k_value'] = pd.to_numeric(filtered_data['k_value'])
    
    # Get unique k values
    k_values = sorted(filtered_data['k_value'].unique())
    print(f"Unique k values: {k_values}")
    
    if len(k_values) == 0:
        print("No data available for the selected theta values")
        return
    
    # Split data into three groups based on provided thresholds
    small_k_data = filtered_data[filtered_data['k_value'] <= thresholds[0]]
    medium_k_data = filtered_data[(filtered_data['k_value'] > thresholds[0]) & (filtered_data['k_value'] <= thresholds[1])]
    large_k_data = filtered_data[filtered_data['k_value'] > thresholds[1]]
    
    # Get unique k values for each group
    small_k_values = sorted(small_k_data['k_value'].unique())
    medium_k_values = sorted(medium_k_data['k_value'].unique())
    large_k_values = sorted(large_k_data['k_value'].unique())
    
    print(f"Small k values (≤ {thresholds[0]}): {small_k_values}")
    print(f"Medium k values ({thresholds[0]} < k ≤ {thresholds[1]}): {medium_k_values}")
    print(f"Large k values (> {thresholds[1]}): {large_k_values}")
    
    # Create melted data for plotting
    melted_data = []
    
    # Process each theta value in the specified order
    # For each k value, add r_strong first, then r_sketch with theta values in the specified order
    for k in k_values:
        # Add r_strong data (only once per k)
        k_data = filtered_data[(filtered_data['k_value'] == k) & (filtered_data['theta'] == selected_theta_values[0])]
        for _, row in k_data.iterrows():
            melted_data.append({
                'k_value': k,
                'estimated_r': row['r_strong'],
                'method': 'r_strong'
            })
        
        # Add r_sketch data for each theta value in the specified order
        for theta in selected_theta_values:
            k_theta_data = filtered_data[(filtered_data['k_value'] == k) & (filtered_data['theta'] == theta)]
            for _, row in k_theta_data.iterrows():
                melted_data.append({
                    'k_value': k,
                    'estimated_r': row['r_sketch'],
                    'method': f'r_sketch (θ={theta})'
                })
    
    melted_df = pd.DataFrame(melted_data)
    
    # Custom color palette - Use the order r_strong, r_sketch (θ=value1), r_sketch (θ=value2)
    methods = ['r_strong', f'r_sketch (θ={selected_theta_values[0]})', f'r_sketch (θ={selected_theta_values[1]})']
    colors = ['#FF0000', '#0000FF', '#FFD700']  # red, blue, yellow (gold)
    custom_palette = {method: color for method, color in zip(methods, colors)}
    
    # Reorder the 'method' column to ensure the desired order in the legend
    method_order = ['r_strong', f'r_sketch (θ={selected_theta_values[0]})', f'r_sketch (θ={selected_theta_values[1]})']
    melted_df['method'] = pd.Categorical(melted_df['method'], categories=method_order, ordered=True)
    
    # Create figure with three subplots
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(24, 10), 
                                       gridspec_kw={'width_ratios': [1, 2, 1]})
    
    # Function to plot on a specific axis
    def plot_on_axis(ax, data, k_values, title, ylim=None):
        if len(data) == 0 or len(k_values) == 0:
            ax.text(0.5, 0.5, 'No data available', 
                    horizontalalignment='center', verticalalignment='center')
            return None
        
        # Create the boxplot
        sns.set(style="whitegrid", font_scale=1.2)
        sns_plot = sns.boxplot(x='k_value', y='estimated_r', hue='method', 
                         data=data,
                         palette=custom_palette,
                         hue_order=method_order,  # Specify the order of methods in legend
                         width=0.8,
                         linewidth=1.0,
                         fliersize=5,
                         ax=ax)
        
        # Add a line showing the true r value
        ax.axhline(y=fixed_r, color='gray', alpha=0.8, linestyle='--', 
                  linewidth=2.0, label=f'True r={fixed_r}')
        
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
        
        return sns_plot
    
    # Plot on each subplot with their respective ylimits
    plot1 = plot_on_axis(ax1, 
                        melted_df[melted_df['k_value'] <= thresholds[0]], 
                        small_k_values, 
                        f'k ≤ {thresholds[0]}',
                        ylim=ylimits[0])
    
    plot2 = plot_on_axis(ax2, 
                        melted_df[(melted_df['k_value'] > thresholds[0]) & (melted_df['k_value'] <= thresholds[1])], 
                        medium_k_values, 
                        f'{thresholds[0]} < k ≤ {thresholds[1]}',
                        ylim=ylimits[1])
    
    plot3 = plot_on_axis(ax3, 
                        melted_df[melted_df['k_value'] > thresholds[1]], 
                        large_k_values, 
                        f'k > {thresholds[1]}',
                        ylim=ylimits[2])
    
    # Handle legends
    for ax in [ax2, ax3]:
        if ax.get_legend():
            ax.get_legend().remove()
    
    if ax1.get_legend():
        # Ensure legend order is correct
        handles, labels = ax1.get_legend_handles_labels()
        # Reorder to match method_order plus the true_r line
        ordered_indices = []
        for method in method_order:
            if method in labels:
                ordered_indices.append(labels.index(method))
        # Add the true_r line at the end
        if f'True r={fixed_r}' in labels:
            ordered_indices.append(labels.index(f'True r={fixed_r}'))
        
        # Create new legend with correct order
        if ordered_indices:
            ordered_handles = [handles[i] for i in ordered_indices]
            ordered_labels = [labels[i] for i in ordered_indices]
            ax1.legend(ordered_handles, ordered_labels, fontsize=14, frameon=True, framealpha=0.9, loc='best')
    
    # Add y-axis label only to the first subplot
    ax1.set_ylabel('Estimated r Value', fontsize=18)
    
    # Add overall title
    plt.suptitle(f'Comparison of r_strong and r_sketch with different θ values (r={fixed_r})', fontsize=24)
    
    # Adjust layout
    plt.tight_layout()
    plt.subplots_adjust(top=0.9)  # Make room for the suptitle
    
    # Save the plot
    plt.savefig(f"boxplot_triple_k_comparison_r{fixed_r}.png", dpi=300)
    print(f"Plot saved as 'boxplot_triple_k_comparison_r{fixed_r}.png'")

# Main execution
if __name__ == "__main__":
    main()