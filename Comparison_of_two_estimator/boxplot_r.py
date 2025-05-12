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
        create_boxplots(
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
            df = pd.read_csv(file_path, header=None, names=['r_hat', 'r_Mash'],
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
                df = pd.DataFrame(rows, columns=['r_hat', 'r_Mash'])
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
        return pd.DataFrame(columns=['r_hat', 'r_Mash', 'true_r'])

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

# Function to plot boxplots on a specific axis
def plot_on_axis(ax, data, r_values, is_right_panel=False):
    # Melt the dataframe to prepare for boxplot
    melted_data = pd.melt(data, 
                        id_vars=['true_r'], 
                        value_vars=['r_hat', 'r_Mash'],
                        var_name='method', 
                        value_name='estimated_r')
    
    # Custom color palette with high contrast colors
    custom_palette = {"r_hat": "#FF0000", "r_Mash": "#FFD700"}
    
    # Create the boxplot with increased linewidth
    sns.set(style="whitegrid", font_scale=1.5)  # Increased font scale for all text
    sns_plot = sns.boxplot(x='true_r', y='estimated_r', hue='method', data=melted_data,
                     palette=custom_palette,
                     width=0.9,
                     linewidth=2.0,  # Increased linewidth from 1.0 to 2.0
                     fliersize=5,
                     ax=ax)
    
    # Add a line showing the true r values (baseline)
    for i, r in enumerate(r_values):
        ax.plot(i, r, 'o', color='gray', alpha=0.8, markersize=8)
    
    # Connect the true r values with a line
    ax.plot(range(len(r_values)), r_values, color='gray', alpha=0.6, linestyle='--', 
            linewidth=2.0, label='True r')
    
    # Set x-axis ticks to show only every third r value
    ax.set_xticks(range(len(r_values)))
    
    # Only label every 3rd tick
    x_labels = []
    for i, r in enumerate(r_values):
        if i % 3 == 0:  # Change this to control spacing (every 3rd tick)
            x_labels.append(f"{r:.3f}")
        else:
            x_labels.append("")
    
    ax.set_xticklabels(x_labels, rotation=45, fontsize=60)  # Added bold
    
    # Set labels with larger font sizes and bold
    ax.set_xlabel('r', fontsize=60)  # Added bold
    
    # Only set y-label for left panel
    if not is_right_panel:
        ax.set_ylabel('Estimated r', fontsize=60)  # Added bold
    else:
        # For right panel, remove the y-label completely
        ax.set_ylabel('')
        # Move y-ticks to the right side
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
    
    # Add grid for easier reading - made grid lines bolder
    ax.grid(True, axis='y', linestyle='--', alpha=0.7, linewidth=1.5)
    
    # Handle the case where we need to adjust y-axis limits
    if len(data) > 0:
        ymin = data[['r_hat', 'r_Mash']].min().min() * 0.9
        ymax = data[['r_hat', 'r_Mash']].max().max() * 1.1
        
        # Ensure the true r values are visible
        ymin = min(ymin, min(r_values) * 0.9)
        ymax = max(ymax, max(r_values) * 1.1)
        
        # Set y-axis limits
        ax.set_ylim(max(0, ymin), ymax)
        
        # Set y-ticks to use clean, round numbers instead of exact values
        # Calculate a nice round step size based on the data range
        data_range = ymax - max(0, ymin)
        if data_range <= 0.1:
            # For very small ranges, use steps of 0.01 or 0.02
            step = 0.01 if data_range <= 0.05 else 0.02
        elif data_range <= 0.5:
            # For small ranges, use steps of 0.1
            step = 0.1
        elif data_range <= 1.0:
            # For medium ranges, use steps of 0.2
            step = 0.2
        else:
            # For larger ranges, use steps of 0.5 or 1.0
            step = 0.5 if data_range <= 3.0 else 1.0
        
        # Create nice round ticks
        start = np.floor(max(0, ymin) / step) * step
        end = np.ceil(ymax / step) * step
        
        # Generate approximately 6 ticks with nice round numbers
        num_steps = min(6, int((end - start) / step) + 1)
        if num_steps < 3:  # Ensure we have at least 3 ticks
            num_steps = 3
            step = (end - start) / (num_steps - 1)
        
        yticks = np.linspace(start, end, num_steps)
        ax.set_yticks(yticks)
        
        # Format y-tick labels to show fewer decimal places for cleaner display
        if step >= 0.5:
            # For larger steps, show integers or .5
            ytick_labels = [f"{y:.1f}" if y % 1 != 0 else f"{int(y)}" for y in yticks]
        elif step >= 0.1:
            # For medium steps, show one decimal place
            ytick_labels = [f"{y:.1f}" for y in yticks]
        else:
            # For small steps, show two decimal places
            ytick_labels = [f"{y:.2f}" for y in yticks]
            
        ax.set_yticklabels(ytick_labels, fontsize=60)
    
    # Increase tick font size and make x-ticks bold
    ax.tick_params(axis='both', which='major', labelsize=40)
    # for tick in ax.xaxis.get_major_ticks():
    #     tick.label1.set_fontweight('bold')
    
    # Create custom legend with bold text
    custom_lines = [
        plt.Line2D([0], [0], color=custom_palette["r_hat"], lw=4),
        plt.Line2D([0], [0], color=custom_palette["r_Mash"], lw=4),
        plt.Line2D([0], [0], color="gray", linestyle='--', lw=2)
    ]
    
    legend = ax.legend(custom_lines, [r'$\hat{r}$', 'r_Mash', 'True r'], 
              fontsize=60, frameon=True, framealpha=0.9, loc='upper left')
    
    # Make legend text bold
    # for text in legend.get_texts():
    #     text.set_fontweight('bold')
    
    return sns_plot

# Function to create boxplots (dual or single panel based on data)
def create_boxplots(data, r_threshold=0.25):
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
    
    # Check if either group is empty
    if len(small_r_values) == 0 or len(large_r_values) == 0:
        # Only one panel needed - use the non-empty group
        if len(small_r_values) == 0:
            print(f"Small r group is empty. Creating single panel for large r values.")
            r_values = large_r_values
            plot_data = large_r_data
        else:
            print(f"Large r group is empty. Creating single panel for small r values.")
            r_values = small_r_values
            plot_data = small_r_data
        
        # Create figure with a single panel
        fig, ax = plt.subplots(figsize=(24, 16))
        
        # Plot on the single axis
        plot_on_axis(ax, plot_data, r_values, is_right_panel=False)
        
    else:
        # Both groups have data - create dual panel
        # Create figure with two subplots
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(24, 16), gridspec_kw={'width_ratios': [3, 1]})
        
        # Plot on first subplot (r < threshold)
        plot1 = plot_on_axis(ax1, small_r_data, small_r_values, is_right_panel=False)
        
        # Plot on second subplot (r >= threshold)
        plot2 = plot_on_axis(ax2, large_r_data, large_r_values, is_right_panel=True)
        
        # Remove duplicate legends
        if ax2.get_legend():
            ax2.get_legend().remove()
    
    # Adjust layout
    plt.tight_layout()
    
    # Save the plot with higher DPI
    output_filename = f"dual_boxplot_k{k}.png"
    plt.savefig(output_filename, dpi=150)  # Increased DPI from 100 to 150
    print(f"Plot saved as '{output_filename}'")

# Main execution
if __name__ == "__main__":
    main()