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
    parser.add_argument('-k', '--kmer', type=int, default=30, 
                        help='k-mer size to analyze')
    args = parser.parse_args()

    global k
    k = args.kmer

    all_data = process_all_files(args.directory)
    
    if all_data is not None:
        create_difference_boxplot(all_data)

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
        
        # Calculate relative differences from true r: (estimate - true) / true
        df['r_hat_rel_diff'] = (df['r_hat'] - true_r) / true_r
        df['r_Mash_rel_diff'] = (df['r_Mash'] - true_r) / true_r
        
        # Print a sample of the processed data
        print(f"Sample of processed data: {df.head(2)}")
        
        return df
        
    except Exception as e:
        print(f"Error processing file {file_path}: {e}")
        # Return an empty DataFrame with the correct structure
        return pd.DataFrame(columns=['r_hat', 'r_Mash', 'true_r', 'r_hat_rel_diff', 'r_Mash_rel_diff'])

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

# Function to create a single boxplot showing relative differences from true r values
def create_difference_boxplot(data):
    # Sort the dataframe by true_r
    data = data.sort_values('true_r')
    
    # Get unique r values (sorted)
    r_values = sorted(data['true_r'].unique())
    
    print(f"All r values: {r_values}")
    
    # Custom color palette with high contrast colors
    custom_palette = {"r_hat_rel_diff": "#FF0000", "r_Mash_rel_diff": "#FFD700"}
    
    # Create figure with a single plot
    fig, ax = plt.subplots(figsize=(18, 12))
    
    # Melt the dataframe to prepare for boxplot
    melted_data = pd.melt(data, 
                        id_vars=['true_r'], 
                        value_vars=['r_hat_rel_diff', 'r_Mash_rel_diff'],
                        var_name='method', 
                        value_name='relative_difference')
    
    # Create the boxplot
    sns.set(style="whitegrid", font_scale=1.2)
    sns_plot = sns.boxplot(x='true_r', y='relative_difference', hue='method', data=melted_data,
                     palette=custom_palette,
                     width=0.9,
                     linewidth=1.0,
                     fliersize=5,
                     ax=ax)
    
    # Add a horizontal line at y=0 (no difference from true r)
    ax.axhline(y=0, color='gray', linestyle='--', linewidth=2.0, alpha=0.8)
    
    # Set x-axis ticks - only show every 3rd tick
    visible_ticks = []
    tick_labels = []
    
    for i, r in enumerate(r_values):
        if i % 3 == 0:  # Every 3rd tick
            visible_ticks.append(i)
            tick_labels.append(f"{r:.3f}")
    
    ax.set_xticks(visible_ticks)
    ax.set_xticklabels(tick_labels, rotation=45, fontsize=55)
    
    # Set labels
    ax.set_xlabel('r', fontsize=60)
    ax.set_ylabel('Relative error', fontsize=60)
    
    # Increase y-axis tick font size
    ax.tick_params(axis='y', labelsize=55)
    
    # Add grid for easier reading
    ax.grid(True, axis='y', linestyle='--', alpha=0.7)
    
    # Set custom legend labels
    handles, labels = ax.get_legend_handles_labels()
    
    if ax.get_legend():
        ax.get_legend().remove()
    
    custom_lines = [
        plt.Line2D([0], [0], color=custom_palette["r_hat_rel_diff"], lw=4),
        plt.Line2D([0], [0], color=custom_palette["r_Mash_rel_diff"], lw=4),
    ]
    
    ax.legend(custom_lines, [r'$(\hat{r} - r)/r$', '(r_Mash - r)/r'], 
               fontsize=55, frameon=True, framealpha=0.9, loc='upper right')
    
    # Adjust layout
    plt.tight_layout()
    
    # Save the plot
    output_filename = f"diff_boxplot_r_k{k}.png"
    plt.savefig(output_filename, dpi=300)
    print(f"Plot saved as '{output_filename}'")

# Main execution
if __name__ == "__main__":
    main()