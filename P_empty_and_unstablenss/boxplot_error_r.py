#!/usr/bin/env python3
import os
import re
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from glob import glob

def parse_args():
    parser = argparse.ArgumentParser(description='Q-value Variation Boxplot')
    parser.add_argument('-d', '--directory', type=str, help='Input directory', required=True)
    parser.add_argument('-t', '--threshold', type=float, default=0.25, 
                        help='Threshold for splitting r values')
    parser.add_argument('-k', '--kmer', type=int, default=30, 
                        help='k-mer size to analyze')
    return parser.parse_args()

def calculate_real_q(r, k):
    """Calculate real q from r and k: q = 1-(1-r)^k"""
    return 1 - (1 - r) ** k

def extract_r_from_filename(filename):
    """Extract r value from filename like r0.1_k30.output"""
    match = re.search(r'r([\d\.]+)_k\d+\.output', os.path.basename(filename))
    if match:
        return float(match.group(1))
    return None

def read_output_file(filepath):
    """Read the output file and extract P_empty and q_hat values"""
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
        
        # First line contains P_empty value
        p_empty = float(lines[0].strip())
        
        # Subsequent lines contain q_hat values
        q_hats = [float(line.strip()) for line in lines[1:] if line.strip()]
        
        return p_empty, q_hats
    except Exception as e:
        print(f"Error reading file {filepath}: {e}")
        return None, []

def main():
    # Parse arguments
    args = parse_args()
    
    # Find all r*_k30.output files in the directory
    pattern = os.path.join(args.directory, f"r*_k{args.kmer}.output")
    files = glob(pattern)
    
    if not files:
        print(f"No matching files found in {args.directory}")
        return
    
    print(f"Found {len(files)} files matching the pattern")
    
    # Store data for plotting
    df_rows = []
    
    # Process each file
    for filepath in files:
        r = extract_r_from_filename(filepath)
        if r is None:
            print(f"Could not extract r value from {filepath}. Skipping.")
            continue
        
        print(f"Processing file: {filepath}, extracted r value: {r}")
        
        p_empty, q_hats = read_output_file(filepath)
        real_q = calculate_real_q(r, args.kmer)
        
        # Skip if there was an error reading the file
        if p_empty is None:
            continue
        
        # Add each q_hat as a row in the dataframe
        for q_hat in q_hats:
            df_rows.append({
                'r_value': r,
                'q_hat': q_hat,
                'p_empty': p_empty,
                'real_q': real_q
            })
    
    if not df_rows:
        print("No valid data found in any of the files")
        return
    
    # Create DataFrame from all data
    all_data = pd.DataFrame(df_rows)
    
    # Sort data by r_value
    all_data = all_data.sort_values('r_value')
    
    # Get unique r values (sorted)
    all_r_values = sorted(all_data['r_value'].unique())
    
    print(f"All r values: {all_r_values}")
    
    # Create figure with a single panel
    fig, ax = plt.subplots(figsize=(20, 12))
    
    # Define colors for visualization elements
    boxplot_color = "#FF0000"  # Red for boxplot edges
    p_empty_color = "#0000FF"  # Blue for p_empty
    
    # Create a mapping of r values to x-coordinates
    x_ticks = np.arange(len(all_r_values))
    r_to_x = {r: x for x, r in zip(x_ticks, all_r_values)}
    
    # For each r value, collect all q_hat values and bounds
    r_boxplot_data = {}
    r_p_empty = {}
    r_real_q = {}
    
    for r in all_r_values:
        r_subset = all_data[all_data['r_value'] == r]
        r_boxplot_data[r] = r_subset['q_hat'].values
        
        # Get the p_empty for this r value (should be the same for all rows with this r)
        if not r_subset.empty:
            r_p_empty[r] = r_subset['p_empty'].iloc[0]
            r_real_q[r] = r_subset['real_q'].iloc[0]
    
    # Create boxplots for each r value
    boxplots = []
    for r in all_r_values:
        box_data = r_boxplot_data[r]
        x_pos = r_to_x[r]
        
        # Create a standard boxplot with median
        bp = ax.boxplot(
            box_data, 
            positions=[x_pos],
            widths=0.6,
            patch_artist=False,  # No fill color
            boxprops={'color': boxplot_color, 'linewidth': 1.5},
            medianprops={'color': 'black', 'linewidth': 1.5},  # Show the median line
            whiskerprops={'color': boxplot_color, 'linewidth': 1.2},
            capprops={'color': boxplot_color, 'linewidth': 1.2},
            flierprops={'marker': 'o', 'markerfacecolor': boxplot_color, 'markersize': 4, 'alpha': 0.5},
            showfliers=True
        )
        
        boxplots.append(bp)
    
    # Plot the P_empty values using triangles
    x_positions = []
    p_empty_values = []
    
    for r in all_r_values:
        x_pos = r_to_x[r]
        x_positions.append(x_pos)
        p_empty_values.append(r_p_empty[r])
    
    # Plot P_empty with down triangles
    ax.scatter(x_positions, p_empty_values, marker='v', color=p_empty_color, 
               s=100, alpha=0.8, label='P_empty (▼)')
    ax.plot(x_positions, p_empty_values, color=p_empty_color, 
            linestyle=':', linewidth=2, alpha=0.7)
    
    # Set x-axis ticks - only show every 3rd tick
    visible_ticks = []
    tick_labels = []
    
    for i, r in enumerate(all_r_values):
        if i % 3 == 0:  # Every 3rd tick
            visible_ticks.append(i)
            tick_labels.append(f"{r:.3f}")
    
    ax.set_xticks(visible_ticks)
    ax.set_xticklabels(tick_labels, rotation=45, fontsize=40)
    
    # Set fixed y-axis limits from 0 to 1.1
    ax.set_ylim(0, 1.1)
    
    # Set custom y-ticks to be sparse and go up to 1.0
    y_ticks = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
    ax.set_yticks(y_ticks)
    
    # Increase y-axis tick font size
    ax.tick_params(axis='y', labelsize=40)
    
    # Set labels
    ax.set_xlabel('r', fontsize=40)
    # 如代码显示，y轴标签已被注释掉
    #ax.set_ylabel('r', fontsize=40)
    
    # Add grid for easier reading
    ax.grid(True, axis='y', linestyle='--', alpha=0.7)
    
    # Add legend with custom elements and updated labels
    handles = [
        plt.Rectangle((0,0), 1, 1, fc="white", ec=boxplot_color, linewidth=1.5),
        plt.Line2D([0], [0], marker='v', color=p_empty_color, linestyle='', markersize=10)
    ]
    labels = [
        r'$\hat{r}$', 
        'P_empty'
    ]
    ax.legend(handles, labels, fontsize=40, frameon=True, framealpha=0.9, loc='upper left')
    
    # Adjust layout
    plt.tight_layout()
    
    # Save the plot
    output_file = f"boxplot_error_r_k{args.kmer}.png"
    plt.savefig(output_file, dpi=300)
    print(f"Boxplot saved to {output_file}")
    
    # Show plot
    plt.show()

if __name__ == "__main__":
    main()