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
    parser = argparse.ArgumentParser(description='Q-value Relative Error Boxplot')
    parser.add_argument('-d', '--directory', type=str, help='Input directory', required=True)
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
    """Read the output file and extract bounds and q_hat values"""
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
        
        # First line contains bounds a, b separated by comma
        bounds = list(map(float, lines[0].strip().split(',')))
        
        # Subsequent lines contain q_hat values
        q_hats = [float(line.strip()) for line in lines[1:] if line.strip()]
        
        return bounds, q_hats
    except Exception as e:
        print(f"Error reading file {filepath}: {e}")
        return [None, None], []

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
        
        bounds, q_hats = read_output_file(filepath)
        real_q = calculate_real_q(r, args.kmer)
        
        # Calculate relative errors: (value - real_q) / real_q
        # Handle case where real_q is very small or zero to avoid division by zero
        if real_q < 1e-10:
            # If real_q is effectively zero, use absolute error instead
            q_hat_rel_errors = [q_hat - real_q for q_hat in q_hats]
            lower_bound_rel_error = bounds[0] - real_q if bounds[0] is not None else None
            upper_bound_rel_error = bounds[1] - real_q if bounds[1] is not None else None
        else:
            # Calculate relative errors
            q_hat_rel_errors = [(q_hat - real_q) / real_q for q_hat in q_hats]
            lower_bound_rel_error = (bounds[0] - real_q) / real_q if bounds[0] is not None else None
            upper_bound_rel_error = (bounds[1] - real_q) / real_q if bounds[1] is not None else None
        
        # Add each q_hat relative error as a row in the dataframe
        for q_hat_rel_error in q_hat_rel_errors:
            df_rows.append({
                'r_value': r,
                'q_hat_rel_error': q_hat_rel_error,
                'lower_bound_rel_error': lower_bound_rel_error,
                'upper_bound_rel_error': upper_bound_rel_error,
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
    r_values = sorted(all_data['r_value'].unique())
    
    print(f"All r values: {r_values}")
    
    # Create figure with a single plot
    fig, ax = plt.subplots(figsize=(20, 12))
    
    # Define colors for visualization elements
    boxplot_color = "#FF0000"  # Red for boxplot edges
    lower_bound_color = "#0000FF"  # Blue for lower bound
    upper_bound_color = "#00AA00"  # Green for upper bound
    zero_line_color = "#777777"    # Gray for zero line (no error)
    
    # Create a mapping of r values to x-coordinates
    x_ticks = np.arange(len(r_values))
    r_to_x = {r: x for x, r in zip(x_ticks, r_values)}
    
    # For each r value, collect all q_hat_rel_error values and bounds
    r_boxplot_data = {}
    r_lower_bounds = {}
    r_upper_bounds = {}
    
    # To calculate the max quartile 3 value (top of box without whiskers)
    all_q3_values = []
    all_q1_values = []
    
    for r in r_values:
        r_subset = all_data[all_data['r_value'] == r]
        r_boxplot_data[r] = r_subset['q_hat_rel_error'].values
        
        # Calculate quartiles to determine box heights
        q1, q3 = np.percentile(r_boxplot_data[r], [25, 75])
        all_q3_values.append(q3)
        all_q1_values.append(q1)
        
        # Get the bounds for this r value (they should be the same for all rows with this r)
        if not r_subset.empty:
            r_lower_bounds[r] = r_subset['lower_bound_rel_error'].iloc[0]
            r_upper_bounds[r] = r_subset['upper_bound_rel_error'].iloc[0]
    
    # Create boxplots for each r value
    boxplots = []
    for r in r_values:
        box_data = r_boxplot_data[r]
        x_pos = r_to_x[r]
        
        # Calculate mean for the mean line
        mean_val = np.mean(box_data)
        
        # Create a boxplot with no whiskers
        bp = ax.boxplot(
            box_data, 
            positions=[x_pos],
            widths=0.6,
            patch_artist=False,  # No fill color
            boxprops={'color': boxplot_color, 'linewidth': 1.5},
            medianprops={'color': 'black', 'linewidth': 0},  # Hide the median line
            whiskerprops={'color': boxplot_color, 'linewidth': 0},  # Hide whiskers
            capprops={'color': boxplot_color, 'linewidth': 0},  # Hide caps
            showfliers=False,  # No outliers
            whis=0  # No whiskers
        )
        
        # Add a horizontal line for the mean
        ax.plot([x_pos-0.3, x_pos+0.3], [mean_val, mean_val], 
                color='black', linewidth=2.0, solid_capstyle='butt')
        
        boxplots.append(bp)
    
    # Plot the reference line at y=0 (no error) WITHOUT adding it to the legend
    ax.axhline(y=0, color=zero_line_color, linestyle='--', linewidth=2, alpha=0.8)
    
    # Plot the theoretical bound errors
    x_positions = []
    lower_bound_errors = []
    upper_bound_errors = []
    
    for r in r_values:
        x_pos = r_to_x[r]
        x_positions.append(x_pos)
        
        # Only append non-None values
        if r_lower_bounds[r] is not None:
            lower_bound_errors.append((x_pos, r_lower_bounds[r]))
        
        if r_upper_bounds[r] is not None:
            upper_bound_errors.append((x_pos, r_upper_bounds[r]))
    
    # Plot lower bounds with down triangles
    if lower_bound_errors:
        x_pos_lower, y_lower = zip(*lower_bound_errors)
        ax.scatter(x_pos_lower, y_lower, marker='v', color=lower_bound_color, 
                s=100, alpha=0.8, label='Lower bound error (▼)')
        ax.plot(x_pos_lower, y_lower, color=lower_bound_color, 
                linestyle=':', linewidth=2, alpha=0.7)
    
    # Plot upper bounds with up triangles
    if upper_bound_errors:
        x_pos_upper, y_upper = zip(*upper_bound_errors)
        ax.scatter(x_pos_upper, y_upper, marker='^', color=upper_bound_color, 
                s=100, alpha=0.8, label='Upper bound error (▲)')
        ax.plot(x_pos_upper, y_upper, color=upper_bound_color, 
                linestyle='-.', linewidth=2, alpha=0.7)
    
    # Set x-axis ticks - show only every 3rd tick
    visible_ticks = []
    tick_labels = []
    
    for i, r in enumerate(r_values):
        if i % 3 == 0:  # Every 3rd tick
            visible_ticks.append(x_ticks[i])
            tick_labels.append(f"{r:.3f}")
    
    ax.set_xticks(visible_ticks)
    ax.set_xticklabels(tick_labels, rotation=45, fontsize=50)
    
    # Set labels with simplified axis names
    ax.set_xlabel('r', fontsize=60)
    ax.set_ylabel('Relative error of q', fontsize=60)
    
    # Add grid for easier reading
    ax.grid(True, axis='y', linestyle='--', alpha=0.7)
    
    # Set y-axis limits based on the data and boxes only (no whiskers)
    max_y = max(max(all_q3_values), 0.05)  # Use the maximum Q3 value or at least 0.05
    min_y = min(min(all_q1_values), -0.05)  # Use the minimum Q1 value or at most -0.05
    
    # Check bounds errors too
    if lower_bound_errors:
        min_y = min(min_y, min(y_lower))
    if upper_bound_errors:
        max_y = max(max_y, max(y_upper))
    
    # Add some padding
    max_y = max_y * 1.1
    min_y = min_y * 1.1 if min_y < 0 else min_y * 0.9
    
    # Ensure zero is visible
    min_y = min(min_y, -0.05)
    
    # Round to simple clean values - using only 0.2, 0.4, 0.6, etc.
    # Find the closest multiple of 0.2 that's greater than max_y
    max_y_rounded = np.ceil(max_y / 0.2) * 0.2
    # Find the closest multiple of 0.2 that's less than min_y
    min_y_rounded = np.floor(min_y / 0.2) * 0.2
    
    # Calculate how many 0.2 steps we need between min and max
    steps_needed = int((max_y_rounded - min_y_rounded) / 0.2) + 1
    
    # If we need too many steps, try using 0.4 as the step
    if steps_needed > 6:
        max_y_rounded = np.ceil(max_y / 0.4) * 0.4
        min_y_rounded = np.floor(min_y / 0.4) * 0.4
        steps_needed = int((max_y_rounded - min_y_rounded) / 0.4) + 1
        yticks = np.linspace(min_y_rounded, max_y_rounded, steps_needed)
    else:
        yticks = np.linspace(min_y_rounded, max_y_rounded, steps_needed)
    
    # Set the limits and ticks
    ax.set_ylim(min_y_rounded, max_y_rounded)
    ax.set_yticks(yticks)
    
    # Format y-tick labels without trailing zeros (0.6 instead of 0.60)
    y_tick_labels = []
    for y in yticks:
        if y == int(y):
            # For whole numbers, show as integers
            y_tick_labels.append(f"{int(y)}")
        else:
            # Remove trailing zeros
            y_str = f"{y:.1f}".rstrip('0').rstrip('.') if y == round(y, 1) else f"{y:.2f}".rstrip('0').rstrip('.')
            y_tick_labels.append(y_str)
    
    ax.set_yticklabels(y_tick_labels, fontsize=50)
    
    # Add legend with custom elements - REMOVED THE "NO ERROR" ITEM
    handles = [
        plt.Rectangle((0,0), 1, 1, fc="white", ec=boxplot_color, linewidth=1.5),
        plt.Line2D([0], [0], color='black', linewidth=2),
        plt.Line2D([0], [0], marker='v', color=lower_bound_color, linestyle='', markersize=10),
        plt.Line2D([0], [0], marker='^', color=upper_bound_color, linestyle='', markersize=10)
    ]
    labels = [
        'q_hat Error (Box)', 
        'Mean Error',
        'Lower bound Error', 
        'Upper bound Error'
    ]
    ax.legend(handles, labels, fontsize=30, frameon=True, framealpha=0.9, loc='upper right')
    
    # Adjust layout
    plt.tight_layout()
    
    # Save the plot
    output_file = f"rel_error_q_k{args.kmer}.png"
    plt.savefig(output_file, dpi=100)
    print(f"Boxplot saved to {output_file}")
    
    # Show plot
    plt.show()

if __name__ == "__main__":
    main()