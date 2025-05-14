import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os

def parse_args():
    parser = argparse.ArgumentParser(description='Generate heatmap from deviation data')
    parser.add_argument('-i', '--input', required=True, help='Input CSV file')
    parser.add_argument('-o', '--output', required=True, help='Output PNG file')
    parser.add_argument('--title', default='Deviation Heatmap', help='Plot title')
    parser.add_argument('--cmap', default='coolwarm', help='Colormap for heatmap')
    parser.add_argument('--vmin', type=float, default=0, help='Minimum value for colorbar (default: 0)')
    parser.add_argument('--vmax', type=float, default=10, help='Maximum value for colorbar (default: 10)')
    parser.add_argument('--figsize', type=str, default='10,16', help='Figure size as width,height')
    parser.add_argument('--annotate', action='store_true', help='Show values in cells')
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Load the data
    df = pd.read_csv(args.input)
    print(f"Data loaded: {len(df)} rows")
    print(f"Column names: {df.columns.tolist()}")
    
    # Get the range of deviation values (for information only)
    min_dev = df['deviation'].min()
    max_dev = df['deviation'].max()
    print(f"Deviation range: {min_dev} to {max_dev}")
    print(f"Using fixed color range: 0 to 10")
    
    # Create pivot table
    try:
        # Convert r to numeric to ensure proper sorting
        df['r'] = pd.to_numeric(df['r'])
        
        # Create pivot table
        pivot_table = df.pivot(index='r', columns='k', values='deviation')
        
        # Sort the pivot table by r in ascending order
        pivot_table = pivot_table.sort_index(ascending=True)
        print(f"Pivot table created with shape: {pivot_table.shape}")
    except Exception as e:
        print(f"Error creating pivot table: {e}")
        print("Trying with pivot_table function...")
        df['r'] = pd.to_numeric(df['r'])
        pivot_table = pd.pivot_table(df, values='deviation', index='r', columns='k', aggfunc=np.mean)
        pivot_table = pivot_table.sort_index(ascending=True)
    
    # Fixed colorbar range
    vmin = 0
    vmax = 1  # 10
    
    # Create the figure
    width, height = map(float, args.figsize.split(','))
    plt.figure(figsize=(width, height))

    # Create the heatmap
    # Remove the origin parameter which is causing the error
    ax = sns.heatmap(pivot_table, 
                     cmap=args.cmap, 
                     vmin=vmin, 
                     vmax=vmax,
                     annot=False,
                     cbar_kws={'label': 'NRMSE'})
    
    # Invert the y-axis to have small values at the bottom
    ax.invert_yaxis()
    
    # Reduce x-axis ticks (k values)
    k_values = pivot_table.columns.tolist()
    if len(k_values) > 4:
        # Select approximately 6 evenly spaced k values
        k_indices = np.linspace(0, len(k_values) - 1, 4, dtype=int)
        k_ticks = [k_values[i] for i in k_indices]
        
        # Get current tick positions (which are cell-centered by default in seaborn)
        current_positions = ax.get_xticks()
        # Select only the positions we want to keep
        new_positions = [current_positions[i] for i in k_indices]
        
        ax.set_xticks(new_positions)
        ax.set_xticklabels([str(tick) for tick in k_ticks])
    
    # Select only 4 y-axis ticks (r values)
    r_values = pivot_table.index.tolist()
    if len(r_values) > 4:
        # Get 4 evenly spaced indices
        r_indices = np.linspace(0, len(r_values) - 1, 4, dtype=int)
        r_ticks = [r_values[i] for i in r_indices]
        
        # Get current tick positions (which are cell-centered by default in seaborn)
        current_positions = ax.get_yticks()
        # Select only the positions we want to keep
        new_positions = [current_positions[i] for i in r_indices]
        
        ax.set_yticks(new_positions)
        ax.set_yticklabels([str(tick) for tick in r_ticks])
    
    # Increase font size for both x and y tick labels
    plt.xticks(fontsize=40)
    plt.yticks(fontsize=40)
    
    plt.xlabel('k', fontsize=50)
    plt.ylabel('r', fontsize=50)
    
    # Format the colorbar label
    cbar = ax.collections[0].colorbar
    cbar.set_label('P_empty', fontsize=50)
    
    # Increase font size for colorbar ticks
    cbar.ax.tick_params(labelsize=40)
    
    plt.tight_layout()
    plt.savefig(args.output, dpi=300)
    plt.close()
    
    print(f"Heatmap saved to {args.output}")

if __name__ == "__main__":
    main()