import os
import argparse
import pandas as pd
import subprocess
import numpy as np
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from Bio import Phylo

def parse_args():
    parser = argparse.ArgumentParser(description="Run mash triangle + NJ/UPGMA tree (fasta input)")
    parser.add_argument("--csv", required=True, help="Input CSV with label,filepath")
    parser.add_argument("-k", type=int, default=31, help="k-mer size (default: 31)")
    parser.add_argument("-s", "--sketch_size", type=int, default=1000, help="Sketch size (default: 1000)")
    parser.add_argument("--temp_dir", default="temp", help="Directory to store triangle output")
    parser.add_argument("-o","--results_dir", default="results", help="Directory to store final outputs")
    parser.add_argument("--prefix", default="tree", help="Prefix for output files in results dir")
    parser.add_argument("--cores", type=int, default=4, help="Number of threads for mash triangle")
    parser.add_argument("--force", action="store_true", help="Force re-run even if output files exist")
    parser.add_argument("--tree-method", choices=["nj", "upgma", "both"], default="nj", 
                       help="Tree construction method: nj (Neighbor Joining), upgma (UPGMA), or both")
    return parser.parse_args()

def build_tree(labels, distances, out_file, method="nj"):
    """
    Build phylogenetic tree using specified method (NJ or UPGMA)
    """
    n = len(labels)
    
    method_name = "NJ" if method == "nj" else "UPGMA"
    print(f"[INFO] Building {method_name} tree for {n} species")
    
    full_matrix = np.zeros((n, n))
    
    for i in range(n):
        for j in range(n):
            if i == j:
                full_matrix[i][j] = 0.0  # set diagonal as 0
            else:
                label_i, label_j = labels[i], labels[j]
                dist = distances.get((label_i, label_j), distances.get((label_j, label_i), 0.0))
                full_matrix[i][j] = dist
    
    dm = DistanceMatrix(names=labels)
    for i in range(n):
        for j in range(i):
            dm[labels[i], labels[j]] = full_matrix[i][j]
    
    constructor = DistanceTreeConstructor()
    
    if method == "nj":
        tree = constructor.nj(dm)
    elif method == "upgma":
        tree = constructor.upgma(dm)
    else:
        raise ValueError(f"Unknown tree construction method: {method}")
    
    Phylo.write(tree, out_file, "newick")
    print(f"[INFO] {method_name} tree successfully written to {out_file}")

def main():
    args = parse_args()
    os.makedirs(args.temp_dir, exist_ok=True)
    os.makedirs(args.results_dir, exist_ok=True)

    triangle_path = os.path.join(args.temp_dir, f"{args.prefix}_triangle.tab")
    dist_matrix_txt = os.path.join(args.results_dir, f"{args.prefix}_distance_matrix.txt")
    pairwise_csv = os.path.join(args.results_dir, f"{args.prefix}_pairwise_distances.csv")
    
    # Define output file paths based on tree method
    if args.tree_method == "nj":
        nj_path = os.path.join(args.results_dir, f"{args.prefix}_nj.nwk")
        output_files = [nj_path, pairwise_csv]
    elif args.tree_method == "upgma":
        upgma_path = os.path.join(args.results_dir, f"{args.prefix}_upgma.nwk")
        output_files = [upgma_path, pairwise_csv]
    else:  # both
        nj_path = os.path.join(args.results_dir, f"{args.prefix}_nj.nwk")
        upgma_path = os.path.join(args.results_dir, f"{args.prefix}_upgma.nwk")
        output_files = [nj_path, upgma_path, pairwise_csv]

    # Check if output files already exist
    if not args.force and all(os.path.exists(f) for f in output_files):
        print(f"Output files already exist:")
        for f in output_files:
            print(f"  - {f}")
        print("Use --force to overwrite existing files")
        return

    # Step 1: Read CSV
    df = pd.read_csv(args.csv)
    labels = df["label"].tolist()
    filepaths = df["filepath"].tolist()
    
    print(f"Found {len(labels)} samples in CSV")

    # Verify all input files exist
    missing_files = []
    for i, fp in enumerate(filepaths):
        if not os.path.exists(fp):
            missing_files.append(f"{labels[i]}: {fp}")
    
    if missing_files:
        print("ERROR: The following input files are missing:")
        for mf in missing_files:
            print(f"  - {mf}")
        return

    # Step 2: Write fasta file list
    list_path = os.path.join(args.temp_dir, f"{args.prefix}_input_list.txt")
    with open(list_path, "w") as f:
        for fp in filepaths:
            f.write(fp.strip() + "\n")

    # Step 3: mash triangle on fasta files (skip if triangle file exists and not forcing)
    if not args.force and os.path.exists(triangle_path):
        print(f"[1/4] Triangle file already exists: {triangle_path}")
    else:
        print(f"[1/4] Running mash triangle directly on FASTA files using {args.cores} threads...")
        cmd = [
            "mash", "triangle",
            "-k", str(args.k),
            "-s", str(args.sketch_size),
            "-p", str(args.cores),
            "-l", list_path
        ]
        with open(triangle_path, "w") as out:
            subprocess.run(cmd, stdout=out, check=True)

    # Step 4: parse triangle result
    print("[2/4] Parsing triangle output...")
    with open(triangle_path) as f:
        lines = [line.strip() for line in f if line.strip()]

    print(f"Triangle file has {len(lines)} lines")
    
    # First line is usually just the number of samples, skip it
    if len(lines) > 0 and lines[0].isdigit():
        print(f"Skipping header line with sample count: {lines[0]}")
        lines = lines[1:]
    
    print(f"Processing {len(lines)} data lines")
    
    if len(lines) != len(labels):
        print(f"WARNING: Triangle data has {len(lines)} lines but CSV has {len(labels)} samples")
        print("This might cause indexing issues")

    matrix = []
    mash_labels = []  # Extract labels from triangle file
    
    for idx, line in enumerate(lines):
        parts = line.split("\t")
        if len(parts) > 0:
            # First part is the file path, extract just the filename/label
            file_path = parts[0]
            mash_labels.append(os.path.basename(file_path))
            
            # Rest are the distance values
            values = list(map(float, parts[1:])) if len(parts) > 1 else []
            while len(values) < idx:
                values.insert(0, 0.0)
            matrix.append(values)
    
    print(f"Extracted {len(mash_labels)} labels from triangle file")

    n = len(matrix)
    print(f"Building {n}x{n} distance matrix")
    
    # Use the original labels from CSV for consistency
    if n != len(labels):
        print(f"ERROR: Matrix size ({n}) doesn't match number of labels ({len(labels)})")
        print("This usually means the triangle output doesn't match the input CSV")
        print("Mash labels found:", mash_labels[:5], "..." if len(mash_labels) > 5 else "")
        print("CSV labels found:", labels[:5], "..." if len(labels) > 5 else "")
        return

    full = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i):
            if j < len(matrix[i]):
                full[i][j] = full[j][i] = matrix[i][j]

    # Step 3: Save outputs
    print("[3/4] Saving distance matrix and pairwise distances...")
    
    # Save matrix (optional - for compatibility)
    print(f"Saving distance matrix to: {dist_matrix_txt}")
    with open(dist_matrix_txt, "w") as f:
        f.write("label\t" + "\t".join(labels) + "\n")
        for i, row in enumerate(full):
            f.write(labels[i] + "\t" + "\t".join(f"{v:.6f}" for v in row) + "\n")

    # Save pairwise distances as CSV
    print(f"Saving pairwise distances to: {pairwise_csv}")
    pairwise_data = []
    for i in range(n):
        for j in range(i+1, n):  # Only upper triangle to avoid duplicates
            pairwise_data.append({
                'label1': labels[i],
                'label2': labels[j], 
                'distance': full[i][j]
            })
    
    pairwise_df = pd.DataFrame(pairwise_data)
    pairwise_df.to_csv(pairwise_csv, index=False)

    # Step 4: Build phylogenetic tree(s)
    print("[4/4] Building phylogenetic tree(s)...")

    # Convert the full matrix to a distances dictionary for the build_tree function
    distances = {}
    for i in range(n):
        for j in range(i+1, n):
            distances[(labels[i], labels[j])] = full[i][j]
    
    # Build tree(s) based on selected method
    if args.tree_method == "nj":
        build_tree(labels, distances, nj_path, "nj")
        print(f"NJ tree written to: {nj_path}")
    elif args.tree_method == "upgma":
        build_tree(labels, distances, upgma_path, "upgma")
        print(f"UPGMA tree written to: {upgma_path}")
    elif args.tree_method == "both":
        # Build both NJ and UPGMA trees
        build_tree(labels, distances, nj_path, "nj")
        build_tree(labels, distances, upgma_path, "upgma")
        print(f"NJ tree written to: {nj_path}")
        print(f"UPGMA tree written to: {upgma_path}")

    print(f"Distance matrix saved to: {dist_matrix_txt}")
    print(f"Pairwise distances saved to: {pairwise_csv}")
    print(f"Temp triangle and list files kept in: {args.temp_dir}")

if __name__ == "__main__":
    main()