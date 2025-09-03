#!/usr/bin/env python3

"""
Author: Syafiq Samsolnizam
Date: 2025-01-02

Description:
    This script reads a TPM table (TSV) with the following format:
      - First column: Gene names
      - Subsequent columns: Sample TPM values

    Optional arguments:
      - A file containing genes of interest (one gene per line) to subset the data.
      - A "tissue file" mapping each tissue to comma-separated sample names for
        color-grouping in the plot.
      - A custom title for the plot.
      - The size of the stripplot points and boxplot outlier markers.
      - A flag to remove outliers from the boxplot display.
      - The box width (to control horizontal spacing) for each category.
      - Unmapped samples (i.e., those not found in the tissue mapping) will be
        excluded from the plot if a tissue file is provided.
      - A user-defined order for the tissues (--tissue_order), so they appear
        in that order in both the plot and legend.

    The script will generate a box plot of TPM values. If a tissue file is
    provided, samples are grouped/color-coded by tissue. The legend
    will show sample counts for each tissue and is displayed
    inside the plot at the top-right corner.
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

def main():
    parser = argparse.ArgumentParser(
        description="Create a box plot from a TPM table with optional gene query, tissue mapping, custom aesthetics, and tissue order."
    )
    parser.add_argument("-i", "--input", required=True,
        help="Path to the input TPM table (TSV).")
    parser.add_argument("-q", "--query", required=False, default=None,
        help="Path to a file containing genes of interest (one per line).")
    parser.add_argument("-t", "--tissue", required=False, default=None,
        help="Path to a tissue assignment file. Format: [tissue]\\t[sample1,sample2,...]")
    parser.add_argument("--tissue_order", required=False, default=None,
        help="Comma-separated list of tissues in the desired order, e.g. 'Leaf,Stem,Root'.")
    parser.add_argument("--title", required=False, default=None,
        help="Custom title for the plot.")
    parser.add_argument("--point_size", required=False, type=float, default=5.0,
        help="Size of the individual points in the strip plot.")
    parser.add_argument("--outlier_size", required=False, type=float, default=5.0,
        help="Size of the outliers in the box plot.")
    parser.add_argument("--remove_outliers", action="store_true", default=False,
        help="If provided, outliers won't be displayed in the boxplot and individual data points will be hidden.")
    parser.add_argument("--box_width", required=False, type=float, default=0.6,
        help="Width of each box in the boxplot to control spacing. Smaller means more spacing.")
    parser.add_argument("-o", "--output", required=False, default="boxplot.png",
        help="Output filename for the plot (PNG/SVG/etc.).")

    args = parser.parse_args()

    # --- 1. Load the TPM table ---
    try:
        df = pd.read_csv(args.input, sep="\t", header=0)
    except Exception as e:
        print(f"[ERROR] Could not read the TPM table: {e}")
        sys.exit(1)

    if df.shape[1] < 2:
        print("[ERROR] The TPM table must have at least 2 columns: Gene + at least one sample.")
        sys.exit(1)
    
    gene_col = df.columns[0]

    # --- 2. Handle query genes ---
    query_genes = None
    if args.query:
        try:
            with open(args.query, "r") as f:
                query_genes = [line.strip() for line in f if line.strip()]
        except Exception as e:
            print(f"[ERROR] Could not read query genes file: {e}")
            sys.exit(1)

        original_count = df.shape[0]
        df = df[df[gene_col].isin(query_genes)]
        filtered_count = df.shape[0]
        print(f"Filtered genes from {original_count} to {filtered_count} based on query file.")

        # enforce gene order from query file
        df[gene_col] = pd.Categorical(df[gene_col], categories=query_genes, ordered=True)
    else:
        # keep original TPM table order
        df[gene_col] = pd.Categorical(df[gene_col], categories=df[gene_col].tolist(), ordered=True)

    # --- 3. Load tissue assignment (optional) ---
    sample_to_tissue = {}
    if args.tissue:
        try:
            tissue_df = pd.read_csv(args.tissue, sep="\t", header=None)
            tissue_df.columns = ["Tissue", "Samples"]
        except Exception as e:
            print(f"[ERROR] Could not read tissue file: {e}")
            sys.exit(1)

        for _, row in tissue_df.iterrows():
            tissue = row["Tissue"]
            samples = row["Samples"].split(",")
            for s in samples:
                sample_to_tissue[s.strip()] = tissue

    # --- 4. Reshape for plotting ---
    df_melted = df.melt(id_vars=gene_col, var_name="Sample", value_name="TPM")

    # --- 5. Apply tissue mapping ---
    if args.tissue:
        df_melted["Tissue"] = df_melted["Sample"].map(sample_to_tissue)
        unmapped = df_melted["Tissue"].isna().sum()
        if unmapped > 0:
            print(f"[INFO] Dropped {unmapped} unmapped samples from plot.")
        df_melted = df_melted.dropna(subset=["Tissue"])
    else:
        df_melted["Tissue"] = "NA"

    if df_melted.empty:
        print("[WARNING] After filtering/unmapped removal, no data remains for plotting.")

    # --- 6. Tissue order ---
    if args.tissue_order:
        desired_tissue_order = [t.strip() for t in args.tissue_order.split(",")]
    else:
        desired_tissue_order = sorted(df_melted["Tissue"].unique())

    df_melted["Tissue"] = pd.Categorical(
        df_melted["Tissue"],
        categories=desired_tissue_order,
        ordered=True
    )

    # --- 7. Seaborn settings ---
    sns.set(style="whitegrid")
    plt.figure(figsize=(8, 6))

    plot_title = args.title if args.title else (
        "TPM Box Plot of Genes" if not args.query else "TPM Box Plot of Queried Genes"
    )

    # --- 8. Boxplot ---
    ax = sns.boxplot(
        data=df_melted,
        x=gene_col,
        y="TPM",
        hue="Tissue",
        palette="Set2",
        dodge=True,
        fliersize=args.outlier_size,
        showfliers=not args.remove_outliers,
        width=args.box_width,
        hue_order=desired_tissue_order,
        order=query_genes if args.query else df[gene_col].cat.categories.tolist()
    )

    # --- 9. Stripplot (skip if removing outliers) ---
    if not args.remove_outliers:
        sns.stripplot(
            data=df_melted,
            x=gene_col,
            y="TPM",
            hue="Tissue",
            palette="Set2",
            dodge=True,
            alpha=0.5,
            linewidth=0.5,
            edgecolor="gray",
            size=args.point_size,
            hue_order=desired_tissue_order,
            order=query_genes if args.query else df[gene_col].cat.categories.tolist(),
            legend=False  # avoid duplicate legend
        )

    # --- 10. Legend ---
    # recompute counts from plotted data
    tissue_counts = df_melted.groupby("Tissue", observed=True)["Sample"].nunique().to_dict()
    handles, labels = ax.get_legend_handles_labels()
    new_handles, new_labels, used = [], [], set()

    for handle, label in zip(handles, labels):
        if label not in used and label in tissue_counts:
            used.add(label)
            label_text = f"{label} (n={tissue_counts[label]})"
            new_handles.append(handle)
            new_labels.append(label_text)

    ax.legend(
        new_handles, new_labels, 
        title="Tissue",
        loc="upper right",
        frameon=True
    )

    # --- 11. Final styling ---
    plt.xticks(rotation=45, ha="right")
    plt.title(plot_title, wrap=True)
    plt.tight_layout()

    plt.savefig(args.output, dpi=600)
    print(f"Box plot saved to {args.output}")

if __name__ == "__main__":
    main()
