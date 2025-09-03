#!/usr/bin/env python3
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict, OrderedDict
import os

def parse_plantpan(file_path, tfbs_filter=None):
    """Parse PlantPAN result file into {gene: [motif dicts]} preserving input order."""
    df = pd.read_csv(file_path, sep="\t")

    # Filter by TFBS ID if requested
    if tfbs_filter:
        df = df[df["TFBS ID"].isin(tfbs_filter)]

    # Preserve gene order from input
    gene_order = list(OrderedDict.fromkeys(df["Sequence ID"]))

    motifs = defaultdict(list)
    for _, row in df.iterrows():
        # Handle missing TFBS Name
        name = str(row["TFBS Name"]).strip()
        if not name or name.lower() == "nan":
            family = str(row.get("TF family", "")).strip()
            if family and family.lower() != "nan":
                name = f"{row['TFBS ID']} ({family})"
            else:
                name = row["TFBS ID"]

        motifs[row["Sequence ID"]].append({
            "motif_id": row["TFBS ID"],
            "motif_name": name,
            "position": int(row["Posistion"]),
            "seq": row["Binding sequence"]
        })

    return motifs, gene_order

def plot_motifs(motifs, gene_order, output_dir, left_margin=20):
    """Plot motifs across promoter regions ignoring strand orientation and using real promoter length."""
    # Get all unique motifs for colors
    unique_motifs = sorted({m["motif_name"] for v in motifs.values() for m in v})
    colors = plt.colormaps["tab20"](np.linspace(0, 1, len(unique_motifs)))
    motif_to_color = {name: colors[i] for i, name in enumerate(unique_motifs)}

    # representative sequences for legend
    motif_to_seq = {}
    for mlist in motifs.values():
        for m in mlist:
            if m["motif_name"] not in motif_to_seq:
                motif_to_seq[m["motif_name"]] = m["seq"]

    plt.figure(figsize=(14, max(6, len(gene_order)*1.2)))
    box_height = 0.6

    # Determine promoter length per gene based on last motif
    promoter_lengths = {gene: max(m["position"] for m in mlist) if mlist else 0
                        for gene, mlist in motifs.items()}

    max_promoter_length = max(promoter_lengths.values())

    for i, gene in enumerate(gene_order, start=1):
        promoter_length = promoter_lengths.get(gene, max_promoter_length)

        # shift everything to the left by left_margin
        shifted_start = -promoter_length - left_margin

        # draw promoter region rectangle
        plt.gca().add_patch(plt.Rectangle(
            (shifted_start, i - box_height/2),
            promoter_length + left_margin, box_height,
            linewidth=1, edgecolor="black", facecolor="lightgrey"
        ))

        # draw motifs
        for m in motifs.get(gene, []):
            x = shifted_start + m["position"]
            color = motif_to_color[m["motif_name"]]
            plt.vlines(x, i - 0.3, i + 0.3, color=color, linewidth=4)

    # axis setup
    #plt.axvline(0, color="black")  # transcription start site
    plt.title("Motif Distribution Across Promoter Regions", fontsize=16)
    plt.xlabel("Position (bp)", fontsize=14)
    plt.ylabel("Genes", fontsize=14)
    plt.yticks(range(1, len(gene_order)+1), gene_order)
    plt.xlim(-max_promoter_length - left_margin, left_margin)
    plt.tight_layout()

    # legend
    legend_handles = [
        plt.Line2D([0], [0], color=motif_to_color[m], lw=4,
                   label=f"{m} ({motif_to_seq[m]})")
        for m in unique_motifs
    ]
    plt.legend(handles=legend_handles, bbox_to_anchor=(0.5, -0.15),
               loc="upper center", ncol=2, fontsize="small", frameon=True)

    # Save PDF + SVG
    os.makedirs(output_dir, exist_ok=True)
    pdf_path = os.path.join(output_dir, "motif_plot.pdf")
    svg_path = os.path.join(output_dir, "motif_plot.svg")
    plt.savefig(pdf_path, dpi=300, bbox_inches="tight")
    plt.savefig(svg_path, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Plots saved to:\n {pdf_path}\n {svg_path}")

def main():
    parser = argparse.ArgumentParser(description="Visualize PlantPAN motif distributions")
    parser.add_argument("-i", "--input", required=True, help="PlantPAN result file (TSV)")
    parser.add_argument("-o", "--output_dir", default="motif_plots", help="Output folder for plots")
    parser.add_argument("-t", "--tfbs", nargs="+",
                        help="One or more TFBS IDs to include (default: all)")
    args = parser.parse_args()

    motifs, gene_order = parse_plantpan(args.input, tfbs_filter=args.tfbs)
    if not motifs:
        print("No motifs found in input file.")
        return

    plot_motifs(motifs, gene_order, args.output_dir)

if __name__ == "__main__":
    main()
