import sys, os, subprocess
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

__usage__ = """
python3 pairwise_heatmap.py \
    --in <FULL_PATH_TO_INPUT_FASTA> \
    --out <FULL_PATH_TO_OUTPUT_FOLDER>
"""

def load_sequences(fasta_file):
    """Load sequences from FASTA into dict (keep input order)."""
    sequences = {}
    order = []
    with open(fasta_file) as f:
        header = None
        seq = []
        for line in f:
            if line.startswith(">"):
                if header:
                    sequences[header] = "".join(seq).upper()
                header = line.strip()[1:]
                order.append(header)
                seq = []
            else:
                seq.append(line.strip())
        if header:
            sequences[header] = "".join(seq).upper()
    return sequences, order

def run_mafft(input_fasta, output_fasta):
    """Run global multiple sequence alignment with MAFFT."""
    cmd = f"mafft --maxiterate 1000 --genafpair {input_fasta} > {output_fasta}"
    subprocess.run(cmd, shell=True, check=True)

def compute_identity_matrix(aln_file, order):
    """Compute pairwise identity from an alignment file."""
    seqs, _ = load_sequences(aln_file)
    aligned = [seqs[name] for name in order]

    n = len(order)
    matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                matrix[i, j] = 100.0
            else:
                matches = sum(a == b for a, b in zip(aligned[i], aligned[j]) if a != "-" and b != "-")
                length = sum((a != "-" and b != "-") for a, b in zip(aligned[i], aligned[j]))
                identity = 100.0 * matches / length if length > 0 else 0
                matrix[i, j] = round(identity, 2)
    return matrix

def plot_heatmap(matrix, order, output_file):
    """Plot and save heatmap of identities."""
    plt.figure(figsize=(10, 8))
    sns.heatmap(matrix, xticklabels=order, yticklabels=order,
                cmap="coolwarm", annot=True, fmt=".1f", vmin=60, vmax=100,
                cbar_kws={'label': 'Identity Percentage'})
    plt.xticks(rotation=90)
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.savefig(output_file, dpi=600)
    plt.close()

def main(arguments):
    input_fasta = arguments[arguments.index("--in")+1]
    output_folder = arguments[arguments.index("--out")+1]
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    aln_file = os.path.join(output_folder, "all_sequences.aln.fasta")

    # run MSA once
    run_mafft(input_fasta, aln_file)

    # load sequences & compute matrix
    seqs, order = load_sequences(input_fasta)
    matrix = compute_identity_matrix(aln_file, order)

    # write text summary
    txt_file = os.path.join(output_folder, "summary.txt")
    with open(txt_file, "w") as out:
        out.write("\t" + "\t".join(order) + "\n")
        for i, name in enumerate(order):
            row = [name] + [f"{matrix[i,j]:.1f}%" for j in range(len(order))]
            out.write("\t".join(row) + "\n")

    # plot heatmap
    fig_file = os.path.join(output_folder, "heatmap.png")
    plot_heatmap(matrix, order, os.path.join(output_folder, "heatmap.pdf"))
    plot_heatmap(matrix, order, os.path.join(output_folder, "heatmap.svg"))

if __name__ == "__main__":
    if "--in" in sys.argv and "--out" in sys.argv:
        main(sys.argv)
    else:
        sys.exit(__usage__)
