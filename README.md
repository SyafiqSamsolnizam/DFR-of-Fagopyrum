# DFR-of-Fagopyrum
This repository contains the scripts and data used in the study of "Expansion of DFR Gene Family in Fagopyrum".

##  TPM Boxplot Script

This script generates **boxplots of gene expression (TPM values)** from a transcript quantification table by Kallisto.  
It is useful for visualizing expression patterns across samples and tissues, with options for customization.  

###  Features
- Reads a **TPM table (TSV)** where:  
  - First column = Gene names  
  - Other columns = TPM values for each sample  
- Optional gene query file to subset genes of interest  
- Optional tissue mapping file to group samples by tissue and color them in the plot  
- Customizable aesthetics:  
  - Title  
  - Point size  
  - Outlier display and size  
  - Box width (spacing)  
  - Tissue display order  
- Legend shows **tissue sample counts**  
- Supports multiple output formats (`.png`, `.svg`, etc.)  

###  Input File Formats
**TPM table** (TSV):

```
Gene Sample1 Sample2 Sample3
FeDFR1 12.3 8.1 5.6
FeDFR2 30.2 25.7 19.4
```

**Gene query file** (optional, one gene ID per line):

```
FeDFR1
FeDFR2
```


**Tissue mapping file** (optional, tab-delimited):

```
Leaf Sample1,Sample2
Root Sample3,Sample4
```


###  Usage
```bash
python3 boxplot_script.py \
  -i TPM_table.tsv \
  -q genes.txt \
  -t tissues.tsv \
  --tissue_order "Leaf,Stem,Root" \
  --title "Expression of DFR Genes" \
  --point_size 6 \
  --outlier_size 4 \
  --box_width 0.5 \
  -o expression_boxplot.png
```

##  PlantPAN Motif Visualizer

This script visualizes **transcription factor binding sites (TFBS)** from **PlantPAN output files** across promoter regions.  
It plots promoter regions as horizontal bars and marks motif positions with colored lines, producing both **PDF and SVG** outputs.  

###  Features
- Reads PlantPAN result tables (TSV format).  
- Preserves **gene order** from the input file.  
- Handles missing TFBS names by replacing with TF family or TFBS ID.  
- Optional filtering to include only specific TFBS IDs.  
- Color-codes motifs consistently across genes.  
- Automatically scales promoter lengths to fit all motifs.  
- Generates **publication-ready motif distribution plots** (PDF + SVG).  

###  Input File Format
The script expects a **PlantPAN result file (TSV)** with at least these columns:  
- `Sequence ID` → Gene identifier  
- `TFBS ID` → Motif identifier  
- `TFBS Name` (or uses TF family if missing)  
- `Posistion` → Position in promoter (bp)  
- `Binding sequence` → Matched sequence  

Example (simplified):

```
Sequence ID	TFBS ID	TFBS Name	Posistion	Strand	Similar Score	Binding sequence	TF family	TF ID
FeDFR1	TFmatrixID_0002		299	+	0.959	caTAATTttt	AT-Hook	AT1G63480
FeDFR1	TFmatrixID_0002		518	-	0.949	aaaAATTAat	AT-Hook	AT1G63480
FeDFR1	TFmatrixID_0002		1101	+	0.951	gaTAATTtat	AT-Hook	AT1G63480
```

###  Usage
```bash
python3 plantpan_visualiser.py \
  -i PlantPAN_results.tsv \
  -o motif_plots \
  -t TF001 TF002
```
###  Output

The script saves two files in the specified output folder (motif_plots/ by default):  
motif_plot.pdf  
motif_plot.svg  

Both show:  
Promoter regions (grey boxes)  
Motif positions (colored vertical bars)  
A legend with motif names and representative sequences  



##  Pairwise Sequence Identity Heatmap

This script computes pairwise **sequence identity percentages** from a set of FASTA sequences and visualizes them in a **heatmap**.  
It uses **MAFFT** for multiple sequence alignment, calculates pairwise identities, and exports both numeric results and plots.  

###  Features
- Loads sequences directly from a FASTA file.  
- Performs **global multiple sequence alignment** using MAFFT (`--genafpair --maxiterate 1000`).  
- Computes pairwise **percentage identity** ignoring gaps.  
- Outputs both a **summary table (TXT)** and a **heatmap (PDF + SVG)**.  
- Annotated heatmap shows identity values and uses a **coolwarm color gradient**.  

###  Input
A standard FASTA file containing nucleotide or protein sequences.  


###  Usage
```bash
python3 pairwise_heatmap.py \
    --in input_sequences.fasta \
    --out results/
```

### Output

The script generates the following files inside the chosen output folder:  

all_sequences.aln.fasta → aligned sequences (via MAFFT)  

summary.txt → tab-separated matrix of pairwise identities  

heatmap.pdf → publication-quality heatmap  

heatmap.svg → vector-format heatmap  
