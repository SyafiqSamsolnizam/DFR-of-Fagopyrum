# DFR-of-Fagopyrum
This repository contains the scripts used in the study of "Expansion of DFR Gene Family in Fagopyrum".

##  TPM Boxplot Script

This script generates **boxplots of gene expression (TPM values)** from a transcript quantification table.  
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

Gene Sample1 Sample2 Sample3
FeDFR1 12.3 8.1 5.6
FeDFR2 30.2 25.7 19.4


**Gene query file** (optional, one gene ID per line):

FeDFR1
FeDFR2


**Tissue mapping file** (optional, tab-delimited):

Leaf Sample1,Sample2
Root Sample3


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
