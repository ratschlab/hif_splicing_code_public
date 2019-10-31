# Characterisation of HIF-dependent alternative isoforms in pancreatic cancer
Splicing analysis scripts for the paper "Characterisation of HIF-dependent alternative isoforms in pancreatic cancer".

# Author
Natalie R. Davidson <natalie.davidson@inf.ethz.ch>

# Requirements
This analysis script was tested on R 3.5.1. External libraries used are: data.table, RColorBrewer, DESeq2, ggplot2, MASS, pheatmap, plyr.

# Usage
The script `splicing_analysis/differential_splicing.R` has 4 command-line parameters. 
- `indir`: The directory containing all splice files 
- `RESULTS_FOLDER`: The directory where you would like all results files to be output.
- `run_splice_only`: TRUE if you would like the analysis to identify HIF dependent splice events that are also independent of expression. FALSE if you would like the events to significantly change in both splicing and expression.
- `translation_file`: A file which translates which contains two columns `ensembl_gene_id` to `hgnc_symbol`.
