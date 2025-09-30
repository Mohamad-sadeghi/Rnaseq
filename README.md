# Rnaseq — Codes for *Transcriptome-based identification of bottlenecks in the paclitaxel biosynthesis pathway*

A compact package for gene expression analysis of KEGG pathways with visualization and statistics.

- Differential Expression Analysis (DEG) with DESeq2  
- KEGG enrichment and results visualization  
- Module activity across samples  
- Precursor pathway analysis with bubble plots  
- Publication-quality plots (Volcano, PCA, Heatmap)  
- Supports both Counts and TPM (TPM optional for plots)


This repository contains the **code** used in the study "Transcriptome-based identification of bottlenecks in the paclitaxel biosynthesis Pathway".
> **Code only** is hosted here. **No data** are included in this repository. Reviewers can run the code by placing required inputs locally and editing the **USER CONFIG** block at the top of `Scripts/R/kegg_pipeline.R`.


Data availability

Data will be provided with the journal’s supplementary materials.
Placeholder: add DOI/URL here after acceptance.


## How to run (quick)
### R
1) Install deps:
```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
cran <- c("readxl","openxlsx","dplyr","tibble","ggplot2","ggrepel","pheatmap","ggpubr","reshape2","tidyr")
bioc <- c("DESeq2","clusterProfiler","enrichplot")
install.packages(setdiff(cran, rownames(installed.packages())))
BiocManager::install(setdiff(bioc, rownames(installed.packages())))
