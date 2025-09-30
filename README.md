# Rnaseq — Reproducible Code for "Transcriptome-based identification of bottlenecks in the paclitaxel biosynthesis Pathway" Article

A complete package for gene expression analysis of KEGG pathways with advanced visualization and statistical analysis capabilities.

Differential Expression Analysis (DEG) with DESeq2

KEGG enrichment and results visualization

Module activity in different samples

Analysis of precursor pathways with bubble plots

Diffusion-quality plots (Volcano, PCA, Heatmap)

Support for Counts and TPM data



This repository contains the **code** used in the study "Transcriptome-based identification of bottlenecks in the paclitaxel biosynthesis Pathway".
**No data are included in this repository.** Reviewers can run the code by placing the required inputs locally and updating the USER CONFIG block in `scripts/R/Pathway Analysis Pipeline.R`.



## How to run (quick)
### R
1) Install deps:
```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
cran <- c("readxl","openxlsx","dplyr","tibble","ggplot2","ggrepel","pheatmap","ggpubr","reshape2","tidyr")
bioc <- c("DESeq2","clusterProfiler","enrichplot")
install.packages(setdiff(cran, rownames(installed.packages())))
BiocManager::install(setdiff(bioc, rownames(installed.packages())))
