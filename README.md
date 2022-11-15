# CSlncR
Inferring cell‑specific lncRNA regulation with single-cell RNA-sequencing data in the developing human neocortex
![](https://github.com/linxi159/CSlncR/blob/main/figures/Figure_1.tif) 

## Description of each directory and each file
data: The preprocessed data from real scRNA-seq data in GEO and other data.

figures: The plot for CSlncR.

Exp_247_lncRNAs_10208_mRNAs_276_single_cells_GSE71315.RData: Matched lncRNA and mRNA expression data across 276 single cells in the human neocortex, Putative lncRNA-target binding information.

CSlncR_network_GSE71315.RData: 276 cell-specific lncRNA regulatory networks.

CSlncR.R: Utility functions for exploring cell-specific lncRNA regulation.

step1_data_preprocessing_dividing_mRNA_lncRNA.py: Dividing mRNAs and lncRNAs in the scRNA-seq data.

step2_case_study.R: Running scripts for exploring cell-specific lncRNA regulation.

## The usage of CSlncR
Paste all files into a single folder (set the folder as the directory of R environment), the workflow of CSlncR is implemented in CSlncR.R. The users can simply run the scripts as follows.
```
source("step2_case_study.R")
```

## Quick example to use CSlncR
For identifying cell-specific lncRNA regulation, users should prepare lncRNA and mRNA single-cell co-expression data. Paste the datasets and our source file (CSlncR.R) into a single folder (set the folder as the directory of R environment), users can use the following scripts to identify cell-specific lncRNA regulation. For convenience, the datasets prepared for users are from our datasets (Exp_247_lncRNAs_10208_mRNAs_276_single_cells_GSE71315.RData).
```
## Load required R packages, please firstly install the following R packages before running scripts
library(pracma)
library(WGCNA)
library(igraph)
library(miRspongeR)
library(biclique)
library(corrplot)
library(dendextend)
library(ggplot2)
library(ggsignif)
library(cowplot)
library(clusterProfiler)
library(msigdbr)
library(vroom)
library(doParallel)

## Load utility functions
source("CSlncR.R")

##############################################################################################################################
############################################## <<1>> Loading data  ###########################################################
##############################################################################################################################
## Load prepared datasets
ASD <- read.csv("./data/ASD_gene_lncRNAs_mRNAs.csv")
ASD <- as.matrix(ASD)

## Load the preprocessing single-cell sequencing data including filtering genes with expression <= 10 cells, and log e(x+1).
lncRNA_scRNA_raw  <- read.csv("./data/GSE71315/GSE71315_scell_ncounts_genes_thresh_lncRNA.csv")
rownames(lncRNA_scRNA_raw ) <- lncRNA_scRNA_raw[,1]
lncRNA_scRNA_raw  <- lncRNA_scRNA_raw[,-1]
lncRNA_scRNA_raw   <- t(lncRNA_scRNA_raw)
lncRNA_scRNA_raw   <- as.data.frame(lncRNA_scRNA_raw)

mRNA_scRNA_raw <- read.csv("./data/GSE71315/GSE71315_scell_ncounts_genes_thresh_mRNA.csv")
rownames(mRNA_scRNA_raw) <- mRNA_scRNA_raw[,1]
mRNA_scRNA_raw <- mRNA_scRNA_raw[,-1]
mRNA_scRNA_raw  <- t(mRNA_scRNA_raw)
mRNA_scRNA_raw <- as.data.frame(mRNA_scRNA_raw)

## compute the average expression values of duplicate genes and remove genes with constant expression values in all cells
# Transformation using log2(x+1)
lncRNA_scRNA_norm <- log2(lncRNA_scRNA_raw+1)
mRNA_scRNA_norm <- log2(mRNA_scRNA_raw+1) 

# Compute the average expression values of duplicate genes
lncRNA_scRNA_norm_average <- Averg_Duplicate(lncRNA_scRNA_norm)
mRNA_scRNA_norm_average <- Averg_Duplicate(mRNA_scRNA_norm)

# Remove genes with constant expression values in all cells
lncRNA_scRNA_norm_sd <- unlist(lapply(seq(dim(lncRNA_scRNA_norm_average)[2]), function(i) sd(lncRNA_scRNA_norm_average[, i])))
lncRNA_scRNA_norm_filter <- lncRNA_scRNA_norm_average[, which(lncRNA_scRNA_norm_sd > 0)]
mRNA_scRNA_norm_sd <- unlist(lapply(seq(dim(mRNA_scRNA_norm_average)[2]), function(i) sd(mRNA_scRNA_norm_average[, i])))
mRNA_scRNA_norm_filter <- mRNA_scRNA_norm_average[, which(mRNA_scRNA_norm_sd > 0)]

# loading prediction lcnRNA-target to filter to human brain lncRNA-mRNA 
lncRNA_mRNA_pre_1 <- unique(read.csv("./data/lncRNA_mRNA-pre.csv"))
lncRNA_mRNA_pre_lncRNA <- as.vector(as.matrix(lncRNA_mRNA_pre_1[1]))
lncRNA_name <- colnames(lncRNA_scRNA_norm_filter)
lncRNA_mRNA_pre_2 <- lncRNA_mRNA_pre_1[which(lncRNA_mRNA_pre_lncRNA %in% lncRNA_name),]
lncRNA_mRNA_pre_mRNA <- as.vector(as.matrix(lncRNA_mRNA_pre_2[2]))
mRNA_name <- colnames(mRNA_scRNA_norm_filter)
lncRNA_mRNA_pre_3 <- lncRNA_mRNA_pre_2[which(lncRNA_mRNA_pre_mRNA %in% mRNA_name),]
lncRNA_mRNA_pre <- lncRNA_mRNA_pre_3 
# loading experimental lcnRNA-target
lncRNA_mRNA_exp <- read.csv("./data/lncRNA_mRNA-exp.csv")
lncRNA_mRNA_exp <- as.matrix(lncRNA_mRNA_exp)
lncRNA_mRNA_exp <- unique(lncRNA_mRNA_exp)

##############################################################################################################################
####################################### <<2>> Exploring cell-specific lncRNA regulation ######################################
##############################################################################################################################
## Discovering cell-specific lncRNA-mRNA regulatory network   
CSlncR_network_null <- CSlncR_net(lncRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter)

```
