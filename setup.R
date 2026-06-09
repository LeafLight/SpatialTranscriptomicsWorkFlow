# ====================================================================
# One-Click Deployment Script for Spatial Transcriptomics Workflow
# ====================================================================

# 1. Install standard package managers
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

message("Checking and installing CRAN packages...")
cran_pkgs <- c("Seurat", "SeuratObject", "tidyverse", "patchwork", "ggplot2", "sp")
new_cran <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(new_cran)) install.packages(new_cran)

message("Checking and installing Bioconductor packages...")
bioc_pkgs <- c("SPOTlight", "STdeconvolve", "scran", "scater", "SingleCellExperiment")
new_bioc <- setdiff(bioc_pkgs, rownames(installed.packages()))
if (length(new_bioc)) BiocManager::install(new_bioc, update = FALSE)

message("Checking and installing Giotto Suite...")
if (!requireNamespace("Giotto", quietly = TRUE)) {
  remotes::install_github("drieslab/Giotto@suite")
}

message("All R dependencies installed successfully!")