#!/usr/bin/env Rscript

# Basic data processing and visualization packages
install.packages(c(
  "dplyr", 
  "tidyr", 
  "ggplot2", 
  "reshape2", 
  "data.table", 
  "RColorBrewer"
))

# Bioinformatics related packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "GenomicRanges",
  "rtracklayer",
  "DESeq2",
  "edgeR",
  "diffHic",
  "HiCcompare",
  "topGO",
  "ComplexHeatmap",
  "chromVAR"
))

# Chromatin structure analysis related
BiocManager::install(c(
  "HiTC",
  "GOTHiC",
  "HiCseg",
  "TADCompare"
))

# Packages for specific functionalities
install.packages(c(
  "parallel", 
  "doParallel", 
  "optparse"
)) 