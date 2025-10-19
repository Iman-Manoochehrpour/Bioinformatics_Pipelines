############################################################
# 🧬 install_packages.R
# Author: Iman Manoochehrpour
# Project: Bioinformatics_Pipelines
# Purpose: Automatically install and load all required packages
############################################################

# ---- Helper function ----
install_if_missing <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("📦 Installing", pkg, "..."))
      install.packages(pkg, dependencies = TRUE)
    } else {
      message(paste("✅", pkg, "is already installed."))
    }
  }
}

# ---- Bioconductor installer ----
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# ---- CRAN packages ----
cran_packages <- c(
  "ggplot2",
  "tidyverse",
  "magrittr",
  "reshape2",
  "gridExtra",
  "viridis",
  "cowplot",
  "ComplexHeatmap"
)

# ---- Bioconductor packages ----
bioc_packages <- c(
  "TCGAbiolinks",
  "SummarizedExperiment",
  "DESeq2",
  "WGCNA",
  "limma",
  "EnhancedVolcano",
  "clusterProfiler",
  "org.Hs.eg.db",
  "AnnotationHub",
  "biomaRt"
)

# ---- Optional utilities ----
optional_packages <- c(
  "enrichR",
  "msigdbr",
  "maftools",
  "bestNormalize",
  "ggpubr"
)

# ---- Install all ----
message("🔹 Installing CRAN packages ...")
install_if_missing(cran_packages)

message("🔹 Installing Bioconductor packages ...")
BiocManager::install(bioc_packages, ask = FALSE, update = FALSE)

message("🔹 Installing optional analysis tools ...")
install_if_missing(optional_packages)

# ---- Load all packages ----
all_packages <- c(cran_packages, bioc_packages, optional_packages)
for (pkg in all_packages) {
  suppressMessages(library(pkg, character.only = TRUE))
}

message("\n✅ All packages successfully installed and loaded!")
############################################################
