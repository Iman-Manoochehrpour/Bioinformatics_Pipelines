# 🧬 Bioinformatics_Pipelines

**Author:** Iman Manoochehrpour  
**Institution:** Department of Genetics, Kharazmi University  
**Project:** Identification of candidate microRNAs for early detection of bladder cancer using bioinformatics and qRT-PCR validation  

---

## 📖 Overview

This repository contains reproducible bioinformatics workflows developed as part of my M.Sc. thesis project.  
The analyses focus on **transcriptomic and regulatory network profiling** in bladder cancer, including:

- Retrieval and preprocessing of TCGA data  
- Differential expression analysis using `DESeq2`  
- Network-based co-expression analysis using `WGCNA`  
- Pathway and enrichment visualization  

The goal of this repository is to demonstrate **structured, modular, and reproducible analysis pipelines** in R for RNA-seq and miRNA-seq datasets.

---

## 📁 Project Structure

Bioinformatics_Pipelines/
│
├── README.md <- Project documentation (this file)
│
├── scripts/ <- Contains analysis scripts (DESeq2, WGCNA, TCGA queries)
│
├── example_data/ <- Simulated or example datasets for demonstration
│
└── visualization/ <- Example plots (e.g., volcano plots, heatmaps)


---

## ⚙️ Planned Scripts

| Script | Description |
|---------|-------------|
| `scripts/TCGA_download.R` | Template for querying TCGA data via `TCGAbiolinks` |
| `scripts/DESeq2_workflow.R` | Example workflow for differential expression analysis |
| `scripts/WGCNA_pipeline.R` | Example workflow for network co-expression analysis |

Each script will include simulated data and annotated R code blocks for clarity and reusability.

---

## 🧪 Example Workflows (Coming Soon)

Workflows will include:

1. **TCGA Data Retrieval** – using `GDCquery()` and `GDCprepare()`  
2. **DESeq2 Differential Expression** – normalization, dispersion estimation, and result extraction  
3. **WGCNA Co-expression Network** – module detection and trait correlation  
4. **Visualization** – generating volcano plots, module heatmaps, and enrichment summaries  

Example:

```R
# Volcano Plot Example
library(ggplot2)

res$significant <- res$padj < 0.05 & abs(res$log2FoldChange) > 1
ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.6) +
  theme_classic() +
  labs(title = "Volcano Plot of Differentially Expressed Genes",
       x = "log2(Fold Change)",
       y = "-log10(Adjusted p-value)") +
  scale_color_manual(values = c("grey60", "red"))
```
##  🧰 Dependencies and Installation

All workflows are implemented in R (≥ 4.2).
The following R packages are required and can be installed via Bioconductor and CRAN:
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Core Bioconductor packages
BiocManager::install(c("TCGAbiolinks", "DESeq2", "WGCNA", "EnhancedVolcano", "clusterProfiler"))

# CRAN packages
install.packages(c("ggplot2", "dplyr", "tibble", "readr"))

💡 Each script will include specific package loading instructions and version details to ensure reproducibility.

##  🧬 Future Directions

Planned extensions for this repository include:

🔹 Integration of DNA methylation and miRNA–mRNA interaction analyses

🔹 Inclusion of target prediction and enrichment analysis pipelines

🔹 Automated workflow documentation using RMarkdown

🔹 Exportable tabular and graphical outputs for downstream reporting

##  🧠 Author Notes

This repository is developed for educational and research demonstration purposes.
All included data are simulated or publicly available, ensuring that no unpublished or confidential thesis data are shared.

For inquiries, collaboration, or citation requests, please contact:
📧 iman.manoochehrpour @ gmail . com

##  📜 License

This repository is distributed under the MIT License.
You are free to reuse, modify, and share the included scripts for non-commercial research and educational purposes, provided proper attribution is given.
