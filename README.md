# 🧬 Bioinformatics_Pipelines  
**Author:** Iman Manoochehrpour  
**Institution:** Department of Genetics, Kharazmi University  
**Project:** Identification of candidate microRNAs for early detection of bladder cancer using bioinformatics and qRT-PCR validation  

---

## 📖 Overview  

This repository contains reproducible **bioinformatics workflows** developed as part of my M.Sc. thesis project.  
The analyses focus on transcriptomic and regulatory network profiling in **bladder cancer (TCGA-BLCA)**, including:  

- Retrieval and preprocessing of TCGA data  
- Differential expression analysis using DESeq2  
- Network-based co-expression analysis using WGCNA  
- Pathway and enrichment visualization  

The goal of this repository is to demonstrate **structured, modular, and reproducible** analysis pipelines in R for RNA-seq and miRNA-seq datasets.

---

## 📁 Project Structure  

Bioinformatics_Pipelines/
│
├── README.md # Project documentation (this file)
├── install_packages.R # Automatic package installation script
│
├── scripts/ # Contains main analysis scripts
│ ├── README.md
│ ├── TCGA_download.R
│ ├── DESeq2_workflow.R
│ └── WGCNA_pipeline.R
│
├── example_data/ # Example or simulated datasets for demonstration
│ └── README.md
│
├── visualization/ # Visualization outputs (plots, tables)
│ └── README.md
│
└── LICENSE # Optional: project license file


---

## ⚙️ Getting Started  

### 🔹 Step 1 — Install All Dependencies  
Before running any scripts, install all required packages by executing the following command in R:  
```R
source("install_packages.R")
```

🔹 Step 2 — Run Example Pipelines

Each workflow can be executed independently using the provided simulated datasets:
```R
source("scripts/TCGA_download.R")     # Retrieves or simulates TCGA-like data
source("scripts/DESeq2_workflow.R")   # Runs DESeq2 differential expression analysis
source("scripts/WGCNA_pipeline.R")    # Performs co-expression network analysis
```
All generated results and plots will be automatically saved in the visualization/ folder.


🧪 Example Workflows

The following analyses are demonstrated using example data:

TCGA Data Retrieval – using GDCquery() and GDCprepare()

DESeq2 Differential Expression – normalization, dispersion estimation, and statistical testing

WGCNA Co-expression Network – module detection and module–trait correlations

Visualization – generation of volcano plots, module heatmaps, and enrichment summaries

Example Volcano Plot:
```R
library(ggplot2)

res$significant <- res$padj < 0.05 & abs(res$log2FoldChange) > 1
ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.6) +
  theme_classic() +
  labs(title = "Volcano Plot of Differentially Expressed Genes",
       x = "log2(Fold Change)",
       y = "-log10(Adjusted p-value)")
```

📦 Dependencies and Installation

All workflows are implemented in R (≥ 4.2) and rely on the following core packages:

Category	Packages
Core	ggplot2, tidyverse, magrittr, reshape2, cowplot
Bioinformatics	TCGAbiolinks, SummarizedExperiment, DESeq2, WGCNA, limma, EnhancedVolcano, clusterProfiler
Annotation / Enrichment	org.Hs.eg.db, biomaRt, AnnotationHub, enrichR, msigdbr
Visualization	ComplexHeatmap, viridis, ggpubr

All dependencies are automatically handled via:
```R
source("install_packages.R")
```
🧬 Future Directions

Planned extensions include:

Integration of DNA methylation and miRNA–mRNA interaction analyses

Automated pipeline documentation using RMarkdown

Exportable results in tabular and graphical formats

🧠 Author Notes

This repository is designed for educational and research purposes.
All data and scripts use synthetic or publicly available resources, ensuring that no unpublished or confidential thesis data are exposed.

For questions or collaboration requests, please contact:
📧 iman.manoochehrpour [at] gmail [dot] com

📜 License

This repository is shared under the MIT License for educational and academic demonstration purposes.
All scripts may be reused with appropriate citation.
