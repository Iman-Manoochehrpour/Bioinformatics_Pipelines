# 🧩 Scripts Documentation

This folder contains modular and reproducible R scripts used in the **Bioinformatics_Pipelines** project.  
Each script represents a distinct stage in the bioinformatics workflow designed for the thesis project:

> **Title:** Identification of candidate microRNAs for early detection of bladder cancer using bioinformatics and qRT-PCR validation  
> **Author:** *Iman Manoochehrpour, Department of Genetics, Kharazmi University*

---

## ⚙️ Script Overview

| Script Name | Description | Key Packages |
|--------------|--------------|---------------|
| **TCGA_download.R** | Downloads, preprocesses, and organizes simulated or example TCGA-BLCA data (mRNA, miRNA, metadata). Generates example datasets stored in `example_data/`. | `TCGAbiolinks`, `SummarizedExperiment`, `tidyverse` |
| **DESeq2_workflow.R** | Demonstrates a clean DESeq2 differential expression pipeline using simulated data. Exports normalized counts, statistical results, and a volcano plot to `visualization/`. | `DESeq2`, `ggplot2`, `EnhancedVolcano` |
| **WGCNA_pipeline.R** | Builds a weighted co-expression network (WGCNA) using simulated or filtered expression data. Produces module detection results and module–trait correlation plots. | `WGCNA`, `limma`, `ComplexHeatmap` |
| **install_packages.R** | Automatically installs and loads all required CRAN and Bioconductor packages for running the project. | `BiocManager`, `DESeq2`, `WGCNA`, etc. |

---

## 🧬 Pipeline Logic

The workflow is **modular**, allowing independent or sequential execution of scripts:

Bioinformatics_Pipelines/ 
│
├── scripts/
│ ├── install_packages.R # Install dependencies
│ ├── TCGA_download.R # Download or simulate expression data
│ ├── DESeq2_workflow.R # Perform differential expression analysis
│ └── WGCNA_pipeline.R # Construct co-expression networks
│
├── example_data/ # Example input data (auto-generated)
├── visualization/ # Generated plots and results
└── README.md # Project documentation


---

## 📊 Output Summary

| Step | Output Folder | Example Files |
|------|----------------|---------------|
| **TCGA_download.R** | `example_data/` | `simulated_counts.csv`, `sample_metadata.csv` |
| **DESeq2_workflow.R** | `visualization/` | `results_deseq2.csv`, `volcano_plot.png` |
| **WGCNA_pipeline.R** | `visualization/` | `module_colors.csv`, `module_trait_heatmap.png` |

---

## 🧠 Notes

- All example data in this repository are **synthetic** or **publicly available** to ensure that **no confidential thesis data** are shared.
- Each script is written in a **teaching-oriented** style with extensive inline comments.
- You can run each script independently via:
  ```R
  source("scripts/TCGA_download.R")
  source("scripts/DESeq2_workflow.R")
  source("scripts/WGCNA_pipeline.R")

Before running any script, make sure to execute:

```R
source("install_packages.R")

```

🧾 Citation

If you use or adapt these scripts, please cite the following thesis:

Manoochehrpour, I. (2025). Identification of candidate microRNAs for early detection of bladder cancer using bioinformatics and qRT-PCR validation.
Department of Genetics, Kharazmi University.

🧩 License

This folder and its contents are shared under the MIT License for educational and academic use.
Reproduction and modification are permitted with proper attribution.
