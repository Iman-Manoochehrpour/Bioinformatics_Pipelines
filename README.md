# 🧬 Bioinformatics_Pipelines

**Reproducible R pipelines for RNA-seq and miRNA differential expression analysis using DESeq2 and edgeR**

---

### 📘 Overview  
This repository contains example R scripts demonstrating common bioinformatics workflows for RNA-seq and miRNA expression data analysis.  
All data used here are **simulated or public demo datasets** for educational and reproducibility purposes.  

The pipelines are based on analytical approaches developed in my master’s thesis:  
> *“Bioinformatics-based identification of candidate microRNAs and validation via qRT-PCR in bladder cancer patients for early detection.”*  

---

### ⚙️ Repository Structure  

| Folder | Description |
|--------|--------------|
| `DESeq2_demo/` | Differential expression workflow using DESeq2 |
| `edgeR_demo/` | Alternative workflow using edgeR |
| `Visualization/` | ggplot2-based data visualization scripts |
| `SampleData/` | Simulated expression and metadata files |
| `Utils/` | Helper functions for normalization and data cleaning |

---

### 🧩 Key Features  
- End-to-end RNA-seq and miRNA analysis pipeline in R  
- Normalization and statistical testing using **DESeq2** and **edgeR**  
- Volcano plot, PCA, and heatmap visualizations  
- Fully annotated and reproducible scripts for educational purposes  
- Modular structure for easy adaptation to custom datasets  

---

### 🧪 Example Workflows  

#### ▶️ 1. Differential Expression with DESeq2
```R
library(DESeq2)

# Load simulated data
counts <- read.csv("SampleData/simulated_counts.csv", row.names = 1)
meta <- read.csv("SampleData/sample_metadata.csv")

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = meta,
                              design = ~ Condition)

# Run differential expression analysis
dds <- DESeq(dds)
res <- results(dds)

# Save results
write.csv(as.data.frame(res), "DESeq2_demo/DESeq2_results.csv")

# Basic MA plot
plotMA(res, ylim = c(-4, 4))
