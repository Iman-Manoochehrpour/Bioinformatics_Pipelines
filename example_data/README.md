# 🧬 Example Data

This folder contains **synthetic datasets** created for demonstration and reproducibility testing of the pipelines developed in the project  
**“Identification of candidate microRNAs for early detection of bladder cancer using bioinformatics and qRT-PCR validation.”**

All data here are **artificially generated** using R scripts and do **not** correspond to any real biological samples or patient-derived data.

---

## 📁 Files Overview

| File Name | Description | Example Use |
|------------|-------------|--------------|
| `simulated_counts.csv` | Randomly generated RNA-seq–like gene expression matrix (100 genes × 20 samples). Values range between 10–500. | Input for **DESeq2_workflow.R** |
| `sample_metadata.csv` | Metadata table for each sample, including `Condition`, `Gender`, and `Age`. 10 Tumor + 10 Normal samples. | Used for grouping and design matrix generation |
| `gene_annotation.csv` | Synthetic gene annotations with gene `Symbol` and `Biotype` (`protein_coding`, `lncRNA`, `miRNA`). | Used for labeling plots and functional filtering |

---

## 🧪 How These Files Were Generated

All datasets were automatically generated using the script  
[`generate_example_data.R`](../generate_example_data.R)  
located in the main project directory.

Example snippet:

```R
source("generate_example_data.R")
```

This script produces all three .csv files inside this folder.
The data are randomized but structured to mimic realistic transcriptomic patterns, including:

Mixed gene types (protein_coding, lncRNA, miRNA)

Balanced Tumor and Normal samples

Heterogeneous expression values

📊 Example Preview
GeneID	Symbol	Biotype
Gene_001	A452	protein_coding
Gene_002	MIR324	miRNA
Gene_003	LNC2043	lncRNA
Gene_004	K889	protein_coding
Gene_005	MIR7890	miRNA
🧠 Usage Notes

These datasets are meant only for demonstration and testing of the analysis workflows.

You may replace these synthetic files with your own real data while maintaining the same structure (column names and formatting).

Each analysis script automatically reads from this folder, so ensure file names remain identical.

⚠️ Ethical and Data Privacy Statement

No confidential or unpublished data from the M.Sc. thesis have been shared in this repository.
All examples are computationally generated and publicly reproducible.

🧩 Related Scripts
Script	Purpose
scripts/DESeq2_workflow.R
	Uses simulated_counts.csv and sample_metadata.csv for differential expression analysis
scripts/WGCNA_pipeline.R
	Demonstrates module detection and trait correlation using the same data
generate_example_data.R
	Creates all files in this folder
📜 License

This folder and its contents are released under the MIT License for educational and research demonstration purposes.
You may freely reuse or adapt these datasets with proper attribution.
