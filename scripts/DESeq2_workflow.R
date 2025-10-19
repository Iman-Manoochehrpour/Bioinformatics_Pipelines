# ================================================================
# Title: Differential Expression Analysis (DESeq2)
# Author: Iman Manoochehrpour
# Institution: Kharazmi University, Department of Genetics
# Description: Template workflow for differential gene expression analysis
# Input: example_data/simulated_counts.csv, example_data/sample_metadata.csv
# Output: results_deseq2.csv + volcano_plot.png
# ================================================================

# ----------------------------
# Step 1. Load required packages
# ----------------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install DESeq2 if not installed
if (!requireNamespace("DESeq2", quietly = TRUE))
    BiocManager::install("DESeq2")

# Load dependencies
library(DESeq2)
library(ggplot2)
library(dplyr)

# ----------------------------
# Step 2. Import example data
# ----------------------------

counts <- read.csv("../example_data/simulated_counts.csv", row.names = 1)
meta <- read.csv("../example_data/sample_metadata.csv")

# Confirm dimensions
cat("\n📊 Data summary:\n")
cat("Genes:", nrow(counts), "\nSamples:", ncol(counts), "\n\n")

# ----------------------------
# Step 3. Prepare DESeq2 dataset
# ----------------------------

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = meta,
  design = ~ Condition
)

# Filter low-count genes (keep genes with >10 total counts)
dds <- dds[rowSums(counts(dds)) > 10, ]

# ----------------------------
# Step 4. Run DESeq2
# ----------------------------

dds <- DESeq(dds)
res <- results(dds)
res <- lfcShrink(dds, coef = 2, type = "apeglm")  # shrink log2FC for stability

# Summary
cat("✅ DESeq2 analysis complete.\n")
cat("Significant genes (padj < 0.05):", sum(res$padj < 0.05, na.rm = TRUE), "\n\n")

# ----------------------------
# Step 5. Save output
# ----------------------------

# Create results directory if not exists
if (!dir.exists("../visualization")) dir.create("../visualization")

# Save results
write.csv(as.data.frame(res), "../visualization/results_deseq2.csv", row.names = TRUE)

# ----------------------------
# Step 6. Visualization: Volcano Plot
# ----------------------------

res_df <- as.data.frame(res)
res_df$significant <- res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1

p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.7, size = 2) +
  theme_classic() +
  labs(
    title = "Volcano Plot of Differentially Expressed Genes",
    x = "log2(Fold Change)",
    y = "-log10(Adjusted p-value)"
  ) +
  scale_color_manual(values = c("grey60", "red")) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("../visualization/volcano_plot.png", plot = p, width = 6, height = 5, dpi = 300)

cat("📁 Outputs generated:\n")
cat("- visualization/results_deseq2.csv\n")
cat("- visualization/volcano_plot.png\n")

# ================================================================
# End of Script
# ================================================================
