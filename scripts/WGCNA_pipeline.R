# ================================================================
# Title: Weighted Gene Co-expression Network Analysis (WGCNA)
# Author: Iman Manoochehrpour
# Institution: Kharazmi University, Department of Genetics
# Description: Template workflow for gene co-expression network construction
# Input: example_data/simulated_counts.csv + sample_metadata.csv
# Output: module_colors.csv + network_dendrogram.png + module_trait_heatmap.png
# ================================================================

# ----------------------------
# Step 1. Load required packages
# ----------------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install WGCNA if missing
if (!requireNamespace("WGCNA", quietly = TRUE))
    BiocManager::install("WGCNA")

# Load libraries
library(WGCNA)
library(ggplot2)
library(reshape2)

# Disable multithreading for reproducibility
allowWGCNAThreads(nThreads = 1)

# ----------------------------
# Step 2. Import and preprocess data
# ----------------------------

# Read simulated expression and metadata
data_expr <- read.csv("../example_data/simulated_counts.csv", row.names = 1)
meta <- read.csv("../example_data/sample_metadata.csv")

# Transpose expression data: genes as columns, samples as rows
datExpr <- as.data.frame(t(data_expr))
rownames(meta) <- meta$Sample

# Check for missing values
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

cat("\n✅ Expression data loaded successfully.\n")
cat("Samples:", nrow(datExpr), " | Genes:", ncol(datExpr), "\n\n")

# ----------------------------
# Step 3. Choose soft-thresholding power
# ----------------------------

powers <- c(1:10)
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = "signed")

# Choose optimal power
bestPower <- sft$powerEstimate
if (is.na(bestPower)) bestPower <- 6
cat("Selected soft-thresholding power:", bestPower, "\n")

# ----------------------------
# Step 4. Build the network and identify modules
# ----------------------------

net <- blockwiseModules(
  datExpr,
  power = bestPower,
  TOMType = "signed",
  minModuleSize = 20,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 3
)

# ----------------------------
# Step 5. Visualize dendrogram and module colors
# ----------------------------

module_colors <- labels2colors(net$colors)

png("../visualization/network_dendrogram.png", width = 800, height = 600)
plotDendroAndColors(
  net$dendrograms[[1]],
  module_colors[net$blockGenes[[1]]],
  "Module colors",
  main = "Gene Co-expression Network Dendrogram"
)
dev.off()

# Save module color table
write.csv(data.frame(Gene = colnames(datExpr), ModuleColor = module_colors),
          "../visualization/module_colors.csv", row.names = FALSE)

cat("📁 Network and module outputs generated.\n")

# ----------------------------
# Step 6. Module–Trait Correlation Analysis
# ----------------------------

# Convert trait (Condition) to numeric
traitData <- meta
traitData$Condition <- ifelse(traitData$Condition == "Tumor", 1, 0)

# Compute module eigengenes
MEs <- moduleEigengenes(datExpr, colors = module_colors)$eigengenes
MEs <- orderMEs(MEs)

# Calculate correlations
moduleTraitCor <- cor(MEs, traitData$Condition, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr))

# Combine correlation and p-value into one data frame
cor_df <- data.frame(
  Module = rownames(moduleTraitCor),
  Correlation = as.vector(moduleTraitCor),
  Pvalue = as.vector(moduleTraitPvalue)
)

write.csv(cor_df, "../visualization/module_trait_correlations.csv", row.names = FALSE)

# ----------------------------
# Step 7. Visualize module–trait relationships
# ----------------------------

# Create a heatmap of module–trait correlations
textMatrix <- paste0(
  signif(moduleTraitCor, 2), "\n(",
  signif(moduleTraitPvalue, 1), ")"
)
dim(textMatrix) <- dim(moduleTraitCor)

png("../visualization/module_trait_heatmap.png", width = 800, height = 600)
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = "Tumor vs Normal",
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.8,
  zlim = c(-1,1),
  main = "Module–Trait Relationships"
)
dev.off()

cat("📈 Module–Trait correlation heatmap generated.\n")

# ----------------------------
# Step 8. Summary
# ----------------------------

cat("\n✅ WGCNA pipeline completed successfully!\n")
cat("Outputs saved in: ../visualization/\n")
cat("- network_dendrogram.png\n")
cat("- module_colors.csv\n")
cat("- module_trait_correlations.csv\n")
cat("- module_trait_heatmap.png\n")

# ================================================================
# End of Script
# ================================================================
