# ================================================================
# Title: TCGA Data Retrieval and Preprocessing (Template)
# Author: Iman Manoochehrpour
# Institution: Kharazmi University, Department of Genetics
# Description: Example workflow for querying and preparing TCGA data
# Note: This script uses publicly available data via TCGAbiolinks
# ================================================================

# ----------------------------
# Step 1. Load required packages
# ----------------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Load or install essential Bioconductor packages
BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment"), ask = FALSE)

# Load packages
library(TCGAbiolinks)
library(SummarizedExperiment)

# ----------------------------
# Step 2. Define project and parameters
# ----------------------------

# Example project: TCGA-BLCA (Bladder Urothelial Carcinoma)
project_id <- "TCGA-BLCA"

# Define query parameters for mRNA expression data (HTSeq-Counts)
query <- GDCquery(
    project = project_id,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "Star - Counts",
    sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

# ----------------------------
# Step 3. Download and prepare data
# ----------------------------

# Download data (this step fetches metadata + count files)
# ⚠️ Note: Uncomment the next line to execute (requires active internet)
GDCdownload(query)

# Prepare data as a SummarizedExperiment object
# ⚠️ Note: Uncomment when actual download is done
data_se <- GDCprepare(query)

# ----------------------------
# Step 4. Preview and structure of the dataset
# ----------------------------

# Example simulated dataset structure
# (Used for demonstration when actual TCGA data is not downloaded)

set.seed(123)
simulated_counts <- matrix(
    rpois(1000, lambda = 10),
    nrow = 100,
    ncol = 10,
    dimnames = list(
        paste0("Gene_", 1:100),
        paste0("Sample_", 1:10)
    )
)

simulated_metadata <- data.frame(
    Sample = paste0("Sample_", 1:10),
    Condition = rep(c("Tumor", "Normal"), each = 5)
)

# ----------------------------
# Step 5. Save example files
# ----------------------------

# Create output directory if not exists
if (!dir.exists("../example_data")) dir.create("../example_data")

# Save simulated count data
write.csv(simulated_counts, "../example_data/simulated_counts.csv", row.names = TRUE)
write.csv(simulated_metadata, "../example_data/sample_metadata.csv", row.names = FALSE)

# ----------------------------
# Step 6. Summary
# ----------------------------

cat("\n✅ TCGA data template completed successfully!\n")
cat("Files generated:\n")
cat("- example_data/simulated_counts.csv\n")
cat("- example_data/sample_metadata.csv\n")
cat("\nYou can now use these files in DESeq2 or WGCNA pipelines.\n")

# ================================================================
# End of script
# ================================================================
