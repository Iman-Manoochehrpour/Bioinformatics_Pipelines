############################################################
# 🧬 generate_example_data.R
# Author: Iman Manoochehrpour
# Project: Bioinformatics_Pipelines
# Purpose: Generate synthetic example datasets for demonstration
############################################################

# ---- Create folder if it doesn't exist ----
if (!dir.exists("example_data")) {
  dir.create("example_data")
  message("📁 Folder 'example_data/' created.")
}

set.seed(123)

############################################################
# 1️⃣ Simulated gene expression counts (RNA-seq-like)
############################################################
genes <- paste0("Gene_", sprintf("%03d", 1:100))
samples <- paste0("Sample_", sprintf("%02d", 1:20))

counts <- matrix(
  sample(10:500, 100 * 20, replace = TRUE),
  nrow = 100, ncol = 20,
  dimnames = list(genes, samples)
)

write.csv(counts, "example_data/simulated_counts.csv", row.names = TRUE)
message("✅ simulated_counts.csv created successfully.")

############################################################
# 2️⃣ Sample metadata (Tumor / Normal, Gender, Age)
############################################################
condition <- rep(c("Tumor", "Normal"), each = 10)
gender <- sample(c("Male", "Female"), 20, replace = TRUE)
age <- round(rnorm(20, mean = 60, sd = 8))

metadata <- data.frame(
  SampleID = samples,
  Condition = condition,
  Gender = gender,
  Age = age
)

write.csv(metadata, "example_data/sample_metadata.csv", row.names = FALSE)
message("✅ sample_metadata.csv created successfully.")

############################################################
# 3️⃣ Gene annotation (Symbol, Biotype)
############################################################
gene_symbol <- paste0("MIR", sample(1000:9999, 100))
biotype <- sample(c("protein_coding", "lncRNA", "miRNA"), 100, replace = TRUE)

annotation <- data.frame(
  GeneID = genes,
  Symbol = gene_symbol,
  Biotype = biotype
)

write.csv(annotation, "example_data/gene_annotation.csv", row.names = FALSE)
message("✅ gene_annotation.csv created successfully.")

############################################################
# 🧩 Summary message
############################################################
message("\n🎉 Example data generation completed!")
message("Files saved in: example_data/")
message("- simulated_counts.csv")
message("- sample_metadata.csv")
message("- gene_annotation.csv")
############################################################
