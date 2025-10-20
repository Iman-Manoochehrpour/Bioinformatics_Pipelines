############################################################
# 🎨 generate_visualization_examples.R
# Author: Iman Manoochehrpour
# Purpose: Generate example visualization outputs
############################################################

if (!dir.exists("visualization")) {
  dir.create("visualization")
  message("📁 Folder 'visualization/' created.")
}

library(ggplot2)
library(pheatmap)

# ---- Load example data ----
counts <- read.csv("example_data/simulated_counts.csv", row.names = 1)
metadata <- read.csv("example_data/sample_metadata.csv")

# ---- Create a dummy DESeq2-like result ----
set.seed(123)
res <- data.frame(
  log2FoldChange = rnorm(nrow(counts), mean = 0, sd = 2),
  padj = runif(nrow(counts), min = 0, max = 1)
)
res$significant <- res$padj < 0.05 & abs(res$log2FoldChange) > 1

# ---- Volcano Plot ----
volcano <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.7) +
  theme_classic() +
  labs(title = "Example Volcano Plot",
       x = "log2(Fold Change)",
       y = "-log10(Adjusted p-value)") +
  scale_color_manual(values = c("grey70", "firebrick"))

ggsave("visualization/example_volcano_plot.png", plot = volcano, width = 7, height = 6, dpi = 300)
message("✅ example_volcano_plot.png created successfully.")

# ---- WGCNA-like Module-Trait Heatmap ----
# Create random module-trait correlation matrix
modules <- paste0("Module_", LETTERS[1:5])
traits <- c("Tumor", "Normal", "Age", "Gender")
cor_mat <- matrix(runif(length(modules)*length(traits), -1, 1),
                  nrow = length(modules),
                  dimnames = list(modules, traits))

# Draw heatmap
pheatmap(cor_mat,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         cluster_rows = FALSE, cluster_cols = FALSE,
         main = "Example Module–Trait Correlation Heatmap",
         filename = "visualization/example_wgcna_heatmap.png",
         width = 6, height = 5)
message("✅ example_wgcna_heatmap.png created successfully.")

message("\n🎉 Example visualizations generated successfully!")
