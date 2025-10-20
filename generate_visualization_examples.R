# ---- Simulate realistic DESeq2-like results ----
set.seed(123)

# بیشتر log2FoldChange‌ها نزدیک به 0 هستند، با چند تغییر بزرگ
log2fc <- c(rnorm(90, 0, 0.7), rnorm(5, -3, 0.5), rnorm(5, 3, 0.5))

# مقادیر p کوچک‌تر برای ژن‌هایی با FoldChange زیاد
pvals <- runif(100, min = 0.001, max = 0.5)
padj <- p.adjust(pvals, method = "BH")

res <- data.frame(
  log2FoldChange = log2fc,
  padj = padj
)
res$significant <- res$padj < 0.05 & abs(res$log2FoldChange) > 1

# ---- Volcano Plot (improved look) ----
library(ggplot2)
volcano <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = c("grey70", "firebrick")) +
  theme_classic(base_size = 13) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", size = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "black", size = 0.4) +
  labs(
    title = "Example Volcano Plot (Simulated Data)",
    x = "log2(Fold Change)",
    y = "-log10(Adjusted p-value)"
  )

ggsave("visualization/example_volcano_plot.png", plot = volcano, width = 7, height = 6, dpi = 300)
message("✅ Improved volcano plot created successfully.")
