library(tidyverse)
library(ggrepel)

df <- read_csv("raw data/5hr-2hr.csv")
colnames(df)[which(names(df) == "p-value")] <- "pval"
colnames(df)[which(names(df) == "p-adj")] <- "padj"
df$diffexpressed <- ifelse(df$padj < 0.05 & df$log2FC > 0.6, "2hr",
                           ifelse(df$padj < 0.05 & df$log2FC < -0.6, "5hr", "NO"))

df$diffexpressed <- factor(df$diffexpressed, levels = c("2hr", "5hr", "NO"))

#filter out lowly expressed genes
df_filtered <- subset(df, base_mean > 0)

#volcano plot of DEGs
ggplot(df_filtered, aes(x = log2FC, y = -log10(padj), color = diffexpressed)) +
  geom_point(size = 2) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "gray", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = "dashed") +
  geom_text_repel(data = subset(df_filtered, diffexpressed %in% c("2hr","5hr")),
                  aes(label = gene_id), size = 3) +
  scale_color_manual(
    values = c("2hr" = "darkorange", "5hr" = "forestgreen", "NO" = "grey"),
    labels = c("Up in 2hr", "Up in 5hr", "Not significant"),
    na.translate = FALSE
  ) +
  labs(title = "Volcano Plot: 2hr vs 5hr Embryo",
       x = "log2 Fold Change",
       y = "-log10(p-adjusted value)",
       color = "Expression change") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#dotplot of GO terms
gprofiler <- read_csv("raw data/df2_output.csv")

ggplot(gprofiler, aes(x = negative_log10_of_adjusted_p_value, y = reorder(term_name, negative_log10_of_adjusted_p_value))) +
  geom_point(aes(size = intersection_size), color = "forestgreen") +
  labs(
    title = "Dotplot of Predicted Biological Processes",
    x = "-log10(p-value)",
    y = "GO Biological Process",
    size = "Intersection Size"
  ) +
  theme_basic()
