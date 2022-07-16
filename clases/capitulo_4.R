# visualizacion

# heatmap
# este heatmap solo se realiza con los genes significativos

sig_norm_counts_wt <- normalized_wt_counts[wt_res_sig$ensgene,]

heat_colors <- brewer.pal(6, "YlOrRd")

pheatmap(sig_norm_counts_wt,
         color = heat_colors,
         cluster_rows = T,
         show_rownames = F,
         annotation = select(wt_metadata, condition),
         scale = "row")


# volcano plot

# priemero se genera una columna de valores logicos con los genes
# que tienen un valor p adjustado

wt_res_all <- wt_res_all %>%
  mutate(threshold = padj < 0.05)

# a continuacion con ggplot2 grafica los valores log2 fold change vs
# los valores log10 de los p values ajustados

ggplot(wt_res_all, aes(x = log2FoldChange, y = -log10(padj),
                       color = threshold)) +
  geom_point()+
  xlab("log2 fold change") +
  ylab(" -log10 adjusted p value") +
  theme(legend.position = "none")
