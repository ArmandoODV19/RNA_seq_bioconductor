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

# podemos hacer un zoom para observar los genes significaivos

ggplot(wt_res_all, aes(x = log2FoldChange, y = -log10(padj),
                       color = threshold)) +
  geom_point()+
  xlab("log2 fold change") +
  ylab(" -log10 adjusted p value") +
  ylim(c(0,15))+
  theme(legend.position = "none")

# otra accion es visualizar la expresion de los 20 genes mas
# significativos, para esto se utilizaran las cuentas normalizadas
# para los genes significativos ordenados por valor p ajustado
# seleccionaremos los 20 genes mas significativos

top_20 <- data.frame(sig_norm_counts_wt)[1:20,] %>%
  rownames_to_column(var = "ensgene")

top_20 <- gather(top_20, key = "samplename", value = "normalized_counts", 2:8)

# para graficar la expresion se necesita unir el metadata para que el color
# del grafico sea acorde a cada grupo

top_20 <- inner_join(top_20,
                     rownames_to_column(wt_metadata, var = "samplename"),
                     by = "samplename")

# garficar

ggplot(top_20, aes(x = ensgene, y = normalized_counts, color = condition))+
  geom_point()+
  scale_y_log10()+
  xlab("Genes")+
  ylab("Normalized counts")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


