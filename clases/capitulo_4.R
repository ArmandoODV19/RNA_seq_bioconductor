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
