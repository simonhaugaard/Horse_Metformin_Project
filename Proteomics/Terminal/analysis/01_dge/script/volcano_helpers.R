# Simplified function to create custom volcano plots with additional features
create_custom_volcano_plot <- function(df, logFC_col = "logFC", pvalue_col = "P.Value", adj_pvalue_col = "adj.P.Val", 
                                       contrast_name = "", fc_cutoff = 0, pvalue_cutoff = 0.05, 
                                       upregulated_color = "#CD2626",  # firebrick3 fill for upregulated genes
                                       downregulated_color = "#000080",  # navy fill for downregulated genes
                                       nonsig_color = "grey70",  # Light grey for non-significant genes
                                       upregulated_rim = "#8B1A1A",  # Darker firebrick for rim
                                       downregulated_rim = "#00004D",  # Darker navy for rim
                                       nonsig_rim = "grey50",  # Darker grey for non-significant genes
                                       plot_width = 6, plot_height = 6, plot_resolution = 600, 
                                       save_plot = FALSE, output_path = "", show_labels = TRUE) {
  
  # Create a new column for significant status based on the specified cutoff values
  df$Significance <- ifelse(df[[adj_pvalue_col]] < pvalue_cutoff & df[[logFC_col]] > fc_cutoff, "Upregulated",
                            ifelse(df[[adj_pvalue_col]] < pvalue_cutoff & df[[logFC_col]] < -fc_cutoff, "Downregulated", 
                                   "Not significant"))
  
  # Create the base volcano plot with custom colors and styles
  base_plot <- ggplot(df, aes(x = !!sym(logFC_col), y = -log10(!!sym(adj_pvalue_col)))) + 
    geom_point(aes(fill = Significance), color = nonsig_rim, size = 3, alpha = 0.8, shape = 21, stroke = 1.5) +
    scale_fill_manual(values = c("Upregulated" = upregulated_color, 
                                 "Downregulated" = downregulated_color, 
                                 "Not significant" = nonsig_color)) + 
    geom_point(data = df %>% filter(Significance == "Upregulated"), 
               aes(x = !!sym(logFC_col), y = -log10(!!sym(adj_pvalue_col))), 
               color = upregulated_rim, size = 3, shape = 21, stroke = 1.5, alpha = 0.8) +
    geom_point(data = df %>% filter(Significance == "Downregulated"), 
               aes(x = !!sym(logFC_col), y = -log10(!!sym(adj_pvalue_col))), 
               color = downregulated_rim, size = 3, shape = 21, stroke = 1.5, alpha = 0.8) +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "grey") + 
    geom_hline(yintercept = -log10(pvalue_cutoff), linetype = "dashed", color = "grey") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "black") +  # Unadjusted p-value line at 0.05
    labs(title = contrast_name, x = "Log2 Fold Change", y = "-Log10 FDR") + 
    theme_minimal(base_size = 14) + 
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      legend.position = "none", 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 14, face = "bold")
    )
  
  # Create the version of the plot without labels
  volcano_plot_no_labels <- base_plot
  
  # Create the version of the plot with gene labels, if required
  if (show_labels) {
    volcano_plot_with_labels <- base_plot +
      geom_text_repel(data = df %>% filter(Significance != "Not significant"), 
                      aes(x = !!sym(logFC_col), y = -log10(!!sym(adj_pvalue_col)), label = Genes),
                      size = 3, max.overlaps = 20, segment.color = "grey50", segment.size = 0.3, 
                      max.iter = 2000,  # More iterations for better placement
                      box.padding = 0.2, point.padding = 0.2)  # Reduced padding for shorter connection lines
  } else {
    volcano_plot_with_labels <- volcano_plot_no_labels
  }
  
  # Save the plots if requested
  if (save_plot) {
    ggsave(filename = paste0(output_path, "volcano_plot_", contrast_name, "_No_Labels.png"), 
           plot = volcano_plot_no_labels, 
           width = plot_width, 
           height = plot_height, 
           dpi = plot_resolution)
    ggsave(filename = paste0(output_path, "volcano_plot_", contrast_name, "_With_Labels.png"), 
           plot = volcano_plot_with_labels, 
           width = plot_width, 
           height = plot_height, 
           dpi = plot_resolution)
  }
  
  # Return a list with both plots for further use
  return(list(No_Labels = volcano_plot_no_labels, With_Labels = volcano_plot_with_labels))
}
