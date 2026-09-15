# --- plot_gsea_dotplots.R ---
suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]

# Find all significant GSEA tables in the target folder
csv_files <- list.files(input_dir, pattern = "01_GSEA_Table_.*\\.csv$", full.names = TRUE)

if (length(csv_files) == 0) {
  message(paste("No significant GSEA tables found in", input_dir))
  q()
}

for (file in csv_files) {
  # Extract the comparison name from the filename
  comparison_name <- gsub("01_GSEA_Table_", "", basename(file))
  comparison_name <- gsub("\\.csv$", "", comparison_name)
  
  df <- read.csv(file)
  
  if (nrow(df) == 0) next
  
  # Calculate Gene Ratio and format variables
  # GSEA outputs driving genes separated by "/" in the 'core_enrichment' column
  df <- df %>%
    mutate(
      Count = str_count(core_enrichment, "/") + 1,
      GeneRatio = Count / setSize,
      Log10P = -log10(p.adjust)
    ) %>%
    # Sort by NES so activated pathways group at the top and suppressed at the bottom
    arrange(desc(NES)) 
  
  # Lock the factor levels to maintain this order in the plot
  df$ID <- factor(df$ID, levels = df$ID)
  
  # Base plotting function
  create_dotplot <- function(plot_data, title) {
    ggplot(plot_data, aes(x = Log10P, y = ID)) +
      geom_point(aes(size = GeneRatio, color = NES)) +
      # Blue for suppressed (negative), Red for activated (positive)
      scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "NES\n(Activation)") +
      theme_bw() +
      labs(title = title,
           x = "-Log10(FDR)",
           y = "",
           size = "Gene Ratio") +
      theme(axis.text.y = element_text(size = 11, face = "bold"),
            axis.text.x = element_text(size = 12),
            plot.title = element_text(size = 14, face = "bold"))
  }
  
  # --- Plot 1: ALL Significant Pathways ---
  # Adjust height dynamically based on how many pathways there are
  plot_height_all <- max(6, nrow(df) * 0.3) 
  p_all <- create_dotplot(df, paste("GSEA Dot Plot (All):", gsub("_", " ", comparison_name)))
  ggsave(file.path(input_dir, paste0("04_GSEA_DotPlot_All_", comparison_name, ".png")), plot = p_all, width = 10, height = plot_height_all)
  
  # --- Plot 2: TOP 10 Significant Pathways ---
  if (nrow(df) > 10) {
    # Isolate the top 10 by highest statistical significance (lowest FDR)
    df_top10 <- df %>% 
      arrange(p.adjust) %>% 
      head(10) %>%
      arrange(desc(NES)) # Re-sort by NES to keep the red/blue separation visually clean
    
    df_top10$ID <- factor(df_top10$ID, levels = df_top10$ID)
    
    p_top10 <- create_dotplot(df_top10, paste("GSEA Dot Plot (Top 10):", gsub("_", " ", comparison_name)))
    ggsave(file.path(input_dir, paste0("05_GSEA_DotPlot_Top10_", comparison_name, ".png")), plot = p_top10, width = 10, height = 6)
  } else {
    # If there are 10 or fewer pathways, just save the "All" plot as the "Top 10" plot too
    ggsave(file.path(input_dir, paste0("05_GSEA_DotPlot_Top10_", comparison_name, ".png")), plot = p_all, width = 10, height = plot_height_all)
  }
  
  message(paste("  -> Generated Dot Plots for:", comparison_name))
}

