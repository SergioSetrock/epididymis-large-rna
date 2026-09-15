# --- run_deseq2_final.R ---
args <- commandArgs(trailingOnly = TRUE)
counts_file <- args[1]
meta_file <- args[2]
outdir <- args[3]

suppressMessages({
  library(DESeq2)
  library(apeglm)
  library(ggplot2)
  library(pheatmap)
  library(dplyr)
  library(tibble)
})

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# --- 1. Load Data ---
counts <- read.table(counts_file, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
meta <- read.csv(meta_file, row.names=1)

common_samples <- intersect(colnames(counts), rownames(meta))
counts <- counts[, common_samples]
meta <- meta[common_samples, ]

# Create a combined Group factor for F1 and F2 specific effects
meta$Group <- factor(paste(meta$Treatment, meta$Generation, sep="_"))
meta$Treatment <- factor(meta$Treatment, levels = c("Control", "iAs"))
meta$Generation <- factor(meta$Generation, levels = c("F1", "F2"))

# --- 2. MODEL A: Specific Effects per Generation (~ Group) ---
message("Fitting Model A for Generation-specific effects...")
dds_group <- DESeqDataSetFromMatrix(countData = round(counts), colData = meta, design = ~ Group)

# Contrast 1: iAs vs Control in F1
dds_group$Group <- relevel(dds_group$Group, ref = "Control_F1")
dds_group <- DESeq(dds_group)
res_F1 <- lfcShrink(dds_group, coef = "Group_iAs_F1_vs_Control_F1", type = "apeglm")

# Contrast 2: iAs vs Control in F2
dds_group$Group <- relevel(dds_group$Group, ref = "Control_F2")
dds_group <- nbinomWaldTest(dds_group)
res_F2 <- lfcShrink(dds_group, coef = "Group_iAs_F2_vs_Control_F2", type = "apeglm")

# --- 3. MODEL B: Overall Effect Across Generations (~ Generation + Treatment) ---
message("Fitting Model B for Overall Treatment effect across generations...")
dds_main <- DESeqDataSetFromMatrix(countData = round(counts), colData = meta, design = ~ Generation + Treatment)
dds_main <- DESeq(dds_main)
res_overall <- lfcShrink(dds_main, coef = "Treatment_iAs_vs_Control", type = "apeglm")

# --- 4. Process Results and Filter ---
process_results <- function(res_shrunk, factor_name) {
  res_df <- as.data.frame(res_shrunk) %>%
    rownames_to_column("GeneID") %>%
    mutate(
      Regulation = case_when(
        !is.na(padj) & padj < 0.1 & log2FoldChange > 0.58 ~ "Upregulated",
        !is.na(padj) & padj < 0.1 & log2FoldChange < -0.58 ~ "Downregulated",
        TRUE ~ "Not_Significant"
      )
    )
  
  sig_genes <- res_df %>% filter(Regulation != "Not_Significant")
  write.csv(sig_genes, file.path(outdir, paste0("01_Significant_DEGs_", factor_name, ".csv")), row.names=FALSE)
  return(res_df)
}

res_df_F1 <- process_results(res_F1, "iAs_vs_Control_in_F1")
res_df_F2 <- process_results(res_F2, "iAs_vs_Control_in_F2")
res_df_overall <- process_results(res_overall, "iAs_vs_Control_Overall_Across_Generations")

# --- 5. Summary Plot ---
summary_data <- data.frame(
  Contrast = c(rep("iAs vs Ctrl (F1)", 2), rep("iAs vs Ctrl (F2)", 2), rep("iAs vs Ctrl (Overall)", 2)),
  Regulation = c("Upregulated", "Downregulated", "Upregulated", "Downregulated", "Upregulated", "Downregulated"),
  Count = c(
    sum(res_df_F1$Regulation == "Upregulated", na.rm=TRUE),
    sum(res_df_F1$Regulation == "Downregulated", na.rm=TRUE),
    sum(res_df_F2$Regulation == "Upregulated", na.rm=TRUE),
    sum(res_df_F2$Regulation == "Downregulated", na.rm=TRUE),
    sum(res_df_overall$Regulation == "Upregulated", na.rm=TRUE),
    sum(res_df_overall$Regulation == "Downregulated", na.rm=TRUE)
  )
)

summary_data$Contrast <- factor(summary_data$Contrast, levels=c("iAs vs Ctrl (F1)", "iAs vs Ctrl (F2)", "iAs vs Ctrl (Overall)"))

p_summary <- ggplot(summary_data, aes(x=Contrast, y=Count, fill=Regulation)) +
  geom_bar(stat="identity", position="dodge", color="black") +
  scale_fill_manual(values=c("Downregulated"="blue", "Upregulated"="red")) +
  geom_text(aes(label=Count), position=position_dodge(width=0.9), vjust=-0.5, fontface="bold") +
  theme_bw() +
  labs(title="DE Genes (apeglm: FDR < 0.1, |LFC| > 0.58)", y="Number of Genes", x="") +
  theme(text = element_text(size=14, face="bold"))

ggsave(file.path(outdir, "02_DE_Gene_Counts_Summary.png"), plot=p_summary, width=10, height=6)

# --- 6. Plotting Functions (Volcano & Heatmap) ---

plot_volcano <- function(res_df, title, filename) {
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Regulation)) +
    geom_point(alpha = 0.8, size = 1.5) +
    scale_color_manual(values = c("Downregulated" = "blue", "Upregulated" = "red", "Not_Significant" = "grey80")) +
    geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "black") +
    theme_bw() +
    labs(title = title, x = "Log2 Fold Change", y = "-Log10(FDR)") +
    theme(text = element_text(size = 14, face = "bold"), legend.title=element_blank())
  
  ggsave(filename, plot = p, width = 8, height = 6)
}

plot_heatmap <- function(sig_genes, vsd_obj, meta_subset, title, filename) {
  if (length(sig_genes) > 1) {
    # Subset matrix to significant genes and ONLY the relevant samples
    mat <- assay(vsd_obj)[sig_genes, rownames(meta_subset)]
    mat_scaled <- t(scale(t(mat)))
    
    ann_col <- data.frame(Treatment = meta_subset$Treatment, Generation = meta_subset$Generation, row.names = rownames(meta_subset))
    ann_colors <- list(
      Treatment = c(Control = "green", iAs = "salmon"),
      Generation = c(F1 = "grey80", F2 = "grey30")
    )
    heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(50)
    
    png(filename, width = 800, height = 800, res = 150)
    pheatmap(mat_scaled, 
             annotation_col = ann_col, 
             annotation_colors = ann_colors,
             color = heatmap_colors,
             show_rownames = FALSE, 
             cluster_cols = TRUE,
             main = paste(title, "-", length(sig_genes), "genes"))
    dev.off()
  } else {
    writeLines("Not enough significant genes to draw a heatmap.", gsub(".png", ".txt", filename))
  }
}

# --- 7. Generate Individual Plots ---

# Obtain VST transformed data for heatmaps
vsd_group <- vst(dds_group, blind=FALSE)
vsd_main <- vst(dds_main, blind=FALSE)

# Subset metadata for specific heatmaps
meta_F1 <- meta[meta$Generation == "F1", ]
meta_F2 <- meta[meta$Generation == "F2", ]

# A. iAs vs Control in F1
sig_F1 <- res_df_F1 %>% filter(Regulation != "Not_Significant") %>% pull(GeneID)
plot_volcano(res_df_F1, "Volcano: iAs vs Control (F1)", file.path(outdir, "03a_Volcano_F1.png"))
plot_heatmap(sig_F1, vsd_group, meta_F1, "DEGs: iAs vs Control (F1)", file.path(outdir, "04a_Heatmap_F1.png"))

# B. iAs vs Control in F2
sig_F2 <- res_df_F2 %>% filter(Regulation != "Not_Significant") %>% pull(GeneID)
plot_volcano(res_df_F2, "Volcano: iAs vs Control (F2)", file.path(outdir, "03b_Volcano_F2.png"))
plot_heatmap(sig_F2, vsd_group, meta_F2, "DEGs: iAs vs Control (F2)", file.path(outdir, "04b_Heatmap_F2.png"))

# C. Overall iAs vs Control
sig_overall <- res_df_overall %>% filter(Regulation != "Not_Significant") %>% pull(GeneID)
plot_volcano(res_df_overall, "Volcano: iAs vs Control (Overall)", file.path(outdir, "03c_Volcano_Overall.png"))
plot_heatmap(sig_overall, vsd_main, meta, "DEGs: iAs vs Control (Overall)", file.path(outdir, "04c_Heatmap_Overall.png"))

message("Final DESeq2 analysis (with individual plots) completed successfully!")
