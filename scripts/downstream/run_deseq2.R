# --- run_deseq2.R ---
args <- commandArgs(trailingOnly = TRUE)
counts_file <- args[1]
meta_file <- args[2]
outdir <- args[3]

# Load required libraries
suppressMessages({
  library(DESeq2)
  library(apeglm)
  library(ashr)
  library(ggplot2)
  library(pheatmap)
  library(dplyr)
  library(tibble)
})

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# --- 1. Load Data ---
counts <- read.table(counts_file, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
meta <- read.csv(meta_file, row.names=1)

# Ensure matching order
common_samples <- intersect(colnames(counts), rownames(meta))
counts <- counts[, common_samples]
meta <- meta[common_samples, ]

# Set reference levels (Control and F1 are the baselines)
meta$Treatment <- factor(meta$Treatment, levels = c("Control", "iAs"))
meta$Generation <- factor(meta$Generation, levels = c("F1", "F2"))

# --- 2. Run DESeq2 (Interaction Model) ---
dds <- DESeqDataSetFromMatrix(countData = round(counts), 
                              colData = meta, 
                              design = ~ Treatment + Generation + Treatment:Generation)
dds <- DESeq(dds)

# Print the names for the log file so you can verify what DESeq2 named them
message("DESeq2 calculated the following coefficients:")
print(resultsNames(dds))

# Use absolute positional indices to avoid any string matching errors!
# 1 = Intercept
# 2 = Main Effect of Treatment (iAs vs Control in F1)
# 3 = Main Effect of Generation (F2 vs F1 in Control)
# 4 = Interaction Effect
coef_treatment_F1 <- 2
coef_generation_Ctrl <- 3
coef_interaction <- 4

# --- 3. QC STEP: Shrinkage Method Comparison (MA Plots) ---
message("Running shrinkage comparison (apeglm vs ashr)...")

res_apeglm_QC <- lfcShrink(dds, coef = coef_treatment_F1, type = "apeglm")
res_ashr_QC <- lfcShrink(dds, coef = coef_treatment_F1, type = "ashr")

png(file.path(outdir, "01_Shrinkage_Comparison_MA_Plot.png"), width=1200, height=600, res=150)
par(mfrow=c(1,2))
plotMA(res_apeglm_QC, main="apeglm Shrinkage (Treatment in F1)", ylim=c(-3,3))
plotMA(res_ashr_QC, main="ashr Shrinkage (Treatment in F1)", ylim=c(-3,3))
dev.off()

# --- 4. Final Shrinkage Application ---
# Change this variable to "ashr" if the MA plots look better for your data
chosen_shrinkage <- "apeglm" 
message(paste("Proceeding with downstream analysis using:", chosen_shrinkage))

res_treat_F1 <- lfcShrink(dds, coef = coef_treatment_F1, type = chosen_shrinkage)
res_gen_Ctrl <- lfcShrink(dds, coef = coef_generation_Ctrl, type = chosen_shrinkage)
res_interaction <- lfcShrink(dds, coef = coef_interaction, type = chosen_shrinkage)

# --- 5. Filter for Significance and Save Tables ---
process_results <- function(res_shrunk, factor_name) {
  res_df <- as.data.frame(res_shrunk) %>%
    rownames_to_column("GeneID") %>%
    mutate(
      Regulation = case_when(
        padj < 0.1 & log2FoldChange > 0.58 ~ "Upregulated",
        padj < 0.1 & log2FoldChange < -0.58 ~ "Downregulated",
        TRUE ~ "Not_Significant"
      )
    )
  
  # Filter significant only
  sig_genes <- res_df %>% filter(Regulation != "Not_Significant")
  
  # Save to CSV
  write.csv(sig_genes, file.path(outdir, paste0("02_Significant_DEGs_", factor_name, ".csv")), row.names=FALSE)
  return(res_df)
}

res_df_treat_F1 <- process_results(res_treat_F1, "Treatment_in_F1")
res_df_gen_Ctrl <- process_results(res_gen_Ctrl, "Generation_in_Control")
res_df_interaction <- process_results(res_interaction, "Interaction_Treatment_x_Generation")

# --- 6. Summary Plot (Up/Down Counts) ---
summary_data <- data.frame(
  Factor = c(rep("Treatment in F1", 2), rep("Generation in Ctrl", 2), rep("Interaction Effect", 2)),
  Regulation = c("Upregulated", "Downregulated", "Upregulated", "Downregulated", "Upregulated", "Downregulated"),
  Count = c(
    sum(res_df_treat_F1$Regulation == "Upregulated", na.rm=TRUE),
    sum(res_df_treat_F1$Regulation == "Downregulated", na.rm=TRUE),
    sum(res_df_gen_Ctrl$Regulation == "Upregulated", na.rm=TRUE),
    sum(res_df_gen_Ctrl$Regulation == "Downregulated", na.rm=TRUE),
    sum(res_df_interaction$Regulation == "Upregulated", na.rm=TRUE),
    sum(res_df_interaction$Regulation == "Downregulated", na.rm=TRUE)
  )
)

summary_data$Factor <- factor(summary_data$Factor, levels=c("Treatment in F1", "Generation in Ctrl", "Interaction Effect"))

p_summary <- ggplot(summary_data, aes(x=Factor, y=Count, fill=Regulation)) +
  geom_bar(stat="identity", position="dodge", color="black") +
  scale_fill_manual(values=c("Downregulated"="blue", "Upregulated"="red")) +
  geom_text(aes(label=Count), position=position_dodge(width=0.9), vjust=-0.5, fontface="bold") +
  theme_bw() +
  labs(title="Count of Significant DE Genes (Interaction Model)", y="Number of Genes", x="") +
  theme(text = element_text(size=14, face="bold"))

ggsave(file.path(outdir, "03_DE_Gene_Counts_Summary.png"), plot=p_summary, width=10, height=6)

# --- 7. Heatmap of ALL Significant DE Genes across the whole model ---
sig_genes_all <- unique(c(
  res_df_treat_F1 %>% filter(Regulation != "Not_Significant") %>% pull(GeneID),
  res_df_gen_Ctrl %>% filter(Regulation != "Not_Significant") %>% pull(GeneID),
  res_df_interaction %>% filter(Regulation != "Not_Significant") %>% pull(GeneID)
))

if (length(sig_genes_all) > 1) {
  vsd <- vst(dds, blind=FALSE)
  mat <- assay(vsd)[sig_genes_all, ]
  
  mat_scaled <- t(scale(t(mat)))
  
  ann_col <- data.frame(Treatment = meta$Treatment, Generation = meta$Generation, row.names = rownames(meta))
  ann_colors <- list(
    Treatment = c(Control = "green", iAs = "salmon"),
    Generation = c(F1 = "grey80", F2 = "grey30")
  )
  heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(50)
  
  png(file.path(outdir, "04_Significant_DEGs_Heatmap_Overall.png"), width=1000, height=800, res=150)
  pheatmap(mat_scaled, 
           annotation_col = ann_col, 
           annotation_colors = ann_colors,
           color = heatmap_colors,
           show_rownames = FALSE, 
           cluster_cols = TRUE,
           main = paste("Significant DE Genes (Model Overview) -", length(sig_genes_all), "genes"))
  dev.off()
} else {
  writeLines("Not enough significant genes to draw a heatmap.", file.path(outdir, "04_NO_SIG_GENES_HEATMAP.txt"))
}

message("DESeq2 interaction analysis completed successfully!")
