# --- run_gsea_from_counts.R ---
suppressMessages({
  library(DESeq2)
  library(clusterProfiler)
  library(msigdbr)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
  library(ggridges)
  library(org.Mm.eg.db)
})

args <- commandArgs(trailingOnly = TRUE)
counts_file <- args[1]
meta_file <- args[2]
out_dir <- args[3]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Prepare MSigDB Hallmark Gene Sets for Mouse
message("Loading MSigDB Hallmark pathways for Mouse (Using ortholog mapping)...")

# We use collection = "H" and let it use the default ortholog mapping
m_df <- msigdbr(species = "Mus musculus", collection = "H") %>% 
  dplyr::select(gs_name, gene_symbol)

# --- QC CHECK ---
message(paste("QC: Successfully loaded", length(unique(m_df$gs_name)), "Hallmark pathways."))
message("QC: Top 5 genes in the MSigDB database:")
print(head(m_df$gene_symbol, 5))
# ----------------

# 2. Load Data
counts <- read.table(counts_file, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
meta <- read.csv(meta_file, row.names=1)

common_samples <- intersect(colnames(counts), rownames(meta))
counts <- counts[, common_samples]
meta <- meta[common_samples, ]

meta$Group <- factor(paste(meta$Treatment, meta$Generation, sep="_"))
meta$Treatment <- factor(meta$Treatment, levels = c("Control", "iAs"))
meta$Generation <- factor(meta$Generation, levels = c("F1", "F2"))

# 3. Run DESeq2 to get robust rankings 
message("\nRunning DESeq2 to generate ranked gene lists...")

dds_group <- DESeqDataSetFromMatrix(countData = round(counts), colData = meta, design = ~ Group)
dds_group$Group <- relevel(dds_group$Group, ref = "Control_F1")
dds_group <- DESeq(dds_group)
res_F1 <- results(dds_group, contrast=c("Group", "iAs_F1", "Control_F1"))

dds_group$Group <- relevel(dds_group$Group, ref = "Control_F2")
dds_group <- nbinomWaldTest(dds_group)
res_F2 <- results(dds_group, contrast=c("Group", "iAs_F2", "Control_F2"))

dds_main <- DESeqDataSetFromMatrix(countData = round(counts), colData = meta, design = ~ Generation + Treatment)
dds_main <- DESeq(dds_main)
res_overall <- results(dds_main, contrast=c("Treatment", "iAs", "Control"))

# 4. GSEA Function
run_gsea <- function(res_obj, comparison_name) {
  message(paste("\n========================================="))
  message(paste("Performing GSEA for:", comparison_name))
  
  df <- as.data.frame(res_obj)
  df$GeneID <- rownames(df)
  df <- df %>% filter(!is.na(stat))
  
  # Map Ensembl IDs to Gene Symbols
  clean_ids <- gsub("\\..*$", "", df$GeneID)
  symbols <- mapIds(org.Mm.eg.db, keys = clean_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
  df$Symbol <- symbols
  
  # Remove missing symbols and handle duplicates
  df <- df %>% 
    filter(!is.na(Symbol) & Symbol != "") %>% 
    arrange(desc(abs(stat))) %>% 
    distinct(Symbol, .keep_all = TRUE)
  
  # Create the ranked named vector
  gene_list <- df$stat
  names(gene_list) <- df$Symbol
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # --- QC CHECK ---
  message("QC: Top 5 genes in our Ranked List (Should look similar to the MSigDB genes above):")
  print(head(names(gene_list), 5))
  # ----------------
  
  # Run GSEA
  set.seed(42)
  gsea_res <- GSEA(geneList = gene_list, 
                   TERM2GENE = m_df, 
                   pvalueCutoff = 0.1, 
                   minGSSize = 10, 
                   maxGSSize = 500,
                   eps = 0)
                   
  if (is.null(gsea_res) || nrow(as.data.frame(gsea_res)) == 0) {
    message("  -> No significant pathways found (FDR < 0.1).")
    writeLines("No significant pathways found.", file.path(out_dir, paste0("01_GSEA_Table_", comparison_name, "_NO_SIG.txt")))
    return(NULL)
  }
  
  # Output 1: Table
  res_table <- as.data.frame(gsea_res)
  res_table$ID <- gsub("HALLMARK_", "", res_table$ID)
  write.csv(res_table, file.path(out_dir, paste0("01_GSEA_Table_", comparison_name, ".csv")), row.names = FALSE)
  
  # Output 2: Ridgeplot (Top 10)
  p_ridge <- ridgeplot(gsea_res, showCategory = 10) + 
    theme_bw() +
    labs(title = paste("GSEA Ridgeplot:", gsub("_", " ", comparison_name)))
  ggsave(file.path(out_dir, paste0("02_GSEA_Ridgeplot_", comparison_name, ".png")), plot = p_ridge, width = 10, height = 8)
  
  # Output 3: Enrichment Plot (Top 3)
  num_pathways <- min(3, nrow(res_table))
  p_enrich <- gseaplot2(gsea_res, geneSetID = 1:num_pathways, 
                        title = paste("Top", num_pathways, "Enriched Pathways"),
                        pvalue_table = TRUE)
  ggsave(file.path(out_dir, paste0("03_GSEA_EnrichmentPlot_", comparison_name, ".png")), plot = p_enrich, width = 10, height = 8)
  
  message("  -> Outputs generated successfully!")
}

# Run for all three comparisons
run_gsea(res_F1, "iAs_vs_Control_in_F1")
run_gsea(res_F2, "iAs_vs_Control_in_F2")
run_gsea(res_overall, "iAs_vs_Control_Overall_Across_Generations")

message("\nAll GSEA analyses completed!")
