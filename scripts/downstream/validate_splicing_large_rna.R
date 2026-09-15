# --- validate_splicing_large_rna.R ---
suppressMessages({
  library(DESeq2)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(dplyr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
counts_file <- args[1]
meta_file   <- args[2]

message("Loading Large RNA Counts and Metadata...")
counts <- read.table(counts_file, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
meta <- read.csv(meta_file, row.names=1)

common_samples <- intersect(colnames(counts), rownames(meta))
counts <- counts[, common_samples]
meta <- meta[common_samples, ]

meta$Group <- factor(paste(meta$Treatment, meta$Generation, sep="_"))
meta$Treatment <- factor(meta$Treatment, levels = c("Control", "iAs"))
meta$Generation <- factor(meta$Generation, levels = c("F1", "F2"))

message("Running DESeq2 for F2 (iAs vs Control)...")
dds <- DESeqDataSetFromMatrix(countData = round(counts), colData = meta, design = ~ Group)
dds$Group <- relevel(dds$Group, ref = "Control_F2")
dds <- DESeq(dds, quiet = TRUE)
res_F2 <- results(dds, contrast=c("Group", "iAs_F2", "Control_F2"))

# Extract significant Differentially Expressed Genes (FDR < 0.05)
sig_genes <- rownames(res_F2)[which(res_F2$padj < 0.05)]

if (length(sig_genes) == 0) {
  stop("No significant DEGs found in F2 Large RNA to run enrichment.")
}

message(paste("Found", length(sig_genes), "significant Large RNA DEGs in F2. Running GO Enrichment..."))

# Clean Ensembl IDs (remove version numbers)
clean_ids <- gsub("\\..*$", "", sig_genes)

# Run GO Enrichment
ego <- enrichGO(gene          = clean_ids,
                OrgDb         = org.Mm.eg.db,
                keyType       = "ENSEMBL",
                ont           = "BP", # Biological Process
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.1,
                readable      = TRUE)

if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
  message("No significant GO pathways were found.")
  q()
}

# --- SEARCH ENGINE: Look for our keywords ---
df_ego <- as.data.frame(ego)

# Use regex to search for splicing, spliceosome, RNA processing, or ribosome
validation_hits <- df_ego %>%
  filter(str_detect(Description, "(?i)splic|ribosom|rna processing|translation")) %>%
  select(ID, Description, p.adjust, Count) %>%
  arrange(p.adjust)

message("\n=========================================================")
message("             MULTI-OMICS VALIDATION RESULTS              ")
message("=========================================================\n")

if (nrow(validation_hits) > 0) {
  message("BINGO! We found significant Splicing/Ribosomal pathways in the Large RNA data!")
  message("This perfectly validates the snRNA/snoRNA findings from the small RNA data.\n")
  print(validation_hits)
} else {
  message("We did not find splicing or ribosomal pathways directly enriched in the F2 DEGs.")
  message("This could mean the changes are happening at the post-transcriptional level,")
  message("or we might need to run a full ranked GSEA on GO terms instead of just DEGs.")
}
message("\n=========================================================")
