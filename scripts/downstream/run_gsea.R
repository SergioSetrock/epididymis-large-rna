# --- run_gsea.R ---
# GSEA Analysis for RNA-seq using clusterProfiler

suppressMessages({
  library(clusterProfiler)
  library(msigdbr)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
  library(ggridges)
})

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
out_dir <- args[2]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- 1. Seleccionar la Base de Datos (MSigDB - Ratón) ---
message("Cargando la base de datos MSigDB Hallmark para Mus musculus...")
m_df <- msigdbr(species = "Mus musculus", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)

# Alternativa: Si quieres Gene Ontology (Procesos Biológicos) en el futuro, comenta la línea de arriba y descomenta esta:
# m_df <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP") %>% dplyr::select(gs_name, gene_symbol)

# --- 2. Función Principal de GSEA ---
run_gsea_for_comparison <- function(csv_file, comparison_name) {
  message(paste("\nCorriendo GSEA para:", comparison_name))
  
  # Cargar datos
  df <- read.csv(csv_file)
  
  # Limpiar datos: Quitar genes sin símbolo y ordenar por Log2FoldChange (apeglm shrunk)
  df <- df %>%
    filter(!is.na(Symbol) & Symbol != "") %>%
    # Si hay símbolos duplicados (isoformas), nos quedamos con el de mayor cambio absoluto
    arrange(desc(abs(log2FoldChange))) %>%
    distinct(Symbol, .keep_all = TRUE)
  
  # Crear la lista ranqueada (Requisito estricto de GSEA: vector nombrado, ordenado de mayor a menor)
  gene_list <- df$log2FoldChange
  names(gene_list) <- df$Symbol
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # Ejecutar GSEA
  # pvalueCutoff = 0.1 para coincidir con tu rigor de DESeq2 (FDR < 0.1)
  set.seed(42) # Para reproducibilidad
  gsea_res <- GSEA(geneList = gene_list, 
                   TERM2GENE = m_df, 
                   pvalueCutoff = 0.1, 
                   minGSSize = 10, 
                   maxGSSize = 500,
                   eps = 0) # eps=0 permite calcular p-valores más exactos
  
  # Verificar si hubo resultados significativos
  if (nrow(as.data.frame(gsea_res)) == 0) {
    message("  -> No se encontraron vías significativas (FDR < 0.1).")
    writeLines("No significant pathways found.", file.path(out_dir, paste0("01_GSEA_Table_", comparison_name, "_NO_SIG.txt")))
    return(NULL)
  }
  
  # --- 3. Generar Outputs ---
  
  # Output 1: Tabla de Vías Significativas (incluye 'core_enrichment' que son tus genes en la vía)
  res_table <- as.data.frame(gsea_res)
  # Limpiar el nombre de la vía para la tabla (quitar "HALLMARK_")
  res_table$ID <- gsub("HALLMARK_", "", res_table$ID)
  write.csv(res_table, file.path(out_dir, paste0("01_GSEA_Table_", comparison_name, ".csv")), row.names = FALSE)
  
  # Output 2: Ridgeplot (Top 10 vías)
  # Muestra la distribución del fold change de los genes dentro de las vías significativas
  p_ridge <- ridgeplot(gsea_res, showCategory = 10) + 
    theme_bw() +
    labs(title = paste("GSEA Ridgeplot:", comparison_name)) +
    theme(axis.text.y = element_text(size = 10, face = "bold"))
  
  ggsave(file.path(out_dir, paste0("02_GSEA_Ridgeplot_", comparison_name, ".png")), plot = p_ridge, width = 10, height = 8)
  
  # Output 3: Enrichment Plot Clásico (Top 3 vías combinadas)
  # Genera la clásica curva de GSEA con el "código de barras" de los genes
  num_pathways_to_plot <- min(3, nrow(res_table))
  p_enrich <- gseaplot2(gsea_res, geneSetID = 1:num_pathways_to_plot, 
                        title = paste("Top", num_pathways_to_plot, "Enriched Pathways"),
                        pvalue_table = TRUE)
  
  ggsave(file.path(out_dir, paste0("03_GSEA_EnrichmentPlot_", comparison_name, ".png")), plot = p_enrich, width = 10, height = 8)
  
  message("  -> ¡Gráficos y tablas generados con éxito!")
}

# --- 4. Iterar sobre los 3 contrastes ---
# Buscamos los archivos "00_ALL_Genes_...csv" que contienen el transcriptoma completo ranqueado
comparisons <- c("iAs_vs_Control_in_F1", "iAs_vs_Control_in_F2", "iAs_vs_Control_Overall_Across_Generations")

for (comp in comparisons) {
  file_path <- file.path(input_dir, paste0("00_ALL_Genes_", comp, "_annotated.csv"))
  
  if (file.exists(file_path)) {
    run_gsea_for_comparison(file_path, comp)
  } else {
    # Si por alguna razón el anotado no existe, intenta con el original
    file_path_unannotated <- file.path(input_dir, paste0("00_ALL_Genes_", comp, ".csv"))
    if(file.exists(file_path_unannotated)){
      run_gsea_for_comparison(file_path_unannotated, comp)
    } else {
      warning(paste("No se encontró el archivo completo para:", comp))
    }
  }
}

message("\n¡Análisis GSEA finalizado!")
