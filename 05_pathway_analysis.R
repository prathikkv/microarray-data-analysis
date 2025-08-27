# Analyze biological pathways for CAMK genes - BEGINNER VERSION
# This script finds what biological processes CAMK genes are involved in

library(clusterProfiler)  # Package for pathway analysis
library(org.Hs.eg.db)     # Human gene database
library(ggplot2)          # Package for plotting

# Set up folders
analysis_results_folder <- "/Users/macbookpro/Desktop/Ananylum/new_meta_analysis/results"
plots_folder <- file.path(analysis_results_folder, "condition_specific_analysis")

print("Loading gene expression analysis results...")
load(file.path(analysis_results_folder, "03_differential_gene_analysis_results.RData"))

# CAMK genes we want to study
important_calcium_genes <- c("CAMK1", "CAMK1D", "CAMK1G", "CAMK2A", "CAMK2B", "CAMK2D", 
                            "CAMK2G", "CAMK4", "CAMKK1", "CAMKK2", "PNCK")

print("Checking which CAMK genes are present in our heart disease data...")

# Find CAMK genes that are actually measured in our datasets
camk_genes_in_af_data <- intersect(important_calcium_genes, rownames(atrial_fib_all_gene_results))
camk_genes_in_hf_data <- intersect(important_calcium_genes, rownames(heart_failure_all_gene_results))

print(paste("CAMK genes found in atrial fibrillation data:", length(camk_genes_in_af_data)))
if(length(camk_genes_in_af_data) > 0) {
  print(paste("  Genes:", paste(camk_genes_in_af_data, collapse = ", ")))
}

print(paste("CAMK genes found in heart failure data:", length(camk_genes_in_hf_data)))
if(length(camk_genes_in_hf_data) > 0) {
  print(paste("  Genes:", paste(camk_genes_in_hf_data, collapse = ", ")))
}

# Combine all CAMK genes found in either dataset
all_camk_genes_found <- unique(c(camk_genes_in_af_data, camk_genes_in_hf_data))
print(paste("Total unique CAMK genes in our data:", length(all_camk_genes_found)))

# Convert gene symbols to Entrez IDs (required for pathway analysis)
print("Converting gene names to database IDs for pathway analysis...")

camk_gene_id_conversion <- bitr(all_camk_genes_found, 
                               fromType = "SYMBOL", 
                               toType = "ENTREZID", 
                               OrgDb = org.Hs.eg.db)

print("CAMK genes with database IDs:")
print(camk_gene_id_conversion)

if(nrow(camk_gene_id_conversion) == 0) {
  print("ERROR: No CAMK genes could be converted to database IDs")
  stop("Cannot proceed with pathway analysis")
}

print("Analyzing biological pathways that involve CAMK genes...")

# Analyze KEGG pathways (well-known biological pathways database)
camk_pathway_analysis <- enrichKEGG(gene = camk_gene_id_conversion$ENTREZID,
                                   organism = 'hsa',  # Human
                                   pAdjustMethod = "BH",  # Multiple testing correction
                                   pvalueCutoff = 1.0,    # Show all results initially
                                   qvalueCutoff = 1.0)

# Check if we found any pathways
if(is.null(camk_pathway_analysis) || nrow(camk_pathway_analysis@result) == 0) {
  print("No significant pathways found for CAMK genes")
  print("This is expected since we only have a small number of CAMK genes")
  
  # Create a summary anyway
  camk_pathway_summary <- data.frame(
    analysis_type = "KEGG_pathways",
    total_camk_genes_analyzed = length(all_camk_genes_found),
    camk_genes_with_database_ids = nrow(camk_gene_id_conversion),
    pathways_found = 0,
    significant_pathways = 0,
    top_pathway = "None found",
    top_pathway_pvalue = NA,
    notes = "Small gene set - limited pathway enrichment expected"
  )
  
} else {
  
  # We found some pathways!
  pathway_results <- camk_pathway_analysis@result
  
  print(paste("Total pathways involving CAMK genes:", nrow(pathway_results)))
  print("")
  
  # Show top pathways
  significant_pathways <- pathway_results[pathway_results$p.adjust < 0.05, ]
  
  if(nrow(significant_pathways) > 0) {
    print(paste("Statistically significant pathways (p < 0.05):", nrow(significant_pathways)))
    print("")
    print("Top significant pathways:")
    top_significant <- head(significant_pathways, 10)
    for(i in 1:nrow(top_significant)) {
      pathway_info <- top_significant[i, ]
      print(paste("", i, ".", pathway_info$Description))
      print(paste("     P-value:", format(pathway_info$pvalue, scientific = TRUE)))
      print(paste("     Adjusted p-value:", format(pathway_info$p.adjust, scientific = TRUE)))
      print(paste("     CAMK genes involved:", pathway_info$Count, "out of", pathway_info$BgRatio))
      print("")
    }
  } else {
    print("No statistically significant pathways found (p < 0.05)")
    print("Showing top pathways regardless of significance:")
    top_pathways <- head(pathway_results, 10)
    for(i in 1:nrow(top_pathways)) {
      pathway_info <- top_pathways[i, ]
      print(paste("", i, ".", pathway_info$Description))
      print(paste("     P-value:", format(pathway_info$pvalue, scientific = TRUE)))
      print(paste("     Adjusted p-value:", format(pathway_info$p.adjust, scientific = TRUE)))
      print(paste("     CAMK genes involved:", pathway_info$Count))
      print("")
    }
  }
  
  # Create summary
  camk_pathway_summary <- data.frame(
    analysis_type = "KEGG_pathways",
    total_camk_genes_analyzed = length(all_camk_genes_found),
    camk_genes_with_database_ids = nrow(camk_gene_id_conversion),
    pathways_found = nrow(pathway_results),
    significant_pathways = nrow(significant_pathways),
    top_pathway = pathway_results$Description[1],
    top_pathway_pvalue = pathway_results$pvalue[1],
    notes = if(nrow(significant_pathways) > 0) "Significant pathways found" else "No significant pathways"
  )
  
  # Save detailed results
  write.csv(pathway_results, 
            file.path(analysis_results_folder, "05_CAMK_KEGG_pathways_detailed.csv"), 
            row.names = FALSE)
  
}

print("Analyzing Gene Ontology (GO) biological processes...")

# Analyze GO Biological Process terms
camk_go_analysis <- enrichGO(gene = camk_gene_id_conversion$ENTREZID,
                            OrgDb = org.Hs.eg.db,
                            ont = "BP",  # Biological Process
                            pAdjustMethod = "BH",
                            pvalueCutoff = 1.0,
                            qvalueCutoff = 1.0,
                            readable = TRUE)

# Check GO results
if(is.null(camk_go_analysis) || nrow(camk_go_analysis@result) == 0) {
  print("No GO biological processes found for CAMK genes")
  
  go_summary <- data.frame(
    analysis_type = "GO_biological_process",
    total_camk_genes_analyzed = length(all_camk_genes_found),
    camk_genes_with_database_ids = nrow(camk_gene_id_conversion),
    processes_found = 0,
    significant_processes = 0,
    top_process = "None found",
    top_process_pvalue = NA,
    notes = "Small gene set - limited GO enrichment expected"
  )
  
} else {
  
  go_results <- camk_go_analysis@result
  
  print(paste("Total GO biological processes involving CAMK genes:", nrow(go_results)))
  
  # Show top GO processes
  significant_go <- go_results[go_results$p.adjust < 0.05, ]
  
  if(nrow(significant_go) > 0) {
    print(paste("Significant GO processes (p < 0.05):", nrow(significant_go)))
    print("")
    print("Top significant GO biological processes:")
    top_go <- head(significant_go, 5)
    for(i in 1:nrow(top_go)) {
      go_info <- top_go[i, ]
      print(paste("", i, ".", go_info$Description))
      print(paste("     P-value:", format(go_info$pvalue, scientific = TRUE)))
      print(paste("     Genes:", go_info$geneID))
      print("")
    }
  } else {
    print("No statistically significant GO processes found")
    print("Top GO processes (regardless of significance):")
    top_go <- head(go_results, 5)
    for(i in 1:nrow(top_go)) {
      go_info <- top_go[i, ]
      print(paste("", i, ".", go_info$Description))
      print(paste("     P-value:", format(go_info$pvalue, scientific = TRUE)))
      print(paste("     Genes:", go_info$geneID))
      print("")
    }
  }
  
  go_summary <- data.frame(
    analysis_type = "GO_biological_process", 
    total_camk_genes_analyzed = length(all_camk_genes_found),
    camk_genes_with_database_ids = nrow(camk_gene_id_conversion),
    processes_found = nrow(go_results),
    significant_processes = nrow(significant_go),
    top_process = go_results$Description[1],
    top_process_pvalue = go_results$pvalue[1],
    notes = if(nrow(significant_go) > 0) "Significant processes found" else "No significant processes"
  )
  
  # Save GO results
  write.csv(go_results,
            file.path(analysis_results_folder, "05_CAMK_GO_biological_processes.csv"),
            row.names = FALSE)
  
}

print("Creating CAMK gene expression summary...")

# Summarize CAMK gene expression patterns
camk_expression_summary <- data.frame()

for(camk_gene in all_camk_genes_found) {
  
  # Get AF data if available
  if(camk_gene %in% rownames(atrial_fib_all_gene_results)) {
    af_data <- atrial_fib_all_gene_results[camk_gene, ]
    af_row <- data.frame(
      gene_name = camk_gene,
      condition = "AF",
      fold_change = af_data$logFC,
      p_value = af_data$P.Value,
      adjusted_p_value = af_data$adj.P.Val,
      average_expression = af_data$AveExpr,
      is_significant_default = (af_data$adj.P.Val < 0.05) & (abs(af_data$logFC) > 0.5),
      stringsAsFactors = FALSE
    )
    camk_expression_summary <- rbind(camk_expression_summary, af_row)
  }
  
  # Get HF data if available
  if(camk_gene %in% rownames(heart_failure_all_gene_results)) {
    hf_data <- heart_failure_all_gene_results[camk_gene, ]
    hf_row <- data.frame(
      gene_name = camk_gene,
      condition = "HF",
      fold_change = hf_data$logFC,
      p_value = hf_data$P.Value,
      adjusted_p_value = hf_data$adj.P.Val,
      average_expression = hf_data$AveExpr,
      is_significant_default = (hf_data$adj.P.Val < 0.05) & (abs(hf_data$logFC) > 0.5),
      stringsAsFactors = FALSE
    )
    camk_expression_summary <- rbind(camk_expression_summary, hf_row)
  }
}

print("Saving pathway analysis results...")

# Combine analysis summaries
complete_pathway_summary <- rbind(camk_pathway_summary, go_summary)

# Save all results
write.csv(complete_pathway_summary,
          file.path(analysis_results_folder, "05_pathway_analysis_summary.csv"),
          row.names = FALSE)

write.csv(camk_expression_summary,
          file.path(analysis_results_folder, "05_CAMK_expression_patterns_summary.csv"),
          row.names = FALSE)

write.csv(camk_gene_id_conversion,
          file.path(analysis_results_folder, "05_CAMK_gene_ID_conversions.csv"),
          row.names = FALSE)

# Save complete analysis
save(camk_pathway_analysis, camk_go_analysis, camk_expression_summary, complete_pathway_summary,
     file = file.path(analysis_results_folder, "05_complete_pathway_analysis_results.RData"))

print("Pathway analysis complete!")
print("")
print("=== FINAL CAMK PATHWAY ANALYSIS SUMMARY ===")
print("")
print("CAMK Genes Analyzed:")
print(paste("  Total CAMK genes in data:", length(all_camk_genes_found)))
print(paste("  Genes:", paste(all_camk_genes_found, collapse = ", ")))
print("")

print("KEGG Pathway Analysis:")
kegg_summary <- complete_pathway_summary[complete_pathway_summary$analysis_type == "KEGG_pathways", ]
print(paste("  Pathways found:", kegg_summary$pathways_found))
print(paste("  Significant pathways:", kegg_summary$significant_pathways))
if(!is.na(kegg_summary$top_pathway_pvalue)) {
  print(paste("  Top pathway:", kegg_summary$top_pathway))
  print(paste("  Top pathway p-value:", format(kegg_summary$top_pathway_pvalue, scientific = TRUE)))
}
print("")

print("GO Biological Process Analysis:")
go_summary_final <- complete_pathway_summary[complete_pathway_summary$analysis_type == "GO_biological_process", ]
print(paste("  Processes found:", go_summary_final$processes_found))
print(paste("  Significant processes:", go_summary_final$significant_processes))
if(!is.na(go_summary_final$top_process_pvalue)) {
  print(paste("  Top process:", go_summary_final$top_process))
  print(paste("  Top process p-value:", format(go_summary_final$top_process_pvalue, scientific = TRUE)))
}
print("")

print("CAMK Gene Expression in Heart Disease:")
significant_camk <- camk_expression_summary[camk_expression_summary$is_significant_default, ]
if(nrow(significant_camk) > 0) {
  print("  Significantly changed CAMK genes:")
  for(i in 1:nrow(significant_camk)) {
    gene_info <- significant_camk[i, ]
    direction <- if(gene_info$fold_change > 0) "increased" else "decreased"
    print(paste("   ", gene_info$gene_name, "in", gene_info$condition, ":", direction, 
                "by", round(abs(gene_info$fold_change), 2), "fold (p =", 
                format(gene_info$adjusted_p_value, scientific = TRUE), ")"))
  }
} else {
  print("  No CAMK genes significantly changed at default thresholds")
}
print("")

print("Files created:")
print("- Detailed KEGG pathway results (if any found)")
print("- GO biological process results (if any found)")
print("- CAMK gene expression pattern summary")
print("- Complete pathway analysis summary")
print("- All results saved for further investigation")
print("")
print("ANALYSIS PIPELINE COMPLETE!")