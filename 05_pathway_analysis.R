# Comprehensive Pathway Analysis for Heart Disease Transcriptomics
# This script performs hierarchical pathway analysis: All DEGs -> Kinases -> CAMK family

library(clusterProfiler)  # Package for pathway analysis
library(org.Hs.eg.db)     # Human gene database
library(ggplot2)          # Package for plotting
library(ReactomePA)       # Reactome pathway analysis
library(enrichplot)       # Enhanced plotting for enrichment
library(DOSE)             # Disease ontology analysis
library(dplyr)            # Data manipulation

# Set up folders
analysis_results_folder <- "/Users/macbookpro/Desktop/Ananylum/new_meta_analysis/results"
plots_folder <- file.path(analysis_results_folder, "condition_specific_analysis")
dir.create(plots_folder, showWarnings = FALSE)

print("Loading differential gene expression results...")
load(file.path(analysis_results_folder, "03_differential_gene_analysis_results.RData"))

print("Setting up comprehensive pathway analysis...")
print("Analysis strategy: Unfiltered DEGs -> Kinase genes -> CAMK family")

# Define gene sets for hierarchical analysis
important_calcium_genes <- c("CAMK1", "CAMK1D", "CAMK1G", "CAMK2A", "CAMK2B", "CAMK2D", 
                            "CAMK2G", "CAMK4", "CAMKK1", "CAMKK2", "PNCK")

# Expanded kinase gene list (major kinase families)
kinase_genes <- c(
  # CAMK family (our main interest)
  "CAMK1", "CAMK1D", "CAMK1G", "CAMK2A", "CAMK2B", "CAMK2D", "CAMK2G", "CAMK4", "CAMKK1", "CAMKK2", "PNCK",
  # Protein kinase A (PKA)
  "PRKACA", "PRKACB", "PRKAR1A", "PRKAR1B", "PRKAR2A", "PRKAR2B",
  # Protein kinase C (PKC)  
  "PRKCA", "PRKCB", "PRKCD", "PRKCE", "PRKCG", "PRKCH", "PRKCI", "PRKCQ", "PRKCZ",
  # MAP kinases
  "MAPK1", "MAPK3", "MAPK8", "MAPK9", "MAPK10", "MAPK11", "MAPK12", "MAPK13", "MAPK14",
  # Cyclin-dependent kinases
  "CDK1", "CDK2", "CDK4", "CDK6", "CDK7", "CDK8", "CDK9",
  # Cardiac-relevant kinases
  "GSK3A", "GSK3B", "AKT1", "AKT2", "AKT3", "ROCK1", "ROCK2", "PIK3CA", "PIK3CB", "PIK3CD"
)

print("Preparing gene lists for hierarchical analysis...")

# Function to prepare gene lists for pathway analysis
prepare_gene_lists_for_pathway_analysis <- function(deg_results, condition_name, significance_levels) {
  
  print(paste("Processing", condition_name, "differential expression results..."))
  
  # Convert gene symbols to Entrez IDs
  all_genes_with_ids <- bitr(rownames(deg_results), 
                            fromType = "SYMBOL", 
                            toType = "ENTREZID", 
                            OrgDb = org.Hs.eg.db)
  
  # Merge with DEG results
  deg_with_ids <- merge(deg_results, all_genes_with_ids, by.x = "row.names", by.y = "SYMBOL")
  names(deg_with_ids)[1] <- "SYMBOL"
  
  gene_lists <- list()
  
  # Create ranked gene list for GSEA (all genes)
  gsea_genes <- deg_with_ids$logFC
  names(gsea_genes) <- deg_with_ids$ENTREZID
  gsea_genes <- sort(gsea_genes, decreasing = TRUE)
  gene_lists$gsea_ranked <- gsea_genes
  
  # Create different significance-based gene sets
  for(level in names(significance_levels)) {
    criteria <- significance_levels[[level]]
    
    significant_genes <- deg_with_ids[
      deg_with_ids$adj.P.Val <= criteria$p_cutoff & 
      abs(deg_with_ids$logFC) >= criteria$fold_cutoff, 
    ]
    
    if(nrow(significant_genes) > 0) {
      gene_lists[[paste0("significant_", level)]] <- significant_genes$ENTREZID
    }
  }
  
  # Filter for kinase genes
  kinase_in_data <- intersect(kinase_genes, deg_with_ids$SYMBOL)
  if(length(kinase_in_data) > 0) {
    kinase_deg_data <- deg_with_ids[deg_with_ids$SYMBOL %in% kinase_in_data, ]
    gene_lists$all_kinases <- kinase_deg_data$ENTREZID
    
    # Significant kinases
    sig_kinases <- kinase_deg_data[
      kinase_deg_data$adj.P.Val <= 0.20 & abs(kinase_deg_data$logFC) >= 0.1, 
    ]
    if(nrow(sig_kinases) > 0) {
      gene_lists$significant_kinases <- sig_kinases$ENTREZID
    }
  }
  
  # Filter for CAMK genes
  camk_in_data <- intersect(important_calcium_genes, deg_with_ids$SYMBOL)
  if(length(camk_in_data) > 0) {
    camk_deg_data <- deg_with_ids[deg_with_ids$SYMBOL %in% camk_in_data, ]
    gene_lists$camk_genes <- camk_deg_data$ENTREZID
  }
  
  return(list(
    gene_lists = gene_lists,
    deg_with_ids = deg_with_ids,
    kinase_genes_found = kinase_in_data,
    camk_genes_found = camk_in_data
  ))
}

# Define significance levels for analysis
significance_levels <- list(
  "lenient" = list(p_cutoff = 0.20, fold_cutoff = 0.1),
  "default" = list(p_cutoff = 0.05, fold_cutoff = 0.5),
  "no_filter" = list(p_cutoff = 1.0, fold_cutoff = 0.0)
)

# Function to perform comprehensive pathway enrichment analysis
perform_pathway_enrichment <- function(gene_list, analysis_name, condition_name) {
  
  if(length(gene_list) == 0) {
    print(paste("No genes provided for", analysis_name, "- skipping"))
    return(NULL)
  }
  
  print(paste("Running pathway enrichment for", analysis_name, "in", condition_name, "( n =", length(gene_list), "genes )"))
  
  results <- list()
  
  # GO Biological Process
  try({
    go_bp <- enrichGO(gene = gene_list,
                      OrgDb = org.Hs.eg.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.20,
                      readable = TRUE)
    if(!is.null(go_bp) && nrow(go_bp@result) > 0) {
      results$GO_BP <- go_bp
    }
  }, silent = TRUE)
  
  # GO Molecular Function
  try({
    go_mf <- enrichGO(gene = gene_list,
                      OrgDb = org.Hs.eg.db,
                      ont = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.20,
                      readable = TRUE)
    if(!is.null(go_mf) && nrow(go_mf@result) > 0) {
      results$GO_MF <- go_mf
    }
  }, silent = TRUE)
  
  # KEGG Pathways
  try({
    kegg <- enrichKEGG(gene = gene_list,
                       organism = 'hsa',
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.20)
    if(!is.null(kegg) && nrow(kegg@result) > 0) {
      results$KEGG <- kegg
    }
  }, silent = TRUE)
  
  # Reactome Pathways
  try({
    reactome <- enrichPathway(gene = gene_list,
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.20,
                             readable = TRUE)
    if(!is.null(reactome) && nrow(reactome@result) > 0) {
      results$Reactome <- reactome
    }
  }, silent = TRUE)
  
  return(results)
}

# Function to perform GSEA
perform_gsea_analysis <- function(ranked_genes, condition_name) {
  
  print(paste("Running GSEA for", condition_name, "with", length(ranked_genes), "ranked genes"))
  
  results <- list()
  
  # GSEA GO Biological Process
  try({
    gsea_go <- gseGO(geneList = ranked_genes,
                     OrgDb = org.Hs.eg.db,
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     verbose = FALSE)
    if(!is.null(gsea_go) && nrow(gsea_go@result) > 0) {
      results$GSEA_GO <- gsea_go
    }
  }, silent = TRUE)
  
  # GSEA KEGG
  try({
    gsea_kegg <- gseKEGG(geneList = ranked_genes,
                         organism = 'hsa',
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         verbose = FALSE)
    if(!is.null(gsea_kegg) && nrow(gsea_kegg@result) > 0) {
      results$GSEA_KEGG <- gsea_kegg
    }
  }, silent = TRUE)
  
  return(results)
}

# ======= MAIN ANALYSIS EXECUTION =======

print("=== PROCESSING ATRIAL FIBRILLATION (AF vs Control) ===")
af_processed <- prepare_gene_lists_for_pathway_analysis(atrial_fib_all_gene_results, "AF", significance_levels)

print("=== PROCESSING HEART FAILURE (HF vs Control) ===") 
hf_processed <- prepare_gene_lists_for_pathway_analysis(heart_failure_all_gene_results, "HF", significance_levels)

# Store all results
all_pathway_results <- list()
all_gsea_results <- list()

print("=== RUNNING COMPREHENSIVE PATHWAY ANALYSIS ===")

# Analyze AF condition
print("ANALYZING AF PATHWAYS...")
for(gene_set_name in names(af_processed$gene_lists)) {
  
  if(gene_set_name == "gsea_ranked") {
    # GSEA analysis
    gsea_result <- perform_gsea_analysis(af_processed$gene_lists[[gene_set_name]], "AF")
    if(length(gsea_result) > 0) {
      all_gsea_results[paste0("AF_", gene_set_name)] <- list(gsea_result)
    }
  } else {
    # Regular enrichment analysis
    enrichment_result <- perform_pathway_enrichment(
      af_processed$gene_lists[[gene_set_name]], 
      gene_set_name, 
      "AF"
    )
    if(!is.null(enrichment_result) && length(enrichment_result) > 0) {
      all_pathway_results[paste0("AF_", gene_set_name)] <- list(enrichment_result)
    }
  }
}

# Analyze HF condition  
print("ANALYZING HF PATHWAYS...")
for(gene_set_name in names(hf_processed$gene_lists)) {
  
  if(gene_set_name == "gsea_ranked") {
    # GSEA analysis
    gsea_result <- perform_gsea_analysis(hf_processed$gene_lists[[gene_set_name]], "HF")
    if(length(gsea_result) > 0) {
      all_gsea_results[paste0("HF_", gene_set_name)] <- list(gsea_result)
    }
  } else {
    # Regular enrichment analysis
    enrichment_result <- perform_pathway_enrichment(
      hf_processed$gene_lists[[gene_set_name]], 
      gene_set_name, 
      "HF"
    )
    if(!is.null(enrichment_result) && length(enrichment_result) > 0) {
      all_pathway_results[paste0("HF_", gene_set_name)] <- list(enrichment_result)
    }
  }
}

print("=== GENERATING COMPREHENSIVE SUMMARY ===")

# Function to summarize pathway results
summarize_pathway_results <- function(results_list, analysis_type) {
  
  summary_data <- data.frame()
  
  for(analysis_name in names(results_list)) {
    result_set <- results_list[[analysis_name]]
    
    for(db_name in names(result_set)) {
      db_result <- result_set[[db_name]]
      
      if(!is.null(db_result) && nrow(db_result@result) > 0) {
        significant_count <- sum(db_result@result$p.adjust < 0.05)
        
        summary_row <- data.frame(
          analysis_name = analysis_name,
          database = db_name,
          analysis_type = analysis_type,
          total_terms = nrow(db_result@result),
          significant_terms = significant_count,
          top_term = if(nrow(db_result@result) > 0) db_result@result$Description[1] else "None",
          top_pvalue = if(nrow(db_result@result) > 0) db_result@result$pvalue[1] else NA,
          top_qvalue = if(nrow(db_result@result) > 0) db_result@result$p.adjust[1] else NA,
          stringsAsFactors = FALSE
        )
        
        summary_data <- rbind(summary_data, summary_row)
      }
    }
  }
  
  return(summary_data)
}

# Create comprehensive summaries
enrichment_summary <- summarize_pathway_results(all_pathway_results, "Enrichment")
gsea_summary <- summarize_pathway_results(all_gsea_results, "GSEA")
complete_analysis_summary <- rbind(enrichment_summary, gsea_summary)

# Create gene-level summaries
print("Creating comprehensive gene summaries...")

# Kinase genes summary
af_kinase_summary <- data.frame()
hf_kinase_summary <- data.frame()

if(length(af_processed$kinase_genes_found) > 0) {
  for(gene in af_processed$kinase_genes_found) {
    if(gene %in% rownames(atrial_fib_all_gene_results)) {
      gene_data <- atrial_fib_all_gene_results[gene, ]
      summary_row <- data.frame(
        gene = gene,
        condition = "AF",
        gene_family = ifelse(gene %in% important_calcium_genes, "CAMK", "Other_Kinase"),
        fold_change = gene_data$logFC,
        p_value = gene_data$P.Value,
        adj_p_value = gene_data$adj.P.Val,
        avg_expression = gene_data$AveExpr,
        is_significant_default = (gene_data$adj.P.Val < 0.05) & (abs(gene_data$logFC) > 0.5),
        is_significant_lenient = (gene_data$adj.P.Val < 0.20) & (abs(gene_data$logFC) > 0.1),
        stringsAsFactors = FALSE
      )
      af_kinase_summary <- rbind(af_kinase_summary, summary_row)
    }
  }
}

if(length(hf_processed$kinase_genes_found) > 0) {
  for(gene in hf_processed$kinase_genes_found) {
    if(gene %in% rownames(heart_failure_all_gene_results)) {
      gene_data <- heart_failure_all_gene_results[gene, ]
      summary_row <- data.frame(
        gene = gene,
        condition = "HF",
        gene_family = ifelse(gene %in% important_calcium_genes, "CAMK", "Other_Kinase"),
        fold_change = gene_data$logFC,
        p_value = gene_data$P.Value,
        adj_p_value = gene_data$adj.P.Val,
        avg_expression = gene_data$AveExpr,
        is_significant_default = (gene_data$adj.P.Val < 0.05) & (abs(gene_data$logFC) > 0.5),
        is_significant_lenient = (gene_data$adj.P.Val < 0.20) & (abs(gene_data$logFC) > 0.1),
        stringsAsFactors = FALSE
      )
      hf_kinase_summary <- rbind(hf_kinase_summary, summary_row)
    }
  }
}

all_kinase_summary <- rbind(af_kinase_summary, hf_kinase_summary)

# Analysis overview summary
analysis_overview <- data.frame(
  condition = c("AF", "HF"),
  total_genes_analyzed = c(nrow(atrial_fib_all_gene_results), nrow(heart_failure_all_gene_results)),
  kinase_genes_found = c(length(af_processed$kinase_genes_found), length(hf_processed$kinase_genes_found)),
  camk_genes_found = c(length(af_processed$camk_genes_found), length(hf_processed$camk_genes_found)),
  significant_genes_default = c(
    nrow(atrial_fib_all_gene_results[atrial_fib_all_gene_results$adj.P.Val < 0.05 & abs(atrial_fib_all_gene_results$logFC) > 0.5, ]),
    nrow(heart_failure_all_gene_results[heart_failure_all_gene_results$adj.P.Val < 0.05 & abs(heart_failure_all_gene_results$logFC) > 0.5, ])
  ),
  significant_genes_lenient = c(
    nrow(atrial_fib_all_gene_results[atrial_fib_all_gene_results$adj.P.Val < 0.20 & abs(atrial_fib_all_gene_results$logFC) > 0.1, ]),
    nrow(heart_failure_all_gene_results[heart_failure_all_gene_results$adj.P.Val < 0.20 & abs(heart_failure_all_gene_results$logFC) > 0.1, ])
  ),
  stringsAsFactors = FALSE
)

print("=== SAVING COMPREHENSIVE RESULTS ===")

# Save pathway analysis results to CSV files
if(nrow(complete_analysis_summary) > 0) {
  write.csv(complete_analysis_summary,
            file.path(analysis_results_folder, "05_comprehensive_pathway_analysis_summary.csv"),
            row.names = FALSE)
}

# Save gene summaries
write.csv(analysis_overview,
          file.path(analysis_results_folder, "05_analysis_overview.csv"),
          row.names = FALSE)

if(nrow(all_kinase_summary) > 0) {
  write.csv(all_kinase_summary,
            file.path(analysis_results_folder, "05_kinase_genes_expression_summary.csv"),
            row.names = FALSE)
}

# Save detailed pathway results
save_detailed_pathway_results <- function(results_list, prefix) {
  for(analysis_name in names(results_list)) {
    result_set <- results_list[[analysis_name]]
    
    for(db_name in names(result_set)) {
      db_result <- result_set[[db_name]]
      
      if(!is.null(db_result) && nrow(db_result@result) > 0) {
        filename <- paste0("05_", prefix, "_", analysis_name, "_", db_name, "_detailed.csv")
        write.csv(db_result@result, 
                  file.path(analysis_results_folder, filename),
                  row.names = FALSE)
      }
    }
  }
}

# Save all detailed results
save_detailed_pathway_results(all_pathway_results, "enrichment")
save_detailed_pathway_results(all_gsea_results, "gsea")

# Save complete analysis objects
save(all_pathway_results, all_gsea_results, complete_analysis_summary, 
     analysis_overview, all_kinase_summary, af_processed, hf_processed,
     file = file.path(analysis_results_folder, "05_comprehensive_pathway_analysis_results.RData"))

print("=== COMPREHENSIVE PATHWAY ANALYSIS COMPLETE! ===")
print("")

# Display comprehensive summary
print("ANALYSIS OVERVIEW:")
for(i in 1:nrow(analysis_overview)) {
  condition_info <- analysis_overview[i, ]
  print(paste("", condition_info$condition, "condition:"))
  print(paste("  Total genes analyzed:", condition_info$total_genes_analyzed))
  print(paste("  Kinase genes found:", condition_info$kinase_genes_found))
  print(paste("  CAMK genes found:", condition_info$camk_genes_found))
  print(paste("  Significant genes (default):", condition_info$significant_genes_default))
  print(paste("  Significant genes (lenient):", condition_info$significant_genes_lenient))
  print("")
}

print("PATHWAY ENRICHMENT RESULTS:")
if(nrow(complete_analysis_summary) > 0) {
  # Group by condition and show summary
  af_results <- complete_analysis_summary[grep("AF_", complete_analysis_summary$analysis_name), ]
  hf_results <- complete_analysis_summary[grep("HF_", complete_analysis_summary$analysis_name), ]
  
  if(nrow(af_results) > 0) {
    print("  AF (Atrial Fibrillation) Results:")
    for(i in 1:nrow(af_results)) {
      result_info <- af_results[i, ]
      print(paste("   ", result_info$analysis_name, "-", result_info$database, ":", 
                  result_info$significant_terms, "significant /", result_info$total_terms, "total terms"))
    }
    print("")
  }
  
  if(nrow(hf_results) > 0) {
    print("  HF (Heart Failure) Results:")
    for(i in 1:nrow(hf_results)) {
      result_info <- hf_results[i, ]
      print(paste("   ", result_info$analysis_name, "-", result_info$database, ":", 
                  result_info$significant_terms, "significant /", result_info$total_terms, "total terms"))
    }
    print("")
  }
} else {
  print("  No significant pathway enrichments found")
  print("  This may be expected for focused gene sets")
  print("")
}

print("KINASE GENES EXPRESSION:")
if(nrow(all_kinase_summary) > 0) {
  # Show significant kinase genes
  significant_kinases <- all_kinase_summary[all_kinase_summary$is_significant_lenient, ]
  if(nrow(significant_kinases) > 0) {
    print("  Significantly changed kinase genes (lenient threshold):")
    for(i in 1:nrow(significant_kinases)) {
      gene_info <- significant_kinases[i, ]
      direction <- if(gene_info$fold_change > 0) "upregulated" else "downregulated"
      print(paste("   ", gene_info$gene, "(", gene_info$gene_family, ") in", gene_info$condition, ":", 
                  direction, "by", round(abs(gene_info$fold_change), 2), "fold"))
    }
  } else {
    print("  No kinase genes significantly changed at lenient threshold")
  }
  print("")
  
  # Show CAMK-specific results
  camk_specific <- all_kinase_summary[all_kinase_summary$gene_family == "CAMK", ]
  if(nrow(camk_specific) > 0) {
    print("  CAMK family gene expression:")
    for(i in 1:nrow(camk_specific)) {
      gene_info <- camk_specific[i, ]
      significance <- ifelse(gene_info$is_significant_lenient, "significant", "not significant")
      print(paste("   ", gene_info$gene, "in", gene_info$condition, ":", 
                  round(gene_info$fold_change, 3), "fold change,", significance))
    }
  }
  print("")
}

print("FILES CREATED:")
print("- 05_comprehensive_pathway_analysis_summary.csv: Overall pathway analysis results")
print("- 05_analysis_overview.csv: High-level analysis statistics") 
print("- 05_kinase_genes_expression_summary.csv: Detailed kinase gene expression")
print("- Multiple detailed CSV files for each enrichment analysis")
print("- 05_comprehensive_pathway_analysis_results.RData: Complete analysis objects")
print("")

print("SUCCESS: Comprehensive hierarchical pathway analysis complete!")
print("Analysis pipeline: Unfiltered DEGs -> Kinase genes -> CAMK family")
print("Methods used: GO, KEGG, GSEA, Reactome for both AF and HF conditions")
print("")
print("=== ANALYSIS PIPELINE COMPLETE! ===")