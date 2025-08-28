# Find genes that are different between heart conditions - BEGINNER VERSION
# This script finds genes that change in atrial fibrillation and heart failure

library(limma)     # Package for gene expression analysis
library(ggplot2)   # Package for making plots
library(pheatmap)  # Package for heatmaps
library(RColorBrewer) # Package for color palettes
library(reshape2)  # Package for data reshaping
library(ComplexHeatmap) # Package for complex heatmaps
library(dplyr)     # Package for data manipulation

# Set up folders
analysis_results_folder <- "results"
dir.create(analysis_results_folder, showWarnings = FALSE)

# List of CAMK genes we're specially interested in
important_calcium_genes <- c("CAMK1", "CAMK1D", "CAMK1G", "CAMK2A", "CAMK2B", "CAMK2D", 
                            "CAMK2G", "CAMK4", "CAMKK1", "CAMKK2", "PNCK")

print("Loading batch-corrected heart disease data...")
load(file.path(analysis_results_folder, "02_condition_specific_corrected.RData"))

print("Finding genes that are different in atrial fibrillation...")

# Set up comparison for atrial fibrillation vs control
atrial_fib_comparison_design <- model.matrix(~0 + group, data = af_metadata)
colnames(atrial_fib_comparison_design) <- levels(factor(af_metadata$group))

# Fit statistical model for atrial fibrillation data
atrial_fib_statistical_model <- lmFit(af_combat_expr, atrial_fib_comparison_design)

# Create comparison: AF samples vs Control samples
atrial_fib_contrast <- makeContrasts(AF_vs_Control = AF - Control, levels = atrial_fib_comparison_design)
atrial_fib_model_with_contrast <- contrasts.fit(atrial_fib_statistical_model, atrial_fib_contrast)
atrial_fib_final_model <- eBayes(atrial_fib_model_with_contrast)

# Get results for all genes
atrial_fib_all_gene_results <- topTable(atrial_fib_final_model, coef = "AF_vs_Control", number = Inf, adjust.method = "BH")

print("Finding genes that are different in heart failure...")

# Set up comparison for heart failure vs control
heart_failure_comparison_design <- model.matrix(~0 + group, data = hf_metadata)
colnames(heart_failure_comparison_design) <- levels(factor(hf_metadata$group))

# Fit statistical model for heart failure data
heart_failure_statistical_model <- lmFit(hf_combat_expr, heart_failure_comparison_design)

# Create comparison: HF samples vs Control samples
heart_failure_contrast <- makeContrasts(HF_vs_Control = HF - Control, levels = heart_failure_comparison_design)
heart_failure_model_with_contrast <- contrasts.fit(heart_failure_statistical_model, heart_failure_contrast)
heart_failure_final_model <- eBayes(heart_failure_model_with_contrast)

# Get results for all genes
heart_failure_all_gene_results <- topTable(heart_failure_final_model, coef = "HF_vs_Control", number = Inf, adjust.method = "BH")

print("Testing different significance thresholds...")

# Define 3 levels of filtering criteria for finding significant genes
significance_thresholds <- data.frame(
  threshold_name = c("default", "lenient", "no_filtering"),
  adjusted_p_value_cutoff = c(0.05, 0.20, 1.0),
  fold_change_cutoff = c(0.5, 0.1, 0.0),
  stringsAsFactors = FALSE
)

# Function to find significant genes at different thresholds
find_significant_genes_at_threshold <- function(gene_results, condition_name, p_cutoff, fold_cutoff, threshold_name) {
  
  # Find genes that meet our criteria
  genes_meet_p_value <- gene_results$adj.P.Val <= p_cutoff
  genes_meet_fold_change <- abs(gene_results$logFC) >= fold_cutoff
  significant_genes <- genes_meet_p_value & genes_meet_fold_change
  
  # Get the significant genes
  significant_gene_results <- gene_results[significant_genes, ]
  significant_gene_results <- significant_gene_results[order(significant_gene_results$adj.P.Val), ]
  
  print(paste("Found", nrow(significant_gene_results), "significant genes for", condition_name, "at", threshold_name, "threshold"))
  
  # Check which CAMK genes are significant
  camk_genes_in_results <- important_calcium_genes[important_calcium_genes %in% rownames(significant_gene_results)]
  
  if(length(camk_genes_in_results) > 0) {
    print(paste("  CAMK genes found:", paste(camk_genes_in_results, collapse = ", ")))
    camk_gene_info <- significant_gene_results[camk_genes_in_results, ]
    print("  CAMK gene details:")
    for(i in 1:nrow(camk_gene_info)) {
      print(paste("   ", rownames(camk_gene_info)[i], ": p =", 
                  round(camk_gene_info$adj.P.Val[i], 6), ", fold change =", 
                  round(camk_gene_info$logFC[i], 3)))
    }
  } else {
    print("  No CAMK genes found at this threshold")
  }
  
  # Save results to file
  output_filename <- file.path(analysis_results_folder, 
                              paste0("03_", condition_name, "_significant_genes_", threshold_name, ".csv"))
  write.csv(significant_gene_results, output_filename, row.names = TRUE)
  
  # Return summary information
  return(data.frame(
    condition = condition_name,
    threshold = threshold_name,
    p_cutoff = p_cutoff,
    fold_cutoff = fold_cutoff,
    total_significant_genes = nrow(significant_gene_results),
    camk_genes_found = length(camk_genes_in_results),
    camk_gene_names = paste(camk_genes_in_results, collapse = "; "),
    stringsAsFactors = FALSE
  ))
}

# Test all threshold combinations
all_threshold_results <- data.frame()

for(i in 1:nrow(significance_thresholds)) {
  current_threshold <- significance_thresholds[i, ]
  
  print(paste("Testing", current_threshold$threshold_name, "threshold: p <", current_threshold$adjusted_p_value_cutoff, 
              "and fold change >", current_threshold$fold_change_cutoff))
  
  # Test atrial fibrillation
  af_threshold_result <- find_significant_genes_at_threshold(
    atrial_fib_all_gene_results, "AF", 
    current_threshold$adjusted_p_value_cutoff, 
    current_threshold$fold_change_cutoff,
    current_threshold$threshold_name
  )
  
  # Test heart failure
  hf_threshold_result <- find_significant_genes_at_threshold(
    heart_failure_all_gene_results, "HF", 
    current_threshold$adjusted_p_value_cutoff, 
    current_threshold$fold_change_cutoff,
    current_threshold$threshold_name
  )
  
  all_threshold_results <- rbind(all_threshold_results, af_threshold_result, hf_threshold_result)
  print("")
}

print("Creating summary of CAMK gene significance...")

# Extract CAMK gene information across all thresholds
camk_gene_tracking_af <- data.frame()
camk_gene_tracking_hf <- data.frame()

for(camk_gene in important_calcium_genes) {
  # Check if gene exists in our data
  if(camk_gene %in% rownames(atrial_fib_all_gene_results)) {
    af_gene_data <- atrial_fib_all_gene_results[camk_gene, ]
    camk_gene_tracking_af <- rbind(camk_gene_tracking_af, data.frame(
      gene_name = camk_gene,
      p_value = af_gene_data$P.Value,
      adjusted_p_value = af_gene_data$adj.P.Val,
      fold_change = af_gene_data$logFC,
      average_expression = af_gene_data$AveExpr,
      condition = "AF",
      stringsAsFactors = FALSE
    ))
  }
  
  if(camk_gene %in% rownames(heart_failure_all_gene_results)) {
    hf_gene_data <- heart_failure_all_gene_results[camk_gene, ]
    camk_gene_tracking_hf <- rbind(camk_gene_tracking_hf, data.frame(
      gene_name = camk_gene,
      p_value = hf_gene_data$P.Value,
      adjusted_p_value = hf_gene_data$adj.P.Val,
      fold_change = hf_gene_data$logFC,
      average_expression = hf_gene_data$AveExpr,
      condition = "HF",
      stringsAsFactors = FALSE
    ))
  }
}

# Combine CAMK tracking data
all_camk_tracking <- rbind(camk_gene_tracking_af, camk_gene_tracking_hf)

print("Creating volcano plots...")

# Function to create volcano plot
create_volcano_plot <- function(gene_results, condition_name, p_cutoff = 0.05, fold_cutoff = 0.5) {
  
  # Prepare data for plotting
  plot_data <- data.frame(
    gene_name = rownames(gene_results),
    fold_change = gene_results$logFC,
    p_value = -log10(gene_results$adj.P.Val),
    stringsAsFactors = FALSE
  )
  
  # Color genes based on significance
  plot_data$significance <- "Not Significant"
  plot_data$significance[abs(plot_data$fold_change) >= fold_cutoff & gene_results$adj.P.Val <= p_cutoff] <- "Significant"
  
  # Highlight CAMK genes
  plot_data$gene_type <- "Other"
  plot_data$gene_type[plot_data$gene_name %in% important_calcium_genes] <- "CAMK Gene"
  
  # Create plot
  volcano_plot <- ggplot(plot_data, aes(x = fold_change, y = p_value, color = significance)) +
    geom_point(alpha = 0.6) +
    geom_point(data = subset(plot_data, gene_type == "CAMK Gene"), 
               aes(x = fold_change, y = p_value), 
               color = "red", size = 3) +
    geom_text(data = subset(plot_data, gene_type == "CAMK Gene"),
              aes(x = fold_change, y = p_value, label = gene_name),
              color = "darkred", vjust = -0.5) +
    scale_color_manual(values = c("Not Significant" = "gray", "Significant" = "blue")) +
    geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(-fold_cutoff, fold_cutoff), linetype = "dashed", color = "red") +
    labs(title = paste(condition_name, "Volcano Plot - Significant Gene Changes"),
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-Value") +
    theme_minimal()
  
  return(volcano_plot)
}

# Create volcano plots
af_volcano_plot <- create_volcano_plot(atrial_fib_all_gene_results, "Atrial Fibrillation")
hf_volcano_plot <- create_volcano_plot(heart_failure_all_gene_results, "Heart Failure")

# Create comprehensive differential expression visualization
pdf(file.path(analysis_results_folder, "03_comprehensive_differential_expression_analysis.pdf"), width = 16, height = 12)

print("Creating volcano plots...")
print(af_volcano_plot)
print(hf_volcano_plot)

print("Creating MA plots...")
# 1. MA plots (log fold change vs average expression)
create_ma_plot <- function(gene_results, condition_name, p_cutoff = 0.05, fold_cutoff = 0.5) {
  
  plot_data <- data.frame(
    gene_name = rownames(gene_results),
    avg_expression = gene_results$AveExpr,
    fold_change = gene_results$logFC,
    p_value = gene_results$adj.P.Val,
    stringsAsFactors = FALSE
  )
  
  # Color genes based on significance
  plot_data$significance <- "Not Significant"
  plot_data$significance[abs(plot_data$fold_change) >= fold_cutoff & plot_data$p_value <= p_cutoff] <- "Significant"
  
  # Highlight CAMK genes
  plot_data$gene_type <- "Other"
  plot_data$gene_type[plot_data$gene_name %in% important_calcium_genes] <- "CAMK Gene"
  
  ma_plot <- ggplot(plot_data, aes(x = avg_expression, y = fold_change, color = significance)) +
    geom_point(alpha = 0.6) +
    geom_point(data = subset(plot_data, gene_type == "CAMK Gene"), 
               aes(x = avg_expression, y = fold_change), 
               color = "red", size = 3) +
    geom_text(data = subset(plot_data, gene_type == "CAMK Gene"),
              aes(x = avg_expression, y = fold_change, label = gene_name),
              color = "darkred", vjust = -0.5, size = 3) +
    scale_color_manual(values = c("Not Significant" = "gray", "Significant" = "blue")) +
    geom_hline(yintercept = c(-fold_cutoff, fold_cutoff), linetype = "dashed", color = "red") +
    labs(title = paste(condition_name, "MA Plot - Average Expression vs Fold Change"),
         x = "Average Log2 Expression",
         y = "Log2 Fold Change") +
    theme_minimal()
  
  return(ma_plot)
}

af_ma_plot <- create_ma_plot(atrial_fib_all_gene_results, "Atrial Fibrillation")
hf_ma_plot <- create_ma_plot(heart_failure_all_gene_results, "Heart Failure")

print(af_ma_plot)
print(hf_ma_plot)

print("Creating threshold comparison plots...")
# 2. Threshold comparison bar plot
threshold_counts <- all_threshold_results %>%
  group_by(condition, threshold) %>%
  summarise(total_genes = sum(total_significant_genes), .groups = 'drop')

threshold_plot <- ggplot(threshold_counts, aes(x = threshold, y = total_genes, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("AF" = "#FF6B6B", "HF" = "#4ECDC4")) +
  labs(title = "Significant Genes Across Different Thresholds",
       subtitle = "Comparison of gene counts at default, lenient, and no_filtering thresholds",
       x = "Significance Threshold", y = "Number of Significant Genes",
       fill = "Condition") +
  theme_minimal() +
  geom_text(aes(label = total_genes), vjust = -0.3, position = position_dodge(width = 0.9))

print(threshold_plot)

print("Creating CAMK gene expression plots...")
# 3. CAMK gene expression boxplots
load(file.path(analysis_results_folder, "02_condition_specific_corrected.RData"))

# Get CAMK genes that are present in both datasets
camk_in_af <- intersect(important_calcium_genes, rownames(af_combat_expr))
camk_in_hf <- intersect(important_calcium_genes, rownames(hf_combat_expr))

if(length(camk_in_af) > 0) {
  # AF CAMK expression
  af_camk_expr <- af_combat_expr[camk_in_af, , drop = FALSE]
  af_camk_data <- melt(as.matrix(af_camk_expr))
  names(af_camk_data) <- c("Gene", "Sample", "Expression")
  af_camk_data$Condition <- af_metadata$group[match(af_camk_data$Sample, rownames(af_metadata))]
  
  af_camk_plot <- ggplot(af_camk_data, aes(x = Gene, y = Expression, fill = Condition)) +
    geom_boxplot() +
    scale_fill_manual(values = c("AF" = "#FF6B6B", "Control" = "#95E4D3")) +
    labs(title = "CAMK Gene Expression in Atrial Fibrillation",
         subtitle = "Expression levels of CAMK genes: AF vs Control",
         x = "CAMK Gene", y = "Log2 Expression Level") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(af_camk_plot)
}

if(length(camk_in_hf) > 0) {
  # HF CAMK expression
  hf_camk_expr <- hf_combat_expr[camk_in_hf, , drop = FALSE]
  hf_camk_data <- melt(as.matrix(hf_camk_expr))
  names(hf_camk_data) <- c("Gene", "Sample", "Expression")
  hf_camk_data$Condition <- hf_metadata$group[match(hf_camk_data$Sample, rownames(hf_metadata))]
  
  hf_camk_plot <- ggplot(hf_camk_data, aes(x = Gene, y = Expression, fill = Condition)) +
    geom_boxplot() +
    scale_fill_manual(values = c("HF" = "#4ECDC4", "Control" = "#95E4D3")) +
    labs(title = "CAMK Gene Expression in Heart Failure",
         subtitle = "Expression levels of CAMK genes: HF vs Control",
         x = "CAMK Gene", y = "Log2 Expression Level") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(hf_camk_plot)
}

print("Creating heatmaps of top significant genes...")
# 4. Heatmaps of top significant genes
# AF top genes heatmap
af_top_genes <- head(atrial_fib_all_gene_results[order(atrial_fib_all_gene_results$adj.P.Val), ], 50)
if(nrow(af_top_genes) > 0) {
  af_top_gene_names <- rownames(af_top_genes)
  af_top_gene_names <- af_top_gene_names[af_top_gene_names %in% rownames(af_combat_expr)]
  
  if(length(af_top_gene_names) > 0) {
    af_heatmap_data <- af_combat_expr[af_top_gene_names, ]
    af_annotation <- data.frame(
      Condition = af_metadata$group,
      Dataset = af_metadata$dataset
    )
    rownames(af_annotation) <- rownames(af_metadata)
    
    # Create annotation colors
    ann_colors <- list(
      Condition = c("AF" = "#FF6B6B", "Control" = "#95E4D3"),
      Dataset = c("GSE115574" = "#FF9999", "GSE79768" = "#FFB366")
    )
    
    pheatmap(af_heatmap_data,
             main = "Top 50 Significant Genes in Atrial Fibrillation",
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             annotation_col = af_annotation,
             annotation_colors = ann_colors,
             scale = "row",
             show_colnames = FALSE,
             fontsize = 8)
  }
}

# HF top genes heatmap
hf_top_genes <- head(heart_failure_all_gene_results[order(heart_failure_all_gene_results$adj.P.Val), ], 50)
if(nrow(hf_top_genes) > 0) {
  hf_top_gene_names <- rownames(hf_top_genes)
  hf_top_gene_names <- hf_top_gene_names[hf_top_gene_names %in% rownames(hf_combat_expr)]
  
  if(length(hf_top_gene_names) > 0) {
    hf_heatmap_data <- hf_combat_expr[hf_top_gene_names, ]
    hf_annotation <- data.frame(
      Condition = hf_metadata$group,
      Dataset = hf_metadata$dataset
    )
    rownames(hf_annotation) <- rownames(hf_metadata)
    
    # Create annotation colors
    ann_colors <- list(
      Condition = c("HF" = "#4ECDC4", "Control" = "#95E4D3"),
      Dataset = c("GSE57338" = "#66D9CC", "GSE76701" = "#7FCCB8")
    )
    
    pheatmap(hf_heatmap_data,
             main = "Top 50 Significant Genes in Heart Failure",
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             annotation_col = hf_annotation,
             annotation_colors = ann_colors,
             scale = "row",
             show_colnames = FALSE,
             fontsize = 8)
  }
}

print("Creating p-value distribution plots...")
# 5. P-value distribution plots
af_pval_data <- data.frame(
  p_value = atrial_fib_all_gene_results$P.Value,
  adj_p_value = atrial_fib_all_gene_results$adj.P.Val,
  condition = "AF"
)

hf_pval_data <- data.frame(
  p_value = heart_failure_all_gene_results$P.Value,
  adj_p_value = heart_failure_all_gene_results$adj.P.Val,
  condition = "HF"
)

all_pval_data <- rbind(af_pval_data, hf_pval_data)

pval_plot <- ggplot(all_pval_data, aes(x = p_value, fill = condition)) +
  geom_histogram(alpha = 0.7, bins = 50, position = "identity") +
  scale_fill_manual(values = c("AF" = "#FF6B6B", "HF" = "#4ECDC4")) +
  labs(title = "P-value Distribution Across Conditions",
       subtitle = "Distribution of raw p-values from differential expression analysis",
       x = "P-value", y = "Count", fill = "Condition") +
  theme_minimal()

print(pval_plot)

dev.off()

print("Comprehensive differential expression plots saved!")
print("- Volcano plots for both conditions")
print("- MA plots showing fold change vs average expression")
print("- Threshold comparison showing gene counts")
print("- CAMK gene expression boxplots")
print("- Heatmaps of top 50 significant genes")
print("- P-value distribution analysis")
print("")

# Also save individual volcano plots (legacy)
pdf(file.path(analysis_results_folder, "03_volcano_plots_gene_changes.pdf"), width = 12, height = 8)
print(af_volcano_plot)
print(hf_volcano_plot)
dev.off()

print("Saving all analysis results...")

# Save summary results
write.csv(all_threshold_results, 
          file.path(analysis_results_folder, "03_threshold_analysis_summary.csv"), 
          row.names = FALSE)

write.csv(all_camk_tracking, 
          file.path(analysis_results_folder, "03_CAMK_gene_detailed_results.csv"), 
          row.names = FALSE)

# Save complete gene results
write.csv(atrial_fib_all_gene_results, 
          file.path(analysis_results_folder, "03_AF_all_genes_statistical_results.csv"), 
          row.names = TRUE)

write.csv(heart_failure_all_gene_results, 
          file.path(analysis_results_folder, "03_HF_all_genes_statistical_results.csv"), 
          row.names = TRUE)

# Save processed data for next steps
save(atrial_fib_all_gene_results, heart_failure_all_gene_results, 
     all_threshold_results, all_camk_tracking,
     file = file.path(analysis_results_folder, "03_differential_gene_analysis_results.RData"))

print("Differential gene expression analysis complete!")
print("")
print("SUMMARY OF FINDINGS:")
print("Key Results:")
for(condition in c("AF", "HF")) {
  condition_results <- all_threshold_results[all_threshold_results$condition == condition, ]
  print(paste("", condition, "condition:"))
  for(i in 1:nrow(condition_results)) {
    result <- condition_results[i, ]
    print(paste("  ", result$threshold, "threshold:", result$total_significant_genes, "genes,", 
                result$camk_genes_found, "CAMK genes"))
    if(result$camk_genes_found > 0) {
      print(paste("    CAMK genes:", result$camk_gene_names))
    }
  }
  print("")
}

print("Files created:")
print("- Volcano plots showing gene changes")
print("- CSV files with significant genes at each threshold")
print("- Complete statistical results for all genes")
print("- CAMK gene tracking across all thresholds")
print("- Ready for machine learning analysis")