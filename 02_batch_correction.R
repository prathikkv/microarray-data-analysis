# Remove batch effects between heart studies - BEGINNER VERSION  
# This script fixes technical differences between the different heart disease studies

library(sva)        # Package for batch effect correction
library(ggplot2)    # Package for making plots
library(pheatmap)   # Package for heatmaps
library(Rtsne)      # Package for tSNE plots
library(gridExtra)  # Package to arrange multiple plots
library(cluster)    # Package for clustering analysis

# Set up our folders
analysis_results_folder <- "results"
dir.create(analysis_results_folder, showWarnings = FALSE)

# Load the combined dataset we created in step 1
load(file.path(analysis_results_folder, "01_combined_4datasets.RData"))
print("Loaded combined heart datasets for batch correction...")

# Set up colors for our plots
dataset_colors_for_plots <- c("GSE115574" = "#E41A1C", "GSE79768" = "#377EB8", "GSE76701" = "#4DAF4A", "GSE57338" = "#984EA3")
condition_colors_for_plots <- c("AF" = "#FF6B6B", "HF" = "#4ECDC4", "Control" = "#95A5A6")

print("Processing atrial fibrillation studies...")

# Step 1: Get only the atrial fibrillation samples
atrial_fib_sample_positions <- which(combined_metadata$study_type == "AF")
atrial_fib_expression_data <- combined_expr[, atrial_fib_sample_positions]
atrial_fib_sample_info <- combined_metadata[atrial_fib_sample_positions, ]

# Step 2: Remove genes with too many missing values (more than 50% missing)
missing_values_per_gene_AF <- rowSums(is.na(atrial_fib_expression_data)) / ncol(atrial_fib_expression_data)
genes_to_keep_AF <- missing_values_per_gene_AF < 0.5
atrial_fib_expression_clean <- atrial_fib_expression_data[genes_to_keep_AF, ]

# Step 3: Remove any genes that are completely missing
genes_with_some_data_AF <- rowSums(!is.na(atrial_fib_expression_clean)) > 0
atrial_fib_expression_clean <- atrial_fib_expression_clean[genes_with_some_data_AF, ]

print(paste("After cleaning: kept", nrow(atrial_fib_expression_clean), "genes x", ncol(atrial_fib_expression_clean), "AF samples"))

# Step 4: Fill in remaining missing values with gene averages
for(gene_row in 1:nrow(atrial_fib_expression_clean)) {
  missing_sample_positions <- is.na(atrial_fib_expression_clean[gene_row, ])
  if(sum(missing_sample_positions) > 0) {
    gene_average <- mean(atrial_fib_expression_clean[gene_row, !missing_sample_positions], na.rm = TRUE)
    atrial_fib_expression_clean[gene_row, missing_sample_positions] <- gene_average
  }
}

# Step 5: Set up batch correction for AF studies
atrial_fib_batch_labels <- as.factor(atrial_fib_sample_info$dataset)
atrial_fib_condition_design <- model.matrix(~group, data = atrial_fib_sample_info)

# Step 6: Run ComBat batch correction for AF studies
atrial_fib_corrected_expression <- ComBat(dat = atrial_fib_expression_clean, 
                                         batch = atrial_fib_batch_labels, 
                                         mod = atrial_fib_condition_design, 
                                         par.prior = TRUE)

print("Processing heart failure studies...")

# Step 1: Get only the heart failure samples
heart_failure_sample_positions <- which(combined_metadata$study_type == "HF")
heart_failure_expression_data <- combined_expr[, heart_failure_sample_positions]
heart_failure_sample_info <- combined_metadata[heart_failure_sample_positions, ]

# Step 2: Remove genes with too many missing values (more than 50% missing)
missing_values_per_gene_HF <- rowSums(is.na(heart_failure_expression_data)) / ncol(heart_failure_expression_data)
genes_to_keep_HF <- missing_values_per_gene_HF < 0.5
heart_failure_expression_clean <- heart_failure_expression_data[genes_to_keep_HF, ]

# Step 3: Remove any genes that are completely missing
genes_with_some_data_HF <- rowSums(!is.na(heart_failure_expression_clean)) > 0
heart_failure_expression_clean <- heart_failure_expression_clean[genes_with_some_data_HF, ]

print(paste("After cleaning: kept", nrow(heart_failure_expression_clean), "genes x", ncol(heart_failure_expression_clean), "HF samples"))

# Step 4: Fill in remaining missing values with gene averages
for(gene_row in 1:nrow(heart_failure_expression_clean)) {
  missing_sample_positions <- is.na(heart_failure_expression_clean[gene_row, ])
  if(sum(missing_sample_positions) > 0) {
    gene_average <- mean(heart_failure_expression_clean[gene_row, !missing_sample_positions], na.rm = TRUE)
    heart_failure_expression_clean[gene_row, missing_sample_positions] <- gene_average
  }
}

# Step 5: Set up batch correction for HF studies
heart_failure_batch_labels <- as.factor(heart_failure_sample_info$dataset)
heart_failure_condition_design <- model.matrix(~group, data = heart_failure_sample_info)

# Step 6: Run ComBat batch correction for HF studies
heart_failure_corrected_expression <- ComBat(dat = heart_failure_expression_clean, 
                                           batch = heart_failure_batch_labels, 
                                           mod = heart_failure_condition_design, 
                                           par.prior = TRUE)

print("Checking how well batch correction worked...")

# Function to calculate how similar samples are within batches (silhouette score)
calculate_batch_mixing_score <- function(original_data, corrected_data, batch_labels, condition_name) {
  
  # Calculate clustering before correction
  sample_distances_before <- dist(t(original_data))
  silhouette_scores_before <- silhouette(as.numeric(batch_labels), sample_distances_before)
  average_silhouette_before <- mean(silhouette_scores_before[,3])
  
  # Calculate clustering after correction
  sample_distances_after <- dist(t(corrected_data))
  silhouette_scores_after <- silhouette(as.numeric(batch_labels), sample_distances_after)
  average_silhouette_after <- mean(silhouette_scores_after[,3])
  
  # Calculate improvement (lower silhouette = better batch mixing)
  improvement_percentage <- round((average_silhouette_before - average_silhouette_after) / average_silhouette_before * 100, 1)
  
  return(data.frame(
    Heart_Condition = condition_name,
    Number_Samples = ncol(original_data),
    Number_Genes = nrow(original_data),
    Datasets_Included = paste(levels(batch_labels), collapse = ", "),
    Batch_Separation_Before = round(average_silhouette_before, 4),
    Batch_Separation_After = round(average_silhouette_after, 4),
    Improvement_Percent = improvement_percentage
  ))
}

# Calculate batch correction effectiveness
atrial_fib_batch_stats <- calculate_batch_mixing_score(atrial_fib_expression_clean, atrial_fib_corrected_expression, atrial_fib_batch_labels, "AF")
heart_failure_batch_stats <- calculate_batch_mixing_score(heart_failure_expression_clean, heart_failure_corrected_expression, heart_failure_batch_labels, "HF")

# Combine statistics
batch_correction_summary <- rbind(atrial_fib_batch_stats, heart_failure_batch_stats)

print("BATCH CORRECTION EFFECTIVENESS:")
print(batch_correction_summary)

print("Creating visualization plots...")

# Function to create before/after comparison plots
create_batch_correction_plots <- function(original_data, corrected_data, sample_info, batch_labels, condition_name, output_filename) {
  
  pdf(output_filename, width = 12, height = 10)
  
  # PCA plots - before correction
  pca_analysis_before <- prcomp(t(original_data), scale = TRUE)
  variance_explained_before <- round(summary(pca_analysis_before)$importance[2, 1:2] * 100, 1)
  
  pca_plot_data_before <- data.frame(
    PC1 = pca_analysis_before$x[,1], 
    PC2 = pca_analysis_before$x[,2], 
    Dataset = sample_info$dataset, 
    Condition = sample_info$group
  )
  
  plot_before <- ggplot(pca_plot_data_before, aes(x = PC1, y = PC2, color = Dataset, shape = Condition)) +
    geom_point(size = 4, alpha = 0.8) + 
    scale_color_manual(values = dataset_colors_for_plots) +
    labs(title = paste(condition_name, "Studies - Before Batch Correction"), 
         x = paste0("PC1 (", variance_explained_before[1], "% variance)"), 
         y = paste0("PC2 (", variance_explained_before[2], "% variance)")) + 
    theme_minimal()
  
  # PCA plots - after correction
  pca_analysis_after <- prcomp(t(corrected_data), scale = TRUE)
  variance_explained_after <- round(summary(pca_analysis_after)$importance[2, 1:2] * 100, 1)
  
  pca_plot_data_after <- data.frame(
    PC1 = pca_analysis_after$x[,1], 
    PC2 = pca_analysis_after$x[,2],
    Dataset = sample_info$dataset, 
    Condition = sample_info$group
  )
  
  plot_after <- ggplot(pca_plot_data_after, aes(x = PC1, y = PC2, color = Dataset, shape = Condition)) +
    geom_point(size = 4, alpha = 0.8) + 
    scale_color_manual(values = dataset_colors_for_plots) +
    labs(title = paste(condition_name, "Studies - After Batch Correction"), 
         x = paste0("PC1 (", variance_explained_after[1], "% variance)"), 
         y = paste0("PC2 (", variance_explained_after[2], "% variance)")) + 
    theme_minimal()
  
  # Show both plots side by side
  grid.arrange(plot_before, plot_after, ncol = 2)
  
  # tSNE plots - before correction
  set.seed(42)
  sample_count <- nrow(sample_info)
  tsne_perplexity <- min(15, floor((sample_count-1)/3))
  
  tsne_analysis_before <- Rtsne(t(original_data), dims = 2, perplexity = tsne_perplexity)
  tsne_plot_data_before <- data.frame(
    tSNE1 = tsne_analysis_before$Y[,1], 
    tSNE2 = tsne_analysis_before$Y[,2],
    Dataset = sample_info$dataset, 
    Condition = sample_info$group
  )
  
  tsne_plot_before <- ggplot(tsne_plot_data_before, aes(x = tSNE1, y = tSNE2, color = Dataset, shape = Condition)) +
    geom_point(size = 4, alpha = 0.8) + 
    scale_color_manual(values = dataset_colors_for_plots) +
    labs(title = paste(condition_name, "Studies tSNE - Before Batch Correction")) + 
    theme_minimal()
  
  # tSNE plots - after correction
  tsne_analysis_after <- Rtsne(t(corrected_data), dims = 2, perplexity = tsne_perplexity)
  tsne_plot_data_after <- data.frame(
    tSNE1 = tsne_analysis_after$Y[,1], 
    tSNE2 = tsne_analysis_after$Y[,2],
    Dataset = sample_info$dataset, 
    Condition = sample_info$group
  )
  
  tsne_plot_after <- ggplot(tsne_plot_data_after, aes(x = tSNE1, y = tSNE2, color = Dataset, shape = Condition)) +
    geom_point(size = 4, alpha = 0.8) + 
    scale_color_manual(values = dataset_colors_for_plots) +
    labs(title = paste(condition_name, "Studies tSNE - After Batch Correction")) + 
    theme_minimal()
  
  # Show tSNE plots side by side
  grid.arrange(tsne_plot_before, tsne_plot_after, ncol = 2)
  
  # Sample correlation heatmaps
  sample_annotation <- data.frame(
    Dataset = sample_info$dataset, 
    Condition = sample_info$group, 
    row.names = colnames(original_data)
  )
  
  # Correlation heatmap - before correction
  pheatmap(cor(original_data), 
          show_rownames = FALSE, 
          show_colnames = FALSE, 
          annotation_col = sample_annotation,
          annotation_colors = list(Dataset = dataset_colors_for_plots, Condition = condition_colors_for_plots),
          main = paste(condition_name, "Studies - Sample Correlation Before Batch Correction"))
  
  # Correlation heatmap - after correction
  pheatmap(cor(corrected_data), 
          show_rownames = FALSE, 
          show_colnames = FALSE, 
          annotation_col = sample_annotation,
          annotation_colors = list(Dataset = dataset_colors_for_plots, Condition = condition_colors_for_plots),
          main = paste(condition_name, "Studies - Sample Correlation After Batch Correction"))
  
  dev.off()
  
  # Return variance information for summary
  return(data.frame(
    Heart_Condition = condition_name,
    PC1_Variance_Before = variance_explained_before[1],
    PC2_Variance_Before = variance_explained_before[2],
    PC1_Variance_After = variance_explained_after[1],
    PC2_Variance_After = variance_explained_after[2],
    PC1_Change = variance_explained_after[1] - variance_explained_before[1],
    PC2_Change = variance_explained_after[2] - variance_explained_before[2]
  ))
}

# Create plots and get variance information
atrial_fib_variance_info <- create_batch_correction_plots(
  atrial_fib_expression_clean, atrial_fib_corrected_expression, atrial_fib_sample_info, atrial_fib_batch_labels, "AF",
  file.path(analysis_results_folder, "AF_batch_correction_report.pdf")
)

heart_failure_variance_info <- create_batch_correction_plots(
  heart_failure_expression_clean, heart_failure_corrected_expression, heart_failure_sample_info, heart_failure_batch_labels, "HF",
  file.path(analysis_results_folder, "HF_batch_correction_report.pdf")
)

# Combine variance information
variance_change_summary <- rbind(atrial_fib_variance_info, heart_failure_variance_info)

print("Saving batch-corrected data and analysis results...")

# Create variables with original names for downstream compatibility
af_combat_expr <- atrial_fib_corrected_expression
af_metadata <- atrial_fib_sample_info
hf_combat_expr <- heart_failure_corrected_expression  
hf_metadata <- heart_failure_sample_info

# Save with original filename that downstream scripts expect
save(af_combat_expr, af_metadata, hf_combat_expr, hf_metadata,
     file = file.path(analysis_results_folder, "02_condition_specific_corrected.RData"))

# Also save individual condition data with beginner-friendly names
save(atrial_fib_corrected_expression, atrial_fib_sample_info, 
     file = file.path(analysis_results_folder, "02_atrial_fib_corrected_data.RData"))
save(heart_failure_corrected_expression, heart_failure_sample_info, 
     file = file.path(analysis_results_folder, "02_heart_failure_corrected_data.RData"))

# Save analysis summaries
write.csv(batch_correction_summary, 
          file.path(analysis_results_folder, "02_batch_correction_effectiveness_summary.csv"), 
          row.names = FALSE)
write.csv(atrial_fib_batch_stats, 
          file.path(analysis_results_folder, "02_atrial_fib_batch_correction_stats.csv"), 
          row.names = FALSE)
write.csv(heart_failure_batch_stats, 
          file.path(analysis_results_folder, "02_heart_failure_batch_correction_stats.csv"), 
          row.names = FALSE)
write.csv(variance_change_summary, 
          file.path(analysis_results_folder, "02_pca_variance_changes.csv"), 
          row.names = FALSE)

print("Batch correction complete!")
print("")
print("SUMMARY:")
print("- Corrected batch effects in atrial fibrillation studies")
print("- Corrected batch effects in heart failure studies") 
print("- Created visualization plots showing before/after comparison")
print("- Saved batch-corrected data for next analysis steps")
print("- Data is now ready for differential gene expression analysis")