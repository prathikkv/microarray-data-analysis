# Load and combine cardiac datasets - BEGINNER VERSION
# This script converts microarray probe data to gene data and combines all datasets

library(GEOquery)       # Package to work with GEO datasets
library(hgu133plus2.db) # Database to convert probes to gene names
library(ggplot2)        # Package for making plots
library(pheatmap)       # Package for heatmaps
library(VennDiagram)    # Package for Venn diagrams
library(RColorBrewer)   # Package for color palettes

# Set up our folder paths
main_data_folder <- "/Users/macbookpro/Desktop/Ananylum/data"
analysis_results_folder <- "results"
# Make sure folders exist
dir.create(analysis_results_folder, showWarnings = FALSE)

print("Loading cardiac datasets from downloaded files...")

# Load each dataset from the downloaded files
atrial_fibrillation_dataset_1 <- getGEO(filename = file.path(main_data_folder, "GSE115574/GSE115574_series_matrix.txt.gz"))
atrial_fibrillation_dataset_2 <- getGEO(filename = file.path(main_data_folder, "GSE79768/GSE79768_series_matrix.txt.gz"))
heart_failure_dataset_1 <- getGEO(filename = file.path(main_data_folder, "GSE76701/GSE76701_series_matrix.txt.gz"))
heart_failure_dataset_2 <- getGEO(filename = file.path(main_data_folder, "GSE57338/GSE57338_series_matrix.txt.gz"))

# Extract gene expression numbers from each dataset
atrial_fib_1_expression_data <- exprs(atrial_fibrillation_dataset_1)
atrial_fib_2_expression_data <- exprs(atrial_fibrillation_dataset_2)
heart_fail_1_expression_data <- exprs(heart_failure_dataset_1)
heart_fail_2_expression_data <- exprs(heart_failure_dataset_2)

print("Determining which samples have which heart conditions...")

# Figure out which samples have atrial fibrillation vs normal heart rhythm
atrial_fib_1_sample_conditions <- ifelse(grepl("atrial fibrillation", pData(atrial_fibrillation_dataset_1)$characteristics_ch1), "AF", "SR")
atrial_fib_2_sample_conditions <- ifelse(grepl("Atrial Fibrillation", pData(atrial_fibrillation_dataset_2)$title), "AF", "SR")

# Figure out which samples have heart failure vs healthy hearts
heart_fail_1_sample_conditions <- ifelse(grepl("Non-failing", pData(heart_failure_dataset_1)$title), "Control", "HF")
heart_fail_2_sample_conditions <- ifelse(grepl("Non-failing", pData(heart_failure_dataset_2)$title), "Control", "HF")

print("Converting microarray probes to actual gene names...")

# Get all probe IDs from first dataset (they use same chip type)
all_probe_ids_on_chip <- rownames(atrial_fib_1_expression_data)

# Look up what genes these probes measure
probe_to_gene_lookup_table <- AnnotationDbi::select(hgu133plus2.db, keys = all_probe_ids_on_chip, columns = "SYMBOL", keytype = "PROBEID")

# Keep only probes that actually have known gene names
good_probe_to_gene_mappings <- probe_to_gene_lookup_table[!is.na(probe_to_gene_lookup_table$SYMBOL) & probe_to_gene_lookup_table$SYMBOL != "", ]

# Initialize detailed tracking data structures
probe_annotation_summary <- data.frame()
gene_mapping_details <- data.frame()
data_integration_metrics <- list()

# Track HGU133Plus2 chip annotation details
hgu133_total_probes_in_annotation <- length(all_probe_ids_on_chip)
hgu133_probes_with_genes <- nrow(good_probe_to_gene_mappings)
hgu133_probes_without_genes <- nrow(probe_to_gene_lookup_table) - hgu133_probes_with_genes
hgu133_mapping_efficiency <- round(hgu133_probes_with_genes / hgu133_total_probes_in_annotation * 100, 2)

print("Converting Dataset 1 (Atrial Fibrillation) probes to genes...")
# Track dataset-specific probe statistics
dataset_1_total_probes <- nrow(atrial_fib_1_expression_data)
dataset_1_mapped_probes <- sum(rownames(atrial_fib_1_expression_data) %in% good_probe_to_gene_mappings$PROBEID)
dataset_1_unmapped_probes <- dataset_1_total_probes - dataset_1_mapped_probes

# Use efficient method: for each gene, if multiple probes measure it, take the median value
atrial_fib_1_gene_data_list <- tapply(1:nrow(good_probe_to_gene_mappings), good_probe_to_gene_mappings$SYMBOL, function(probe_indices) {
  probe_ids_for_gene <- good_probe_to_gene_mappings$PROBEID[probe_indices]
  
  if(length(probe_ids_for_gene) == 1) {
    # Only one probe for this gene - use it directly
    atrial_fib_1_expression_data[probe_ids_for_gene, ]
  } else {
    # Multiple probes for this gene - take median across probes for each sample
    apply(atrial_fib_1_expression_data[probe_ids_for_gene, , drop = FALSE], 2, median)
  }
})
# Convert list to matrix
atrial_fib_1_gene_data <- do.call(rbind, atrial_fib_1_gene_data_list)
dataset_1_genes_measured <- nrow(atrial_fib_1_gene_data)

print("Converting Dataset 2 (Atrial Fibrillation) probes to genes...")
# Track dataset-specific probe statistics
dataset_2_total_probes <- nrow(atrial_fib_2_expression_data)
dataset_2_mapped_probes <- sum(rownames(atrial_fib_2_expression_data) %in% good_probe_to_gene_mappings$PROBEID)
dataset_2_unmapped_probes <- dataset_2_total_probes - dataset_2_mapped_probes

# Same efficient process for second atrial fibrillation dataset
atrial_fib_2_gene_data_list <- tapply(1:nrow(good_probe_to_gene_mappings), good_probe_to_gene_mappings$SYMBOL, function(probe_indices) {
  probe_ids_for_gene <- good_probe_to_gene_mappings$PROBEID[probe_indices]
  
  if(length(probe_ids_for_gene) == 1) {
    atrial_fib_2_expression_data[probe_ids_for_gene, ]
  } else {
    apply(atrial_fib_2_expression_data[probe_ids_for_gene, , drop = FALSE], 2, median)
  }
})
atrial_fib_2_gene_data <- do.call(rbind, atrial_fib_2_gene_data_list)
dataset_2_genes_measured <- nrow(atrial_fib_2_gene_data)

print("Converting Dataset 3 (Heart Failure) probes to genes...")
# Track dataset-specific probe statistics
dataset_3_total_probes <- nrow(heart_fail_1_expression_data)
dataset_3_mapped_probes <- sum(rownames(heart_fail_1_expression_data) %in% good_probe_to_gene_mappings$PROBEID)
dataset_3_unmapped_probes <- dataset_3_total_probes - dataset_3_mapped_probes

# Same efficient process for first heart failure dataset
heart_fail_1_gene_data_list <- tapply(1:nrow(good_probe_to_gene_mappings), good_probe_to_gene_mappings$SYMBOL, function(probe_indices) {
  probe_ids_for_gene <- good_probe_to_gene_mappings$PROBEID[probe_indices]
  
  if(length(probe_ids_for_gene) == 1) {
    heart_fail_1_expression_data[probe_ids_for_gene, ]
  } else {
    apply(heart_fail_1_expression_data[probe_ids_for_gene, , drop = FALSE], 2, median)
  }
})
heart_fail_1_gene_data <- do.call(rbind, heart_fail_1_gene_data_list)
dataset_3_genes_measured <- nrow(heart_fail_1_gene_data)

print("Converting Dataset 4 (Heart Failure) - uses different chip type...")
# Fourth dataset uses different microarray chip, need different conversion method
different_chip_annotation <- getGEO("GPL11532", destdir = tempdir())
different_chip_lookup_table <- Table(different_chip_annotation)

# Track different chip annotation details
gpl11532_total_probes_in_annotation <- nrow(different_chip_lookup_table)

# Extract gene names from the chip annotation
gene_names_from_different_chip <- c()
for(i in 1:nrow(different_chip_lookup_table)) {
  gene_assignment_text <- different_chip_lookup_table$gene_assignment[i]
  
  if(is.na(gene_assignment_text) || gene_assignment_text == "---") {
    gene_names_from_different_chip[i] <- NA
  } else {
    # Gene name is in the second part after splitting by "//"
    assignment_parts <- strsplit(gene_assignment_text, "//")[[1]]
    if(length(assignment_parts) >= 2) {
      gene_names_from_different_chip[i] <- trimws(assignment_parts[2])
    } else {
      gene_names_from_different_chip[i] <- NA
    }
  }
}

# Create lookup table for this different chip
different_chip_probe_to_gene <- data.frame(
  probe = different_chip_lookup_table$ID, 
  gene = gene_names_from_different_chip, 
  stringsAsFactors = FALSE
)
different_chip_probe_to_gene <- different_chip_probe_to_gene[!is.na(different_chip_probe_to_gene$gene) & different_chip_probe_to_gene$gene != "", ]

# Track GPL11532 chip mapping efficiency
gpl11532_probes_with_genes <- nrow(different_chip_probe_to_gene)
gpl11532_probes_without_genes <- gpl11532_total_probes_in_annotation - gpl11532_probes_with_genes
gpl11532_mapping_efficiency <- round(gpl11532_probes_with_genes / gpl11532_total_probes_in_annotation * 100, 2)

# Find which probes from dataset 4 we can convert
dataset_4_probe_ids <- rownames(heart_fail_2_expression_data)
convertible_probes <- intersect(dataset_4_probe_ids, different_chip_probe_to_gene$probe)
probe_gene_mapping_dataset_4 <- different_chip_probe_to_gene[match(convertible_probes, different_chip_probe_to_gene$probe), ]
heart_fail_2_expression_subset <- heart_fail_2_expression_data[convertible_probes, ]

# Track dataset 4 specific probe statistics
dataset_4_total_probes <- length(dataset_4_probe_ids)
dataset_4_mapped_probes <- length(convertible_probes)
dataset_4_unmapped_probes <- dataset_4_total_probes - dataset_4_mapped_probes

# Convert dataset 4 probes to genes using efficient method
heart_fail_2_gene_data_list <- tapply(1:nrow(heart_fail_2_expression_subset), probe_gene_mapping_dataset_4$gene, function(probe_indices) {
  if(length(probe_indices) == 1) {
    heart_fail_2_expression_subset[probe_indices, ]
  } else {
    apply(heart_fail_2_expression_subset[probe_indices, , drop = FALSE], 2, median)
  }
})
heart_fail_2_gene_data <- do.call(rbind, heart_fail_2_gene_data_list)
dataset_4_genes_measured <- nrow(heart_fail_2_gene_data)

print("Combining all datasets together...")

# Find ALL genes measured across any of the 4 studies
all_unique_gene_names <- unique(c(
  rownames(atrial_fib_1_gene_data), 
  rownames(atrial_fib_2_gene_data), 
  rownames(heart_fail_1_gene_data), 
  rownames(heart_fail_2_gene_data)
))

# Calculate total number of samples across all studies
total_samples_count <- (
  ncol(atrial_fib_1_gene_data) + 
  ncol(atrial_fib_2_gene_data) + 
  ncol(heart_fail_1_gene_data) + 
  ncol(heart_fail_2_gene_data)
)

# Create big combined expression matrix
combined_expression_matrix <- matrix(NA, nrow = length(all_unique_gene_names), ncol = total_samples_count)
rownames(combined_expression_matrix) <- all_unique_gene_names

# Create combined sample names
combined_sample_names <- c(
  colnames(atrial_fib_1_gene_data), 
  colnames(atrial_fib_2_gene_data), 
  colnames(heart_fail_1_gene_data), 
  colnames(heart_fail_2_gene_data)
)
colnames(combined_expression_matrix) <- combined_sample_names

print("Filling in gene expression data for each study...")

# Fill in data from each study
for(current_gene in all_unique_gene_names) {
  sample_position <- 1
  
  # Add data from atrial fibrillation study 1
  if(current_gene %in% rownames(atrial_fib_1_gene_data)) {
    end_position <- sample_position + ncol(atrial_fib_1_gene_data) - 1
    combined_expression_matrix[current_gene, sample_position:end_position] <- atrial_fib_1_gene_data[current_gene, ]
  }
  sample_position <- sample_position + ncol(atrial_fib_1_gene_data)
  
  # Add data from atrial fibrillation study 2
  if(current_gene %in% rownames(atrial_fib_2_gene_data)) {
    end_position <- sample_position + ncol(atrial_fib_2_gene_data) - 1
    combined_expression_matrix[current_gene, sample_position:end_position] <- atrial_fib_2_gene_data[current_gene, ]
  }
  sample_position <- sample_position + ncol(atrial_fib_2_gene_data)
  
  # Add data from heart failure study 1
  if(current_gene %in% rownames(heart_fail_1_gene_data)) {
    end_position <- sample_position + ncol(heart_fail_1_gene_data) - 1
    combined_expression_matrix[current_gene, sample_position:end_position] <- heart_fail_1_gene_data[current_gene, ]
  }
  sample_position <- sample_position + ncol(heart_fail_1_gene_data)
  
  # Add data from heart failure study 2
  if(current_gene %in% rownames(heart_fail_2_gene_data)) {
    end_position <- sample_position + ncol(heart_fail_2_gene_data) - 1
    combined_expression_matrix[current_gene, sample_position:end_position] <- heart_fail_2_gene_data[current_gene, ]
  }
}

print("Creating sample information table...")

# Create table describing each sample (using original variable names for compatibility)
combined_metadata <- data.frame(
  sample_id = combined_sample_names,
  dataset = c(
    rep("GSE115574", ncol(atrial_fib_1_gene_data)), 
    rep("GSE79768", ncol(atrial_fib_2_gene_data)),
    rep("GSE76701", ncol(heart_fail_1_gene_data)), 
    rep("GSE57338", ncol(heart_fail_2_gene_data))
  ),
  condition = c(
    atrial_fib_1_sample_conditions, 
    atrial_fib_2_sample_conditions, 
    heart_fail_1_sample_conditions, 
    heart_fail_2_sample_conditions
  ),
  study_type = c(
    rep("AF", ncol(atrial_fib_1_gene_data) + ncol(atrial_fib_2_gene_data)), 
    rep("HF", ncol(heart_fail_1_gene_data) + ncol(heart_fail_2_gene_data))
  ),
  stringsAsFactors = FALSE
)

# Add standardized group column for downstream analysis (make "SR" into "Control" for consistency)
combined_metadata$group <- combined_metadata$condition
combined_metadata$group[combined_metadata$group == "SR"] <- "Control"

print("Checking for important cardiac genes...")

# Check that important cardiac and CAMK genes are present in our final dataset
important_cardiac_genes <- c("NPPA", "NPPB", "MYH6", "MYH7", "CAMK2A", "CAMK2B", 
                            "CAMK2D", "CAMK2G", "CAMK1", "CAMK1D", "CAMK4", "TNNT2", "TNNI3")
genes_found_in_data <- important_cardiac_genes[important_cardiac_genes %in% all_unique_gene_names]

print(paste("Important cardiac genes found:", length(genes_found_in_data), "out of", length(important_cardiac_genes)))
print(paste("Present genes:", paste(genes_found_in_data, collapse = ", ")))
if(length(genes_found_in_data) < length(important_cardiac_genes)) {
  missing_genes <- setdiff(important_cardiac_genes, genes_found_in_data)
  print(paste("Missing genes:", paste(missing_genes, collapse = ", ")))
}

# Collect comprehensive gene mapping and integration statistics
print("Collecting detailed statistics for CSV reports...")

# Calculate probe-to-gene mapping ratios for HGU133Plus2 chip
gene_probe_counts_hgu133 <- table(good_probe_to_gene_mappings$SYMBOL)
single_probe_genes_hgu133 <- sum(gene_probe_counts_hgu133 == 1)
multi_probe_genes_hgu133 <- sum(gene_probe_counts_hgu133 > 1)
max_probes_per_gene_hgu133 <- max(gene_probe_counts_hgu133)
avg_probes_per_gene_hgu133 <- round(mean(gene_probe_counts_hgu133), 2)

# Calculate probe-to-gene mapping ratios for GPL11532 chip
gene_probe_counts_gpl11532 <- table(different_chip_probe_to_gene$gene)
single_probe_genes_gpl11532 <- sum(gene_probe_counts_gpl11532 == 1)
multi_probe_genes_gpl11532 <- sum(gene_probe_counts_gpl11532 > 1)
max_probes_per_gene_gpl11532 <- max(gene_probe_counts_gpl11532)
avg_probes_per_gene_gpl11532 <- round(mean(gene_probe_counts_gpl11532), 2)

# Gene coverage across datasets
genes_dataset_1 <- rownames(atrial_fib_1_gene_data)
genes_dataset_2 <- rownames(atrial_fib_2_gene_data)
genes_dataset_3 <- rownames(heart_fail_1_gene_data)
genes_dataset_4 <- rownames(heart_fail_2_gene_data)

# Gene overlap analysis
genes_in_all_4_datasets <- Reduce(intersect, list(genes_dataset_1, genes_dataset_2, genes_dataset_3, genes_dataset_4))
genes_in_3_datasets <- length(all_unique_gene_names[rowSums(cbind(
  all_unique_gene_names %in% genes_dataset_1,
  all_unique_gene_names %in% genes_dataset_2,
  all_unique_gene_names %in% genes_dataset_3,
  all_unique_gene_names %in% genes_dataset_4
)) >= 3])
genes_in_2_datasets <- length(all_unique_gene_names[rowSums(cbind(
  all_unique_gene_names %in% genes_dataset_1,
  all_unique_gene_names %in% genes_dataset_2,
  all_unique_gene_names %in% genes_dataset_3,
  all_unique_gene_names %in% genes_dataset_4
)) >= 2])
genes_in_1_dataset_only <- length(all_unique_gene_names[rowSums(cbind(
  all_unique_gene_names %in% genes_dataset_1,
  all_unique_gene_names %in% genes_dataset_2,
  all_unique_gene_names %in% genes_dataset_3,
  all_unique_gene_names %in% genes_dataset_4
)) == 1])

# Missing data analysis in final combined matrix
total_possible_values <- nrow(combined_expression_matrix) * ncol(combined_expression_matrix)
missing_values <- sum(is.na(combined_expression_matrix))
data_completeness_percentage <- round((total_possible_values - missing_values) / total_possible_values * 100, 2)

print("Final dataset summary:")
print(paste("- Total samples:", ncol(combined_expression_matrix)))
print(paste("- Total genes:", nrow(combined_expression_matrix)))
print(paste("- AF samples:", sum(combined_metadata$study_type == "AF")))
print(paste("- HF samples:", sum(combined_metadata$study_type == "HF")))
print(paste("- Data completeness:", data_completeness_percentage, "%"))
print("")

print("Creating comprehensive quality control and visualization plots...")

# Create comprehensive plots for data loading QC
pdf(file.path(analysis_results_folder, "01_data_loading_quality_control.pdf"), width = 14, height = 10)

# 1. Probe mapping summary bar plot
probe_mapping_data <- data.frame(
  Dataset = c("GSE115574", "GSE79768", "GSE76701", "GSE57338"),
  Dataset_Name = c("AF_1", "AF_2", "HF_1", "HF_2"),
  Total_Probes = c(dataset_1_total_probes, dataset_2_total_probes, dataset_3_total_probes, dataset_4_total_probes),
  Mapped_Probes = c(dataset_1_mapped_probes, dataset_2_mapped_probes, dataset_3_mapped_probes, dataset_4_mapped_probes),
  Unmapped_Probes = c(dataset_1_unmapped_probes, dataset_2_unmapped_probes, dataset_3_unmapped_probes, dataset_4_unmapped_probes)
)

library(reshape2)
probe_mapping_long <- melt(probe_mapping_data[,c("Dataset_Name", "Mapped_Probes", "Unmapped_Probes")], 
                          id.vars = "Dataset_Name", variable.name = "Mapping_Status", value.name = "Count")

mapping_plot <- ggplot(probe_mapping_long, aes(x = Dataset_Name, y = Count, fill = Mapping_Status)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Mapped_Probes" = "#2E8B57", "Unmapped_Probes" = "#CD5C5C")) +
  labs(title = "Probe Mapping Success Across Datasets",
       subtitle = "Green: Successfully mapped to genes, Red: Failed mapping",
       x = "Dataset", y = "Number of Probes", fill = "Mapping Status") +
  theme_minimal() +
  geom_text(data = probe_mapping_data, aes(x = Dataset_Name, y = Total_Probes + 500, 
                                          label = paste0(round(Mapped_Probes/Total_Probes*100, 1), "%")),
           inherit.aes = FALSE, vjust = 0)

print(mapping_plot)

# 2. Gene coverage Venn diagram (simplified for 4 datasets)
print("Creating gene overlap visualization...")

# For visualization, we'll show the overlap between conditions
af_genes <- unique(c(genes_dataset_1, genes_dataset_2))
hf_genes <- unique(c(genes_dataset_3, genes_dataset_4))
all_genes_set <- unique(c(af_genes, hf_genes))

# Create overlap data
gene_overlap_data <- data.frame(
  Gene_Set = c("AF_only", "HF_only", "Shared", "AF_total", "HF_total", "Combined_total"),
  Count = c(
    length(setdiff(af_genes, hf_genes)),
    length(setdiff(hf_genes, af_genes)),
    length(intersect(af_genes, hf_genes)),
    length(af_genes),
    length(hf_genes),
    length(all_genes_set)
  )
)

overlap_plot <- ggplot(gene_overlap_data[1:3,], aes(x = Gene_Set, y = Count, fill = Gene_Set)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("AF_only" = "#FF6B6B", "HF_only" = "#4ECDC4", "Shared" = "#45B7D1")) +
  labs(title = "Gene Coverage Overlap Between AF and HF Studies",
       subtitle = "Showing genes unique to each condition vs shared genes",
       x = "Gene Set", y = "Number of Genes") +
  theme_minimal() +
  geom_text(aes(label = Count), vjust = -0.3)

print(overlap_plot)

# 3. Expression distribution plots (sample from combined matrix)
set.seed(123)
sample_genes <- sample(rownames(combined_expression_matrix)[!is.na(rowSums(combined_expression_matrix))], 1000)
sample_expression <- combined_expression_matrix[sample_genes, ]

# Box plot of expression distributions by dataset
sample_metadata <- data.frame(
  Sample_ID = colnames(sample_expression),
  Dataset = combined_metadata$dataset,
  Condition = combined_metadata$group
)

expression_long <- melt(as.matrix(sample_expression))
names(expression_long) <- c("Gene", "Sample_ID", "Expression")
expression_long <- merge(expression_long, sample_metadata, by = "Sample_ID")
expression_long <- expression_long[!is.na(expression_long$Expression), ]

expression_boxplot <- ggplot(expression_long, aes(x = Dataset, y = Expression, fill = Condition)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c("AF" = "#FF6B6B", "HF" = "#4ECDC4", "Control" = "#95E4D3", "SR" = "#95E4D3")) +
  labs(title = "Expression Distributions Across Datasets (Sample of 1000 genes)",
       subtitle = "Boxplots show median, quartiles, and outliers",
       x = "Dataset", y = "Log2 Expression Level") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(expression_boxplot)

# 4. Missing data heatmap
missing_data_summary <- data.frame(
  Dataset = c("GSE115574", "GSE79768", "GSE76701", "GSE57338"),
  Total_Genes = c(dataset_1_genes_measured, dataset_2_genes_measured, dataset_3_genes_measured, dataset_4_genes_measured),
  Missing_in_Combined = c(
    sum(!all_unique_gene_names %in% genes_dataset_1),
    sum(!all_unique_gene_names %in% genes_dataset_2),
    sum(!all_unique_gene_names %in% genes_dataset_3),
    sum(!all_unique_gene_names %in% genes_dataset_4)
  )
)

missing_data_summary$Present_in_Combined <- missing_data_summary$Total_Genes
missing_data_long <- melt(missing_data_summary[,c("Dataset", "Present_in_Combined", "Missing_in_Combined")],
                         id.vars = "Dataset", variable.name = "Status", value.name = "Count")

missing_plot <- ggplot(missing_data_long, aes(x = Dataset, y = Count, fill = Status)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Present_in_Combined" = "#2E8B57", "Missing_in_Combined" = "#FF6B6B")) +
  labs(title = "Gene Availability Across Datasets in Combined Matrix",
       subtitle = "Green: Genes available in combined dataset, Red: Genes missing",
       x = "Dataset", y = "Number of Genes") +
  theme_minimal()

print(missing_plot)

# 5. Important genes tracking
important_genes_plot_data <- data.frame(
  Gene = important_cardiac_genes,
  Found = important_cardiac_genes %in% all_unique_gene_names,
  AF_1 = important_cardiac_genes %in% genes_dataset_1,
  AF_2 = important_cardiac_genes %in% genes_dataset_2,
  HF_1 = important_cardiac_genes %in% genes_dataset_3,
  HF_2 = important_cardiac_genes %in% genes_dataset_4
)

# Convert to long format for heatmap
important_genes_long <- melt(important_genes_plot_data, id.vars = "Gene", variable.name = "Dataset", value.name = "Present")
important_genes_long$Present <- as.numeric(important_genes_long$Present)

important_genes_heatmap <- ggplot(important_genes_long, aes(x = Dataset, y = Gene, fill = factor(Present))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("0" = "#FFE5E5", "1" = "#2E8B57"), labels = c("Absent", "Present")) +
  labs(title = "Important Cardiac Genes Availability Across Datasets",
       subtitle = "Tracking key cardiac and CAMK genes across all 4 studies",
       x = "Dataset", y = "Gene", fill = "Gene Present") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(important_genes_heatmap)

dev.off()

print("Quality control plots saved successfully!")
print(paste("- Saved to:", file.path(analysis_results_folder, "01_data_loading_quality_control.pdf")))
print("")

print("VISUALIZATION SUMMARY:")
print("- Probe mapping success rates across datasets")
print("- Gene coverage overlap between AF and HF studies") 
print("- Expression distributions to check data quality")
print("- Missing data patterns in combined matrix")
print("- Important cardiac genes availability tracking")
print("")

print("Generating detailed CSV reports...")

# CSV 1: Probe Annotation Summary
probe_annotation_summary <- data.frame(
  Annotation_Database = c("HGU133Plus2.db", "HGU133Plus2.db", "GPL11532", "GPL11532"),
  Dataset = c("GSE115574,GSE79768,GSE76701", "GSE115574,GSE79768,GSE76701", "GSE57338", "GSE57338"),
  Chip_Type = c("HGU133Plus2", "HGU133Plus2", "GPL11532", "GPL11532"),
  Metric = c("Total_Probes_in_Annotation", "Probes_with_Gene_Mapping", "Total_Probes_in_Annotation", "Probes_with_Gene_Mapping"),
  Count = c(hgu133_total_probes_in_annotation, hgu133_probes_with_genes, gpl11532_total_probes_in_annotation, gpl11532_probes_with_genes),
  Percentage = c(100, hgu133_mapping_efficiency, 100, gpl11532_mapping_efficiency),
  stringsAsFactors = FALSE
)

# Add dataset-specific probe statistics
dataset_probe_stats <- data.frame(
  Dataset = c("GSE115574", "GSE79768", "GSE76701", "GSE57338"),
  Dataset_Name = c("Atrial_Fibrillation_1", "Atrial_Fibrillation_2", "Heart_Failure_1", "Heart_Failure_2"),
  Total_Probes_in_Raw_Data = c(dataset_1_total_probes, dataset_2_total_probes, dataset_3_total_probes, dataset_4_total_probes),
  Mapped_Probes = c(dataset_1_mapped_probes, dataset_2_mapped_probes, dataset_3_mapped_probes, dataset_4_mapped_probes),
  Unmapped_Probes = c(dataset_1_unmapped_probes, dataset_2_unmapped_probes, dataset_3_unmapped_probes, dataset_4_unmapped_probes),
  Mapping_Efficiency_Percent = c(
    round(dataset_1_mapped_probes/dataset_1_total_probes*100, 2),
    round(dataset_2_mapped_probes/dataset_2_total_probes*100, 2),
    round(dataset_3_mapped_probes/dataset_3_total_probes*100, 2),
    round(dataset_4_mapped_probes/dataset_4_total_probes*100, 2)
  ),
  Final_Genes_Measured = c(dataset_1_genes_measured, dataset_2_genes_measured, dataset_3_genes_measured, dataset_4_genes_measured),
  stringsAsFactors = FALSE
)

write.csv(rbind(
  data.frame(Section = "Annotation_Overview", probe_annotation_summary),
  data.frame(Section = "Dataset_Specific_Mapping", 
             Annotation_Database = dataset_probe_stats$Dataset,
             Dataset = dataset_probe_stats$Dataset_Name,
             Chip_Type = c(rep("HGU133Plus2", 3), "GPL11532"),
             Metric = "Dataset_Stats",
             Count = dataset_probe_stats$Total_Probes_in_Raw_Data,
             Percentage = dataset_probe_stats$Mapping_Efficiency_Percent)
), file = file.path(analysis_results_folder, "01_probe_annotation_summary.csv"), row.names = FALSE)

# CSV 2: Gene Mapping Details
gene_coverage_matrix <- data.frame(
  Gene = all_unique_gene_names,
  Present_in_GSE115574 = all_unique_gene_names %in% genes_dataset_1,
  Present_in_GSE79768 = all_unique_gene_names %in% genes_dataset_2,
  Present_in_GSE76701 = all_unique_gene_names %in% genes_dataset_3,
  Present_in_GSE57338 = all_unique_gene_names %in% genes_dataset_4,
  stringsAsFactors = FALSE
)
gene_coverage_matrix$Total_Datasets_Present <- rowSums(gene_coverage_matrix[,2:5])

# Probe-to-gene mapping details
probe_gene_mapping_summary <- data.frame(
  Chip_Type = c("HGU133Plus2", "HGU133Plus2", "HGU133Plus2", "GPL11532", "GPL11532", "GPL11532"),
  Mapping_Type = c("Single_Probe_Genes", "Multi_Probe_Genes", "Max_Probes_Per_Gene", 
                   "Single_Probe_Genes", "Multi_Probe_Genes", "Max_Probes_Per_Gene"),
  Count = c(single_probe_genes_hgu133, multi_probe_genes_hgu133, max_probes_per_gene_hgu133,
            single_probe_genes_gpl11532, multi_probe_genes_gpl11532, max_probes_per_gene_gpl11532),
  Average_Probes_Per_Gene = c(avg_probes_per_gene_hgu133, avg_probes_per_gene_hgu133, avg_probes_per_gene_hgu133,
                             avg_probes_per_gene_gpl11532, avg_probes_per_gene_gpl11532, avg_probes_per_gene_gpl11532),
  stringsAsFactors = FALSE
)

# Important genes tracking
important_genes_status <- data.frame(
  Gene = important_cardiac_genes,
  Found_in_Final_Dataset = important_cardiac_genes %in% all_unique_gene_names,
  Present_in_GSE115574 = important_cardiac_genes %in% genes_dataset_1,
  Present_in_GSE79768 = important_cardiac_genes %in% genes_dataset_2,
  Present_in_GSE76701 = important_cardiac_genes %in% genes_dataset_3,
  Present_in_GSE57338 = important_cardiac_genes %in% genes_dataset_4,
  stringsAsFactors = FALSE
)

# Combine all gene mapping details into one dataframe
gene_coverage_summary <- data.frame(
  Section = "Gene_Coverage_Summary",
  Metric = c("Genes_in_all_4_datasets", "Genes_in_3_or_more_datasets", "Genes_in_2_or_more_datasets", "Genes_in_only_1_dataset"),
  Count = c(length(genes_in_all_4_datasets), genes_in_3_datasets, genes_in_2_datasets, genes_in_1_dataset_only),
  stringsAsFactors = FALSE
)

probe_gene_mapping_summary$Section <- "Probe_Gene_Mapping"
important_genes_status$Section <- "Important_Cardiac_Genes"

# Ensure all data frames have consistent column structure
gene_mapping_data_1 <- data.frame(
  Section = gene_coverage_summary$Section,
  Metric = gene_coverage_summary$Metric,
  Count = gene_coverage_summary$Count,
  Additional_Info = "",
  stringsAsFactors = FALSE
)

gene_mapping_data_2 <- data.frame(
  Section = probe_gene_mapping_summary$Section, 
  Metric = probe_gene_mapping_summary$Mapping_Type,
  Count = probe_gene_mapping_summary$Count,
  Additional_Info = paste("Avg_per_gene:", probe_gene_mapping_summary$Average_Probes_Per_Gene),
  stringsAsFactors = FALSE
)

gene_mapping_data_3 <- data.frame(
  Section = important_genes_status$Section,
  Metric = important_genes_status$Gene,
  Count = as.numeric(important_genes_status$Found_in_Final_Dataset),
  Additional_Info = paste("Datasets:", rowSums(important_genes_status[,3:6])),
  stringsAsFactors = FALSE
)

write.csv(rbind(gene_mapping_data_1, gene_mapping_data_2, gene_mapping_data_3), 
          file = file.path(analysis_results_folder, "01_gene_mapping_details.csv"), row.names = FALSE)

# CSV 3: Data Integration Metrics
sample_summary <- data.frame(
  Dataset = c("GSE115574", "GSE79768", "GSE76701", "GSE57338", "Total"),
  Dataset_Name = c("Atrial_Fibrillation_1", "Atrial_Fibrillation_2", "Heart_Failure_1", "Heart_Failure_2", "Combined"),
  Study_Type = c("AF", "AF", "HF", "HF", "All"),
  Total_Samples = c(ncol(atrial_fib_1_gene_data), ncol(atrial_fib_2_gene_data), ncol(heart_fail_1_gene_data), ncol(heart_fail_2_gene_data), ncol(combined_expression_matrix)),
  AF_Samples = c(sum(atrial_fib_1_sample_conditions == "AF"), sum(atrial_fib_2_sample_conditions == "AF"), 0, 0, sum(combined_metadata$condition == "AF")),
  HF_Samples = c(0, 0, sum(heart_fail_1_sample_conditions == "HF"), sum(heart_fail_2_sample_conditions == "HF"), sum(combined_metadata$condition == "HF")),
  Control_Samples = c(sum(atrial_fib_1_sample_conditions == "SR"), sum(atrial_fib_2_sample_conditions == "SR"), 
                     sum(heart_fail_1_sample_conditions == "Control"), sum(heart_fail_2_sample_conditions == "Control"), 
                     sum(combined_metadata$group == "Control")),
  stringsAsFactors = FALSE
)

quality_metrics <- data.frame(
  Metric = c("Total_Unique_Genes", "Total_Samples", "Total_Possible_Data_Points", "Missing_Data_Points", 
             "Data_Completeness_Percent", "Important_Cardiac_Genes_Found", "Important_Cardiac_Genes_Total"),
  Value = c(nrow(combined_expression_matrix), ncol(combined_expression_matrix), total_possible_values, 
            missing_values, data_completeness_percentage, length(genes_found_in_data), length(important_cardiac_genes)),
  stringsAsFactors = FALSE
)

# Combine all integration metrics into one dataframe
sample_summary$Section <- "Sample_Summary"
quality_metrics$Section <- "Quality_Control_Metrics"
dataset_probe_stats$Section <- "Dataset_Contribution"

write.csv(rbind(
  data.frame(Section = sample_summary$Section, 
             Metric = paste(sample_summary$Dataset, sample_summary$Dataset_Name, sep="_"),
             Value = sample_summary$Total_Samples,
             Additional_Info = paste("AF:", sample_summary$AF_Samples, "HF:", sample_summary$HF_Samples, "Control:", sample_summary$Control_Samples)),
  data.frame(Section = quality_metrics$Section,
             Metric = quality_metrics$Metric,
             Value = quality_metrics$Value,
             Additional_Info = ""),
  data.frame(Section = dataset_probe_stats$Section,
             Metric = paste(dataset_probe_stats$Dataset, "mapping_efficiency", sep="_"),
             Value = dataset_probe_stats$Mapping_Efficiency_Percent,
             Additional_Info = paste("Mapped:", dataset_probe_stats$Mapped_Probes, "of", dataset_probe_stats$Total_Probes_in_Raw_Data))
), file = file.path(analysis_results_folder, "01_data_integration_metrics.csv"), row.names = FALSE)

print("CSV reports generated successfully!")
print("- 01_probe_annotation_summary.csv: Probe and annotation details")
print("- 01_gene_mapping_details.csv: Gene coverage and mapping ratios") 
print("- 01_data_integration_metrics.csv: Integration quality and sample summary")

print("Saving combined dataset...")

# Save the final combined dataset for next analysis step (using original filename)
# Change variable name for compatibility with downstream scripts
combined_expr <- combined_expression_matrix
save(combined_expr, combined_metadata, 
     file = file.path(analysis_results_folder, "01_combined_4datasets.RData"))

print("Data loading and combination complete! Ready for batch correction.")