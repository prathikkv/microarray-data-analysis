# Load and combine cardiac datasets - BEGINNER VERSION
# This script converts microarray probe data to gene data and combines all datasets

library(GEOquery)       # Package to work with GEO datasets
library(hgu133plus2.db) # Database to convert probes to gene names

# Set up our folder paths
main_data_folder <- "/Users/macbookpro/Desktop/Ananylum/data"
analysis_results_folder <- "/Users/macbookpro/Desktop/Ananylum/new_meta_analysis/results"

# Make sure results folder exists
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

print("Converting Dataset 1 (Atrial Fibrillation) probes to genes...")
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

print("Converting Dataset 2 (Atrial Fibrillation) probes to genes...")
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

print("Converting Dataset 3 (Heart Failure) probes to genes...")
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

print("Converting Dataset 4 (Heart Failure) - uses different chip type...")
# Fourth dataset uses different microarray chip, need different conversion method
different_chip_annotation <- getGEO("GPL11532", destdir = tempdir())
different_chip_lookup_table <- Table(different_chip_annotation)

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

# Find which probes from dataset 4 we can convert
dataset_4_probe_ids <- rownames(heart_fail_2_expression_data)
convertible_probes <- intersect(dataset_4_probe_ids, different_chip_probe_to_gene$probe)
probe_gene_mapping_dataset_4 <- different_chip_probe_to_gene[match(convertible_probes, different_chip_probe_to_gene$probe), ]
heart_fail_2_expression_subset <- heart_fail_2_expression_data[convertible_probes, ]

# Convert dataset 4 probes to genes using efficient method
heart_fail_2_gene_data_list <- tapply(1:nrow(heart_fail_2_expression_subset), probe_gene_mapping_dataset_4$gene, function(probe_indices) {
  if(length(probe_indices) == 1) {
    heart_fail_2_expression_subset[probe_indices, ]
  } else {
    apply(heart_fail_2_expression_subset[probe_indices, , drop = FALSE], 2, median)
  }
})
heart_fail_2_gene_data <- do.call(rbind, heart_fail_2_gene_data_list)

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

print("Final dataset summary:")
print(paste("- Total samples:", ncol(combined_expression_matrix)))
print(paste("- Total genes:", nrow(combined_expression_matrix)))
print(paste("- AF samples:", sum(combined_metadata$study_type == "AF")))
print(paste("- HF samples:", sum(combined_metadata$study_type == "HF")))
print("")

print("Saving combined dataset...")

# Save the final combined dataset for next analysis step (using original filename)
# Change variable name for compatibility with downstream scripts
combined_expr <- combined_expression_matrix
save(combined_expr, combined_metadata, 
     file = file.path(analysis_results_folder, "01_combined_4datasets.RData"))

print("Data loading and combination complete! Ready for batch correction.")