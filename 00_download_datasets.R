# Download cardiac datasets from GEO database - BEGINNER VERSION
# This script downloads 4 heart disease datasets for analysis

library(GEOquery)  # Package to download GEO datasets
library(Biobase)   # Package for gene expression data handling

# Create folders to store our data
main_data_folder <- "/Users/macbookpro/Desktop/Ananylum/data"
analysis_results_folder <- "/Users/macbookpro/Desktop/Ananylum/new_meta_analysis/results"

# Make sure these folders exist (create them if they don't)
dir.create(main_data_folder, showWarnings = FALSE, recursive = TRUE)
dir.create(analysis_results_folder, showWarnings = FALSE, recursive = TRUE)

print("Starting to download cardiac datasets...")

# Download first dataset: Atrial Fibrillation study (GSE115574)
print("Downloading atrial fibrillation dataset 1...")
atrial_fibrillation_dataset_1 <- getGEO("GSE115574", destdir = main_data_folder, GSEMatrix = TRUE)
if(is.list(atrial_fibrillation_dataset_1)) {
  atrial_fibrillation_dataset_1 <- atrial_fibrillation_dataset_1[[1]]
}

# Download second dataset: Atrial Fibrillation study (GSE79768)  
print("Downloading atrial fibrillation dataset 2...")
atrial_fibrillation_dataset_2 <- getGEO("GSE79768", destdir = main_data_folder, GSEMatrix = TRUE)
if(is.list(atrial_fibrillation_dataset_2)) {
  atrial_fibrillation_dataset_2 <- atrial_fibrillation_dataset_2[[1]]
}

# Download third dataset: Heart Failure study (GSE76701)
print("Downloading heart failure dataset 1...")
heart_failure_dataset_1 <- getGEO("GSE76701", destdir = main_data_folder, GSEMatrix = TRUE)
if(is.list(heart_failure_dataset_1)) {
  heart_failure_dataset_1 <- heart_failure_dataset_1[[1]]
}

# Download fourth dataset: Heart Failure study (GSE57338)
print("Downloading heart failure dataset 2...")
heart_failure_dataset_2 <- getGEO("GSE57338", destdir = main_data_folder, GSEMatrix = TRUE)
if(is.list(heart_failure_dataset_2)) {
  heart_failure_dataset_2 <- heart_failure_dataset_2[[1]]
}

print("All cardiac datasets downloaded successfully!")

# Save all downloaded datasets in one file for next step
save(atrial_fibrillation_dataset_1, atrial_fibrillation_dataset_2, 
     heart_failure_dataset_1, heart_failure_dataset_2,
     file = file.path(analysis_results_folder, "00_downloaded_cardiac_datasets.RData"))

print("Download complete! Datasets saved and ready for analysis.")