# Predict heart conditions using machine learning - BEGINNER VERSION
# This script uses gene expression to predict atrial fibrillation and heart failure

library(randomForest)  # Package for random forest machine learning
library(e1071)         # Package for support vector machines
library(caret)         # Package for machine learning tools
library(pROC)          # Package for ROC curves and performance metrics
library(ggplot2)       # Package for plotting

# Set up folders
analysis_results_folder <- "/Users/macbookpro/Desktop/Ananylum/new_meta_analysis/results"
plots_folder <- file.path(analysis_results_folder, "condition_specific_analysis")

print("Loading batch-corrected heart disease data and gene analysis results...")
load(file.path(analysis_results_folder, "02_condition_specific_corrected.RData"))
load(file.path(analysis_results_folder, "03_differential_gene_analysis_results.RData"))

# CAMK genes we're especially interested in
important_calcium_genes <- c("CAMK1", "CAMK1D", "CAMK1G", "CAMK2A", "CAMK2B", "CAMK2D", 
                            "CAMK2G", "CAMK4", "CAMKK1", "CAMKK2", "PNCK")

print(paste("Number of AF samples:", ncol(af_combat_expr)))
print(paste("Number of HF samples:", ncol(hf_combat_expr)))

# Function to prepare data and train machine learning models
train_prediction_models_for_condition <- function(gene_expression_data, sample_information, condition_name) {
  
  print(paste("Training prediction models for", condition_name, "condition..."))
  
  # Step 1: Clean up the data
  samples_with_complete_data <- complete.cases(t(gene_expression_data))
  clean_gene_expression <- gene_expression_data[, samples_with_complete_data]
  clean_sample_info <- sample_information[samples_with_complete_data, ]
  
  # Step 2: Remove genes that are missing in too many samples
  genes_present_in_most_samples <- rowSums(!is.na(clean_gene_expression)) / ncol(clean_gene_expression)
  genes_to_keep <- genes_present_in_most_samples > 0.8
  final_gene_expression <- clean_gene_expression[genes_to_keep, ]
  
  print(paste("Using", nrow(final_gene_expression), "genes and", ncol(final_gene_expression), "samples"))
  
  # Step 3: Create labels for machine learning (what we want to predict)
  if(condition_name == "AF") {
    condition_labels <- clean_sample_info$group
    positive_class <- "AF"
    negative_class <- "Control"
  } else if(condition_name == "HF") {
    condition_labels <- clean_sample_info$group
    positive_class <- "HF"
    negative_class <- "Control"
  }
  
  print(paste("Samples with", positive_class, ":", sum(condition_labels == positive_class)))
  print(paste("Samples with", negative_class, ":", sum(condition_labels == negative_class)))
  
  # Step 4: Find the most important genes (top 100) WITHOUT using test data
  
  # This is a simplified feature selection - in practice, this should be done inside cross-validation
  # For beginners, we'll use a simpler approach but note this limitation
  print("Selecting most variable genes as features...")
  
  gene_variance <- apply(final_gene_expression, 1, var, na.rm = TRUE)
  most_variable_genes <- names(sort(gene_variance, decreasing = TRUE))[1:100]
  selected_gene_expression <- final_gene_expression[most_variable_genes, ]
  
  print(paste("Selected", length(most_variable_genes), "most variable genes for prediction"))
  
  # Check if any CAMK genes are in our selected features
  camk_genes_selected <- most_variable_genes[most_variable_genes %in% important_calcium_genes]
  if(length(camk_genes_selected) > 0) {
    print(paste("CAMK genes selected as features:", paste(camk_genes_selected, collapse = ", ")))
  } else {
    print("No CAMK genes were among the most variable genes")
  }
  
  # Step 5: Prepare data for machine learning
  # Transpose so samples are rows and genes are columns
  ml_feature_matrix <- t(selected_gene_expression)
  ml_labels <- factor(condition_labels, levels = c(negative_class, positive_class))
  
  # Step 6: Split data into training and testing sets
  set.seed(42)  # For reproducible results
  total_samples <- length(ml_labels)
  training_sample_indices <- sample(1:total_samples, size = round(0.7 * total_samples))
  testing_sample_indices <- setdiff(1:total_samples, training_sample_indices)
  
  # Training data
  training_features <- ml_feature_matrix[training_sample_indices, ]
  training_labels <- ml_labels[training_sample_indices]
  
  # Testing data
  testing_features <- ml_feature_matrix[testing_sample_indices, ]
  testing_labels <- ml_labels[testing_sample_indices]
  
  print(paste("Training samples:", length(training_labels)))
  print(paste("Testing samples:", length(testing_labels)))
  
  # Step 7: Train different machine learning models
  
  print("Training Random Forest model...")
  random_forest_model <- randomForest(x = training_features, y = training_labels, ntree = 500, importance = TRUE)
  
  print("Training Support Vector Machine model...")
  svm_model <- svm(x = training_features, y = training_labels, probability = TRUE, kernel = "radial")
  
  print("Training Logistic Regression model...")
  # Convert to data frame for glm
  training_data_df <- data.frame(training_features)
  training_data_df$condition <- training_labels
  logistic_model <- glm(condition ~ ., data = training_data_df, family = binomial)
  
  # Step 8: Make predictions on test set
  
  print("Making predictions on test data...")
  
  # Random Forest predictions
  rf_predictions <- predict(random_forest_model, testing_features, type = "prob")[, positive_class]
  rf_class_predictions <- predict(random_forest_model, testing_features)
  
  # SVM predictions
  svm_predictions <- predict(svm_model, testing_features, probability = TRUE)
  svm_probabilities <- attr(svm_predictions, "probabilities")[, positive_class]
  
  # Logistic Regression predictions
  testing_data_df <- data.frame(testing_features)
  lr_probabilities <- predict(logistic_model, testing_data_df, type = "response")
  lr_class_predictions <- ifelse(lr_probabilities > 0.5, positive_class, negative_class)
  lr_class_predictions <- factor(lr_class_predictions, levels = c(negative_class, positive_class))
  
  # Step 9: Calculate performance metrics
  
  print("Calculating model performance...")
  
  # Random Forest performance
  rf_confusion <- confusionMatrix(rf_class_predictions, testing_labels, positive = positive_class)
  rf_roc <- roc(testing_labels, rf_predictions)
  rf_auc <- auc(rf_roc)
  
  # SVM performance  
  svm_class_predictions <- factor(svm_predictions, levels = c(negative_class, positive_class))
  svm_confusion <- confusionMatrix(svm_class_predictions, testing_labels, positive = positive_class)
  svm_roc <- roc(testing_labels, svm_probabilities)
  svm_auc <- auc(svm_roc)
  
  # Logistic Regression performance
  lr_confusion <- confusionMatrix(lr_class_predictions, testing_labels, positive = positive_class)
  lr_roc <- roc(testing_labels, lr_probabilities)
  lr_auc <- auc(lr_roc)
  
  # Step 10: Get feature importance from Random Forest
  feature_importance <- importance(random_forest_model)
  important_features <- rownames(feature_importance)[order(feature_importance[, "MeanDecreaseGini"], decreasing = TRUE)][1:20]
  
  # Check for CAMK genes in top features
  top_camk_features <- important_features[important_features %in% important_calcium_genes]
  
  print("Performance Results:")
  print(paste("Random Forest AUC:", round(rf_auc, 3)))
  print(paste("SVM AUC:", round(svm_auc, 3)))
  print(paste("Logistic Regression AUC:", round(lr_auc, 3)))
  
  if(length(top_camk_features) > 0) {
    print(paste("CAMK genes in top 20 features:", paste(top_camk_features, collapse = ", ")))
  }
  
  # Return comprehensive results
  return(list(
    condition = condition_name,
    sample_info = list(
      total_samples = ncol(final_gene_expression),
      training_samples = length(training_labels),
      testing_samples = length(testing_labels),
      positive_samples = sum(condition_labels == positive_class),
      negative_samples = sum(condition_labels == negative_class)
    ),
    feature_info = list(
      total_genes_used = nrow(final_gene_expression),
      selected_features = length(most_variable_genes),
      camk_genes_selected = camk_genes_selected,
      top_20_features = important_features,
      top_camk_features = top_camk_features
    ),
    model_performance = data.frame(
      model = c("Random_Forest", "SVM", "Logistic_Regression"),
      auc = c(rf_auc, svm_auc, lr_auc),
      accuracy = c(rf_confusion$overall["Accuracy"], svm_confusion$overall["Accuracy"], lr_confusion$overall["Accuracy"]),
      sensitivity = c(rf_confusion$byClass["Sensitivity"], svm_confusion$byClass["Sensitivity"], lr_confusion$byClass["Sensitivity"]),
      specificity = c(rf_confusion$byClass["Specificity"], svm_confusion$byClass["Specificity"], lr_confusion$byClass["Specificity"]),
      stringsAsFactors = FALSE
    ),
    models = list(
      random_forest = random_forest_model,
      svm = svm_model,
      logistic_regression = logistic_model
    ),
    test_data = list(
      features = testing_features,
      labels = testing_labels,
      rf_probabilities = rf_predictions,
      svm_probabilities = svm_probabilities,
      lr_probabilities = lr_probabilities
    )
  ))
}

# Train models for both conditions
print("=== TRAINING MACHINE LEARNING MODELS ===")
print("")

af_ml_results <- train_prediction_models_for_condition(af_combat_expr, af_metadata, "AF")
hf_ml_results <- train_prediction_models_for_condition(hf_combat_expr, hf_metadata, "HF")

print("Creating performance comparison plots...")

# Combine performance results
all_performance_results <- rbind(
  data.frame(condition = "AF", af_ml_results$model_performance),
  data.frame(condition = "HF", hf_ml_results$model_performance)
)

# Create AUC comparison plot
auc_plot <- ggplot(all_performance_results, aes(x = model, y = auc, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("AF" = "#FF6B6B", "HF" = "#4ECDC4")) +
  labs(title = "Machine Learning Model Performance Comparison",
       x = "Machine Learning Model",
       y = "AUC (Area Under Curve)",
       fill = "Heart Condition") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Create accuracy comparison plot
accuracy_plot <- ggplot(all_performance_results, aes(x = model, y = accuracy, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("AF" = "#FF6B6B", "HF" = "#4ECDC4")) +
  labs(title = "Machine Learning Model Accuracy Comparison",
       x = "Machine Learning Model",
       y = "Accuracy",
       fill = "Heart Condition") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save plots
pdf(file.path(plots_folder, "04_machine_learning_performance_comparison.pdf"), width = 12, height = 8)
print(auc_plot)
print(accuracy_plot)
dev.off()

print("Creating CAMK gene biomarker summary...")

# Summarize CAMK gene involvement
camk_biomarker_summary <- data.frame(
  condition = c("AF", "HF"),
  camk_genes_selected_as_features = c(
    paste(af_ml_results$feature_info$camk_genes_selected, collapse = "; "),
    paste(hf_ml_results$feature_info$camk_genes_selected, collapse = "; ")
  ),
  camk_genes_in_top_20_features = c(
    paste(af_ml_results$feature_info$top_camk_features, collapse = "; "),
    paste(hf_ml_results$feature_info$top_camk_features, collapse = "; ")
  ),
  best_model_auc = c(
    max(af_ml_results$model_performance$auc),
    max(hf_ml_results$model_performance$auc)
  ),
  stringsAsFactors = FALSE
)

print("Saving machine learning results...")

# Save performance summary
write.csv(all_performance_results, 
          file.path(analysis_results_folder, "04_ml_model_performance_summary.csv"), 
          row.names = FALSE)

# Save CAMK biomarker summary
write.csv(camk_biomarker_summary, 
          file.path(analysis_results_folder, "04_camk_biomarker_potential_summary.csv"), 
          row.names = FALSE)

# Save feature importance for top models
af_top_features <- data.frame(
  gene = af_ml_results$feature_info$top_20_features,
  rank = 1:20,
  condition = "AF",
  is_camk_gene = af_ml_results$feature_info$top_20_features %in% important_calcium_genes,
  stringsAsFactors = FALSE
)

hf_top_features <- data.frame(
  gene = hf_ml_results$feature_info$top_20_features,
  rank = 1:20,
  condition = "HF", 
  is_camk_gene = hf_ml_results$feature_info$top_20_features %in% important_calcium_genes,
  stringsAsFactors = FALSE
)

all_top_features <- rbind(af_top_features, hf_top_features)
write.csv(all_top_features, 
          file.path(analysis_results_folder, "04_top_predictive_genes.csv"), 
          row.names = FALSE)

# Save complete results for further analysis
save(af_ml_results, hf_ml_results, all_performance_results, camk_biomarker_summary,
     file = file.path(analysis_results_folder, "04_complete_machine_learning_results.RData"))

print("Machine learning analysis complete!")
print("")
print("SUMMARY OF MACHINE LEARNING RESULTS:")
print("")
print("Model Performance:")
for(condition in c("AF", "HF")) {
  condition_results <- all_performance_results[all_performance_results$condition == condition, ]
  print(paste("", condition, "condition:"))
  for(i in 1:nrow(condition_results)) {
    model_result <- condition_results[i, ]
    print(paste("  ", model_result$model, ": AUC =", round(model_result$auc, 3), 
                ", Accuracy =", round(model_result$accuracy, 3)))
  }
  print("")
}

print("CAMK Gene Biomarker Potential:")
for(i in 1:nrow(camk_biomarker_summary)) {
  condition_summary <- camk_biomarker_summary[i, ]
  print(paste("", condition_summary$condition, ":"))
  if(condition_summary$camk_genes_selected_as_features != "") {
    print(paste("  Selected as features:", condition_summary$camk_genes_selected_as_features))
  } else {
    print("  No CAMK genes selected as features")
  }
  if(condition_summary$camk_genes_in_top_20_features != "") {
    print(paste("  In top 20 features:", condition_summary$camk_genes_in_top_20_features))
  } else {
    print("  No CAMK genes in top 20 features")
  }
  print(paste("  Best model AUC:", round(condition_summary$best_model_auc, 3)))
  print("")
}

print("Files created:")
print("- Model performance comparison plots")
print("- Performance metrics for all models")
print("- CAMK gene biomarker potential summary") 
print("- Top predictive genes for each condition")
print("- Complete results saved for pathway analysis")