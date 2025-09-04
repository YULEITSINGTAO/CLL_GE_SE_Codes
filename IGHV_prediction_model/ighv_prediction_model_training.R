#------------------------------------------------------------------------------
# Title: CLL Gene Expression and Splicing Integration for IGHV Prediction
# Author: Lei Yu
# Email: lyi062@ucr.edu
# Script Created: 2024-07-26
# Last Modified: 2025-09-04
# Description: Train LASSO regression models to predict IGHV mutation status
#              using gene expression (GE), splicing events (PSI), and combined
#              features in chronic lymphocytic leukemia samples
#
# Inputs: 
#   - Training_PSI_Profile.rds: Splicing event profiles (PSI values)
#   - Training_Gene_Expression_Profile.rds: Gene expression profiles
#   - Clinical_Outcome.rds: Clinical metadata including IGHV mutation status
#
# Outputs:
#   - Trained LASSO models for IGHV prediction
#   - Venn diagram of selected features (Venn_plot.pdf)
#
# Dependencies: R >= 4.0.0, glmnet, caret, dplyr, ggvenn
#------------------------------------------------------------------------------

# Load required packages --------------------------------------------------
suppressPackageStartupMessages({
  library(glmnet)
  library(caret)
  library(dplyr)
  library(survival)
  library(plotROC)
  library(ggvenn)
  library(data.table)
})

# Function definitions ----------------------------------------------------

#' Train LASSO binomial model with cross-validation
#'
#' @param x_train Training feature matrix
#' @param y_train Training labels (binary factor)
#' @param alpha_value Elastic net mixing parameter (1 = LASSO, 0 = ridge)
#' @param seed Random seed for reproducibility (default: 123)
#' @return Fitted glmnet model object
#'
lasso_binomial_training <- function(x_train, y_train, alpha_value, seed = 123) {
  
  # Validate inputs
  if (!is.matrix(x_train)) {
    stop("x_train must be a matrix")
  }
  if (!is.factor(y_train)) {
    stop("y_train must be a factor")
  }
  if (length(levels(y_train)) != 2) {
    stop("y_train must be a binary factor")
  }
  
  set.seed(seed)
  
  # Cross-validation to find optimal lambda
  cv_model <- cv.glmnet(x_train, y_train, 
                        alpha = alpha_value, 
                        family = "binomial",
                        nfolds = 10,
                        type.measure = "class")
  
  # Fit final model with optimal lambda
  final_model <- glmnet(x_train, y_train, 
                        alpha = alpha_value, 
                        family = "binomial",
                        lambda = cv_model$lambda.min)
  
  # Store CV results in model object for later reference
  final_model$cv_results <- cv_model
  
  return(final_model)
}

#' Extract selected features from LASSO model
#'
#' @param model Fitted glmnet model object
#' @return Character vector of selected feature names
#'
extract_selected_features <- function(model) {
  coef_matrix <- coef(model, s = "lambda.min")
  selected_indices <- which(coef_matrix[, 1] != 0)
  feature_names <- rownames(coef_matrix)[selected_indices]
  # Remove intercept term
  feature_names[feature_names != "(Intercept)"]
}

# Data loading and preprocessing ------------------------------------------

cat("Loading training data...\n")

# Load feature matrices
x_psi <- readRDS("data/Training_PSI_Profile.rds")
x_ge <- readRDS("data/Training_Gene_Expression_Profile.rds")

# Validate data dimensions
cat(sprintf("PSI features: %d samples × %d features\n", 
            nrow(x_psi), ncol(x_psi)))
cat(sprintf("Gene expression features: %d samples × %d features\n", 
            nrow(x_ge), ncol(x_ge)))

# Check sample alignment
if (!identical(rownames(x_psi), rownames(x_ge))) {
  stop("Sample names do not match between PSI and gene expression data")
}

# Combine feature matrices
x_combined <- cbind(x_psi, x_ge)
cat(sprintf("Combined features: %d samples × %d features\n", 
            nrow(x_combined), ncol(x_combined)))

# Load clinical data
training_cli <- readRDS("data/Clinical_Outcome.rds")

# Validate clinical data alignment
if (!identical(rownames(x_psi), rownames(training_cli))) {
  stop("Sample names do not match between features and clinical data")
}

# Prepare response variable
y_raw <- training_cli$IGHV_mut
y_binary <- ifelse(y_raw == "mutated", 1, 0)
y <- as.factor(y_binary)

# Data summary
cat(sprintf("IGHV mutation status: %d mutated, %d unmutated\n", 
            sum(y == 1), sum(y == 0)))

# Model training ----------------------------------------------------------

cat("Training LASSO models...\n")

# Train individual models
cat("  - Gene expression model...\n")
ge_lasso_model <- lasso_binomial_training(x_ge, y, alpha_value = 1)

cat("  - Splicing event model...\n")  
psi_lasso_model <- lasso_binomial_training(x_psi, y, alpha_value = 1)

cat("  - Combined model...\n")
combined_lasso_model <- lasso_binomial_training(x_combined, y, alpha_value = 1)

# Store models in list
ighv_prediction_models <- list(
  gene_expression = ge_lasso_model,
  splicing_events = psi_lasso_model,
  combined = combined_lasso_model
)

cat("Model training complete.\n")

# Feature selection analysis ----------------------------------------------

cat("Analyzing selected features...\n")

# Extract selected features for each model
selected_features <- list(
  GE = extract_selected_features(ge_lasso_model),
  PSI = extract_selected_features(psi_lasso_model),
  Combined = extract_selected_features(combined_lasso_model)
)

# Print feature counts
feature_counts <- sapply(selected_features, length)
cat("Selected features per model:\n")
for (i in seq_along(feature_counts)) {
  cat(sprintf("  %s: %d features\n", names(feature_counts)[i], feature_counts[i]))
}

# Visualization -----------------------------------------------------------

cat("Creating feature selection Venn diagram...\n")

# Ensure output directory exists
output_dir <- "figures"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Create Venn diagram
pdf(file = file.path(output_dir, "IGHV_feature_selection_venn.pdf"),
    width = 6, height = 5)

venn_plot <- ggvenn(
  data = list(
    "Gene Expression" = selected_features$GE,
    "Splicing Events" = selected_features$PSI,
    "Combined Model" = selected_features$Combined
  ),
  fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF"),
  stroke_size = 1, 
  set_name_size = 4,
  text_size = 3
) +
  labs(title = "Selected Features Across IGHV Prediction Models",
       subtitle = paste("Total features: GE =", length(selected_features$GE),
                        "| PSI =", length(selected_features$PSI), 
                        "| Combined =", length(selected_features$Combined)))

print(venn_plot)
dev.off()

cat("Venn diagram saved to:", file.path(output_dir, "IGHV_feature_selection_venn.pdf"), "\n")

# Model performance summary -----------------------------------------------

cat("\nModel Performance Summary:\n")
for (model_name in names(ighv_prediction_models)) {
  model <- ighv_prediction_models[[model_name]]
  cv_results <- model$cv_results
  
  cat(sprintf("  %s model:\n", tools::toTitleCase(gsub("_", " ", model_name))))
  cat(sprintf("    - Optimal lambda: %.4f\n", cv_results$lambda.min))
  cat(sprintf("    - CV error: %.4f ± %.4f\n", 
              min(cv_results$cvm), 
              cv_results$cvsd[which.min(cv_results$cvm)]))
  cat(sprintf("    - Selected features: %d\n", 
              length(extract_selected_features(model))))
}

# Save results ------------------------------------------------------------

cat("\nSaving results...\n")

# Save models
saveRDS(ighv_prediction_models, file = "results/ighv_prediction_models.rds")

# Save feature lists
saveRDS(selected_features, file = "results/selected_features.rds")

# Create session info for reproducibility
session_info <- sessionInfo()
saveRDS(session_info, file = "results/session_info.rds")

cat("Analysis complete! Results saved to 'results/' directory.\n")

# Print session info for manuscript methods section
cat("\n" + rep("=", 60) + "\n")
cat("SESSION INFO FOR METHODS SECTION:\n")
cat(rep("=", 60) + "\n")
print(session_info)