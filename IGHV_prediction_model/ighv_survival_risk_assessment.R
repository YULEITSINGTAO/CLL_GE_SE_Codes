#------------------------------------------------------------------------------
# Title: CLL IGHV Prediction Model Validation and Survival Risk Assessment
# Author: Lei Yu
# Script Created: 2024-07-26
# Last Modified: 2025-09-04
# Description: Apply trained IGHV prediction models to determine optimal risk 
#              stratification cutoffs and evaluate survival outcomes in CLL patients
#
# Inputs:
#   - ighv_prediction_model.rds: Trained LASSO models (GE, PSI, combined)
#   - Clinical data: training_cli.rds, validation_cli.rds, testing_cli.rds
#   - Feature matrices: PSI and gene expression profiles for each cohort
#
# Outputs:
#   - Prediction tables with risk stratification for each cohort
#   - Optimal cutoff values for OS and FFS
#   - Kaplan-Meier survival plots
#   - Cutoff optimization visualization
#
# Dependencies: R >= 4.0.0, survival, survminer, glmnet, dplyr, ggplot2
#------------------------------------------------------------------------------

# Load required packages --------------------------------------------------
suppressPackageStartupMessages({
  library(glmnet)
  library(caret)
  library(dplyr)
  library(survival)
  library(survminer)
  library(plotROC)
  library(ggvenn)
  library(data.table)
  library(survcomp)
  library(ggsurvfit)
  library(mclust)
  library(classInt)
  library(tidyverse)
  library(reshape2)
})

# Configuration ------------------------------------------------------------
set.seed(123)

# Define paths
paths <- list(
  models = "results/ighv_prediction_models.rds",
  
  # Clinical data
  train_cli = "data/training_cli.rds",
  valid_cli = "data/validation_cli.rds", 
  test_cli = "data/testing_cli.rds",
  
  # Feature data
  train_psi = "data/whole_filtered_training_psi_table_var_0_05.rds",
  train_ge = "data/training_gene_count.rds",
  valid_psi = "data/whole_validation_psi_table.rds",
  valid_ge = "data/validation_gene_count.rds",
  test_psi = "data/whole_testing_psi_table.rds",
  test_ge = "data/testing_gene_count.rds",
  
  # Output directories
  predictions_dir = "results/predictions",
  cutoff_dir = "results/cutoff_optimization", 
  figures_dir = "figures/survival_analysis"
)

# Create output directories
for (dir in paths[grepl("_dir$", names(paths))]) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

# Function definitions ----------------------------------------------------

#' Apply IGHV prediction models to feature matrices
#'
#' @param models List of trained glmnet models
#' @param x_ge Gene expression matrix
#' @param x_psi Splicing event matrix  
#' @param x_combined Combined feature matrix
#' @return Data frame with prediction probabilities
#'
apply_ighv_models <- function(models, x_ge, x_psi, x_combined) {
  predictions <- data.frame(
    sample_id = rownames(x_ge),
    GE_pred = as.numeric(predict(models$gene_expression, newx = x_ge, type = "response")),
    PSI_pred = as.numeric(predict(models$splicing_events, newx = x_psi, type = "response")),
    Combined_pred = as.numeric(predict(models$combined, newx = x_combined, type = "response"))
  )
  return(predictions)
}

#' Create prediction table with clinical outcomes
#'
#' @param predictions Data frame from apply_ighv_models
#' @param clinical Clinical data frame
#' @return Combined prediction and clinical data
#'
create_prediction_table <- function(predictions, clinical) {
  
  # Validate sample alignment
  if (!identical(predictions$sample_id, rownames(clinical))) {
    stop("Sample IDs do not match between predictions and clinical data")
  }
  
  combined_table <- predictions %>%
    mutate(
      real_IGHV_binary = ifelse(clinical$IGHV_mut == "unmutated", 0, 1),
      real_IGHV = clinical$IGHV_mut,
      ffs_days = clinical$ffs_days,
      ffs_event = clinical$ffs,
      os_days = clinical$os_days,
      os_event = clinical$os
    )
  
  return(combined_table)
}

#' Calculate survival p-values for risk stratification
#'
#' @param data Data frame with predictions and survival data
#' @param cutoff_value Probability cutoff for risk stratification
#' @param cohort_name Name of the cohort (for output)
#' @return Data frame with p-values for each model and outcome
#'
calculate_survival_pvalues <- function(data, cutoff_value, cohort_name) {
  
  # Apply risk stratification
  data_with_risk <- data %>%
    mutate(
      GE_risk = ifelse(GE_pred <= cutoff_value, "High Risk", "Low Risk"),
      PSI_risk = ifelse(PSI_pred <= cutoff_value, "High Risk", "Low Risk"),
      Combined_risk = ifelse(Combined_pred <= cutoff_value, "High Risk", "Low Risk")
    )
  
  # Calculate survival p-values
  calc_pval <- function(risk_var, time_var, event_var) {
    if (length(unique(data_with_risk[[risk_var]])) == 2) {
      survdiff(Surv(data_with_risk[[time_var]], data_with_risk[[event_var]]) ~ 
                 data_with_risk[[risk_var]])$pvalue
    } else {
      NA
    }
  }
  
  # Overall survival p-values
  os_ge <- calc_pval("GE_risk", "os_days", "os_event")
  os_psi <- calc_pval("PSI_risk", "os_days", "os_event")
  os_combined <- calc_pval("Combined_risk", "os_days", "os_event")
  os_ighv <- survdiff(Surv(os_days, os_event) ~ real_IGHV, data_with_risk)$pvalue
  
  # Failure-free survival p-values
  ffs_ge <- calc_pval("GE_risk", "ffs_days", "ffs_event")
  ffs_psi <- calc_pval("PSI_risk", "ffs_days", "ffs_event")
  ffs_combined <- calc_pval("Combined_risk", "ffs_days", "ffs_event")
  ffs_ighv <- survdiff(Surv(ffs_days, ffs_event) ~ real_IGHV, data_with_risk)$pvalue
  
  # Return results
  pvalue_table <- data.frame(
    model = c("GE", "PSI", "Combined", "IGHV"),
    OS_pvalue = c(os_ge, os_psi, os_combined, os_ighv),
    FFS_pvalue = c(ffs_ge, ffs_psi, ffs_combined, ffs_ighv),
    cutoff = cutoff_value,
    cohort = cohort_name,
    stringsAsFactors = FALSE
  )
  
  return(pvalue_table)
}

#' Create Kaplan-Meier survival plot
#'
#' @param fit Survival fit object
#' @param data Data frame with survival data
#' @param title Plot title
#' @param palette Color palette
#' @return ggsurvplot object
#'
create_km_plot <- function(fit, data, title, palette = c("#98677E", "#679881")) {
  ggsurvplot(
    fit, 
    data = data,
    title = title,
    pval = TRUE,
    pval.size = 4,
    risk.table = FALSE,
    risk.table.col = "strata",
    palette = palette,
    ggtheme = theme_bw(),
    font.main = c(12, "bold"),
    font.x = c(12, "bold"),
    font.y = c(12, "bold"),
    font.tickslab = c(10, "plain"),
    legend.title = "Risk Group",
    legend = "bottom"
  )
}

# Load trained models -----------------------------------------------------
cat("Loading trained IGHV prediction models...\n")
trained_models <- readRDS(paths$models)
cat("Models loaded successfully.\n")

# Load clinical data ------------------------------------------------------
cat("Loading clinical data...\n")
train_cli <- readRDS(paths$train_cli)
valid_cli <- readRDS(paths$valid_cli) 
test_cli <- readRDS(paths$test_cli)

cat(sprintf("Clinical data loaded: %d training, %d validation, %d testing samples\n",
            nrow(train_cli), nrow(valid_cli), nrow(test_cli)))

# Load feature data -------------------------------------------------------
cat("Loading feature matrices...\n")

# Training data
train_psi <- readRDS(paths$train_psi)
train_ge <- readRDS(paths$train_ge)
train_combined <- cbind(train_psi, train_ge)

# Validation data  
valid_psi <- readRDS(paths$valid_psi)[, colnames(train_psi)]
valid_ge <- readRDS(paths$valid_ge)
valid_combined <- cbind(valid_psi, valid_ge)

# Testing data
test_psi <- readRDS(paths$test_psi)[, colnames(train_psi)]
test_ge <- readRDS(paths$test_ge)
test_combined <- cbind(test_psi, test_ge)

cat("Feature matrices loaded and aligned.\n")

# Generate predictions ----------------------------------------------------
cat("Generating IGHV predictions for all cohorts...\n")

# Apply models to each cohort
train_pred <- apply_ighv_models(trained_models, train_ge, train_psi, train_combined)
valid_pred <- apply_ighv_models(trained_models, valid_ge, valid_psi, valid_combined)
test_pred <- apply_ighv_models(trained_models, test_ge, test_psi, test_combined)

# Create prediction tables with clinical data
train_table <- create_prediction_table(train_pred, train_cli)
valid_table <- create_prediction_table(valid_pred, valid_cli)
test_table <- create_prediction_table(test_pred, test_cli)

# Save prediction tables
timestamp <- format(Sys.Date(), "%Y_%m_%d")
fwrite(train_table, file.path(paths$predictions_dir, paste0("training_predictions_", timestamp, ".csv")))
fwrite(valid_table, file.path(paths$predictions_dir, paste0("validation_predictions_", timestamp, ".csv")))
fwrite(test_table, file.path(paths$predictions_dir, paste0("testing_predictions_", timestamp, ".csv")))

cat("Prediction tables saved.\n")

# Optimize cutoff values --------------------------------------------------
cat("Optimizing risk stratification cutoffs...\n")

# Initialize results table
all_pvalues <- data.frame()

# Test cutoff values from 0.01 to 0.75 by 0.05
cutoff_values <- seq(0.01, 0.75, 0.05)

cat(sprintf("Testing %d cutoff values across cohorts...\n", length(cutoff_values)))

# Progress tracking
pb <- txtProgressBar(min = 0, max = length(cutoff_values), style = 3)

for (i in seq_along(cutoff_values)) {
  cutoff <- cutoff_values[i]
  
  # Calculate p-values for each cohort
  train_pvals <- calculate_survival_pvalues(train_table, cutoff, "training")
  valid_pvals <- calculate_survival_pvalues(valid_table, cutoff, "validation") 
  test_pvals <- calculate_survival_pvalues(test_table, cutoff, "testing")
  
  # Combine results
  all_pvalues <- bind_rows(all_pvalues, train_pvals, valid_pvals, test_pvals)
  
  setTxtProgressBar(pb, i)
}
close(pb)

# Add transformed p-values for visualization
all_pvalues <- all_pvalues %>%
  mutate(
    neg_log10_OS = -log10(OS_pvalue),
    neg_log10_FFS = -log10(FFS_pvalue)
  )

cat("Cutoff optimization complete.\n")

# Find optimal cutoffs ----------------------------------------------------
cat("Identifying optimal cutoffs...\n")

# Find cutoffs that minimize p-values in validation set
optimal_cutoffs <- all_pvalues %>%
  filter(cohort == "validation", model != "IGHV") %>%
  group_by(model) %>%
  summarise(
    optimal_OS_cutoff = cutoff[which.min(OS_pvalue)],
    optimal_OS_pvalue = min(OS_pvalue, na.rm = TRUE),
    optimal_FFS_cutoff = cutoff[which.min(FFS_pvalue)],
    optimal_FFS_pvalue = min(FFS_pvalue, na.rm = TRUE),
    .groups = "drop"
  )

print("Optimal cutoffs:")
print(optimal_cutoffs)

# Save optimal cutoffs
fwrite(optimal_cutoffs, file.path(paths$cutoff_dir, "optimal_cutoffs.csv"))

# Visualization -----------------------------------------------------------
cat("Creating cutoff optimization visualization...\n")

# Prepare data for plotting
plot_data <- all_pvalues %>%
  filter(cohort %in% c("training", "validation")) %>%
  select(model, cutoff, cohort, neg_log10_OS, neg_log10_FFS) %>%
  pivot_longer(cols = c(neg_log10_OS, neg_log10_FFS), 
               names_to = "outcome", 
               values_to = "neg_log10_pvalue") %>%
  mutate(
    outcome = ifelse(outcome == "neg_log10_OS", "Overall Survival", "Failure-Free Survival"),
    model = factor(model, levels = c("GE", "PSI", "Combined", "IGHV")),
    cohort = factor(cohort, levels = c("training", "validation"))
  )

# Create cutoff optimization plot
cutoff_plot <- ggplot(plot_data, aes(x = cutoff, y = neg_log10_pvalue, color = model)) +
  geom_line(aes(linetype = model), size = 0.8) +
  geom_point(aes(shape = model), size = 2) +
  scale_color_manual(values = c("#49B69A", "#6449B6", "#B64965", "#9BB649")) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash")) +
  scale_shape_manual(values = c(16, 17, 18, 15)) +
  facet_grid(cohort ~ outcome, scales = "free_y") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10, color = "black"),
    strip.text = element_text(size = 12, face = "bold", color = "black"),
    strip.background = element_rect(fill = "white"),
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Risk Stratification Cutoff",
    y = "-log₁₀(P-value)",
    color = "Model Type",
    linetype = "Model Type", 
    shape = "Model Type",
    title = "Cutoff Optimization for Risk Stratification",
    subtitle = "Lower cutoffs assign more samples to high-risk group"
  )

# Save cutoff optimization plot
ggsave(file.path(paths$figures_dir, "cutoff_optimization.pdf"), 
       cutoff_plot, width = 10, height = 6)

# Apply optimal cutoffs to test set ---------------------------------------
cat("Applying optimal cutoffs to test set...\n")

# Extract optimal cutoffs for easy access
os_cutoffs <- setNames(optimal_cutoffs$optimal_OS_cutoff, optimal_cutoffs$model)
ffs_cutoffs <- setNames(optimal_cutoffs$optimal_FFS_cutoff, optimal_cutoffs$model)

# Apply optimal risk stratification to test set
test_final <- test_table %>%
  mutate(
    # OS-optimized risk groups
    GE_OS_risk = ifelse(GE_pred <= os_cutoffs["GE"], "High Risk", "Low Risk"),
    PSI_OS_risk = ifelse(PSI_pred <= os_cutoffs["PSI"], "High Risk", "Low Risk"),
    Combined_OS_risk = ifelse(Combined_pred <= os_cutoffs["Combined"], "High Risk", "Low Risk"),
    
    # FFS-optimized risk groups  
    GE_FFS_risk = ifelse(GE_pred <= ffs_cutoffs["GE"], "High Risk", "Low Risk"),
    PSI_FFS_risk = ifelse(PSI_pred <= ffs_cutoffs["PSI"], "High Risk", "Low Risk"),
    Combined_FFS_risk = ifelse(Combined_pred <= ffs_cutoffs["Combined"], "High Risk", "Low Risk")
  )

# Generate survival plots -------------------------------------------------
cat("Creating Kaplan-Meier survival plots...\n")

# Create survival fits
survival_fits <- list(
  # Overall survival
  os_ge = survfit(Surv(os_days, os_event) ~ GE_OS_risk, data = test_final),
  os_psi = survfit(Surv(os_days, os_event) ~ PSI_OS_risk, data = test_final),
  os_combined = survfit(Surv(os_days, os_event) ~ Combined_OS_risk, data = test_final),
  os_ighv = survfit(Surv(os_days, os_event) ~ real_IGHV, data = test_final),
  
  # Failure-free survival
  ffs_ge = survfit(Surv(ffs_days, ffs_event) ~ GE_FFS_risk, data = test_final),
  ffs_psi = survfit(Surv(ffs_days, ffs_event) ~ PSI_FFS_risk, data = test_final),
  ffs_combined = survfit(Surv(ffs_days, ffs_event) ~ Combined_FFS_risk, data = test_final),
  ffs_ighv = survfit(Surv(ffs_days, ffs_event) ~ real_IGHV, data = test_final)
)

# Create and save survival plots
plot_configs <- list(
  list(fit = "os_ge", title = "Gene Expression Model\n(Overall Survival)", 
       palette = c("#98677E", "#679881"), file = "GE_OS_testing.pdf"),
  list(fit = "os_psi", title = "Splicing Events Model\n(Overall Survival)", 
       palette = c("#98677E", "#679881"), file = "PSI_OS_testing.pdf"),
  list(fit = "os_combined", title = "Combined Model\n(Overall Survival)", 
       palette = c("#98677E", "#679881"), file = "Combined_OS_testing.pdf"),
  list(fit = "os_ighv", title = "Actual IGHV Status\n(Overall Survival)", 
       palette = c("#AF6950", "#5096AF"), file = "IGHV_OS_testing.pdf"),
  
  list(fit = "ffs_ge", title = "Gene Expression Model\n(Failure-Free Survival)", 
       palette = c("#98677E", "#679881"), file = "GE_FFS_testing.pdf"),
  list(fit = "ffs_psi", title = "Splicing Events Model\n(Failure-Free Survival)", 
       palette = c("#98677E", "#679881"), file = "PSI_FFS_testing.pdf"),
  list(fit = "ffs_combined", title = "Combined Model\n(Failure-Free Survival)", 
       palette = c("#98677E", "#679881"), file = "Combined_FFS_testing.pdf"),
  list(fit = "ffs_ighv", title = "Actual IGHV Status\n(Failure-Free Survival)", 
       palette = c("#AF6950", "#5096AF"), file = "IGHV_FFS_testing.pdf")
)

for (config in plot_configs) {
  km_plot <- create_km_plot(
    survival_fits[[config$fit]], 
    test_final, 
    config$title, 
    config$palette
  )
  
  ggsave(file.path(paths$figures_dir, config$file), 
         km_plot$plot, width = 6, height = 5)
}

# Performance summary -----------------------------------------------------
cat("Generating performance summary...\n")

# Calculate test set performance
test_performance <- data.frame(
  Model = c("Gene Expression", "Splicing Events", "Combined", "Actual IGHV"),
  
  OS_pvalue = c(
    survdiff(Surv(os_days, os_event) ~ GE_OS_risk, test_final)$pvalue,
    survdiff(Surv(os_days, os_event) ~ PSI_OS_risk, test_final)$pvalue, 
    survdiff(Surv(os_days, os_event) ~ Combined_OS_risk, test_final)$pvalue,
    survdiff(Surv(os_days, os_event) ~ real_IGHV, test_final)$pvalue
  ),
  
  FFS_pvalue = c(
    survdiff(Surv(ffs_days, ffs_event) ~ GE_FFS_risk, test_final)$pvalue,
    survdiff(Surv(ffs_days, ffs_event) ~ PSI_FFS_risk, test_final)$pvalue,
    survdiff(Surv(ffs_days, ffs_event) ~ Combined_FFS_risk, test_final)$pvalue,
    survdiff(Surv(ffs_days, ffs_event) ~ real_IGHV, test_final)$pvalue
  )
) %>%
  mutate(
    OS_significant = ifelse(OS_pvalue < 0.05, "Yes", "No"),
    FFS_significant = ifelse(FFS_pvalue < 0.05, "Yes", "No")
  )

print("Test Set Performance Summary:")
print(test_performance)

# Save final results ------------------------------------------------------
fwrite(test_final, file.path(paths$predictions_dir, paste0("final_testing_with_risk_groups_", timestamp, ".csv")))
fwrite(test_performance, file.path(paths$cutoff_dir, "test_performance_summary.csv"))
fwrite(all_pvalues, file.path(paths$cutoff_dir, "all_cutoff_pvalues.csv"))

# Session info for reproducibility
session_info <- sessionInfo()
saveRDS(session_info, file.path(paths$cutoff_dir, "session_info.rds"))

cat("\n" + rep("=", 70) + "\n")
cat("SURVIVAL ANALYSIS COMPLETE\n")
cat(rep("=", 70) + "\n")
cat("Key Results:\n")
cat(sprintf("- Optimal cutoffs determined using validation set (n=%d)\n", nrow(valid_table)))
cat(sprintf("- Test set performance evaluated (n=%d)\n", nrow(test_final)))
cat(sprintf("- Generated %d survival plots\n", length(plot_configs)))
cat(sprintf("- Results saved to: %s\n", paths$cutoff_dir))
cat(sprintf("- Figures saved to: %s\n", paths$figures_dir))

print(session_info)