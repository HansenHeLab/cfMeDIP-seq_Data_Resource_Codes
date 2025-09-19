###############################################################################
# File:      10.1_Process_Yuanchang_Data_Cancer_Type_Classification
# Author:    Dory Abelman
# Date:      2025‑03‑30
#
# Purpose:
#   - Load and merge per‑cancer classification summaries (.sum.csv) and ROC data (.roc.csv)
#     from cfMeDIP-based ML results.
#   - Clean and harmonize cancer labels (e.g., rename “Eye Cancer” to “Uveal Melanoma”;
#     exclude AML due to low sample size).
#   - Compute and export:
#       • Best-performing fragmentomic feature per cancer type (mean ± SD AUC).
#       • Best-performing ML model per cancer type (mean ± SD AUC).
#       • Fold‑level AUC summaries for each cancer/feature/model/iteration.
#       • Overall “best model” selection and its detailed per‑run metrics.
#   - Merge ROC curves with summary run information.
#   - Aggregate and plot:
#       • Mean ROC curves with confidence ribbons, faceted by cancer.
#       • Confusion‑matrix–derived sensitivity at 95% specificity.
#       • Cleveland dotplots of sensitivity and AUC by cancer and by feature/model.
#   - Compare methylation‑only vs. combined‑feature classifiers (AUC and sensitivity gains).
#   - Generate result tables (CSV) and figures (PNG/PDF) for downstream reporting.
###############################################################################



# Load required libraries
library(dplyr)
library(ggplot2)
library(readr)
library(data.table)
library(ggpubr)
library(tidyr)


# Set the path to the directory containing the CSV files
path <- "/Users/dabelman/Documents/Thesis_work/R/M4/Projects/Fragmentation_CfMeDIP_Yong/Final Data Dec 2024/Yuanchang_ML_results/cfMEDIP_clf_res_20250330/results_cancers"

# List all .sum.csv and .roc.csv files in the directory
sum_files <- list.files(path, pattern = "\\.sum\\.csv$")
roc_files <- list.files(path, pattern = "\\.roc\\.csv$")

# ---------------------------
# Load and combine summary files
# ---------------------------
sum_data_list <- lapply(sum_files, function(file) {
  # Use file.path to ensure the correct path is used when reading each file
  data <- read_csv(file.path(path, file), col_types = cols())
  # Record the file source (optional)
  data$source_file <- file
  return(data)
})
sum_data <- bind_rows(sum_data_list)

# ---------------------------
# Load and combine ROC files
# ---------------------------
roc_data_list <- lapply(roc_files, function(file) {
  data <- read_csv(file.path(path, file), col_types = cols())
  data$source_file <- file
  return(data)
})
roc_data <- bind_rows(roc_data_list)

# Remove AML after combining the summary files since only have a few samples:
sum_data <- bind_rows(sum_data_list) %>%
  filter(cancer != "AML")

sum_data <- sum_data %>%
  mutate(cancer = if_else(cancer == "Eye Cancer", "Uveal Melanoma", cancer))

# And similarly for roc_data:
roc_data <- bind_rows(roc_data_list) %>%
  filter(!grepl("^AML", run))

# Preview the loaded data (optional)
print("Summary data head:")
print(head(sum_data))
print("ROC data head:")
print(head(roc_data))

# ---------------------------
# Summarize best fragmentomic feature per cancer type
# (Aggregating test_auc across folds/iterations)
# ---------------------------
best_feature_summary <- sum_data %>%
  group_by(cancer, feature) %>%
  summarise(
    mean_auc = mean(test_auc, na.rm = TRUE),
    sd_auc   = sd(test_auc, na.rm = TRUE),
    min_auc  = min(test_auc, na.rm = TRUE),
    max_auc  = max(test_auc, na.rm = TRUE),
    n        = n(),
    .groups  = "drop"
  ) %>%
  arrange(cancer, desc(mean_auc))
print("Best Feature Summary:")
print(best_feature_summary)

# ---------------------------
# Summarize best models per cancer type
# ---------------------------
best_model_summary <- sum_data %>%
  group_by(cancer, model) %>%
  summarise(
    mean_auc = mean(test_auc, na.rm = TRUE),
    sd_auc   = sd(test_auc, na.rm = TRUE),
    min_auc  = min(test_auc, na.rm = TRUE),
    max_auc  = max(test_auc, na.rm = TRUE),
    n        = n(),
    .groups  = "drop"
  ) %>%
  arrange(cancer, desc(mean_auc))
print("Best Model Summary:")
print(best_model_summary)

## see nuber of iterations run 
# Example: results_df has columns Feature, cancer_type, Model, iteration
iterations_summary <- sum_data %>%
  group_by(feature, cancer, model) %>%
  summarise(
    n_iter = n(),                         # total number of iterations
    mean_iter = mean(n_iter),             # mean (trivial here since n_iter is single)
    median_iter = median(n_iter),
    min_iter = min(n_iter),
    max_iter = max(n_iter),
    range_iter = max(n_iter) - min(n_iter),
    .groups = "drop"
  )

print(iterations_summary)

# ---------------------------
# Summarize fold-level (iteration) data
# ---------------------------
fold_summary <- sum_data %>%
  group_by(cancer, feature, model, iteration) %>%
  summarise(
    mean_fold_auc = mean(test_auc, na.rm = TRUE),
    min_fold_auc  = min(test_auc, na.rm = TRUE),
    max_fold_auc  = max(test_auc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(cancer, feature, model, iteration)
print("Fold Summary:")
print(fold_summary)


# Identify the best model for each cancer type based on average test AUC
best_model_by_cancer <- sum_data %>%
  group_by(cancer, model) %>%
  summarise(mean_auc = mean(test_auc, na.rm = TRUE), .groups = "drop") %>%
  group_by(cancer) %>%
  filter(mean_auc == max(mean_auc, na.rm = TRUE)) %>%
  ungroup()

# If more than one model ties for best in a cancer type, you may take one (e.g., the first)
best_model_by_cancer <- best_model_by_cancer %>%
  group_by(cancer) %>%
  slice(1) %>%
  ungroup()

# Filter sum_data to keep only rows for the best model per cancer type
best_model_data <- inner_join(sum_data, best_model_by_cancer, by = c("cancer", "model"))

# Optionally, remove the extra column (mean_auc) if you don't need it
best_model_data <- best_model_data %>% select(-mean_auc)


# Optionally write the summary tables to CSV files
write_csv(best_feature_summary, "best_feature_summary.csv")
write_csv(best_model_summary, "best_model_summary.csv")
write_csv(fold_summary, "fold_summary.csv")
write_csv(best_model_by_cancer, "best_model_by_cancer.csv")



# ---------------------------
# Merge ROC data with summary/model information
# (Assumes that the 'run' column uniquely identifies the model instance)
# ---------------------------
# Extract unique run info from summary data
run_info <- sum_data %>%
  select(run, cancer, feature, model) %>%
  distinct()

# Merge with ROC data
roc_data_joined <- merge(roc_data, run_info, by = "run", all.x = TRUE)


#### Testing - skip ahead
# ---------------------------
# Plot ROC curves by cancer type
# ---------------------------
roc_plot <- ggplot(roc_data_joined, aes(x = fpr, y = tpr, group = run, color = model)) +
  geom_line(alpha = 0.7) +
  facet_wrap(~ cancer, scales = "free") +
  theme_minimal() +
  labs(title = "ROC Curves by Cancer Type",
       x = "False Positive Rate",
       y = "True Positive Rate")

# Save the plot to a file
ggsave("ROC_curves_by_cancer.png", roc_plot, width = 10, height = 8, dpi = 300)

# Print the plot
print(roc_plot)

# Extract the best run for each combination of cancer type and feature based on test_auc
best_run_info <- sum_data %>%
  group_by(cancer, feature) %>%
  filter(test_auc == max(test_auc, na.rm = TRUE)) %>%
  distinct(run, .keep_all = TRUE) %>%
  ungroup() %>%
  select(run, cancer, feature)

# Merge the best run info with the ROC data
roc_data_best <- merge(roc_data, best_run_info, by = "run", all.x = FALSE)

# Plot ROC curves by cancer type with each color representing a different feature
roc_plot <- ggplot(roc_data_best, aes(x = fpr, y = tpr, group = run, color = feature)) +
  geom_line(alpha = 0.7) +
  facet_wrap(~ cancer, scales = "free") +
  theme_minimal() +
  labs(title = "ROC Curves by Cancer Type (Best Model per Feature)",
       x = "False Positive Rate",
       y = "True Positive Rate")

# Save the plot to a file
ggsave("ROC_curves_by_cancer_best_model.png", roc_plot, width = 10, height = 8, dpi = 300)

# Print the plot
print(roc_plot)





#### Updated 
## Get mean ROC
# Merge run info from the summary with the ROC data.
roc_data_joined <- merge(roc_data, best_model_data %>% select(run, cancer, feature), by = "run", all.x = FALSE)

# Initialize list to store aggregated ROC curves for each cancer-feature combination

# Convert to data.table for efficient processing
roc_dt <- as.data.table(roc_data_joined)

# Ensure fpr and tpr are numeric
roc_dt[, fpr := as.numeric(fpr)]
roc_dt[, tpr := as.numeric(tpr)]

# List to store results
roc_avg_list <- list()

# Get unique combinations of cancer and feature
unique_combos <- unique(roc_dt[, .(cancer, feature)])

# Optionally, reduce grid resolution (e.g., 50 points instead of 100)
grid_points <- 50
fpr_grid <- seq(0, 1, length.out = grid_points)

# Loop over each combination
for (i in seq_len(nrow(unique_combos))) {
  c_val <- unique_combos$cancer[i]
  f_val <- unique_combos$feature[i]
  
  # Subset for the current cancer-feature combination
  subset_dt <- roc_dt[cancer == c_val & feature == f_val & !is.na(fpr) & !is.na(tpr)]
  
  # Get unique run identifiers
  runs <- unique(subset_dt$run)
  
  # Initialize list to store interpolated TPR for each run
  tpr_list <- vector("list", length(runs))
  
  for (j in seq_along(runs)) {
    run_data <- subset_dt[run == runs[j]]
    # Pre-aggregate duplicate FPR values by averaging TPR
    run_data <- run_data[, .(tpr = mean(tpr)), by = fpr]
    # Order by fpr
    setorder(run_data, fpr)
    # Interpolate TPR values on the common grid
    interp_vals <- approx(x = run_data$fpr, y = run_data$tpr, xout = fpr_grid)$y
    tpr_list[[j]] <- interp_vals
  }
  
  # Combine the interpolated TPR values from all runs
  tpr_matrix <- do.call(cbind, tpr_list)
  mean_tpr <- rowMeans(tpr_matrix, na.rm = TRUE)
  sd_tpr <- apply(tpr_matrix, 1, sd, na.rm = TRUE)
  
  # Store results in a data.table
  roc_avg_list[[paste(c_val, f_val, sep = "_")]] <- data.table(
    cancer = c_val,
    feature = f_val,
    fpr = fpr_grid,
    tpr = mean_tpr,
    tpr_lower = pmax(mean_tpr - sd_tpr, 0),
    tpr_upper = pmin(mean_tpr + sd_tpr, 1)
  )
  
  # Clean up temporary objects to free memory
  rm(subset_dt, tpr_list, tpr_matrix)
  gc()
}

# Combine all results
roc_avg_dt <- rbindlist(roc_avg_list)

# Plot using ggplot2
roc_plot <- ggplot(roc_avg_dt, aes(x = fpr, y = tpr, color = feature)) +
  geom_line(size = 1, alpha = 0.8) +
  geom_ribbon(aes(ymin = tpr_lower, ymax = tpr_upper, fill = feature), alpha = 0.2, color = NA) +
  facet_wrap(~ cancer, scales = "free") +
  theme_minimal() +
  labs(title = "Mean ROC Curves by Cancer Type (Averaged Across Folds)",
       x = "False Positive Rate",
       y = "True Positive Rate",
       color = "Feature",
       fill = "Feature")

# Save and display the plot
ggsave("Mean_ROC_curves_by_cancer_with_CI_best_model.png", roc_plot, width = 10, height = 8, dpi = 300)
print(roc_plot)



### Update style 
library(dplyr)
library(ggplot2)
library(ggpubr)

### ---------------------------
### Generate performance (CM) data from sum_data
### ---------------------------
# Here we assume sum_data (from your .sum.csv files) has a column named
# test_sensitivity_95specificity representing the sensitivity at 95% specificity.
# We aggregate this metric by cancer type and feature.
# Aggregate sensitivity (at 95% specificity) across folds for each cancer-feature combination.
data_cm <- best_model_data %>%
  group_by(cancer, feature) %>%
  summarise(
    mean_sensitivity = mean(test_sensitivity_95specificity, na.rm = TRUE),
    sd_sensitivity   = sd(test_sensitivity_95specificity, na.rm = TRUE),
    # Specificity is assumed to be fixed at 0.95 because sensitivity was measured at that threshold.
    mean_specificity = 0.95,
    sd_specificity   = 0,
    .groups = "drop"
  )

# Map raw feature names to descriptive labels so they match the custom palette.
data_cm$Metric_label <- factor(data_cm$feature,
                               levels = c("all", "delfi", "motif", "insert_size", "methyl", "ns_peaks",
                                          "motif+methyl", "motif+methyl+delfi", "delfi+methyl"),
                               labels = c("Combined All Fragment Features",
                                          "Fragment Ratio",
                                          "End Motifs",
                                          "Insert Size",
                                          "Methylation",
                                          "Nucleosome Peak",
                                          "Combined Motif + Methylation",
                                          "Combined Motif + Methylation + Ratios",
                                          "Combined Ratios + Methylation")
)


### ---------------------------
### Define custom color palette (for both ROC and CM plots)
### ---------------------------
roc_color_palette <- c(
  "Combined All Fragment Features"             = "#6A3D9A",
  "Fragment Ratio"                           = "#009E73",
  "End Motifs"                               = "#E69F00",
  "Insert Size"                              = "#56B4E9",
  "Methylation"                              = "#CC79A7",
  "Nucleosome Peak"                          = "#0072B2",
  "Combined Motif + Methylation"        = "#32CD32",
  "Combined Motif + Methylation + Ratios" = "#FFD700",
  "Combined Ratios + Methylation"       = "#FF4500"
)

### ---------------------------
### ROC Plot (using aggregated roc_avg_dt data)
### ---------------------------
# We assume that roc_avg_dt has already been generated and includes columns:
# cancer, fpr, tpr, tpr_lower, tpr_upper, and a raw 'feature' column.
# First, map the raw feature names to descriptive labels.
roc_avg_dt$feature_label <- factor(roc_avg_dt$feature,
                                   levels = c("all", "delfi", "motif", "insert_size", "methyl", "ns_peaks",
                                              "motif+methyl", "motif+methyl+delfi", "delfi+methyl"),
                                   labels = c("Combined All Fragment Features",
                                              "Fragment Ratio",
                                              "End Motifs",
                                              "Insert Size",
                                              "Methylation",
                                              "Nucleosome Peak",
                                              "Combined Motif + Methylation",
                                              "Combined Motif + Methylation + Ratios",
                                              "Combined Ratios + Methylation")
)

roc_plot <- ggplot(roc_avg_dt, aes(x = fpr, y = tpr, color = feature_label, fill = feature_label)) +
  geom_ribbon(aes(ymin = tpr_lower, ymax = tpr_upper), alpha = 0.2, color = NA) +
  geom_line(size = 1, alpha = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50", alpha = 0.5) +
  facet_wrap(~ cancer, scales = "free") +
  scale_color_manual(values = roc_color_palette) +
  scale_fill_manual(values = roc_color_palette) +
  labs(title = "Mean ROC Curves by Cancer Type (Averaged Across Folds)",
       x = "False Positive Rate",
       y = "True Positive Rate",
       color = "Feature",
       fill = "Feature") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8))

ggsave("Mean_ROC_curves_by_cancer_with_CI_custom_best_model.png", roc_plot, width = 10, height = 8, dpi = 500)
print(roc_plot)

### ---------------------------
### Confusion Matrix (CM) Plot
### ---------------------------
CM_plot <- ggplot(data_cm, aes(x = mean_sensitivity, y = mean_specificity, color = Metric_label)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey") +
  geom_errorbarh(aes(xmin = pmax(0, mean_sensitivity - sd_sensitivity),
                     xmax = pmin(1, mean_sensitivity + sd_sensitivity)),
                 height = 0.02, alpha = 0.7) +
  geom_errorbar(aes(ymin = pmax(0, mean_specificity - sd_specificity),
                    ymax = pmin(1, mean_specificity + sd_specificity)),
                width = 0.02, alpha = 0.7) +
  scale_color_manual(values = roc_color_palette) +
  xlab("Sensitivity") +
  ylab("Specificity") +
  ggtitle("Performance") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.title = element_text(size = 10),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 8),
        legend.position = "none") +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1))

print(CM_plot)

### Change to bar plot 

# Create a bar plot for sensitivity at 95% specificity.
CM_bar_plot <- ggplot(data_cm, aes(x = Metric_label, y = mean_sensitivity, fill = Metric_label)) +
  geom_bar(stat = "identity", width = 0.7, color = "black") +
  geom_errorbar(aes(ymin = pmax(0, mean_sensitivity - sd_sensitivity),
                    ymax = pmin(1, mean_sensitivity + sd_sensitivity)),
                width = 0.2) +
  facet_wrap(~ cancer) +
  scale_fill_manual(values = roc_color_palette) +
  labs(title = "Sensitivity at 95% Specificity by Feature",
       x = "Feature",
       y = "Sensitivity at 95% Specificity") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# Display the plot
print(CM_bar_plot)

ggsave("CM_bar_plot_best_model.png", CM_bar_plot, width = 12, height = 16, dpi = 500)



#### Now do a cleavaland dot plot 
# Optionally reorder cancers by overall mean sensitivity (across features)
data_cm <- data_cm %>%
  group_by(cancer) %>%
  mutate(overall_mean = mean(mean_sensitivity, na.rm = TRUE)) %>%
  ungroup() %>%
  # Reorder 'cancer' so those with higher overall_mean appear at the bottom
  mutate(cancer = reorder(cancer, overall_mean))

# Build a Cleveland dot plot
ClevelandPlot <- ggplot(data_cm, aes(x = mean_sensitivity, y = cancer, color = Metric_label)) +
  # Plot the mean sensitivity as points
  geom_point(size = 3, alpha = 0.8, 
             position = position_dodge(width = 0.7)) +
  # Add horizontal error bars for standard deviation
  geom_errorbarh(aes(xmin = pmax(0, mean_sensitivity - sd_sensitivity),
                     xmax = pmin(1, mean_sensitivity + sd_sensitivity)),
                 height = 0.3, alpha = 0.7,
                 position = position_dodge(width = 0.7)) +
  # Use your custom color palette (ensure it matches Metric_label levels)
  scale_color_manual(values = roc_color_palette) +
  labs(title = "Sensitivity at 95% Specificity by Cancer & Feature",
       x = "Mean Sensitivity (95% Specificity)",
       y = "Cancer Type",
       color = "Feature") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.line  = element_line(color = "black"),
    legend.position = "bottom"        
    
  )

# Print and save the plot
print(ClevelandPlot)
ggsave("ClevelandPlot_Sensitivity95_wide_best_model.png", ClevelandPlot, width = 14, height = 10, dpi = 500)
ggsave("ClevelandPlot_Sensitivity95_narrow_best_model.png", ClevelandPlot, width = 7, height = 10, dpi = 500)
ggsave("ClevelandPlot_Sensitivity95_narrow_best_model2.png", ClevelandPlot, width = 7, height = 13, dpi = 500)


### If want larger text 
ClevelandPlot <- ggplot(data_cm, aes(x = mean_sensitivity, y = cancer, color = Metric_label)) +
  # Plot the mean sensitivity as points
  geom_point(size = 3, alpha = 0.8, position = position_dodge(width = 0.7)) +
  # Add horizontal error bars for standard deviation
  geom_errorbarh(
    aes(xmin = pmax(0, mean_sensitivity - sd_sensitivity),
        xmax = pmin(1, mean_sensitivity + sd_sensitivity)),
    height = 0.3, alpha = 0.7, position = position_dodge(width = 0.7)
  ) +
  # Use your custom color palette (ensure it matches Metric_label levels)
  scale_color_manual(values = roc_color_palette) +
  labs(
    title = "Sensitivity at 95% Specificity by Cancer & Feature",
    x = "Mean Sensitivity (95% Specificity)",
    y = "Cancer Type",
    color = "Feature"
  ) +
  theme_classic(base_size = 14) +  # base_size increases overall text size
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.line     = element_line(color = "black"),
    axis.text.y   = element_text(size = 14),  # Increase y-axis text size
    axis.text.x   = element_text(size = 12),
    axis.title.x  = element_text(size = 16),
    axis.title.y  = element_text(size = 16),
    legend.position  = "bottom",
    legend.title  = element_text(size = 14),
    legend.text   = element_text(size = 12)
  )

# Print and save the plot
print(ClevelandPlot)
ggsave("ClevelandPlot_Sensitivity95_LargerText_narrow.png", ClevelandPlot, width = 7, height = 10, dpi = 500)


### ---------------------------
### Combine ROC and CM Plots
### ---------------------------
final_fig <- ggarrange(roc_plot, CM_plot, ncol = 2)
final_fig <- annotate_figure(final_fig,
                             top = text_grob("Classifier Performance", face = "bold", size = 14))

ggsave("Combined_ROC_and_CM_custom.png", final_fig, width = 12, height = 6, dpi = 300)
print(final_fig)



### Now get more info 
# 1. Summarize performance by cancer type and feature using sensitivity at 95% specificity and test AUC (with ranges)
performance_summary <- best_model_data %>%
  group_by(cancer, feature) %>%
  summarise(
    mean_auc  = mean(test_auc, na.rm = TRUE),
    sd_auc    = sd(test_auc, na.rm = TRUE),
    min_auc   = min(test_auc, na.rm = TRUE),
    max_auc   = max(test_auc, na.rm = TRUE),
    mean_sens = mean(test_sensitivity_95specificity, na.rm = TRUE),
    sd_sens   = sd(test_sensitivity_95specificity, na.rm = TRUE),
    min_sens  = min(test_sensitivity_95specificity, na.rm = TRUE),
    max_sens  = max(test_sensitivity_95specificity, na.rm = TRUE),
    n         = n(),
    .groups   = "drop"
  )


# 2. Identify the best feature per cancer type (by highest mean AUC)
best_feature_by_cancer <- performance_summary %>%
  group_by(cancer) %>%
  filter(mean_auc == max(mean_auc, na.rm = TRUE)) %>%
  arrange(cancer) %>%
  ungroup()

# 3. Determine the overall performance by feature (averaging over all cancers, with ranges)
overall_feature_performance <- performance_summary %>%
  group_by(feature) %>%
  summarise(
    overall_mean_auc  = mean(mean_auc, na.rm = TRUE),
    overall_sd_auc    = mean(sd_auc, na.rm = TRUE),
    overall_min_auc   = min(min_auc, na.rm = TRUE),
    overall_max_auc   = max(max_auc, na.rm = TRUE),
    overall_mean_sens = mean(mean_sens, na.rm = TRUE),
    overall_sd_sens   = mean(sd_sens, na.rm = TRUE),
    overall_min_sens  = min(min_sens, na.rm = TRUE),
    overall_max_sens  = max(max_sens, na.rm = TRUE),
    count             = n(),
    .groups           = "drop"
  ) %>%
  arrange(desc(overall_mean_auc))


# Summarize methylation-only performance (feature == "methyl") with ranges
methyl_summary <- performance_summary %>%
  filter(feature == "methyl") %>%
  select(cancer, 
         mean_auc_methyl = mean_auc, 
         sd_auc_methyl   = sd_auc,
         min_auc_methyl  = min_auc,
         max_auc_methyl  = max_auc,
         mean_sens_methyl = mean_sens,
         sd_sens_methyl   = sd_sens,
         min_sens_methyl  = min_sens,
         max_sens_methyl  = max_sens)

# Summarize combined performance (for any feature containing "Combined")
# Adjust the grepl pattern if your "combined" features have different naming.
# Summarize combined performance (for any feature containing "all" or "+") with ranges
combined_summary <- performance_summary %>%
  filter(grepl("all|\\+", feature, ignore.case = TRUE)) %>%
  group_by(cancer) %>%
  summarise(
    mean_auc_combined = mean(mean_auc, na.rm = TRUE),
    sd_auc_combined   = mean(sd_auc, na.rm = TRUE),
    min_auc_combined  = min(min_auc, na.rm = TRUE),
    max_auc_combined  = max(max_auc, na.rm = TRUE),
    mean_sens_combined = mean(mean_sens, na.rm = TRUE),
    sd_sens_combined   = mean(sd_sens, na.rm = TRUE),
    min_sens_combined  = min(min_sens, na.rm = TRUE),
    max_sens_combined  = max(max_sens, na.rm = TRUE),
    .groups = "drop"
  )



# Join methylation and combined summaries by cancer and compute gains with ranges
comparison <- left_join(methyl_summary, combined_summary, by = "cancer") %>%
  mutate(
    auc_gain_mean  = mean_auc_combined - mean_auc_methyl,
    auc_gain_sd    = sd_auc_combined - sd_auc_methyl,
    auc_gain_min   = min_auc_combined - min_auc_methyl,
    auc_gain_max   = max_auc_combined - max_auc_methyl,
    sens_gain_mean = mean_sens_combined - mean_sens_methyl,
    sens_gain_sd   = sd_sens_combined - sd_sens_methyl,
    sens_gain_min  = min_sens_combined - min_sens_methyl,
    sens_gain_max  = max_sens_combined - max_sens_methyl
  )
print("Comparison of Methylation-only vs. Combined Models (with Ranges):")
print(comparison)



# ---------------------------
# Comparison 2: Specific Feature ("motif+methyl+delfi")
# ---------------------------
# Comparison 2: Specific Feature ("motif+methyl+delfi") with ranges
specific_summary <- performance_summary %>%
  filter(tolower(feature) == "motif+methyl+delfi") %>%
  select(cancer, 
         mean_auc_specific = mean_auc, 
         sd_auc_specific   = sd_auc,
         min_auc_specific  = min_auc,
         max_auc_specific  = max_auc,
         mean_sens_specific = mean_sens, 
         sd_sens_specific   = sd_sens,
         min_sens_specific  = min_sens,
         max_sens_specific  = max_sens)

comparison_specific <- left_join(methyl_summary, specific_summary, by = "cancer") %>%
  mutate(
    auc_gain_specific  = mean_auc_specific - mean_auc_methyl,
    sd_gain_specific   = sd_auc_specific - sd_auc_methyl,
    min_gain_specific  = min_auc_specific - min_auc_methyl,
    max_gain_specific  = max_auc_specific - max_auc_methyl,
    sens_gain_specific = mean_sens_specific - mean_sens_methyl,
    sd_sens_gain_specific = sd_sens_specific - sd_sens_methyl,
    min_sens_gain_specific = min_sens_specific - min_sens_methyl,
    max_sens_gain_specific = max_sens_specific - max_sens_methyl
  )
print("Comparison between methylation-only and 'motif+methyl+delfi' (with Ranges):")
print(comparison_specific)

# Optionally, calculate overall average gains across cancers for each comparison:
# Optionally, calculate overall gains across cancers for each comparison (with ranges):
overall_combined <- comparison %>%
  summarise(
    overall_auc_gain_mean  = mean(auc_gain_mean, na.rm = TRUE),
    overall_auc_gain_min   = min(auc_gain_mean, na.rm = TRUE),
    overall_auc_gain_max   = max(auc_gain_max, na.rm = TRUE),
    overall_sens_gain_mean = mean(sens_gain_mean, na.rm = TRUE),
    overall_sens_gain_min  = min(sens_gain_mean, na.rm = TRUE),
    overall_sens_gain_max  = max(sens_gain_max, na.rm = TRUE)
  )
print("Overall Gains (Combined vs. Methylation-only):")
print(overall_combined)

overall_specific <- comparison_specific %>%
  summarise(
    overall_auc_gain_specific  = mean(auc_gain_specific, na.rm = TRUE),
    overall_auc_gain_specific_min = min(auc_gain_specific, na.rm = TRUE),
    overall_auc_gain_specific_max = max(auc_gain_specific, na.rm = TRUE),
    overall_sens_gain_specific = mean(sens_gain_specific, na.rm = TRUE),
    overall_sens_gain_specific_min = min(sens_gain_specific, na.rm = TRUE),
    overall_sens_gain_specific_max = max(sens_gain_specific, na.rm = TRUE)
  )
print("Overall Gains (motif+methyl+delfi vs. Methylation-only):")
print(overall_specific)



# 5. Print the summarized data to inspect
print("Best feature per cancer type:")
print(best_feature_by_cancer)

print("Overall feature performance:")
print(overall_feature_performance)

print("Comparison methylation vs. combined features:")
print(comparison)



### Now get model summaries and make plot: 

# Summarize performance by cancer type and model
model_summary <- sum_data %>%
  group_by(cancer, model) %>%
  summarise(
    mean_auc = mean(test_auc, na.rm = TRUE),
    sd_auc   = sd(test_auc, na.rm = TRUE),
    n        = n(),
    .groups  = "drop"
  )

# Optionally, reorder cancers by the maximum mean AUC (or any desired metric)
cancer_order <- model_summary %>%
  group_by(cancer) %>%
  summarise(max_auc = max(mean_auc, na.rm = TRUE)) %>%
  arrange(max_auc) %>%
  pull(cancer)
model_summary$cancer <- factor(model_summary$cancer, levels = cancer_order)

# Build a Cleveland dot plot: x-axis = mean AUC, y-axis = cancer type; points colored by model
cleveland_plot <- ggplot(model_summary, aes(x = mean_auc, y = cancer, color = model)) +
  geom_point(size = 3, position = position_dodge(width = 0.7)) +
  geom_errorbarh(aes(xmin = pmax(0, mean_auc - sd_auc), xmax = pmin(1, mean_auc + sd_auc)),
                 height = 0.2, position = position_dodge(width = 0.7)) +
  labs(
    title = "Model Performance by Cancer Type",
    x = "Mean Test AUC",
    y = "Cancer Type",
    color = "Model"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.y  = element_text(size = 12),
    axis.text.x  = element_text(size = 12),
    axis.title   = element_text(size = 14),
    legend.position = "bottom",
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12)
  )

# Save and print the Cleveland plot
ggsave("ClevelandPlot_Models_AUC.png", cleveland_plot, width = 10, height = 8, dpi = 300)
print(cleveland_plot)

### Updated version 
# Map short model codes to descriptive names
model_summary <- model_summary %>%
  mutate(ModelName = dplyr::recode_factor(model,
                                          "lr"  = "Logistic Regression",
                                          "svm" = "Support Vector Machine",
                                          "knn" = "K-Nearest Neighbors",
                                          "rf"  = "Random Forest",
                                          "xgb" = "XGBoost",
                                          "gbm" = "Gradient Boosting Machine",
                                          "lda" = "Linear Discriminant Analysis",
                                          .default = "Other"
  ))

# Optionally, reorder cancers by their maximum mean AUC
cancer_order <- model_summary %>%
  group_by(cancer) %>%
  summarise(max_auc = max(mean_auc, na.rm = TRUE)) %>%
  arrange(max_auc) %>%
  pull(cancer)

model_summary$cancer <- factor(model_summary$cancer, levels = cancer_order)

# Build a Cleveland dot plot: x-axis = mean AUC, y-axis = cancer type; points colored by model
cleveland_plot <- ggplot(model_summary, aes(x = mean_auc, y = cancer, color = ModelName)) +
  geom_point(size = 3, position = position_dodge(width = 0.7)) +
  geom_errorbarh(
    aes(xmin = pmax(0, mean_auc - sd_auc), xmax = pmin(1, mean_auc + sd_auc)),
    height = 0.2,
    position = position_dodge(width = 0.7)
  ) +
  labs(
    title = "Model Performance by Cancer Type",
    x = "Mean Test AUC",
    y = "Cancer Type",
    color = "Model"
  ) +
  # Use a neutral ColorBrewer palette, e.g. Set2 or Dark2
  scale_color_brewer(palette = "Set2") +
  theme_classic(base_size = 14) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.y     = element_text(size = 12),
    axis.text.x     = element_text(size = 12),
    axis.title      = element_text(size = 14),
    legend.position = "bottom",
    legend.title    = element_text(size = 14),
    legend.text     = element_text(size = 12)
  )

# Save and print the Cleveland plot
ggsave("ClevelandPlot_Models_AUC_updated.png", cleveland_plot, width = 13, height = 10, dpi = 500)
print(cleveland_plot)

ggsave("ClevelandPlot_Models_AUC_updated_cancer_type_long.png", cleveland_plot, width = 7, height = 10, dpi = 500)




print(paste("Overall AUC gain (combined vs. methylation alone):", round(overall_auc_gain, 3)))
print(paste("Overall Sensitivity gain (combined vs. methylation alone):", round(overall_sens_gain, 3)))



# -------------------------------------------------------------------------
# Create a paragraph for the Results section based on the computed values:
# -------------------------------------------------------------------------

### Summaries of Results ###

# 1) Summaries of the best features
# best_feature_by_cancer has columns: cancer, feature, mean_auc, mean_sens, ...
# overall_feature_performance has columns: feature, overall_mean_auc, overall_mean_sens, count, ...
cat("\n\n--- Best Feature Summaries ---\n\n")

# A short paragraph about which feature was best for each cancer:
best_by_cancer_text <- best_feature_by_cancer %>%
  mutate(txt = paste0(
    "For ", cancer, ", the best feature was '", feature, 
    "' with a mean AUC of ", round(mean_auc, 3),
    " (range: ", round(min_auc, 3), " - ", round(max_auc, 3), "), and a mean sensitivity of ", 
    round(mean_sens, 3), " (range: ", round(min_sens, 3), " - ", round(max_sens, 3), ")."
  )) %>%
  pull(txt) %>%
  paste(collapse = "\n")
cat("Best Feature per Cancer Type:\n")
cat(best_by_cancer_text, "\n\n")


cat("Best Feature per Cancer Type:\n")
cat(best_by_cancer_text, "\n\n")

# Identify the top feature overall by AUC
top_feature <- overall_feature_performance$feature[1]
top_feature_auc <- round(overall_feature_performance$overall_mean_auc[overall_feature_performance$feature == top_feature], 3)
top_feature_sens <- round(overall_feature_performance$overall_mean_sens[overall_feature_performance$feature == top_feature], 3)
overall_features_text <- paste0(
  "Overall, the highest-performing feature (averaged across all cancer types) was '", 
  top_feature, "', which achieved a mean AUC of ", top_feature_auc, 
  " (range: ", round(overall_feature_performance$overall_min_auc[overall_feature_performance$feature == top_feature], 3), 
  " - ", round(overall_feature_performance$overall_max_auc[overall_feature_performance$feature == top_feature], 3), 
  ") and a mean sensitivity of ", top_feature_sens, 
  " (range: ", round(overall_feature_performance$overall_min_sens[overall_feature_performance$feature == top_feature], 3), 
  " - ", round(overall_feature_performance$overall_max_sens[overall_feature_performance$feature == top_feature], 3), ")."
)
cat(overall_features_text, "\n\n")


# 2) Summaries of combined vs. methylation-only
# comparison has columns: cancer, mean_auc_methyl, mean_auc_combined, ...
# overall_combined has columns: overall_auc_gain_mean, overall_auc_gain_max, ...
cat("\n\n--- Comparison: Methylation vs. Combined Features ---\n\n")
mean_gain_auc  <- round(overall_combined$overall_auc_gain_mean, 3)
max_gain_auc   <- round(overall_combined$overall_auc_gain_max, 3)
min_gain_auc   <- round(overall_combined$overall_auc_gain_min, 3)
mean_gain_sens <- round(overall_combined$overall_sens_gain_mean, 3)
max_gain_sens  <- round(overall_combined$overall_sens_gain_max, 3)
min_gain_sens  <- round(overall_combined$overall_sens_gain_min, 3)
combo_text <- paste0(
  "When comparing methylation-only classifiers to those that included combined features (i.e., features containing 'all' or a plus sign), ",
  "we observed an average AUC gain of ", mean_gain_auc, " (range: ", min_gain_auc, " - ", max_gain_auc, 
  "), and an average sensitivity gain of ", mean_gain_sens, " (range: ", min_gain_sens, " - ", max_gain_sens, 
  "). These gains highlight the additional predictive value provided by combining fragmentomic features with methylation."
)
cat(combo_text, "\n\n")


# 3) Summaries of the specific feature "motif+methyl+delfi"
# comparison_specific has columns: mean_auc_specific, mean_auc_methyl, ...
# overall_specific has columns: overall_auc_gain_specific, overall_sens_gain_specific
cat("\n\n--- Comparison: Methylation vs. 'motif+methyl+delfi' ---\n\n")
spec_auc_gain  <- round(overall_specific$overall_auc_gain_specific, 3)
spec_sens_gain <- round(overall_specific$overall_sens_gain_specific, 3)
min_gain_specific <- round(overall_specific$overall_auc_gain_specific_min, 3)
max_gain_specific <- round(overall_specific$overall_auc_gain_specific_max, 3)
min_sens_gain_specific <- round(overall_specific$overall_sens_gain_specific_min, 3)
max_sens_gain_specific <- round(overall_specific$overall_sens_gain_specific_max, 3)
specific_text <- paste0(
  "We also compared methylation-only classifiers with the specific combined feature 'motif+methyl+delfi'. ",
  "Across the cancers examined, this feature yielded an average AUC gain of ", spec_auc_gain, 
  " (range: ", min_gain_specific, " - ", max_gain_specific, 
  ") and an average sensitivity gain of ", spec_sens_gain, 
  " (range: ", min_sens_gain_specific, " - ", max_sens_gain_specific, 
  ") compared to methylation alone. ",
  "In certain cancer types, this improvement was even higher, underscoring the potential synergy of combining motif and fragmentomic signals with methylation."
)
cat(specific_text, "\n\n")


# 4) Summaries of which cancer types performed best overall (by top mean AUC)
#   and how each model performed from model_summary
cat("\n\n--- Cancer Type & Model Performance ---\n\n")

# model_summary has columns: cancer, model, mean_auc, sd_auc, n
# Identify the top 3 best-performing cancers by maximum mean_auc
best_cancers <- model_summary %>%
  group_by(cancer) %>%
  summarise(max_auc_cancer = max(mean_auc, na.rm = TRUE)) %>%
  arrange(desc(max_auc_cancer)) %>%
  slice(1:3)

best_cancers_text <- best_cancers %>%
  mutate(txt = paste0(
    cancer, " (max mean AUC: ", round(max_auc_cancer, 3), ")"
  )) %>%
  pull(txt) %>%
  paste(collapse = ", ")

cancer_perf_text <- paste0(
  "Among all cancer types, the top three in terms of maximum mean AUC were: ", 
  best_cancers_text, ". These results suggest that these particular cancer types ",
  "are more readily distinguishable by the features/models employed."
)

cat(cancer_perf_text, "\n\n")

# Identify the best 2-3 models overall by mean AUC across all cancers
model_rank <- model_summary %>%
  group_by(model) %>%
  summarise(
    global_mean_auc = mean(mean_auc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(global_mean_auc))

top_models <- model_rank %>% slice(1:3)
top_models_text <- top_models %>%
  mutate(txt = paste0(
    model, " (mean AUC: ", round(global_mean_auc, 3), ")"
  )) %>%
  pull(txt) %>%
  paste(collapse = ", ")

model_perf_text <- paste0(
  "In terms of algorithms, the top three models by average AUC across all cancer types were: ",
  top_models_text, ". This indicates that these methods provided the most robust performance overall."
)

cat(model_perf_text, "\n\n")


### Export important plots 
# Define the export directory
export_dir <- "Final Data Dec 2024/ML_Data/Yuanchang_processed_cancer_type"

# Create the directory (and any parent directories) if it doesn't exist
if (!dir.exists(export_dir)) {
  dir.create(export_dir, recursive = TRUE)
}



## See which feature is best for each feature type 
#best_feature_summary <- read.csv("Final Data Dec 2024/ML_data/Yuanchang_processed_cancer_type/best_feature_summary.csv")
#model_summary <- read.csv("Final Data Dec 2024/ML_data/Yuanchang_processed_cancer_type/model_summary.csv")

# Assuming best_feature_summary is already loaded
top_feature_per_cancer <- best_feature_summary %>%
  group_by(cancer) %>%
  slice_max(order_by = mean_auc, n = 1, with_ties = FALSE) %>%
  arrange(desc(mean_auc))


# Calculate performance gains for cancer types where a non-"methyl" model outperformed methylation alone
best_vs_methyl <- performance_summary %>%
  group_by(cancer) %>%
  filter(any(feature != "methyl") & any(feature == "methyl")) %>%
  summarize(
    best_model_auc = max(mean_auc[feature != "methyl"], na.rm = TRUE),
    best_model_sens = mean_sens[which.max(if_else(feature != "methyl", mean_auc, -Inf))],
    methyl_auc = mean_auc[feature == "methyl"][1],
    methyl_sens = mean_sens[feature == "methyl"][1]
  ) %>%
  filter(best_model_auc > methyl_auc) %>%
  mutate(
    auc_gain = best_model_auc - methyl_auc,
    sens_gain = best_model_sens - methyl_sens
  )

# Print the table of gains
print(best_vs_methyl)

# Average gains
avg_auc_gain  <- round(mean(best_vs_methyl$auc_gain), 3)
avg_sens_gain <- round(mean(best_vs_methyl$sens_gain), 3)
num_cases     <- nrow(best_vs_methyl)

# Range of gains + cancer type labels
min_gain_auc  <- round(min(best_vs_methyl$auc_gain), 3)
max_gain_auc  <- round(max(best_vs_methyl$auc_gain), 3)
min_auc_cancer <- best_vs_methyl$cancer[which.min(best_vs_methyl$auc_gain)]
max_auc_cancer <- best_vs_methyl$cancer[which.max(best_vs_methyl$auc_gain)]

min_gain_sens <- round(min(best_vs_methyl$sens_gain), 3)
max_gain_sens <- round(max(best_vs_methyl$sens_gain), 3)
min_sens_cancer <- best_vs_methyl$cancer[which.min(best_vs_methyl$sens_gain)]
max_sens_cancer <- best_vs_methyl$cancer[which.max(best_vs_methyl$sens_gain)]

# Total cancer types evaluated
total_cases <- performance_summary %>%
  group_by(cancer) %>%
  filter(any(feature == "methyl")) %>%
  summarise() %>%
  nrow()

# Generate detailed per-cancer summary
best_vs_methyl_text <- best_vs_methyl %>%
  mutate(txt = paste0("For ", cancer, ", the best non-methyl model improved mean AUC by ", 
                      round(auc_gain, 3), " (", round(methyl_auc, 3), " vs. ", round(best_model_auc, 3),
                      ") and sensitivity by ", round(sens_gain, 3), " (", round(methyl_sens, 3), " vs. ", 
                      round(best_model_sens, 3), ").")) %>%
  pull(txt) %>%
  paste(collapse = "\n")

cat(best_vs_methyl_text, "\n\n")

# Overall summary with range + cancer names
overall_text <- paste0(
  "Across ", num_cases, " of ", total_cases, " cancer types where a non-methyl model outperformed methylation alone, ",
  "the average gain was ", avg_auc_gain, " in mean AUC (range: ", min_gain_auc, " in ", min_auc_cancer, 
  " to ", max_gain_auc, " in ", max_auc_cancer, ") and ", avg_sens_gain, " in sensitivity (range: ", 
  min_gain_sens, " in ", min_sens_cancer, " to ", max_gain_sens, " in ", max_sens_cancer, ")."
)
cat(overall_text, "\n")



# View the result
print(top_feature_per_cancer)
# Export important tables

# 1. Best Feature Summary by Cancer Type
write_csv(best_feature_summary, file.path(export_dir, "best_feature_summary.csv"))

write_csv(top_feature_per_cancer, file.path(export_dir, "top_feature_per_cancer.csv"))

write_csv(sum_data, file.path(export_dir, "sum_data_all_models.csv"))
write_csv(roc_data_joined, file.path(export_dir, "roc_data_joined_all_models.csv"))


# 2. Best Model Summary by Cancer Type
write_csv(best_model_summary, file.path(export_dir, "best_model_summary.csv"))

# 3. Fold-level Summary Data
write_csv(fold_summary, file.path(export_dir, "fold_summary.csv"))

# 4. Best Model by Cancer (the chosen model per cancer type)
write_csv(best_model_by_cancer, file.path(export_dir, "best_model_by_cancer.csv"))

# 5. Performance Summary (based on best model data; by cancer and feature)
write_csv(performance_summary, file.path(export_dir, "performance_summary.csv"))

# 6. Overall Feature Performance (averaged across cancers)
write_csv(overall_feature_performance, file.path(export_dir, "overall_feature_performance.csv"))

# 7. Methylation-only Summary
write_csv(methyl_summary, file.path(export_dir, "methyl_summary.csv"))

# 8. Combined Features Summary (general, i.e. features containing "all" or "+")
write_csv(combined_summary, file.path(export_dir, "combined_summary.csv"))

# 9. Comparison of Methylation-only vs. Combined Models
write_csv(comparison, file.path(export_dir, "comparison_methyl_vs_combined.csv"))

# 10. Comparison of Methylation-only vs. 'motif+methyl+delfi'
write_csv(comparison_specific, file.path(export_dir, "comparison_methyl_vs_motif+methyl+delfi.csv"))

# 11. Overall Gains (Combined vs. Methylation-only)
write_csv(overall_combined, file.path(export_dir, "overall_gains_combined_vs_methyl.csv"))

# 12. Overall Gains (motif+methyl+delfi vs. Methylation-only)
write_csv(overall_specific, file.path(export_dir, "overall_gains_motif+methyl+delfi_vs_methyl.csv"))

# 13. Model Summary (by cancer and model; before mapping to longer names)
write_csv(model_summary, file.path(export_dir, "model_summary.csv"))


# Save the best_vs_methyl dataframe to CSV
write_csv(best_vs_methyl, file.path(export_dir, "best_vs_methyl_summary.csv"))

# Save detailed per-cancer summary text to a text file
writeLines(best_vs_methyl_text, file.path(export_dir, "best_vs_methyl_per_cancer_summary.txt"))

# Save overall summary text to a text file
writeLines(overall_text, file.path(export_dir, "overall_best_vs_methyl_summary.txt"))



cat("Export complete. All summary tables have been saved to:\n", export_dir)

