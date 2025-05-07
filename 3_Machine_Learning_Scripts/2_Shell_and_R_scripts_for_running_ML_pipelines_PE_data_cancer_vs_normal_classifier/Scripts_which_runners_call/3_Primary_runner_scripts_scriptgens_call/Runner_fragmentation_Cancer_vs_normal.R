# runner 
# for cluster
library(caret)
library(dplyr)
library(tidyverse)
library(pROC)
library(doParallel)

### Set variables
args <- commandArgs(trailingOnly = TRUE)
technology_name <- args[1]
input_data <- args[2]
metadata_file <- args[3]
outdir <- args[4]
title <- args[5]
title2 <- args[6]
sample_id_column <- args[7]
type1 <- args[8]
type2 <- args[9]
external_data_file <- args[10] # e.g., "/path/to/external_validation_data.rds"
scripts <- "/cluster/projects/tcge/cell_free_epigenomics/processed_data/fragmentomics_dory/Run2/Machine_Learning/January_2025/Scripts"

### Import data (Starting with the 5Mb ratios), from the delfi ratio script
data <- readRDS(file = input_data)
data <- as.data.frame(data)

print(paste("Data dimensions: ", dim(data)))

## Import the external data file 
if (file.exists(external_data_file)) {
  data.external.validation <- readRDS(external_data_file)
  message("Successfully loaded external validation data.")
} else {
  data.external.validation <- NULL
  message("External validation file not found or invalid path: ", external_data_file)
}

print(paste("External data dimensions: ", dim(data.external.validation)))

## Load metadata
metadata_df <- readRDS(file = metadata_file)
metadata_df <- as.data.frame(metadata_df)

print(paste("Metadata dimensions: ", dim(metadata_df)))

scaleFUN <- function(x) sprintf("%.2f", x)
theme <- theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10), 
               axis.line = element_line(colour = "black"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank(),
               legend.key = element_rect(fill = "white"),
               legend.text = element_text(size = 12),
               legend.position = "none",
               legend.background = element_blank(),
               axis.text = element_text(size = 12),
               axis.title = element_text(size = 12),
               axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

print("Theme set")

# path to fun classifier
source(file.path(scripts, "FuncClassifier.R"))

#path to classifier.R
class_script <- file.path(scripts, "final_classifier.R")

### Set models to run
#classifiers <- c("rf", "glm", "gbm", "svmRadial", "knn")
#names <- c("RF", "GLM", "GBM", "SVM", "KNN")
classifiers <- c("rf", "glm", "gbm", "svmRadial", "knn", "xgbTree", "glmnet", "pls", "bayesglm", "lda")
names <- c("RF", "GLM", "GBM", "SVM", "KNN", "XGB", "LASSO", "PLS", "BayesGLM", "LDA")


# The list of cancer types
#cancer_types <- c("Normal", "Blood", "Brain", "Breast", "Colorectal", "Head_and_Neck", "LTX", "Lung", "Pancreatic", "Prostate")

results <- data.frame(test = c(), kappa = c(), CI = c())

# Get data
samples1 <- metadata_df[metadata_df$CN_classifier == type1,]
samples2 <- metadata_df[metadata_df$CN_classifier == type2,]

# Troubleshooting 
print(paste("Samples1 dimensions: ", dim(samples1)))
print(paste("Samples2 dimensions: ", dim(samples2)))

# Check if samples1 and samples2 are not empty
if (nrow(samples1) <= 1 | nrow(samples2) <= 1) {
  stop("Not enough samples in one or both groups.")
}

# Inner loop through classifiers
for (k in 1:length(classifiers)) {
  # Extract training samples from group "type1" (e.g., healthy)
  normal <- data[rownames(data) %in% samples1[[sample_id_column]], ]
  # Extract samples from group "type2" (e.g., cancer)
  cancer <- data[rownames(data) %in% samples2[[sample_id_column]], ]
  # Combine (here we assume you want to combine normal and the other group)
  Data <- bind_rows(cancer, normal)
  
  # Create the class label vector: type2 for samples from samples1, and type1 for samples from samples2
  y <- c(rep(type2, sum(rownames(Data) %in% samples2[[sample_id_column]])),
         rep(type1, sum(rownames(Data) %in% samples1[[sample_id_column]])))
  
  y <- factor(y, levels = c(type2, type1))
  
  ### Set variables
  class <- classifiers[[k]]
  name <- names[[k]]
  # Rename variable "class" to "algorithm" to avoid conflicts.
  algorithm <- classifiers[[k]]

  ### Run suite of classifiers
  set.seed(123)
  source(file.path(scripts, "classifier_opt.R"))

  ### Get Kappa performance
  mean <- mean(kappa)
  margin <- qt(0.975, df = length(kappa) - 1) * sd(kappa) / sqrt(length(kappa))

  ### Append to dataframe
  a <- data.frame(test = name, kappa = mean, error = margin)
  results <- rbind(results, a)
}

### Plot Kappa values
results <- results[order(results$kappa), ]
results$test <- factor(results$test, levels = results$test)
plot <- ggplot(results) +
  geom_point(aes(x = test, y = kappa)) +
  geom_errorbar(aes(x = test, ymin = kappa - error, ymax = kappa + error), width = 0) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5) +
  xlab("Algorithm") + 
  ylab("Kappa") +
  labs(title = title2, subtitle = paste0(type1, " vs ", type2), color = "") +
  theme +
  scale_y_continuous(labels = scaleFUN) +
  coord_flip()
plot
ggsave(file.path(outdir, paste0(technology_name, "_kappa_", type1, "_vs_", type2, ".pdf")), plot, width = 4, height = 5)
write.table(results, file = file.path(outdir, paste0(technology_name, "_results_", type1, "_vs_", type2, ".txt")), row.names = F, quote = F, sep = "\t")

# Now do classification
name_to_classifier <- c(
  "RF"       = "rf",
  "GLM"      = "glm",
  "GBM"      = "gbm",
  "SVM"      = "svmRadial",
  "KNN"      = "knn",
  "XGB"      = "xgbTree",
  "LASSO"    = "glmnet",
  "PLS"      = "pls",
  "BayesGLM" = "bayesglm",
  "LDA"      = "lda"
)

best_model_name <- as.character(tail(results$test, n = 1))

# Clear any existing definition of 'algorithm'
if (exists("algorithm")) {
  rm(algorithm)
  cat("DEBUG: 'algorithm' has been removed from the environment.\n")
}

algorithm <- name_to_classifier[best_model_name]

cat("DEBUG: New 'algorithm' is set to:", algorithm, "\n")

# Source the classification script to perform predictions
source(class_script)

# Save the internal cross-validation results (All.kFold)
dynamic_var_name <- paste0("outcomes_", type1, "_vs_", type2)
assign(dynamic_var_name, All.kFold)
internal_cv_path <- file.path(outdir, paste0(technology_name, "_outcomes_", type1, "_vs_", type2, ".rds"))
saveRDS(get(dynamic_var_name), internal_cv_path)
message("Internal cross-validation results saved to: ", internal_cv_path)

# Save the external validation results, if available (External.kfold)
if (exists("External.kfold") && !is.null(External.kfold)) {
  validation_var_name <- paste0("validation_outcomes_", type1, "_vs_", type2)
  assign(validation_var_name, External.kfold)
  external_validation_path <- file.path(outdir, paste0(technology_name, "_validation_outcomes_", type1, "_vs_", type2, ".rds"))
  saveRDS(get(validation_var_name), external_validation_path)
  message("External validation predictions saved to: ", external_validation_path)
} else {
  message("No external validation predictions to save (External.kfold is NULL or does not exist).")
}

# Save the best caret cross-validation model, if used
if (exists("caret_cv_model") && !is.null(caret_cv_model)) {
  caret_model_path <- file.path(outdir, paste0(technology_name, "_caret_cv_model_", type1, "_vs_", type2, ".rds"))
  saveRDS(caret_cv_model, caret_model_path)
  message("Caret cross-validation model saved to: ", caret_model_path)
}

# Save the final model trained on all data
if (exists("final_model") && !is.null(final_model)) {
  final_model_path <- file.path(outdir, paste0(technology_name, "_final_model_", type1, "_vs_", type2, ".rds"))
  saveRDS(final_model, final_model_path)
  message("Final model saved to: ", final_model_path)
}

# Save external validation results from the final model, if applicable
if (exists("external_results") && !is.null(external_results)) {
  external_results_path <- file.path(outdir, paste0(technology_name, "_final_external_results_", type1, "_vs_", type2, ".rds"))
  saveRDS(external_results, external_results_path)
  message("Final external validation predictions saved to: ", external_results_path)
} else {
  message("No final external validation predictions to save (external_results is NULL or does not exist).")
}

# Save any additional metrics, if computed
if (exists("final_metrics") && !is.null(final_metrics)) {
  final_metrics_path <- file.path(outdir, paste0(technology_name, "_final_metrics_", type1, "_vs_", type2, ".rds"))
  saveRDS(final_metrics, final_metrics_path)
  message("Final metrics saved to: ", final_metrics_path)
}

# Summary message
message("All relevant outputs have been saved in the specified output directory: ", outdir)