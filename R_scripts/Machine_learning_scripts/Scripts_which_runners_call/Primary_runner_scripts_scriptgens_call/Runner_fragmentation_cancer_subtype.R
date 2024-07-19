# runner 
# for cluster
library(caret)
library(dplyr)
library(tidyverse)
library(pROC)

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
scripts <- "/cluster/projects/tcge/cell_free_epigenomics/processed_data/fragmentomics_dory/Run2/Machine_Learning/Delfi_ratios_5Mb/Scripts"

### Import data (Starting with the 5Mb ratios), from the delfi ratio script
data <- readRDS(file = input_data)
data <- as.data.frame(data)

print(paste("Data dimensions: ", dim(data)))


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
class_script <- file.path(scripts, "classifier_updated.R")

### Set models to run
classifiers <- c("rf", "glm", "gbm", "svmRadial", "knn")
names <- c("RF", "GLM", "GBM", "SVM", "KNN")

# The list of cancer types
#cancer_types <- c("Normal", "Blood", "Brain", "Breast", "Colorectal", "Head_and_Neck", "LTX", "Lung", "Pancreatic", "Prostate")

results <- data.frame(test = c(), kappa = c(), CI = c())

# Get data
samples1 <- metadata_df[metadata_df$cancer_subtype_corrected == type1,]
samples2 <- metadata_df[metadata_df$cancer_subtype_corrected == type2,]

# Troubleshooting 
print(paste("Samples1 dimensions: ", dim(samples1)))
print(paste("Samples2 dimensions: ", dim(samples2)))

# Check if samples1 and samples2 are not empty
if (nrow(samples1) <= 1 | nrow(samples2) <= 1) {
  stop("Not enough samples in one or both groups.")
}

# Inner loop through classifiers
for (k in 1:length(classifiers)) {
  normal <- data[row.names(data) %in% samples1[[sample_id_column]], ]
  Data <- data[row.names(data) %in% samples2[[sample_id_column]], ]
  Data <- bind_rows(Data, normal)

  y <- c(rep(type2, sum(row.names(Data) %in% samples1[[sample_id_column]])), 
         rep(type1, sum(row.names(Data) %in% samples2[[sample_id_column]])))

  y <- factor(y, levels = c(type1, type2))

  ### Set variables
  class <- classifiers[[k]]
  name <- names[[k]]

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
ggsave(file.path(outdir, paste0(technology_name, "_kappa_", type1, "_vs_", type2, ".pdf")), plot, width = 2, height = 2.5)
write.table(results, file = file.path(outdir, paste0(technology_name, "_results_", type1, "_vs_", type2, ".txt")), row.names = F, quote = F, sep = "\t")

# Now do classification
name_to_classifier <- c("RF" = "rf", "GLM" = "glm", "GBM" = "gbm", "SVM" = "svmRadial", "KNN" = "knn")
best_model_name <- as.character(tail(results$test, n = 1))
algorithm <- name_to_classifier[best_model_name]
source(class_script)
dynamic_var_name <- paste0("outcomes_", type1, "_vs_", type2)
assign(dynamic_var_name, All.kFold)
saveRDS(get(dynamic_var_name), file.path(outdir, paste0(technology_name, "_outcomes_", type1, "_vs_", type2, ".rds")))
