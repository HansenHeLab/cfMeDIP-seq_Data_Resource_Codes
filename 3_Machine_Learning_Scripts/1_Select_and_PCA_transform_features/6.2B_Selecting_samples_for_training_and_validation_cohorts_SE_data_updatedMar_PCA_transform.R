###############################################################################
# File:      6.2B_Selecting_samples_for_training_and_validation_cohorts_SE_data_updatedMar_PCA_transform.R
# Author:    Dory Abelman
# Date:      Dec 2024
#
# Purpose:
#   1. Load & harmonize single‑end (SE) fragmentomic and methylation features:
#        - End‑motifs and combat‑seq adjusted DMR methylation counts.
#   2. Filter SE samples and derive binary cancer vs normal (CN) classifier labels 
#      (cancer vs. healthy) from metadata.
#   3. Stratified 90/10 train–validation split by CN_classifier:
#        - Ensures ≥5 samples (or 10%) per class in the hold‑out set.
#   4. Align feature matrices to metadata subsets:
#        - Subset, transpose (methylation), and split into train/validation.
#   5. PCA reduction on each feature set (train & validation):
#        - Compute & export scree plots, cumulative‑variance plots.
#        - Select PCs at multiple variance thresholds (e.g. 80, 85, 90, 97%).
#        - Save PCA objects, PC score tables, and retention summaries.
#   6. Merge PC features across end‑motif & methylation (>1% & >5% variance):
#        - Export combined training/validation matrices for downstream ML.
#   7. Validate outputs & NA diagnostics.
#
# Usage:
#   Source this script to prepare SE datasets for classification modeling.
###############################################################################


## Prepare data for ML classification - Dory Abelman
## Dec 2024
## Splits data into 90/10 for SE data 

# Load necessary libraries
library(ggplot2)
library(factoextra) # For PCA visualizations
library(tibble)
library(dplyr)


#### Load in all raw features
### Split by validation with HBC and then seperately based on 70/30 split 
### Run PCA on each one seperately 
### Export scree plots as needed 
### Then run classifier


### Get raw features 

# Define the directory containing the files
raw_features_dir <- "Final Data Dec 2024/Raw_features/"

# Load each file and assign it a specific variable name based on its filename
end_motifs_se <- readRDS(file.path(raw_features_dir, "Dataframe matrix for end motifs, all frags, May 2024, SE.rds"))

## Load in metadata df 
metadata_df <- readRDS("Metadata df all samples with CN classifier Dec 2024.rds")


#### Now load in methylation and determine which one is best
### First run PCA on the methylation data 
# Define input files
file_PE <- "Final Data Dec 2024/Methylation_data/Updated_Feb3/PanCancer_hyper_DMRs_combat_deseq2_normalized_count_with_blood_cell_age_sex_filtering_for_PE_samples_after_QC.RDS"
file_SE <- "Final Data Dec 2024/Methylation_data/Updated_Feb3/PanCancer_hyper_DMRs_combat_deseq2_normalized_count_with_blood_cell_age_sex_filtering_for_SE_samples_after_QC.RDS"

#file_beta <- "Final Data Dec 2024/Methylation_data/PanCancer_hyper_DMRs_qsea_beta_with_PBL_filtering_for_all_samples_after_QC.RDS"
#file_combat_adjusted <- "Final Data Dec 2024/Methylation_data/PanCancer_hyper_DMRs_combat_adjusted_count_with_PBL_filtering_for_all_samples_after_QC.RDS"

# Define output directory
outdir <- "PCA_results/Run_PCA_on_everything_SE_updated_tables/"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Load the data
data_SE <- readRDS(file_SE)




#### Start from here
#### Now split into training and validation - first for PE 
### 1) Filter metadata to remove TCGE-CFMe-HCC and 'SE' samples
metadata_filtered_SE <- metadata_df %>%
  filter(
    type == "SE"                   # Remove SE samples
  )

## Make a column for cancer vs normal 
## Edit metadata df
metadata_filtered_SE <- metadata_filtered_SE %>%
  dplyr::mutate(CN_classifier = ifelse(cancer_type_corrected_updated == "healthy", "healthy", "cancer"))

table(metadata_filtered_SE$CN_classifier)

# Set a seed for reproducibility
set.seed(123)

# This function will split the metadata into training and validation sets
# while ensuring each cancer type is proportionally represented and that
# we have at least 5 samples (or 10%) in validation for each group.
split_metadata_by_cancer <- function(metadata,
                                     cancer_type_col = "CN_classifier",
                                     proportion_validation = 0.1) {
  
  # Group by the cancer type column
  metadata_split <- metadata %>%
    group_by(!!sym(cancer_type_col)) %>%
    group_split()
  
  splitted_dfs <- lapply(metadata_split, function(df_subset) {
    # Calculate how many samples are 10% of this subset
    n_val_temp <- ceiling(nrow(df_subset) * proportion_validation)
    # Ensure at least 5, but cannot exceed the total number of samples
    n_val <- min(nrow(df_subset), max(5, n_val_temp))
    
    # Randomly pick n_val indices for validation
    val_indices <- sample(seq_len(nrow(df_subset)), size = n_val)
    
    val_subset   <- df_subset[val_indices, ]
    train_subset <- df_subset[-val_indices, ]
    
    list(train = train_subset, validation = val_subset)
  })
  
  # Combine all train subsets into one data frame and validation subsets into another
  training_dfs   <- lapply(splitted_dfs, `[[`, "train")
  validation_dfs <- lapply(splitted_dfs, `[[`, "validation")
  
  metadata_train <- bind_rows(training_dfs)
  metadata_val   <- bind_rows(validation_dfs)
  
  list(train = metadata_train, validation = metadata_val)
}

# Run on metadata
metadata_splits <- split_metadata_by_cancer(
  metadata_filtered_SE,
  cancer_type_col       = "CN_classifier",
  proportion_validation = 0.1
)

metadata_train <- metadata_splits$train
metadata_val   <- metadata_splits$validation

# Check distribution
table(metadata_train$CN_classifier)
table(metadata_val$CN_classifier)

## Check natural cancer type distribution 
table(metadata_filtered_SE$CN_classifier)

### 3) Define a helper function to split datasets
split_feature_data <- function(feature_data, row_id_col) {
  # Subset for training and validation
  train_data <- feature_data[rownames(feature_data) %in% metadata_train[[row_id_col]], , drop = FALSE]
  val_data <- feature_data[rownames(feature_data) %in% metadata_val[[row_id_col]], , drop = FALSE]
  
  # Manually add rownames back
  rownames(train_data) <- rownames(feature_data)[rownames(feature_data) %in% metadata_train[[row_id_col]]]
  rownames(val_data) <- rownames(feature_data)[rownames(feature_data) %in% metadata_val[[row_id_col]]]
  
  # Return both subsets
  list(train = train_data, validation = val_data)
}



### 4) Combine feature matrices and split into training and validation

#### End motifs datasets
end_motifs_split <- split_feature_data(end_motifs_se, "endmotif_name")

# Transpose so samples are rows 
data_SE <- t(data_SE)

#### Methylation 
methylation_split_SE <- split_feature_data(data_SE, "sequencing_id") # to do later


### 5) Final combined datasets for each feature
# Training datasets
training_end_motifs <- end_motifs_split$train
training_methylation <- methylation_split_SE$train

# Validation datasets
validation_end_motifs <- end_motifs_split$validation
validation_methylation <- methylation_split_SE$validation


### Run PCA on everything and export plots

## Set function
run_pca_and_export <- function(data_matrix, feature_name, cohort_type, outdir, thresholds = c(0.90, 0.97)) {
  # Ensure data is numeric
  if (!is.matrix(data_matrix)) {
    data_matrix <- as.matrix(data_matrix)
  }
  
  # Remove constant/zero-variance columns
  col_sd <- apply(data_matrix, 2, sd, na.rm = TRUE)  # Calculate standard deviation of columns
  data_matrix <- data_matrix[, col_sd > 0, drop = FALSE]  # Keep only columns with non-zero SD
  
  if (ncol(data_matrix) == 0) {
    stop(paste("All columns in", feature_name, cohort_type, "are constant or zero-variance."))
  }
  
  # Perform PCA
  pca_result <- prcomp(data_matrix, center = TRUE, scale. = TRUE)
  
  # Extract variance explained
  var_explained <- summary(pca_result)$importance[2, ]  # Proportion of variance explained
  cum_variance <- cumsum(var_explained)                # Cumulative variance explained
  
  # Determine the number of components for each threshold
  num_components <- sapply(thresholds, function(th) which(cum_variance >= th)[1])
  
  # Create directory for saving files
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  
  # If you want it optional, you can do:
  if (tolower(cohort_type) %in% c("training", "validation")) {
    # Drop it from the title
    title_text <- paste("Cumulative Variance -", feature_name)
  } else {
    # Keep the cohort type in the title
    title_text <- paste("Cumulative Variance -", feature_name, "(", cohort_type, ")")
  }
  
  # Plot cumulative explained variance
  cum_var_df <- tibble(
    PC = seq_along(cum_variance),
    CumulativeVariance = cum_variance
  )
  
  for (i in seq_along(thresholds)) {
    th <- thresholds[i]
    n_comp <- num_components[i]
    
    p_cum <- ggplot(cum_var_df, aes(x = PC, y = CumulativeVariance)) +
      geom_line(color = "black") +
      geom_point(color = "black") +
      geom_hline(yintercept = cum_variance[n_comp], color = "red", linetype = "dashed") +
      geom_vline(xintercept = n_comp, color = "red", linetype = "dashed") +
      labs(
        title = title_text,  # use our adjusted title
        subtitle = paste0("Selected PCs: ", n_comp, 
                          " (Cumulative: ", round(cum_variance[n_comp], 2), ")"),
        x = "Number of Principal Components",
        y = "Cumulative Explained Variance"
      ) +
      theme_classic() +
      theme(
        # Make text smaller:
        plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 8),
        axis.title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 8)
      )
    
    # Save plot with smaller dimensions
    plot_file <- file.path(
      outdir,
      paste0("Cumulative_Explained_Variance_",
             feature_name, "_",
             cohort_type, "_", th * 100, "percent.png")
    )
    ggsave(plot_file, plot = p_cum, width = 5, height = 5, dpi = 500)
  }
  
  ## Also save table 
  cumvar_file <- file.path(
    outdir,
    paste0("Cumulative_Explained_Variance_Table_", feature_name, "_", cohort_type, ".csv")
  )
  write.csv(cum_var_df, cumvar_file, row.names = FALSE)
  

  # Scree plot
  scree_plot <- fviz_eig(pca_result, addlabels = TRUE) +
    ggtitle(paste("Scree Plot for", feature_name, "(", cohort_type, ")")) +
    theme_classic() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    ) +  ylim(0, 53)  # Lock y-axis between 0 and 53%

  
  # Save scree plot
  scree_plot_file <- file.path(outdir, paste0("Scree_Plot_", feature_name, "_", cohort_type, ".png"))
  ggsave(scree_plot_file, plot = scree_plot, width = 10, height = 6, dpi = 300)
  
  # Export results for each threshold
  results <- list()
  for (i in seq_along(thresholds)) {
    th <- thresholds[i]
    n_comp <- num_components[i]
    
    # Extract the principal components up to the threshold
    pca_scores <- pca_result$x[, 1:n_comp, drop = FALSE]
    
    # Save PCA scores as CSV
    csv_file <- file.path(outdir, paste0(feature_name, "_", cohort_type, "_PCA_", th * 100, "percent.csv"))
    write.csv(as.data.frame(pca_scores), csv_file, row.names = TRUE)
    
    # Save PCA object
    rds_file <- file.path(outdir, paste0(feature_name, "_", cohort_type, "_PCA_", th * 100, "percent.rds"))
    saveRDS(pca_result, rds_file)
    
    # Store results in the list
    results[[paste0(th * 100, "%")]] <- list(
      num_components = n_comp,
      pca_scores = pca_scores,
      cumulative_variance = cum_variance[n_comp]
    )
  }
  
  # Create summary of number of PCs retained
  pc_summary <- tibble(
    Feature = feature_name,
    Cohort = cohort_type,
    Threshold = paste0(thresholds * 100, "%"),
    Num_Components = num_components,
    Cumulative_Variance = cum_variance[num_components]
  )
  
  # Save PC summary as CSV
  pc_summary_file <- file.path(outdir, paste0("PC_Summary_", feature_name, "_", cohort_type, ".csv"))
  write.csv(pc_summary, pc_summary_file, row.names = FALSE)
  
  # Return a summary
  return(list(
    pca_result = pca_result,
    var_explained = var_explained,
    cum_variance = cum_variance,
    thresholds = results,
    pc_summary = pc_summary
  ))
}


## Run function
# Define output directory
outdir <- "PCA_results/Run_PCA_on_everything_SE_updated_tables/"

# List of datasets
datasets <- list(
  training_end_motifs = list(data = training_end_motifs, cohort = "Training"),
  validation_end_motifs = list(data = validation_end_motifs, cohort = "Validation"),
  training_methylation = list(data = training_methylation, cohort = "Training"),
  validation_methylation = list(data = validation_methylation, cohort = "Validation")
)

# Run PCA for all datasets
pca_results <- list()

for (feature_name in names(datasets)) {
  dataset <- datasets[[feature_name]]
  
  cat("Running PCA for:", feature_name, "\n")
  
  # Run PCA and save results
  pca_results[[feature_name]] <- run_pca_and_export(
    data_matrix = dataset$data,
    feature_name = feature_name,
    cohort_type = dataset$cohort,
    outdir = outdir,
    thresholds = c(0.8, 0.85, 0.90, 0.97)  # Retain 90% and 97% variance
  )
}


#### Now run on train and fit on validation 

# Set function to run PCA on training and apply it to validation
run_training_pca_and_fit_validation <- function(training_data, validation_data, feature_name, outdir, thresholds = c(0.90, 0.97)) {
  # Ensure training data is numeric and has no constant columns
  if (!is.matrix(training_data)) training_data <- as.matrix(training_data)
  col_sd <- apply(training_data, 2, sd, na.rm = TRUE)
  training_data <- training_data[, col_sd > 0, drop = FALSE]
  
  if (ncol(training_data) == 0) {
    stop(paste("All columns in training data for", feature_name, "are constant or zero-variance."))
  }
  
  # Perform PCA on training data
  pca_result <- prcomp(training_data, center = TRUE, scale. = TRUE)
  
  # Calculate variance explained
  var_explained <- summary(pca_result)$importance[2, ]
  cum_variance <- cumsum(var_explained)
  
  # Determine the number of components for each threshold
  num_components <- sapply(thresholds, function(th) which(cum_variance >= th)[1])
  
  # Save training PCA scores for each threshold
  training_results <- list()
  for (i in seq_along(thresholds)) {
    th <- thresholds[i]
    n_comp <- num_components[i]
    
    # Extract the top principal components
    pca_scores_train <- pca_result$x[, 1:n_comp, drop = FALSE]
    
    # Save PCA scores for training data as both CSV and RDS
    write.csv(as.data.frame(pca_scores_train), file.path(outdir, paste0(feature_name, "_Training_PCA_", th * 100, "percent.csv")), row.names = TRUE)
    saveRDS(as.data.frame(pca_scores_train), file.path(outdir, paste0(feature_name, "_Training_PCA_", th * 100, "percent.rds")))
    
    # Save PCA object
    saveRDS(pca_result, file.path(outdir, paste0(feature_name, "_PCA_Object.rds")))
    
    # Store results
    training_results[[paste0(th * 100, "%")]] <- list(
      num_components = n_comp,
      pca_scores = pca_scores_train,
      cumulative_variance = cum_variance[n_comp]
    )
  }
  
  # Apply PCA to validation data
  if (!is.matrix(validation_data)) validation_data <- as.matrix(validation_data)
  validation_scaled <- scale(validation_data, center = pca_result$center, scale = pca_result$scale)
  validation_results <- list()
  
  for (i in seq_along(thresholds)) {
    th <- thresholds[i]
    n_comp <- num_components[i]
    
    # Project validation data onto training PCA space
    pca_scores_validation <- validation_scaled %*% pca_result$rotation[, 1:n_comp, drop = FALSE]
    
    # Save PCA scores for validation data as both CSV and RDS
    write.csv(as.data.frame(pca_scores_validation), file.path(outdir, paste0(feature_name, "_Validation_PCA_", th * 100, "percent.csv")), row.names = TRUE)
    saveRDS(as.data.frame(pca_scores_validation), file.path(outdir, paste0(feature_name, "_Validation_PCA_", th * 100, "percent.rds")))
    
    # Store results
    validation_results[[paste0(th * 100, "%")]] <- list(
      num_components = n_comp,
      pca_scores = pca_scores_validation
    )
  }
  
  # Export summary of retained components
  pc_summary <- tibble(
    Feature = feature_name,
    Threshold = paste0(thresholds * 100, "%"),
    Num_Components = num_components,
    Cumulative_Variance = cum_variance[num_components]
  )
  
  # Save PC summary as both CSV and RDS
  write.csv(pc_summary, file.path(outdir, paste0("PC_Summary_", feature_name, ".csv")), row.names = FALSE)
  saveRDS(pc_summary, file.path(outdir, paste0("PC_Summary_", feature_name, ".rds")))
  
  # Return PCA results for training and validation
  return(list(training = training_results, validation = validation_results, summary = pc_summary))
}



# Define output directory
outdir <- "PCA_results/Run_PCA_on_train_fit_on_validation_SE/"

# Run PCA on training and apply to validation for each dataset
for (feature_name in names(datasets)) {
  dataset <- datasets[[feature_name]]
  
  if (grepl("training_", feature_name)) {
    cat("Running PCA for:", feature_name, "\n")
    
    # Determine corresponding validation data
    validation_feature_name <- gsub("training_", "validation_", feature_name)
    validation_data <- datasets[[validation_feature_name]]$data
    
    # Run PCA on training and fit to validation
    run_training_pca_and_fit_validation(
      training_data = dataset$data,
      validation_data = validation_data,
      feature_name = gsub("training_", "", feature_name),
      outdir = outdir,
      thresholds = c(0.8, 0.85, 0.90, 0.97)  # Retain 90% and 97% variance
    )
  }
}


#### Now integrate the components 
#------------------------------------------------------------------------------
# 1) Define helper function to filter PCs by variance (based on training PCA)
#------------------------------------------------------------------------------
filter_pcs_by_variance <- function(pca_object, pc_scores, threshold = 1) {
  # pca_object: a prcomp object from the training set
  # pc_scores:  the matrix/data.frame of PC scores (training or validation)
  # threshold:  how much % variance a PC must have to be retained
  
  # Calculate the percentage of variance explained by each PC
  explained_variance <- pca_object$sdev^2 / sum(pca_object$sdev^2) * 100
  
  # Select the PCs that explain more than 'threshold'%
  selected_pcs <- which(explained_variance > threshold)
  
  # Subset the pc_scores accordingly
  pcs_sub <- as.data.frame(pc_scores[, selected_pcs, drop = FALSE])
  
  # Add sample_id as a column for subsequent merges
  pcs_sub <- pcs_sub %>%
    rownames_to_column(var = "sample_id")
  
  return(pcs_sub)
}


#------------------------------------------------------------------------------
#  2) Now process. For each feature, load:
#    - The training PCA object (e.g. "end_motifs_PCA_Object.rds")
#    - The training PC scores file (e.g. "end_motifs_Training_PCA_90percent.rds")
#    - The validation PC scores file (e.g. "end_motifs_Validation_PCA_90percent.rds")
#    Then create filtered data frames for both training and validation at >1% & >5%.
#------------------------------------------------------------------------------


## End_motifs
pca_obj_end_motifs <- readRDS("PCA_results/Run_PCA_on_train_fit_on_validation_SE/end_motifs_PCA_Object.rds")
train_scores_end_motifs <- readRDS("PCA_results/Run_PCA_on_train_fit_on_validation_SE/end_motifs_Training_PCA_97percent.rds")
val_scores_end_motifs   <- readRDS("PCA_results/Run_PCA_on_train_fit_on_validation_SE/end_motifs_Validation_PCA_97percent.rds")

training_end_motifs_5pct <- filter_pcs_by_variance(pca_obj_end_motifs, train_scores_end_motifs, threshold = 5)
training_end_motifs_1pct <- filter_pcs_by_variance(pca_obj_end_motifs, train_scores_end_motifs, threshold = 1)

validation_end_motifs_5pct <- filter_pcs_by_variance(pca_obj_end_motifs, val_scores_end_motifs, threshold = 5)
validation_end_motifs_1pct <- filter_pcs_by_variance(pca_obj_end_motifs, val_scores_end_motifs, threshold = 1)

# Fix sample IDs (to facilitate merges)
training_end_motifs_5pct$sample_id <- gsub("_motifs", "", training_end_motifs_5pct$sample_id)
training_end_motifs_1pct$sample_id <- gsub("_motifs", "", training_end_motifs_1pct$sample_id)
validation_end_motifs_5pct$sample_id <- gsub("_motifs", "", validation_end_motifs_5pct$sample_id)
validation_end_motifs_1pct$sample_id <- gsub("_motifs", "", validation_end_motifs_1pct$sample_id)



## Methylation
pca_obj_methylation <- readRDS("PCA_results/Run_PCA_on_train_fit_on_validation_SE/methylation_PCA_Object.rds")
train_scores_methylation <- readRDS("PCA_results/Run_PCA_on_train_fit_on_validation_SE/methylation_Training_PCA_97percent.rds")
val_scores_methylation   <- readRDS("PCA_results/Run_PCA_on_train_fit_on_validation_SE/methylation_Validation_PCA_97percent.rds")

training_methylation_5pct <- filter_pcs_by_variance(pca_obj_methylation, train_scores_methylation, threshold = 5)
training_methylation_1pct <- filter_pcs_by_variance(pca_obj_methylation, train_scores_methylation, threshold = 1)

validation_methylation_5pct <- filter_pcs_by_variance(pca_obj_methylation, val_scores_methylation, threshold = 5)
validation_methylation_1pct <- filter_pcs_by_variance(pca_obj_methylation, val_scores_methylation, threshold = 1)

## Ajust naming 
training_methylation_5pct$sample_id <- paste0(training_methylation_5pct$sample_id, "_dedup")
training_methylation_1pct$sample_id <- paste0(training_methylation_1pct$sample_id, "_dedup")
validation_methylation_5pct$sample_id <- paste0(validation_methylation_5pct$sample_id, "_dedup")
validation_methylation_1pct$sample_id <- paste0(validation_methylation_1pct$sample_id, "_dedup")


## Nucleosome_peak
pca_obj_nucleosome <- readRDS("PCA_results/Run_PCA_on_train_fit_on_validation_SE/nucleosome_peak_PCA_Object.rds")
train_scores_nucleosome <- readRDS("PCA_results/Run_PCA_on_train_fit_on_validation_SE/nucleosome_peak_Training_PCA_97percent.rds")
val_scores_nucleosome   <- readRDS("PCA_results/Run_PCA_on_train_fit_on_validation_SE/nucleosome_peak_Validation_PCA_97percent.rds")

training_nucleosome_5pct <- filter_pcs_by_variance(pca_obj_nucleosome, train_scores_nucleosome, threshold = 5)
training_nucleosome_1pct <- filter_pcs_by_variance(pca_obj_nucleosome, train_scores_nucleosome, threshold = 1)

validation_nucleosome_5pct <- filter_pcs_by_variance(pca_obj_nucleosome, val_scores_nucleosome, threshold = 5)
validation_nucleosome_1pct <- filter_pcs_by_variance(pca_obj_nucleosome, val_scores_nucleosome, threshold = 1)

# Fix sample IDs
training_nucleosome_5pct$sample_id <- gsub("_peak_distance", "", training_nucleosome_5pct$sample_id)
training_nucleosome_1pct$sample_id <- gsub("_peak_distance", "", training_nucleosome_1pct$sample_id)
validation_nucleosome_5pct$sample_id <- gsub("_peak_distance", "", validation_nucleosome_5pct$sample_id)
validation_nucleosome_1pct$sample_id <- gsub("_peak_distance", "", validation_nucleosome_1pct$sample_id)

#------------------------------------------------------------------------------
# 3) Now merge the >5% variance training sets across features
#------------------------------------------------------------------------------
training_combined_5pct <- training_end_motifs_5pct %>%
  full_join(training_methylation_5pct, by = "sample_id")

# Convert sample_id to rownames
rownames(training_combined_5pct) <- training_combined_5pct$sample_id
training_combined_5pct <- training_combined_5pct %>% dplyr::select(-sample_id)

# Save the >5% variance training data
saveRDS(training_combined_5pct, "PCA_results/Run_PCA_on_train_fit_on_validation_SE/Training_combined_5pct_PCs.rds")


#------------------------------------------------------------------------------
# 4) Merge the >1% variance training sets across features
#------------------------------------------------------------------------------
training_combined_1pct <- training_end_motifs_1pct %>%
  full_join(training_methylation_1pct, by = "sample_id")

rownames(training_combined_1pct) <- training_combined_1pct$sample_id
training_combined_1pct <- training_combined_1pct %>% dplyr::select(-sample_id)

# Save the >1% variance training data
saveRDS(training_combined_1pct, "PCA_results/Run_PCA_on_train_fit_on_validation_SE/Training_combined_1pct_PCs.rds")


#------------------------------------------------------------------------------
# 5) Merge the >5% variance validation sets across features
#------------------------------------------------------------------------------
validation_combined_5pct <- validation_end_motifs_5pct %>%
  full_join(validation_methylation_5pct, by = "sample_id")
rownames(validation_combined_5pct) <- validation_combined_5pct$sample_id
validation_combined_5pct <- validation_combined_5pct %>% dplyr::select(-sample_id)

saveRDS(validation_combined_5pct, "PCA_results/Run_PCA_on_train_fit_on_validation_SE/Validation_combined_5pct_PCs.rds")


#------------------------------------------------------------------------------
# 6) Merge the >1% variance validation sets across features
#------------------------------------------------------------------------------
validation_combined_1pct <- validation_end_motifs_1pct %>%
  full_join(validation_methylation_1pct, by = "sample_id") 

rownames(validation_combined_1pct) <- validation_combined_1pct$sample_id
validation_combined_1pct <- validation_combined_1pct %>% dplyr::select(-sample_id)

saveRDS(validation_combined_1pct, "PCA_results/Run_PCA_on_train_fit_on_validation_SE/Validation_combined_1pct_PCs.rds")


#------------------------------------------------------------------------------
# 7) Example: Merging only End Motifs + Methylation (>1% Variance) for Training
#------------------------------------------------------------------------------
combined_end_motifs_methylation_1pct_training <- training_end_motifs_1pct %>%
  full_join(training_methylation_1pct, by = "sample_id")

rownames(combined_end_motifs_methylation_1pct_training) <- combined_end_motifs_methylation_1pct_training$sample_id
combined_end_motifs_methylation_1pct_training <- combined_end_motifs_methylation_1pct_training %>% dplyr::select(-sample_id)

saveRDS(combined_end_motifs_methylation_1pct_training, 
        "PCA_results/Run_PCA_on_train_fit_on_validation_SE/Combined_EndMotif_Methylation_1pct_PCs_training.rds")


#------------------------------------------------------------------------------
# 8) Example: Merging only End Motifs + Methylation (>1% Variance) for Validation
#------------------------------------------------------------------------------
combined_end_motifs_methylation_1pct_validation <- validation_end_motifs_1pct %>%
  full_join(validation_methylation_1pct, by = "sample_id")

rownames(combined_end_motifs_methylation_1pct_validation) <- combined_end_motifs_methylation_1pct_validation$sample_id
combined_end_motifs_methylation_1pct_validation <- combined_end_motifs_methylation_1pct_validation %>% dplyr::select(-sample_id)

saveRDS(combined_end_motifs_methylation_1pct_validation, 
        "PCA_results/Run_PCA_on_train_fit_on_validation_SE/Combined_EndMotif_Methylation_1pct_PCs_validation.rds")




#### Now check for NA rows in any of the exported dataframes 
# Define a helper function to check for NA rows in a dataframe
check_na_rows <- function(df_name, df) {
  if (!is.data.frame(df)) {
    return(NULL)
  }
  na_rows <- sum(apply(df, 1, function(row) any(is.na(row))))
  return(data.frame(
    DataFrameName = df_name,
    TotalRows = nrow(df),
    RowsWithNA = na_rows,
    PercentageNA = round((na_rows / nrow(df)) * 100, 2)
  ))
}

# List of dataframes to check
dataframes_to_check <- list(
  "pca_obj_end_motifs" = pca_obj_end_motifs,
  "train_scores_end_motifs" = train_scores_end_motifs,
  "val_scores_end_motifs" = val_scores_end_motifs,
  "training_end_motifs_5pct" = training_end_motifs_5pct,
  "training_end_motifs_1pct" = training_end_motifs_1pct,
  "validation_end_motifs_5pct" = validation_end_motifs_5pct,
  "validation_end_motifs_1pct" = validation_end_motifs_1pct,
  "pca_obj_methylation" = pca_obj_methylation,
  "train_scores_methylation" = train_scores_methylation,
  "val_scores_methylation" = val_scores_methylation,
  "training_methylation_5pct" = training_methylation_5pct,
  "training_methylation_1pct" = training_methylation_1pct,
  "validation_methylation_5pct" = validation_methylation_5pct,
  "validation_methylation_1pct" = validation_methylation_1pct
)

# Check each dataframe and summarize results
na_summary <- do.call(rbind, lapply(names(dataframes_to_check), function(df_name) {
  df <- dataframes_to_check[[df_name]]
  check_na_rows(df_name, df)
}))

# Filter out NULLs from results
na_summary <- na_summary[!is.na(na_summary$RowsWithNA), ]

# Save the summary as a table
na_summary_path <- "PCA_results/Run_PCA_on_train_fit_on_validation_SE/NA_Summary.csv"
write.csv(na_summary, na_summary_path, row.names = FALSE)

# Print the summary
print(na_summary)
message("NA summary saved to: ", na_summary_path)




### Save
saveRDS(metadata_filtered_SE, file = "Metadata_df_SE_CN_classifier_Dec_2024.rds")





### Now do this for the SE samples 

## Load in methylation and end motifs for the SE 
## Split to leave 10% for validation 





### Lastly redo by cancer type and subtype using the same cohort


###tmp 
tmp <- readRDS("PCA_results/Run_PCA_on_train_fit_on_validation_SE/Combined_EndMotif_Methylation_1pct_PCs_training.rds")
