###############################################################################
# File:      6.2A_Selecting_samples_for_training_and_validation_cohorts_PE_data_updatedJan_PCA_transform.R
# Author:    Dory Abelman
# Date:      Dec 2024, Updated March 2025
#
# Purpose:
#   1. Load & harmonize raw feature matrices for PE (paired-end) and SE (single-end). SE is for exploratory purposes only in this script and not required for script to run:
#        - End‑motifs, nucleosome peak distances, DELFI fragment ratios, insert‑size,
#          and combat‑seq adjusted methylation counts.
#   2. Load, clean & augment sample metadata:
#        - Fill missing sequencing `type`, normalize cancer subtype labels,
#          derive binary cancer vs normal (CN) classifier labels.
#   3. Sanity checks & cohort splitting:
#        - Compare feature distributions between training (90%) and hold‐out (10%) 
#          sets via feature‑wise t‑tests with BH correction.
#        - Aggregate cohort means, plot boxplots with jitter, and run two‐sample tests.
#   4. PCA for dimensionality reduction:
#        - Run PCA independently on each feature set (PE/SE methylation & fragmentomics).
#        - Generate scree & cumulative variance plots, select PCs to explain X% variance.
#        - Export PCA objects, selected PC scores, and summary tables.
#   5. Train/validation split & alignment:
#        - Define training vs. validation samples based on project cohorts (INSPIRE, LFS, HCC).
#        - Align feature matrices by sample IDs, split into matched training/validation.
#   6. Export processed data for downstream ML:
#        - Save raw and PCA‑reduced training/validation matrices as RDS/CSV.
#        - Provide summary counts of retained PCs and NA diagnostics.
#
# Usage:
#   Source this script before running classification pipelines (e.g., random forest,
#   logistic regression) on the PCA‑reduced feature sets.
###############################################################################


## Prepare data for ML classification - Dory Abelman
## Splits data into 90/10 for PE data 

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
metadata_df <- readRDS("Metadata_df_all_with_corrected_cancer_subtype_Jan2025.rds")

## Fix error 
metadata_df <- metadata_df %>%
  mutate(type = ifelse(is.na(type), "PE", type))

saveRDS(metadata_df, file = "Metadata_df_all_samples_with_CN_classifier_Jan_2025.rds")

data_dir <- "Final Data Dec 2024/Raw_features/"

# Load the datasets
matrix_nucleosome_peak <- readRDS(file.path(data_dir, "Matrix for nucleosome peak.rds"))
matrix_nucleosome_peak_validation <- readRDS(file.path(data_dir, "Matrix for nucleosome peak validation.rds"))

delfi_ratios <- readRDS(file.path(data_dir, "DELFI_Ratios_updated_Dec2024.rds"))
delfi_ratios_validation <- readRDS(file.path(data_dir, "DELFI_Ratios_updated_Dec2024_validation_cohort.rds"))

insert_matrix <- readRDS(file.path(data_dir, "insert_matrix_df_for_Althaf_updated.rds"))
insert_matrix_validation <- readRDS(file.path(data_dir, "insert_matrix_df_validation_for_Althaf_updated.rds"))

end_motifs_all_frags <- readRDS(file.path(data_dir, "Dataframe matrix for end motifs, all frags, May 2024.rds"))
end_motifs_validation <- readRDS(file.path(data_dir, "Dataframe matrix for end motifs, all frags, May 2024, validation.rds"))

### Check correct number of rows 
# Check row counts for non-validation datasets
non_val <- list(
  matrix_nucleosome_peak = matrix_nucleosome_peak,
  delfi_ratios = delfi_ratios,
  insert_matrix = insert_matrix,
  end_motifs_all_frags = end_motifs_all_frags
)
non_val_rows <- sapply(non_val, nrow)
print("Non-validation row counts:")
print(non_val_rows)
print("Unique row counts in non-validation datasets:")
print(unique(non_val_rows))

# Check row counts for validation datasets
val <- list(
  matrix_nucleosome_peak_validation = matrix_nucleosome_peak_validation,
  delfi_ratios_validation = delfi_ratios_validation,
  insert_matrix_validation = insert_matrix_validation,
  end_motifs_validation = end_motifs_validation
)
val_rows <- sapply(val, nrow)
print("Validation row counts:")
print(val_rows)
print("Unique row counts in validation datasets:")
print(unique(val_rows))


# Compare each training vs. validation pair using t-tests.
# Define a helper function to compare two matrices (training vs. validation) feature-by-feature.
compare_matrices <- function(train, validation, test = "t.test") {
  # Ensure that the matrices have the same column names
  features <- intersect(colnames(train), colnames(validation))
  results <- data.frame(feature = features, p_value = NA, stringsAsFactors = FALSE)
  for(feat in features){
    # Extract feature values and remove NAs
    train_feat <- train[, feat]
    valid_feat <- validation[, feat]
    train_feat <- train_feat[!is.na(train_feat)]
    valid_feat <- valid_feat[!is.na(valid_feat)]
    # Perform test: t-test or Wilcoxon (if non-normal)
    if(test == "t.test"){
      test_res <- t.test(train_feat, valid_feat)
    } else if(test == "wilcox"){
      test_res <- wilcox.test(train_feat, valid_feat)
    }
    results$p_value[results$feature == feat] <- test_res$p.value
  }
  # Adjust p-values using Benjamini-Hochberg (BH) correction.
  results$adj_p <- p.adjust(results$p_value, method = "BH")
  return(results)
}


nucleosome_diff <- compare_matrices(matrix_nucleosome_peak, matrix_nucleosome_peak_validation, test = "t.test")
delfi_diff <- compare_matrices(delfi_ratios, delfi_ratios_validation, test = "t.test")
insert_diff <- compare_matrices(insert_matrix, insert_matrix_validation, test = "t.test")
end_motifs_diff <- compare_matrices(end_motifs_all_frags, end_motifs_validation, test = "t.test")

# Optionally, print the results
print(nucleosome_diff)
print(delfi_diff)
print(insert_diff)
print(end_motifs_diff)

# Define a helper function to compute aggregate values, plot them, and run a t-test.
compare_aggregate <- function(train_matrix, valid_matrix, dataset_name) {
  # Compute the row means for each sample (aggregate measure)
  train_means <- rowMeans(train_matrix, na.rm = TRUE)
  valid_means <- rowMeans(valid_matrix, na.rm = TRUE)
  
  # Combine into a data frame
  df <- data.frame(
    MeanValue = c(train_means, valid_means),
    Cohort = rep(c("Training", "Validation"), 
                 times = c(length(train_means), length(valid_means)))
  )
  
  # Create a boxplot with jittered points
  p <- ggplot(df, aes(x = Cohort, y = MeanValue, fill = Cohort)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "darkblue") +
    ggtitle(paste("Aggregate Feature Values for", dataset_name)) +
    xlab("") +
    ylab("Row Mean (Aggregate Value)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  print(p)
  ggsave(paste0(dataset_name, "_Aggregate_Boxplot.pdf"), p, width = 6, height = 4)
  
  # Perform a two-sample t-test between training and validation aggregate values
  ttest_res <- t.test(train_means, valid_means)
  print(ttest_res)
  
  # Return both the plot and test result
  return(list(plot = p, ttest = ttest_res))
}

# Run the comparison for each dataset
res_nucleosome <- compare_aggregate(matrix_nucleosome_peak, matrix_nucleosome_peak_validation, "Nucleosome Peak")
res_delfi <- compare_aggregate(delfi_ratios, delfi_ratios_validation, "DELFI Ratios")
res_insert <- compare_aggregate(insert_matrix, insert_matrix_validation, "Insert Matrix")
res_endmotifs <- compare_aggregate(end_motifs_all_frags, end_motifs_validation, "End Motifs")



#### Now load in methylation and determine which one is best
### First run PCA on the methylation data 
# Define input files
file_PE <- "Final Data Dec 2024/Methylation_data/Updated_Feb3/PanCancer_hyper_DMRs_combat_deseq2_normalized_count_with_blood_cell_age_sex_filtering_for_PE_samples_after_QC.RDS"
file_SE <- "Final Data Dec 2024/Methylation_data/Updated_Feb3/PanCancer_hyper_DMRs_combat_deseq2_normalized_count_with_blood_cell_age_sex_filtering_for_SE_samples_after_QC.RDS"

#file_beta <- "Final Data Dec 2024/Methylation_data/PanCancer_hyper_DMRs_qsea_beta_with_PBL_filtering_for_all_samples_after_QC.RDS"
#file_combat_adjusted <- "Final Data Dec 2024/Methylation_data/PanCancer_hyper_DMRs_combat_adjusted_count_with_PBL_filtering_for_all_samples_after_QC.RDS"

# Define output directory
outdir <- "PCA_results/Updated Mar 2025/PE"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Load the data
data_PE <- readRDS(file_PE)
data_SE <- readRDS(file_SE)
#data_beta <- readRDS(file_beta)
#data_combat_adjusted <- readRDS(file_combat_adjusted)

### Load in the combat for the validaiton cohort 
file_validation <- "Final Data Dec 2024/Methylation_data/Updated_Feb3/PanCancer_hyper_DMRs_combat_deseq2_normalized_count_with_blood_cell_age_sex_filtering_for_validation_samples_after_QC_vialA_only.RDS"
data_validation <- readRDS(file_validation)


# Function for PCA analysis and feature selection
run_pca_with_feature_selection <- function(data_matrix, feature_name, target_components = 20:50, outdir = ".") {
  # Ensure data is numeric and scaled
  if (!is.matrix(data_matrix)) data_matrix <- as.matrix(data_matrix)
  pca_result <- prcomp(data_matrix, center = TRUE, scale. = TRUE)
  
  # Variance explained
  var_explained <- summary(pca_result)$importance[2, ]  # Proportion of variance
  cum_variance  <- cumsum(var_explained)                # Cumulative variance
  
  # Determine number of components in the target range explaining significant variance
  num_components <- target_components[
    which(cum_variance[target_components] == min(cum_variance[target_components]))
  ]
  
  # Plot cumulative variance
  cum_var_df <- tibble(
    PC = seq_along(cum_variance),
    CumulativeVariance = cum_variance
  )
  
  p_cum <- ggplot(cum_var_df, aes(x = PC, y = CumulativeVariance)) +
    geom_line(color = "black") +
    geom_point(color = "black") +
    geom_hline(yintercept = cum_variance[num_components], color = "red", linetype = "dashed") +
    geom_vline(xintercept = num_components, color = "red", linetype = "dashed") +
    labs(
      title = paste("Cumulative Variance -", feature_name),
      subtitle = paste0("Selected PCs: ", num_components, " (Cumulative: ", round(cum_variance[num_components], 2), ")"),
      x = "Number of Principal Components",
      y = "Cumulative Explained Variance"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10)
    )
  
  # Save the plot
  plot_file <- file.path(outdir, paste0("Cumulative_Explained_Variance_", feature_name, ".png"))
  ggsave(plot_file, plot = p_cum, width = 7, height = 5, dpi = 300)
  
  # Extract principal components
  selected_pcs <- pca_result$x[, 1:num_components, drop = FALSE]
  selected_pcs_df <- as.data.frame(selected_pcs)
  
  # Save selected PCs to file
  pc_csv <- file.path(outdir, paste0("Selected_PCs_", feature_name, "_", num_components, ".csv"))
  write.csv(selected_pcs_df, pc_csv, row.names = TRUE)
  
  # Save PCA object for reproducibility
  saveRDS(pca_result, file = file.path(outdir, paste0("PCA_Object_", feature_name, ".rds")))
  
  # Return PCA summary
  list(
    pca_result = pca_result,
    num_components = num_components,
    var_explained = var_explained,
    cum_variance = cum_variance,
    selected_pcs = selected_pcs_df
  )
}

## transpose 
data_PE <- t(data_PE)
data_SE <- t(data_SE)

# Run PCA for PE dataset
pca_PE <- run_pca_with_feature_selection(
  data_matrix = data_PE,
  feature_name = "PanCancer_PE",
  target_components = 20:50,  # Retain 20-50 components
  outdir = outdir
)

# Run PCA for SE dataset
pca_SE <- run_pca_with_feature_selection(
  data_matrix = data_SE,
  feature_name = "PanCancer_SE",
  target_components = 20:50,  # Retain 20-50 components
  outdir = outdir
)

# Run PCA for beta dataset

# Function to remove constant or non-numeric columns
remove_constant_columns <- function(data_matrix) {
  # Ensure the input is numeric
  if (!is.matrix(data_matrix)) {
    data_matrix <- data.matrix(data_matrix)  # Convert to a numeric matrix
  }
  
  # Calculate standard deviation for each column
  col_sd <- apply(data_matrix, 2, sd, na.rm = TRUE)
  
  # Keep only columns with non-zero standard deviation
  non_constant_cols <- col_sd > 0
  
  if (sum(non_constant_cols) == 0) {
    stop("All columns are constant or non-numeric!")
  }
  
  # Return cleaned matrix
  return(data_matrix[, non_constant_cols, drop = FALSE])
}

# Apply the function to your dataset
data_beta <- t(data_beta)

# Replace NA with 0
data_beta[is.na(data_beta)] <- 0

data_beta_clean <- remove_constant_columns(data_beta)

# Check if the cleaning worked
cat("Original dimensions:", dim(data_beta), "\n")
cat("Cleaned dimensions:", dim(data_beta_clean), "\n")

# Rerun PCA after cleaning
pca_beta <- run_pca_with_feature_selection(
  data_matrix = data_beta_clean,
  feature_name = "PanCancer_Beta",
  target_components = 20:50,  # Retain 20-50 components
  outdir = outdir
)

# Apply the function to your dataset
data_combat_adjusted <- t(data_combat_adjusted)

# Rerun PCA after cleaning
pca_combat_adjusted <- run_pca_with_feature_selection(
  data_matrix = data_combat_adjusted,
  feature_name = "PanCancer_combat_adjusted",
  target_components = 20:50,  # Retain 20-50 components
  outdir = outdir
)


### The beta values seem to be the best, so use this one for now 
### Yong prefers the combat-seq adjusted one, seperate for PE, SE and validation








#### Now split into training and validation - first for PE 
### 1) Filter metadata to remove TCGE-CFMe-HCC and 'SE' samples
metadata_df$type[is.na(metadata_df$type)] <- "PE"

metadata_filtered <- metadata_df %>%
  filter(
    !(project_id == "TCGE-CFMe-HCC" & cancer_type == "Liver Cancer"),  # Remove HCC samples only if cancer_type is Liver Cancer
    type != "SE"  # Remove SE samples
  )


# 2. Create the validation set from TCGE-CFMe-INSPIRE and TCGE-CFMe-LFS,
#    but exclude LFS Survivor within TCGE-CFMe-LFS.
metadata_val <- metadata_filtered %>%
  filter(project_id %in% c("TCGE-CFMe-INSPIRE", "TCGE-CFMe-LFS", "TCGE-CFMe-HCC"))

# 3. Create the training set from everything else (i.e., not in validation),
#    also excluding any leftover from the cohorts we said not to include (HCC, INSPIRE, LFS).
metadata_train <- metadata_filtered %>%
  filter(!project_id %in% c("TCGE-CFMe-INSPIRE", "TCGE-CFMe-LFS", "TCGE-CFMe-HCC"))

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

#### Nucleosome peak datasets
matrix_nucleosome_peak_combined <- rbind(matrix_nucleosome_peak, matrix_nucleosome_peak_validation)
nucleosome_peak_split <- split_feature_data(matrix_nucleosome_peak_combined, "nucleosome_peak_sampleid")

#### DELFI ratios datasets
delfi_ratios_combined <- rbind(delfi_ratios, delfi_ratios_validation)
delfi_ratios_split <- split_feature_data(delfi_ratios_combined, "sample_id")

#### Insert size datasets
insert_matrix_combined <- rbind(insert_matrix, insert_matrix_validation)
insert_matrix_split <- split_feature_data(insert_matrix_combined, "sample_id")

#### End motifs datasets
end_motifs_combined <- rbind(end_motifs_all_frags, end_motifs_validation)
end_motifs_split <- split_feature_data(end_motifs_combined, "endmotif_name")


##### Methylation dataset 

### First add the validation to the PE 
data_validation <- t(data_validation)
#data_PE <- t(data_PE) # ensure samples are rows

# Check if the column names match
if (!identical(colnames(data_PE), colnames(data_validation))) {
  cat("Columns differ between data_PE and data_validation\n")
  
  # Identify columns in data_PE but not in data_validation
  columns_in_PE_not_in_validation <- setdiff(colnames(data_PE), colnames(data_validation))
  if (length(columns_in_PE_not_in_validation) > 0) {
    cat("Columns in data_PE but not in data_validation:\n")
    print(columns_in_PE_not_in_validation)
  } else {
    cat("No extra columns in data_PE.\n")
  }
  
  # Identify columns in data_validation but not in data_PE
  columns_in_validation_not_in_PE <- setdiff(colnames(data_validation), colnames(data_PE))
  if (length(columns_in_validation_not_in_PE) > 0) {
    cat("Columns in data_validation but not in data_PE:\n")
    print(columns_in_validation_not_in_PE)
  } else {
    cat("No extra columns in data_validation.\n")
  }
} else {
  cat("Column names match between data_PE and data_validation.\n")
}

# Ensure both datasets have the same columns
# Combine the two data frames
# Create a union of column names
all_columns <- union(colnames(data_PE), colnames(data_validation))

# Align both matrices by adding missing columns filled with NA
data_PE_aligned <- matrix(NA, nrow = nrow(data_PE), ncol = length(all_columns),
                          dimnames = list(rownames(data_PE), all_columns))
data_PE_aligned[, colnames(data_PE)] <- data_PE

data_validation_aligned <- matrix(NA, nrow = nrow(data_validation), ncol = length(all_columns),
                                  dimnames = list(rownames(data_validation), all_columns))
data_validation_aligned[, colnames(data_validation)] <- data_validation

# Bind the matrices row-wise
data_methylation_combined <- rbind(data_PE_aligned, data_validation_aligned)
rm(data_validation_aligned)
rm(data_PE_aligned)

methylation_split_PE <- split_feature_data(data_methylation_combined, "sequencing_id")
#methylation_split_SE <- split_feature_data(data_SE, "sample_name") to do later

### 5) Final combined datasets for each feature
# Training datasets
training_nucleosome_peak <- nucleosome_peak_split$train
training_delfi_ratios <- delfi_ratios_split$train
training_insert_matrix <- insert_matrix_split$train
training_end_motifs <- end_motifs_split$train
training_methylation <- methylation_split_PE$train

# Validation datasets
validation_nucleosome_peak <- nucleosome_peak_split$validation
validation_delfi_ratios <- delfi_ratios_split$validation
validation_insert_matrix <- insert_matrix_split$validation
validation_end_motifs <- end_motifs_split$validation
validation_methylation <- methylation_split_PE$validation


### Run PCA on everything and export plots

## Set function
run_pca_and_export <- function(data_matrix, feature_name, cohort_type, outdir, thresholds = c(0.90, 0.97)) {
  # Ensure data is numeric
  if (!is.matrix(data_matrix)) {
    data_matrix <- as.matrix(data_matrix)
  }
  
  # Remove constant/zero-variance columns
  col_sd <- apply(data_matrix, 2, sd, na.rm = TRUE)  # Calculate standard deviation of columns
 # data_matrix <- data_matrix[, col_sd > 0, drop = FALSE]  # Keep only columns with non-zero SD
  data_matrix <- data_matrix[, !is.na(col_sd) & col_sd > 0, drop = FALSE] # Keep only columns with non-zero SD
  
  
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
  
  # Export cumulative variance plot
  cum_var_df <- tibble(
    PC = seq_along(cum_variance),
    CumulativeVariance = cum_variance
  )
  
  # Plot cumulative explained variance
  for (i in seq_along(thresholds)) {
    th <- thresholds[i]
    n_comp <- num_components[i]
    
    # Plot cumulative variance
    p_cum <- ggplot(cum_var_df, aes(x = PC, y = CumulativeVariance)) +
      geom_line(color = "black") +
      geom_point(color = "black") +
      geom_hline(yintercept = cum_variance[n_comp], color = "red", linetype = "dashed") +
      geom_vline(xintercept = n_comp, color = "red", linetype = "dashed") +
      labs(
        title = paste("Cumulative Variance -", feature_name, "(", cohort_type, ")"),
        subtitle = paste0("Selected PCs: ", n_comp, " (Cumulative: ", round(cum_variance[n_comp], 2), ")"),
        x = "Number of Principal Components",
        y = "Cumulative Explained Variance"
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10)
      )
    
    # Save cumulative variance plot
    plot_file <- file.path(outdir, paste0("Cumulative_Explained_Variance_", feature_name, "_", cohort_type, "_", th * 100, "percent.png"))
    ggsave(plot_file, plot = p_cum, width = 7, height = 5, dpi = 300)
  }
  
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
   # write.csv(as.data.frame(pca_scores), csv_file, row.names = TRUE)
    
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



### Updated code 
run_pca_and_export <- function(data_matrix, feature_name, cohort_type, outdir, thresholds = c(0.90, 0.97)) {
  # Ensure data is numeric
  if (!is.matrix(data_matrix)) {
    data_matrix <- as.matrix(data_matrix)
  }
  
  # Remove constant/zero-variance columns
  col_sd <- apply(data_matrix, 2, sd, na.rm = TRUE) 
  data_matrix <- data_matrix[, !is.na(col_sd) & col_sd > 0, drop = FALSE]
  
  if (ncol(data_matrix) == 0) {
    stop(paste("All columns in", feature_name, cohort_type, "are constant or zero-variance."))
  }
  
  # Perform PCA
  pca_result <- prcomp(data_matrix, center = TRUE, scale. = TRUE)
  
  # Extract variance explained
  var_explained <- summary(pca_result)$importance[2, ]  # Proportion of variance explained
  cum_variance <- cumsum(var_explained)                 # Cumulative variance explained
  
  # Determine the number of components for each threshold
  num_components <- sapply(thresholds, function(th) which(cum_variance >= th)[1])
  
  # Create directory for saving files
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  # --- Helper to rename "delfi_ratios" -> "fragment_ratios" in feature_name
  # Apply both substitutions on feature_name:
  pretty_feature_name <- gsub("insert_matrix", "insert_size",
                              gsub("delfi_ratios", "fragment_ratios", feature_name))
  
  # --- Optionally remove references to "training" or "validation" from the title:
  #     For instance, if your cohort_type is exactly "Training" or "Validation":
  #     If you truly don't want it in the title, just remove it:
  # title_text <- paste("Cumulative Variance -", pretty_feature_name)
  #
  # If you want it optional, you can do:
  if (tolower(cohort_type) %in% c("training", "validation")) {
    # Drop it from the title
    title_text <- paste("Cumulative Variance -", pretty_feature_name)
  } else {
    # Keep the cohort type in the title
    title_text <- paste("Cumulative Variance -", pretty_feature_name, "(", cohort_type, ")")
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
             pretty_feature_name, "_",
             cohort_type, "_", th * 100, "percent.png")
    )
    ggsave(plot_file, plot = p_cum, width = 4, height = 5, dpi = 500)
  }
  
  ## Also save table 
  cumvar_file <- file.path(
    outdir,
    paste0("Cumulative_Explained_Variance_Table_", pretty_feature_name, "_", cohort_type, ".csv")
  )
  write.csv(cum_var_df, cumvar_file, row.names = FALSE)
  
  # Scree plot
  scree_plot <- fviz_eig(pca_result, addlabels = TRUE) +
    ggtitle(paste("Scree Plot for", pretty_feature_name)) +  # or remove 'cohort_type'
    theme_classic() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 10, face = "bold"),
      axis.text = element_text(size = 8),
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
    ) +
    ylim(0, 53)  
  
  scree_plot_file <- file.path(
    outdir,
    paste0("Scree_Plot_", pretty_feature_name, "_", cohort_type, ".png")
  )
  ggsave(scree_plot_file, plot = scree_plot, width = 4, height = 5, dpi = 500)
  
  # Export results for each threshold
  results <- list()
  for (i in seq_along(thresholds)) {
    th <- thresholds[i]
    n_comp <- num_components[i]
    
    pca_scores <- pca_result$x[, 1:n_comp, drop = FALSE]
    
    csv_file <- file.path(outdir, paste0(pretty_feature_name, "_", cohort_type,
                                         "_PCA_", th * 100, "percent.csv"))
    # write.csv(as.data.frame(pca_scores), csv_file, row.names = TRUE)
    
    rds_file <- file.path(outdir, paste0(pretty_feature_name, "_", cohort_type,
                                         "_PCA_", th * 100, "percent.rds"))
  #  saveRDS(pca_result, rds_file) ## unblock to save file
    
    results[[paste0(th * 100, "%")]] <- list(
      num_components = n_comp,
      pca_scores = pca_scores,
      cumulative_variance = cum_variance[n_comp]
    )
  }
  
  pc_summary <- tibble(
    Feature = pretty_feature_name,
    Cohort = cohort_type,
    Threshold = paste0(thresholds * 100, "%"),
    Num_Components = num_components,
    Cumulative_Variance = cum_variance[num_components]
  )
  
  pc_summary_file <- file.path(
    outdir,
    paste0("PC_Summary_", pretty_feature_name, "_", cohort_type, ".csv")
  )
  write.csv(pc_summary, pc_summary_file, row.names = FALSE)
  
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
outdir <- "PCA_results/Run_PCA_on_everything_March2025_updated_with_dataframe/"

# Create the directory if it does not exist
if (!dir.exists(outdir)) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
}

# List of datasets
datasets <- list(
  training_nucleosome_peak = list(data = training_nucleosome_peak, cohort = "Training"),
  validation_nucleosome_peak = list(data = validation_nucleosome_peak, cohort = "Validation"),
  training_delfi_ratios = list(data = training_delfi_ratios, cohort = "Training"),
  validation_delfi_ratios = list(data = validation_delfi_ratios, cohort = "Validation"),
  training_insert_matrix = list(data = training_insert_matrix, cohort = "Training"),
  validation_insert_matrix = list(data = validation_insert_matrix, cohort = "Validation"),
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
    thresholds = c(0.5, 0.6, 0.7, 0.8, 0.85, 0.90, 0.97)  # Retain 90% and 97% variance
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
 #   write.csv(as.data.frame(pca_scores_train), file.path(outdir, paste0(feature_name, "_Training_PCA_", th * 100, "percent.csv")), row.names = TRUE)
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
  #  write.csv(as.data.frame(pca_scores_validation), file.path(outdir, paste0(feature_name, "_Validation_PCA_", th * 100, "percent.csv")), row.names = TRUE)
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
outdir <- "PCA_results/Run_PCA_on_train_fit_on_validation_March2025/"

# Create the directory if it does not exist
if (!dir.exists(outdir)) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
}

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
      thresholds = c(0.7, 0.8, 0.85, 0.90, 0.97)  # Retain 90% and 97% variance
    )
  }
}


#### Now integrate the components 
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
# 2) Define helper function to add feature-specific prefixes to column names
#------------------------------------------------------------------------------
add_prefix_to_pcs <- function(df, prefix) {
  df %>%
    rename_with(~ paste0(prefix, .), -sample_id)
}

#------------------------------------------------------------------------------
#  3) Process each feature: load PCA objects and scores, filter PCs by variance,
#     fix sample IDs if necessary, and add prefixes.
#------------------------------------------------------------------------------

## End_motifs
pca_obj_end_motifs <- readRDS("PCA_results/Run_PCA_on_train_fit_on_validation_March2025/end_motifs_PCA_Object.rds")
train_scores_end_motifs <- readRDS("PCA_results/Run_PCA_on_train_fit_on_validation_March2025/end_motifs_Training_PCA_97percent.rds")
val_scores_end_motifs   <- readRDS("PCA_results/Run_PCA_on_train_fit_on_validation_March2025/end_motifs_Validation_PCA_97percent.rds")

training_end_motifs_5pct <- filter_pcs_by_variance(pca_obj_end_motifs, train_scores_end_motifs, threshold = 5)
training_end_motifs_1pct <- filter_pcs_by_variance(pca_obj_end_motifs, train_scores_end_motifs, threshold = 1)
validation_end_motifs_5pct <- filter_pcs_by_variance(pca_obj_end_motifs, val_scores_end_motifs, threshold = 5)
validation_end_motifs_1pct <- filter_pcs_by_variance(pca_obj_end_motifs, val_scores_end_motifs, threshold = 1)

# Fix sample IDs if needed
training_end_motifs_5pct$sample_id <- gsub("_motifs", "", training_end_motifs_5pct$sample_id)
training_end_motifs_1pct$sample_id <- gsub("_motifs", "", training_end_motifs_1pct$sample_id)
validation_end_motifs_5pct$sample_id <- gsub("_motifs", "", validation_end_motifs_5pct$sample_id)
validation_end_motifs_1pct$sample_id <- gsub("_motifs", "", validation_end_motifs_1pct$sample_id)

## Methylation
pca_obj_methylation <- readRDS("PCA_results/Run_PCA_on_train_fit_on_validation_March2025/methylation_PCA_Object.rds")
train_scores_methylation <- readRDS("PCA_results/Run_PCA_on_train_fit_on_validation_March2025/methylation_Training_PCA_97percent.rds")
val_scores_methylation   <- readRDS("PCA_results/Run_PCA_on_train_fit_on_validation_March2025/methylation_Validation_PCA_97percent.rds")

training_methylation_5pct <- filter_pcs_by_variance(pca_obj_methylation, train_scores_methylation, threshold = 5)
training_methylation_1pct <- filter_pcs_by_variance(pca_obj_methylation, train_scores_methylation, threshold = 1)
validation_methylation_5pct <- filter_pcs_by_variance(pca_obj_methylation, val_scores_methylation, threshold = 5)
validation_methylation_1pct <- filter_pcs_by_variance(pca_obj_methylation, val_scores_methylation, threshold = 1)

# Adjust naming by appending "_dedup" for methylation sample IDs
training_methylation_5pct$sample_id <- paste0(training_methylation_5pct$sample_id, "_dedup")
training_methylation_1pct$sample_id <- paste0(training_methylation_1pct$sample_id, "_dedup")
validation_methylation_5pct$sample_id <- paste0(validation_methylation_5pct$sample_id, "_dedup")
validation_methylation_1pct$sample_id <- paste0(validation_methylation_1pct$sample_id, "_dedup")

## Nucleosome_peak
pca_obj_nucleosome <- readRDS("PCA_results/Run_PCA_on_train_fit_on_validation_March2025/nucleosome_peak_PCA_Object.rds")
train_scores_nucleosome <- readRDS("PCA_results/Run_PCA_on_train_fit_on_validation_March2025/nucleosome_peak_Training_PCA_97percent.rds")
val_scores_nucleosome   <- readRDS("PCA_results/Run_PCA_on_train_fit_on_validation_March2025/nucleosome_peak_Validation_PCA_97percent.rds")

training_nucleosome_5pct <- filter_pcs_by_variance(pca_obj_nucleosome, train_scores_nucleosome, threshold = 5)
training_nucleosome_1pct <- filter_pcs_by_variance(pca_obj_nucleosome, train_scores_nucleosome, threshold = 1)
validation_nucleosome_5pct <- filter_pcs_by_variance(pca_obj_nucleosome, val_scores_nucleosome, threshold = 5)
validation_nucleosome_1pct <- filter_pcs_by_variance(pca_obj_nucleosome, val_scores_nucleosome, threshold = 1)

# Fix sample IDs for nucleosome_peak
training_nucleosome_5pct$sample_id <- gsub("_peak_distance", "", training_nucleosome_5pct$sample_id)
training_nucleosome_1pct$sample_id <- gsub("_peak_distance", "", training_nucleosome_1pct$sample_id)
validation_nucleosome_5pct$sample_id <- gsub("_peak_distance", "", validation_nucleosome_5pct$sample_id)
validation_nucleosome_1pct$sample_id <- gsub("_peak_distance", "", validation_nucleosome_1pct$sample_id)

## Delfi_ratios
pca_obj_delfi <- readRDS("PCA_results/Run_PCA_on_train_fit_on_validation_March2025/delfi_ratios_PCA_Object.rds")
train_scores_delfi <- readRDS("PCA_results/Run_PCA_on_train_fit_on_validation_March2025/delfi_ratios_Training_PCA_97percent.rds")
val_scores_delfi   <- readRDS("PCA_results/Run_PCA_on_train_fit_on_validation_March2025/delfi_ratios_Validation_PCA_97percent.rds")

training_delfi_5pct <- filter_pcs_by_variance(pca_obj_delfi, train_scores_delfi, threshold = 5)
training_delfi_1pct <- filter_pcs_by_variance(pca_obj_delfi, train_scores_delfi, threshold = 1)
validation_delfi_5pct <- filter_pcs_by_variance(pca_obj_delfi, val_scores_delfi, threshold = 5)
validation_delfi_1pct <- filter_pcs_by_variance(pca_obj_delfi, val_scores_delfi, threshold = 1)

## Insert size
pca_obj_insert <- readRDS("PCA_results/Run_PCA_on_train_fit_on_validation_March2025/insert_matrix_PCA_Object.rds")
train_scores_insert <- readRDS("PCA_results/Run_PCA_on_train_fit_on_validation_March2025/insert_matrix_Training_PCA_97percent.rds")
val_scores_insert   <- readRDS("PCA_results/Run_PCA_on_train_fit_on_validation_March2025/insert_matrix_Validation_PCA_97percent.rds")

training_insert_sizes_5pct <- filter_pcs_by_variance(pca_obj_insert, train_scores_insert, threshold = 5)
training_insert_sizes_1pct <- filter_pcs_by_variance(pca_obj_insert, train_scores_insert, threshold = 1)
validation_insert_sizes_5pct <- filter_pcs_by_variance(pca_obj_insert, val_scores_insert, threshold = 5)
validation_insert_sizes_1pct <- filter_pcs_by_variance(pca_obj_insert, val_scores_insert, threshold = 1)

#------------------------------------------------------------------------------
# 4) Apply feature-specific prefixes to each training and validation data frame
#------------------------------------------------------------------------------
# For training sets
training_end_motifs_5pct <- add_prefix_to_pcs(training_end_motifs_5pct, "End_motifs_")
training_end_motifs_1pct <- add_prefix_to_pcs(training_end_motifs_1pct, "End_motifs_")

training_methylation_5pct <- add_prefix_to_pcs(training_methylation_5pct, "Methylation_")
training_methylation_1pct <- add_prefix_to_pcs(training_methylation_1pct, "Methylation_")

training_nucleosome_5pct <- add_prefix_to_pcs(training_nucleosome_5pct, "Nucleosome_")
training_nucleosome_1pct <- add_prefix_to_pcs(training_nucleosome_1pct, "Nucleosome_")

training_delfi_5pct <- add_prefix_to_pcs(training_delfi_5pct, "Fragment_Ratios_")
training_delfi_1pct <- add_prefix_to_pcs(training_delfi_1pct, "Fragment_Ratios_")

training_insert_sizes_5pct <- add_prefix_to_pcs(training_insert_sizes_5pct, "Insert_size_")
training_insert_sizes_1pct <- add_prefix_to_pcs(training_insert_sizes_1pct, "Insert_size_")

# For validation sets
validation_end_motifs_5pct <- add_prefix_to_pcs(validation_end_motifs_5pct, "End_motifs_")
validation_end_motifs_1pct <- add_prefix_to_pcs(validation_end_motifs_1pct, "End_motifs_")

validation_methylation_5pct <- add_prefix_to_pcs(validation_methylation_5pct, "Methylation_")
validation_methylation_1pct <- add_prefix_to_pcs(validation_methylation_1pct, "Methylation_")

validation_nucleosome_5pct <- add_prefix_to_pcs(validation_nucleosome_5pct, "Nucleosome_")
validation_nucleosome_1pct <- add_prefix_to_pcs(validation_nucleosome_1pct, "Nucleosome_")

validation_delfi_5pct <- add_prefix_to_pcs(validation_delfi_5pct, "Fragment_Ratios_")
validation_delfi_1pct <- add_prefix_to_pcs(validation_delfi_1pct, "Fragment_Ratios_")

validation_insert_sizes_5pct <- add_prefix_to_pcs(validation_insert_sizes_5pct, "Insert_size_")
validation_insert_sizes_1pct <- add_prefix_to_pcs(validation_insert_sizes_1pct, "Insert_size_")

#------------------------------------------------------------------------------
# 5) Now merge the >5% variance training sets across features
#------------------------------------------------------------------------------
training_combined_5pct <- training_end_motifs_5pct %>%
  full_join(training_methylation_5pct, by = "sample_id") %>%
  full_join(training_nucleosome_5pct,  by = "sample_id") %>%
  full_join(training_delfi_5pct,       by = "sample_id") %>%
  full_join(training_insert_sizes_5pct,by = "sample_id")

# Convert sample_id to rownames and remove the column
rownames(training_combined_5pct) <- training_combined_5pct$sample_id
training_combined_5pct <- training_combined_5pct %>% dplyr::select(-sample_id)

saveRDS(training_combined_5pct, "PCA_results/Run_PCA_on_train_fit_on_validation_March2025/Training_combined_5pct_PCs.rds")

#------------------------------------------------------------------------------
# 6) Merge the >1% variance training sets across features
#------------------------------------------------------------------------------
training_combined_1pct <- training_end_motifs_1pct %>%
  full_join(training_methylation_1pct, by = "sample_id") %>%
  full_join(training_nucleosome_1pct,  by = "sample_id") %>%
  full_join(training_delfi_1pct,       by = "sample_id") %>%
  full_join(training_insert_sizes_1pct,by = "sample_id")

rownames(training_combined_1pct) <- training_combined_1pct$sample_id
training_combined_1pct <- training_combined_1pct %>% dplyr::select(-sample_id)

saveRDS(training_combined_1pct, "PCA_results/Run_PCA_on_train_fit_on_validation_March2025/Training_combined_1pct_PCs.rds")

#------------------------------------------------------------------------------
# 7) Merge the >5% variance validation sets across features
#------------------------------------------------------------------------------
validation_combined_5pct <- validation_end_motifs_5pct %>%
  full_join(validation_methylation_5pct, by = "sample_id") %>%
  full_join(validation_nucleosome_5pct,  by = "sample_id") %>%
  full_join(validation_delfi_5pct,       by = "sample_id") %>%
  full_join(validation_insert_sizes_5pct,by = "sample_id")

rownames(validation_combined_5pct) <- validation_combined_5pct$sample_id
validation_combined_5pct <- validation_combined_5pct %>% dplyr::select(-sample_id)

saveRDS(validation_combined_5pct, "PCA_results/Run_PCA_on_train_fit_on_validation_March2025/Validation_combined_5pct_PCs.rds")

#------------------------------------------------------------------------------
# 8) Merge the >1% variance validation sets across features
#------------------------------------------------------------------------------
validation_combined_1pct <- validation_end_motifs_1pct %>%
  full_join(validation_methylation_1pct, by = "sample_id") %>%
  full_join(validation_nucleosome_1pct,  by = "sample_id") %>%
  full_join(validation_delfi_1pct,       by = "sample_id") %>%
  full_join(validation_insert_sizes_1pct,by = "sample_id")

rownames(validation_combined_1pct) <- validation_combined_1pct$sample_id
validation_combined_1pct <- validation_combined_1pct %>% dplyr::select(-sample_id)

saveRDS(validation_combined_1pct, "PCA_results/Run_PCA_on_train_fit_on_validation_March2025/Validation_combined_1pct_PCs.rds")

#------------------------------------------------------------------------------
# 9) Examples of merging two features (with >1% variance) for training and validation
#------------------------------------------------------------------------------
combined_end_motifs_methylation_1pct_training <- training_end_motifs_1pct %>%
  full_join(training_methylation_1pct, by = "sample_id")
rownames(combined_end_motifs_methylation_1pct_training) <- combined_end_motifs_methylation_1pct_training$sample_id
combined_end_motifs_methylation_1pct_training <- combined_end_motifs_methylation_1pct_training %>% dplyr::select(-sample_id)
saveRDS(combined_end_motifs_methylation_1pct_training, 
        "PCA_results/Run_PCA_on_train_fit_on_validation_March2025/Combined_EndMotif_Methylation_1pct_PCs_training.rds")

combined_end_motifs_methylation_1pct_validation <- validation_end_motifs_1pct %>%
  full_join(validation_methylation_1pct, by = "sample_id")
rownames(combined_end_motifs_methylation_1pct_validation) <- combined_end_motifs_methylation_1pct_validation$sample_id
combined_end_motifs_methylation_1pct_validation <- combined_end_motifs_methylation_1pct_validation %>% dplyr::select(-sample_id)
saveRDS(combined_end_motifs_methylation_1pct_validation, 
        "PCA_results/Run_PCA_on_train_fit_on_validation_March2025/Combined_EndMotif_Methylation_1pct_PCs_validation.rds")

# Also include End Motifs + Methylation + Fragment Ratios (>1% Variance)
combined_end_motifs_methylation_fragment_ratios_1pct_training <- training_end_motifs_1pct %>%
  full_join(training_methylation_1pct, by = "sample_id") %>%
  full_join(training_delfi_1pct, by = "sample_id")
rownames(combined_end_motifs_methylation_fragment_ratios_1pct_training) <- combined_end_motifs_methylation_fragment_ratios_1pct_training$sample_id
combined_end_motifs_methylation_fragment_ratios_1pct_training <- combined_end_motifs_methylation_fragment_ratios_1pct_training %>% dplyr::select(-sample_id)
saveRDS(combined_end_motifs_methylation_fragment_ratios_1pct_training, 
        "PCA_results/Run_PCA_on_train_fit_on_validation_March2025/Combined_EndMotif_Methylation_FragmentRatios_1pct_PCs_training.rds")

combined_end_motifs_methylation_fragment_ratios_1pct_validation <- validation_end_motifs_1pct %>%
  full_join(validation_methylation_1pct, by = "sample_id") %>%
  full_join(validation_delfi_1pct, by = "sample_id")
rownames(combined_end_motifs_methylation_fragment_ratios_1pct_validation) <- combined_end_motifs_methylation_fragment_ratios_1pct_validation$sample_id
combined_end_motifs_methylation_fragment_ratios_1pct_validation <- combined_end_motifs_methylation_fragment_ratios_1pct_validation %>% dplyr::select(-sample_id)
saveRDS(combined_end_motifs_methylation_fragment_ratios_1pct_validation, 
        "PCA_results/Run_PCA_on_train_fit_on_validation_March2025/Combined_EndMotif_Methylation_FragmentRatios_1pct_PCs_validation.rds")

# Another example: Merging only Fragment Ratios + Methylation (>1% Variance) for Training and Validation
combined_ratios_methylation_1pct_training <- training_delfi_1pct %>%
  full_join(training_methylation_1pct, by = "sample_id")
rownames(combined_ratios_methylation_1pct_training) <- combined_ratios_methylation_1pct_training$sample_id
combined_ratios_methylation_1pct_training <- combined_ratios_methylation_1pct_training %>% dplyr::select(-sample_id)
saveRDS(combined_ratios_methylation_1pct_training, 
        "PCA_results/Run_PCA_on_train_fit_on_validation_March2025/Combined_ratios_Methylation_1pct_PCs_training.rds")

combined_ratios_methylation_1pct_validation <- validation_delfi_1pct %>%
  full_join(validation_methylation_1pct, by = "sample_id")
rownames(combined_ratios_methylation_1pct_validation) <- combined_ratios_methylation_1pct_validation$sample_id
combined_ratios_methylation_1pct_validation <- combined_ratios_methylation_1pct_validation %>% dplyr::select(-sample_id)
saveRDS(combined_ratios_methylation_1pct_validation, 
        "PCA_results/Run_PCA_on_train_fit_on_validation_March2025/Combined_ratios_Methylation_1pct_PCs_validation.rds")










##### Now train on just the cancer for Yuanchang 
#### Now split into training and validation - first for PE 
### 1) Filter metadata to remove TCGE-CFMe-HCC and 'SE' samples
metadata_filtered <- metadata_df %>%
  filter(
    !(project_id == "TCGE-CFMe-HCC" & cancer_type == "Liver Cancer"),  # Remove HCC samples only if cancer_type is Liver Cancer
    type != "SE"  # Remove SE samples
  )

metadata_filtered <- metadata_filtered %>%
  filter(
    !(cancer_type == "Normal")
  ) # Remove healthy 

# 2. Create the validation set from TCGE-CFMe-INSPIRE and TCGE-CFMe-LFS,
#    but exclude LFS Survivor within TCGE-CFMe-LFS.
metadata_val <- metadata_filtered %>%
  filter(project_id %in% c("TCGE-CFMe-INSPIRE", "TCGE-CFMe-LFS", "TCGE-CFMe-HCC"))

# 3. Create the training set from everything else (i.e., not in validation),
#    also excluding any leftover from the cohorts we said not to include (HCC, INSPIRE, LFS).
metadata_train <- metadata_filtered %>%
  filter(!project_id %in% c("TCGE-CFMe-INSPIRE", "TCGE-CFMe-LFS", "TCGE-CFMe-HCC"))

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

#### Nucleosome peak datasets
matrix_nucleosome_peak_combined <- rbind(matrix_nucleosome_peak, matrix_nucleosome_peak_validation)
nucleosome_peak_split <- split_feature_data(matrix_nucleosome_peak_combined, "nucleosome_peak_sampleid")

#### DELFI ratios datasets
delfi_ratios_combined <- rbind(delfi_ratios, delfi_ratios_validation)
delfi_ratios_split <- split_feature_data(delfi_ratios_combined, "sample_id")

#### Insert size datasets
insert_matrix_combined <- rbind(insert_matrix, insert_matrix_validation)
insert_matrix_split <- split_feature_data(insert_matrix_combined, "sample_id")

#### End motifs datasets
end_motifs_combined <- rbind(end_motifs_all_frags, end_motifs_validation)
end_motifs_split <- split_feature_data(end_motifs_combined, "endmotif_name")


##### Methylation dataset 

### First add the validation to the PE 
data_validation <- t(data_validation)
#data_PE <- t(data_PE) # ensure samples are rows

# Check if the column names match
if (!identical(colnames(data_PE), colnames(data_validation))) {
  cat("Columns differ between data_PE and data_validation\n")
  
  # Identify columns in data_PE but not in data_validation
  columns_in_PE_not_in_validation <- setdiff(colnames(data_PE), colnames(data_validation))
  if (length(columns_in_PE_not_in_validation) > 0) {
    cat("Columns in data_PE but not in data_validation:\n")
    print(columns_in_PE_not_in_validation)
  } else {
    cat("No extra columns in data_PE.\n")
  }
  
  # Identify columns in data_validation but not in data_PE
  columns_in_validation_not_in_PE <- setdiff(colnames(data_validation), colnames(data_PE))
  if (length(columns_in_validation_not_in_PE) > 0) {
    cat("Columns in data_validation but not in data_PE:\n")
    print(columns_in_validation_not_in_PE)
  } else {
    cat("No extra columns in data_validation.\n")
  }
} else {
  cat("Column names match between data_PE and data_validation.\n")
}

# Ensure both datasets have the same columns
# Combine the two data frames
# Create a union of column names
all_columns <- union(colnames(data_PE), colnames(data_validation))

# Align both matrices by adding missing columns filled with NA
data_PE_aligned <- matrix(NA, nrow = nrow(data_PE), ncol = length(all_columns),
                          dimnames = list(rownames(data_PE), all_columns))
data_PE_aligned[, colnames(data_PE)] <- data_PE

data_validation_aligned <- matrix(NA, nrow = nrow(data_validation), ncol = length(all_columns),
                                  dimnames = list(rownames(data_validation), all_columns))
data_validation_aligned[, colnames(data_validation)] <- data_validation

# Bind the matrices row-wise
data_methylation_combined <- rbind(data_PE_aligned, data_validation_aligned)
rm(data_validation_aligned)
rm(data_PE_aligned)

methylation_split_PE <- split_feature_data(data_methylation_combined, "sequencing_id")
#methylation_split_SE <- split_feature_data(data_SE, "sample_name") to do later

### 5) Final combined datasets for each feature
# Training datasets
training_nucleosome_peak <- nucleosome_peak_split$train
training_delfi_ratios <- delfi_ratios_split$train
training_insert_matrix <- insert_matrix_split$train
training_end_motifs <- end_motifs_split$train
training_methylation <- methylation_split_PE$train

# Validation datasets
validation_nucleosome_peak <- nucleosome_peak_split$validation
validation_delfi_ratios <- delfi_ratios_split$validation
validation_insert_matrix <- insert_matrix_split$validation
validation_end_motifs <- end_motifs_split$validation
validation_methylation <- methylation_split_PE$validation


### Run PCA on everything and export plots

## Set function
run_pca_and_export <- function(data_matrix, feature_name, cohort_type, outdir, thresholds = c(0.90, 0.97)) {
  # Ensure data is numeric
  if (!is.matrix(data_matrix)) {
    data_matrix <- as.matrix(data_matrix)
  }
  
  # Remove constant/zero-variance columns
  col_sd <- apply(data_matrix, 2, sd, na.rm = TRUE)  # Calculate standard deviation of columns
  # data_matrix <- data_matrix[, col_sd > 0, drop = FALSE]  # Keep only columns with non-zero SD
  data_matrix <- data_matrix[, !is.na(col_sd) & col_sd > 0, drop = FALSE] # Keep only columns with non-zero SD
  
  
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
  
  # Export cumulative variance plot
  cum_var_df <- tibble(
    PC = seq_along(cum_variance),
    CumulativeVariance = cum_variance
  )
  
  # Plot cumulative explained variance
  for (i in seq_along(thresholds)) {
    th <- thresholds[i]
    n_comp <- num_components[i]
    
    # Plot cumulative variance
    p_cum <- ggplot(cum_var_df, aes(x = PC, y = CumulativeVariance)) +
      geom_line(color = "black") +
      geom_point(color = "black") +
      geom_hline(yintercept = cum_variance[n_comp], color = "red", linetype = "dashed") +
      geom_vline(xintercept = n_comp, color = "red", linetype = "dashed") +
      labs(
        title = paste("Cumulative Variance -", feature_name, "(", cohort_type, ")"),
        subtitle = paste0("Selected PCs: ", n_comp, " (Cumulative: ", round(cum_variance[n_comp], 2), ")"),
        x = "Number of Principal Components",
        y = "Cumulative Explained Variance"
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10)
      )
    
    # Save cumulative variance plot
    plot_file <- file.path(outdir, paste0("Cumulative_Explained_Variance_", feature_name, "_", cohort_type, "_", th * 100, "percent.png"))
    ggsave(plot_file, plot = p_cum, width = 7, height = 5, dpi = 300)
  }
  
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
    # write.csv(as.data.frame(pca_scores), csv_file, row.names = TRUE)
    
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
outdir <- "PCA_results/Run_PCA_on_everything_March2025_cancer_only/"

# Create the directory if it does not exist
if (!dir.exists(outdir)) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
}

# List of datasets
datasets <- list(
  training_nucleosome_peak = list(data = training_nucleosome_peak, cohort = "Training"),
  validation_nucleosome_peak = list(data = validation_nucleosome_peak, cohort = "Validation"),
  training_delfi_ratios = list(data = training_delfi_ratios, cohort = "Training"),
  validation_delfi_ratios = list(data = validation_delfi_ratios, cohort = "Validation"),
  training_insert_matrix = list(data = training_insert_matrix, cohort = "Training"),
  validation_insert_matrix = list(data = validation_insert_matrix, cohort = "Validation"),
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
    thresholds = c(0.5, 0.6, 0.7, 0.8, 0.85, 0.90, 0.97)  # Retain 90% and 97% variance
  )
}




































########## Old code
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
pca_obj_end_motifs <- readRDS("PCA_results/Run_PCA_on_train_fit_on_validation/end_motifs_PCA_Object.rds")
train_scores_end_motifs <- readRDS("PCA_results/Run_PCA_on_train_fit_on_validation/end_motifs_Training_PCA_97percent.rds")
val_scores_end_motifs   <- readRDS("PCA_results/Run_PCA_on_train_fit_on_validation/end_motifs_Validation_PCA_97percent.rds")

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
pca_obj_methylation <- readRDS("PCA_results/Run_PCA_on_train_fit_on_validation/methylation_PCA_Object.rds")
train_scores_methylation <- readRDS("PCA_results/Run_PCA_on_train_fit_on_validation/methylation_Training_PCA_97percent.rds")
val_scores_methylation   <- readRDS("PCA_results/Run_PCA_on_train_fit_on_validation/methylation_Validation_PCA_97percent.rds")

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
pca_obj_nucleosome <- readRDS("PCA_results/Run_PCA_on_train_fit_on_validation/nucleosome_peak_PCA_Object.rds")
train_scores_nucleosome <- readRDS("PCA_results/Run_PCA_on_train_fit_on_validation/nucleosome_peak_Training_PCA_97percent.rds")
val_scores_nucleosome   <- readRDS("PCA_results/Run_PCA_on_train_fit_on_validation/nucleosome_peak_Validation_PCA_97percent.rds")

training_nucleosome_5pct <- filter_pcs_by_variance(pca_obj_nucleosome, train_scores_nucleosome, threshold = 5)
training_nucleosome_1pct <- filter_pcs_by_variance(pca_obj_nucleosome, train_scores_nucleosome, threshold = 1)

validation_nucleosome_5pct <- filter_pcs_by_variance(pca_obj_nucleosome, val_scores_nucleosome, threshold = 5)
validation_nucleosome_1pct <- filter_pcs_by_variance(pca_obj_nucleosome, val_scores_nucleosome, threshold = 1)

# Fix sample IDs
training_nucleosome_5pct$sample_id <- gsub("_peak_distance", "", training_nucleosome_5pct$sample_id)
training_nucleosome_1pct$sample_id <- gsub("_peak_distance", "", training_nucleosome_1pct$sample_id)
validation_nucleosome_5pct$sample_id <- gsub("_peak_distance", "", validation_nucleosome_5pct$sample_id)
validation_nucleosome_1pct$sample_id <- gsub("_peak_distance", "", validation_nucleosome_1pct$sample_id)


## Delfi_ratios
pca_obj_delfi <- readRDS("PCA_results/Run_PCA_on_train_fit_on_validation/delfi_ratios_PCA_Object.rds")
train_scores_delfi <- readRDS("PCA_results/Run_PCA_on_train_fit_on_validation/delfi_ratios_Training_PCA_97percent.rds")
val_scores_delfi   <- readRDS("PCA_results/Run_PCA_on_train_fit_on_validation/delfi_ratios_Validation_PCA_97percent.rds")

training_delfi_5pct <- filter_pcs_by_variance(pca_obj_delfi, train_scores_delfi, threshold = 5)
training_delfi_1pct <- filter_pcs_by_variance(pca_obj_delfi, train_scores_delfi, threshold = 1)

validation_delfi_5pct <- filter_pcs_by_variance(pca_obj_delfi, val_scores_delfi, threshold = 5)
validation_delfi_1pct <- filter_pcs_by_variance(pca_obj_delfi, val_scores_delfi, threshold = 1)


## Insert size
pca_obj_insert <- readRDS("PCA_results/Run_PCA_on_train_fit_on_validation/insert_matrix_PCA_Object.rds")
train_scores_insert <- readRDS("PCA_results/Run_PCA_on_train_fit_on_validation/insert_matrix_Training_PCA_97percent.rds")
val_scores_insert   <- readRDS("PCA_results/Run_PCA_on_train_fit_on_validation/insert_matrix_Validation_PCA_97percent.rds")

training_insert_sizes_5pct <- filter_pcs_by_variance(pca_obj_insert, train_scores_insert, threshold = 5)
training_insert_sizes_1pct <- filter_pcs_by_variance(pca_obj_insert, train_scores_insert, threshold = 1)

validation_insert_sizes_5pct <- filter_pcs_by_variance(pca_obj_insert, val_scores_insert, threshold = 5)
validation_insert_sizes_1pct <- filter_pcs_by_variance(pca_obj_insert, val_scores_insert, threshold = 1)


### Rename to make easier to know what came from where 
# Apply prefix to training sets
add_prefix_to_pcs <- function(df, prefix) {
  df %>%
    rename_with(~ paste0(prefix, .), -sample_id)  # Add prefix to all columns except 'sample_id'
}


training_end_motifs_5pct <- add_prefix_to_pcs(training_end_motifs_5pct, "EndMotifs_")
training_methylation_5pct <- add_prefix_to_pcs(training_methylation_5pct, "Methylation_")
training_nucleosome_5pct <- add_prefix_to_pcs(training_nucleosome_5pct, "Nucleosome_")
training_delfi_5pct <- add_prefix_to_pcs(training_delfi_5pct, "FragmentRatios_")
training_insert_sizes_5pct <- add_prefix_to_pcs(training_insert_sizes_5pct, "InsertSize_")

training_end_motifs_1pct <- add_prefix_to_pcs(training_end_motifs_1pct, "EndMotifs_")
training_methylation_1pct <- add_prefix_to_pcs(training_methylation_1pct, "Methylation_")
training_nucleosome_1pct <- add_prefix_to_pcs(training_nucleosome_1pct, "Nucleosome_")
training_delfi_1pct <- add_prefix_to_pcs(training_delfi_1pct, "FragmentRatios_")
training_insert_sizes_1pct <- add_prefix_to_pcs(training_insert_sizes_1pct, "InsertSize_")

# Apply prefix to validation sets
validation_end_motifs_5pct <- add_prefix_to_pcs(validation_end_motifs_5pct, "EndMotifs_")
validation_methylation_5pct <- add_prefix_to_pcs(validation_methylation_5pct, "Methylation_")
validation_nucleosome_5pct <- add_prefix_to_pcs(validation_nucleosome_5pct, "Nucleosome_")
validation_delfi_5pct <- add_prefix_to_pcs(validation_delfi_5pct, "FragmentRatios_")
validation_insert_sizes_5pct <- add_prefix_to_pcs(validation_insert_sizes_5pct, "InsertSize_")

validation_end_motifs_1pct <- add_prefix_to_pcs(validation_end_motifs_1pct, "EndMotifs_")
validation_methylation_1pct <- add_prefix_to_pcs(validation_methylation_1pct, "Methylation_")
validation_nucleosome_1pct <- add_prefix_to_pcs(validation_nucleosome_1pct, "Nucleosome_")
validation_delfi_1pct <- add_prefix_to_pcs(validation_delfi_1pct, "FragmentRatios_")
validation_insert_sizes_1pct <- add_prefix_to_pcs(validation_insert_sizes_1pct, "InsertSize_")


#------------------------------------------------------------------------------
# 3) Now merge the >5% variance training sets across features
#------------------------------------------------------------------------------
training_combined_5pct <- training_end_motifs_5pct %>%
  full_join(training_methylation_5pct, by = "sample_id") %>%
  full_join(training_nucleosome_5pct,  by = "sample_id") %>%
  full_join(training_delfi_5pct,       by = "sample_id") %>%
  full_join(training_insert_sizes_5pct,by = "sample_id")

# Convert sample_id to rownames
rownames(training_combined_5pct) <- training_combined_5pct$sample_id
training_combined_5pct <- training_combined_5pct %>% dplyr::select(-sample_id)

# Save the >5% variance training data
saveRDS(training_combined_5pct, "PCA_results/Run_PCA_on_train_fit_on_validation/Training_combined_5pct_PCs.rds")


#------------------------------------------------------------------------------
# 4) Merge the >1% variance training sets across features
#------------------------------------------------------------------------------
training_combined_1pct <- training_end_motifs_1pct %>%
  full_join(training_methylation_1pct, by = "sample_id") %>%
  full_join(training_nucleosome_1pct,  by = "sample_id") %>%
  full_join(training_delfi_1pct,       by = "sample_id") %>%
  full_join(training_insert_sizes_1pct,by = "sample_id")

rownames(training_combined_1pct) <- training_combined_1pct$sample_id
training_combined_1pct <- training_combined_1pct %>% dplyr::select(-sample_id)

# Save the >1% variance training data
saveRDS(training_combined_1pct, "PCA_results/Run_PCA_on_train_fit_on_validation/Training_combined_1pct_PCs.rds")


#------------------------------------------------------------------------------
# 5) Merge the >5% variance validation sets across features
#------------------------------------------------------------------------------
validation_combined_5pct <- validation_end_motifs_5pct %>%
  full_join(validation_methylation_5pct, by = "sample_id") %>%
  full_join(validation_nucleosome_5pct,  by = "sample_id") %>%
  full_join(validation_delfi_5pct,       by = "sample_id") %>%
  full_join(validation_insert_sizes_5pct,by = "sample_id")

rownames(validation_combined_5pct) <- validation_combined_5pct$sample_id
validation_combined_5pct <- validation_combined_5pct %>% dplyr::select(-sample_id)

saveRDS(validation_combined_5pct, "PCA_results/Run_PCA_on_train_fit_on_validation/Validation_combined_5pct_PCs.rds")


#------------------------------------------------------------------------------
# 6) Merge the >1% variance validation sets across features
#------------------------------------------------------------------------------
validation_combined_1pct <- validation_end_motifs_1pct %>%
  full_join(validation_methylation_1pct, by = "sample_id") %>%
  full_join(validation_nucleosome_1pct,  by = "sample_id") %>%
  full_join(validation_delfi_1pct,       by = "sample_id") %>%
  full_join(validation_insert_sizes_1pct,by = "sample_id")

rownames(validation_combined_1pct) <- validation_combined_1pct$sample_id
validation_combined_1pct <- validation_combined_1pct %>% dplyr::select(-sample_id)

saveRDS(validation_combined_1pct, "PCA_results/Run_PCA_on_train_fit_on_validation/Validation_combined_1pct_PCs.rds")


#------------------------------------------------------------------------------
# 7) Example: Merging only End Motifs + Methylation (>1% Variance) for Training
#------------------------------------------------------------------------------
combined_end_motifs_methylation_1pct_training <- training_end_motifs_1pct %>%
  full_join(training_methylation_1pct, by = "sample_id")

rownames(combined_end_motifs_methylation_1pct_training) <- combined_end_motifs_methylation_1pct_training$sample_id
combined_end_motifs_methylation_1pct_training <- combined_end_motifs_methylation_1pct_training %>% dplyr::select(-sample_id)

saveRDS(combined_end_motifs_methylation_1pct_training, 
        "PCA_results/Run_PCA_on_train_fit_on_validation/Combined_EndMotif_Methylation_1pct_PCs_training.rds")


#------------------------------------------------------------------------------
# 8) Example: Merging only End Motifs + Methylation (>1% Variance) for Validation
#------------------------------------------------------------------------------
combined_end_motifs_methylation_1pct_validation <- validation_end_motifs_1pct %>%
  full_join(validation_methylation_1pct, by = "sample_id")

rownames(combined_end_motifs_methylation_1pct_validation) <- combined_end_motifs_methylation_1pct_validation$sample_id
combined_end_motifs_methylation_1pct_validation <- combined_end_motifs_methylation_1pct_validation %>% dplyr::select(-sample_id)

saveRDS(combined_end_motifs_methylation_1pct_validation, 
        "PCA_results/Run_PCA_on_train_fit_on_validation/Combined_EndMotif_Methylation_1pct_PCs_validation.rds")



### Also bnow include end motif with methylation and fragment ratios 
combined_end_motifs_methylation_fragment_ratios_1pct_training <- training_end_motifs_1pct %>%
  full_join(training_methylation_1pct, by = "sample_id") %>%
  full_join(training_delfi_1pct, by = "sample_id")

rownames(combined_end_motifs_methylation_fragment_ratios_1pct_training) <- combined_end_motifs_methylation_fragment_ratios_1pct_training$sample_id
combined_end_motifs_methylation_fragment_ratios_1pct_training <- combined_end_motifs_methylation_fragment_ratios_1pct_training %>% dplyr::select(-sample_id)

saveRDS(combined_end_motifs_methylation_fragment_ratios_1pct_training, 
        "PCA_results/Run_PCA_on_train_fit_on_validation/Combined_EndMotif_Methylation_FragmentRatios_1pct_PCs_training.rds")


combined_end_motifs_methylation_fragment_ratios_1pct_validation <- validation_end_motifs_1pct %>%
  full_join(validation_methylation_1pct, by = "sample_id") %>%
  full_join(validation_delfi_1pct, by = "sample_id")

rownames(combined_end_motifs_methylation_fragment_ratios_1pct_validation) <- combined_end_motifs_methylation_fragment_ratios_1pct_validation$sample_id
combined_end_motifs_methylation_fragment_ratios_1pct_validation <- combined_end_motifs_methylation_fragment_ratios_1pct_validation %>% dplyr::select(-sample_id)

saveRDS(combined_end_motifs_methylation_fragment_ratios_1pct_validation, 
        "PCA_results/Run_PCA_on_train_fit_on_validation/Combined_EndMotif_Methylation_FragmentRatios_1pct_PCs_validation.rds")


## Add for one more - methylation + fragment ratios 
#------------------------------------------------------------------------------
# 7) Example: Merging only Ratios + Methylation (>1% Variance) for Training
#------------------------------------------------------------------------------
combined_ratios_methylation_1pct_training <- training_delfi_1pct %>%
  full_join(training_methylation_1pct, by = "sample_id")

rownames(combined_ratios_methylation_1pct_training) <- combined_ratios_methylation_1pct_training$sample_id
combined_ratios_methylation_1pct_training <- combined_ratios_methylation_1pct_training %>% dplyr::select(-sample_id)

saveRDS(combined_ratios_methylation_1pct_training, 
        "PCA_results/Run_PCA_on_train_fit_on_validation/Combined_ratios_Methylation_1pct_PCs_training.rds")


#------------------------------------------------------------------------------
# 8) Example: Merging only Ratios + Methylation (>1% Variance) for Validation
#------------------------------------------------------------------------------
combined_ratios_methylation_1pct_validation <- validation_delfi_1pct %>%
  full_join(validation_methylation_1pct, by = "sample_id")

rownames(combined_ratios_methylation_1pct_validation) <- combined_ratios_methylation_1pct_validation$sample_id
combined_ratios_methylation_1pct_validation <- combined_ratios_methylation_1pct_validation %>% dplyr::select(-sample_id)

saveRDS(combined_ratios_methylation_1pct_validation, 
        "PCA_results/Run_PCA_on_train_fit_on_validation/Combined_ratios_Methylation_1pct_PCs_validation.rds")









########### Continue here
### Now count the number of PCs 
## Do minus 1 for the non-combined 

# Function to count the number of components (excluding `sample_id` for individual features)
count_pcs <- function(df, is_combined = FALSE) {
  if (is_combined) {
    return(ncol(df))  # Keep all columns for combined datasets
  } else {
    return(ncol(df) - 1)  # Exclude sample_id for individual feature datasets
  }
}

# Create a summary table with the number of components for each dataset
training_component_counts <- data.frame(
  Dataset = c(
    "Training_EndMotifs_1pct", "Training_Methylation_1pct", "Training_Nucleosome_1pct", 
    "Training_FragmentRatios_1pct", "Training_InsertSize_1pct",
    
    "Training_EndMotifs_5pct", "Training_Methylation_5pct", "Training_Nucleosome_5pct", 
    "Training_FragmentRatios_5pct", "Training_InsertSize_5pct",
    
    "Combined_Training_1pct", "Combined_Training_5pct", 
    "Combined_EndMotif_Methylation_1pct", "Combined_EndMotif_Methylation_FragmentRatios_1pct",
    "Combined_Ratios_Methylation_1pct"
  ),
  Components = c(
    count_pcs(training_end_motifs_1pct, FALSE), count_pcs(training_methylation_1pct, FALSE), count_pcs(training_nucleosome_1pct, FALSE), 
    count_pcs(training_delfi_1pct, FALSE), count_pcs(training_insert_sizes_1pct, FALSE),
    
    count_pcs(training_end_motifs_5pct, FALSE), count_pcs(training_methylation_5pct, FALSE), count_pcs(training_nucleosome_5pct, FALSE), 
    count_pcs(training_delfi_5pct, FALSE), count_pcs(training_insert_sizes_5pct, FALSE),
    
    count_pcs(training_combined_1pct, TRUE), count_pcs(training_combined_5pct, TRUE), 
    count_pcs(combined_end_motifs_methylation_1pct_training, TRUE), 
    count_pcs(combined_end_motifs_methylation_fragment_ratios_1pct_training, TRUE),
    count_pcs(combined_ratios_methylation_1pct_training, TRUE)
  )
)

# Display the table
print(training_component_counts)

# Save table as CSV
write.csv(training_component_counts, "PCA_results/Run_PCA_on_train_fit_on_validation_March2025/Training_Component_Counts_cancer_type.csv", row.names = FALSE)







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
  "validation_methylation_1pct" = validation_methylation_1pct,
  "pca_obj_nucleosome" = pca_obj_nucleosome,
  "train_scores_nucleosome" = train_scores_nucleosome,
  "val_scores_nucleosome" = val_scores_nucleosome,
  "training_nucleosome_5pct" = training_nucleosome_5pct,
  "training_nucleosome_1pct" = training_nucleosome_1pct,
  "validation_nucleosome_5pct" = validation_nucleosome_5pct,
  "validation_nucleosome_1pct" = validation_nucleosome_1pct,
  "pca_obj_delfi" = pca_obj_delfi,
  "train_scores_delfi" = train_scores_delfi,
  "val_scores_delfi" = val_scores_delfi,
  "training_delfi_5pct" = training_delfi_5pct,
  "training_delfi_1pct" = training_delfi_1pct,
  "validation_delfi_5pct" = validation_delfi_5pct,
  "validation_delfi_1pct" = validation_delfi_1pct,
  "pca_obj_insert" = pca_obj_insert,
  "train_scores_insert" = train_scores_insert,
  "val_scores_insert" = val_scores_insert,
  "training_insert_sizes_5pct" = training_insert_sizes_5pct,
  "training_insert_sizes_1pct" = training_insert_sizes_1pct,
  "validation_insert_sizes_5pct" = validation_insert_sizes_5pct,
  "validation_insert_sizes_1pct" = validation_insert_sizes_1pct,
  "training_combined_5pct" = training_combined_5pct,
  "training_combined_1pct" = training_combined_1pct,
  "validation_combined_5pct" = validation_combined_5pct,
  "validation_combined_1pct" = validation_combined_1pct,
  "combined_end_motifs_methylation_1pct_training" = combined_end_motifs_methylation_1pct_training,
  "combined_end_motifs_methylation_1pct_validation" = combined_end_motifs_methylation_1pct_validation,
  "combined_end_motifs_methylation_fragment_ratios_1pct_training" = combined_end_motifs_methylation_fragment_ratios_1pct_training,
  "combined_end_motifs_methylation_fragment_ratios_1pct_validation" = combined_end_motifs_methylation_fragment_ratios_1pct_validation
)

# Check each dataframe and summarize results
na_summary <- do.call(rbind, lapply(names(dataframes_to_check), function(df_name) {
  df <- dataframes_to_check[[df_name]]
  check_na_rows(df_name, df)
}))

# Filter out NULLs from results
na_summary <- na_summary[!is.na(na_summary$RowsWithNA), ]

# Save the summary as a table
na_summary_path <- "PCA_results/Run_PCA_on_train_fit_on_validation/NA_Summary.csv"
write.csv(na_summary, na_summary_path, row.names = FALSE)

# Print the summary
print(na_summary)
message("NA summary saved to: ", na_summary_path)









## Edit metadata df
metadata_df <- metadata_df %>%
  dplyr::mutate(CN_classifier = ifelse(cancer_type_corrected_updated == "healthy", "healthy", "cancer"))

saveRDS(metadata_df, file = "Metadata df all samples with CN classifier Jan 2025.rds")


metadata_df <- metadata_df %>%
  dplyr::mutate(cancer_subtype_corrected = tolower(gsub(" ", "_", cancer_subtype)))

## Minor edit for uveal_melanoma
metadata_df <- metadata_df %>%
  mutate(cancer_subtype_corrected = ifelse(cancer_subtype_corrected == "uveal_melanoma", group, cancer_subtype_corrected))

saveRDS(metadata_df, file = "Metadata_df_all_with_corrected_cancer_subtype_Jan2025.rds")




# Define the output directory (change this path as needed)
outdir <- "Final Data Dec 2024/PE/"

# Export combined datasets as RDS files
saveRDS(matrix_nucleosome_peak_combined, file = file.path(outdir, "matrix_nucleosome_peak_combined.rds"))
saveRDS(delfi_ratios_combined, file = file.path(outdir, "delfi_ratios_combined.rds"))
saveRDS(insert_matrix_combined, file = file.path(outdir, "insert_matrix_combined.rds"))
saveRDS(end_motifs_combined, file = file.path(outdir, "end_motifs_combined.rds"))
saveRDS(data_methylation_combined, file = file.path(outdir, "data_methylation_combined.rds"))


### Export the actual classes of the validation samples 
### Now do this for the SE samples 

## Load in methylation and end motifs for the SE 
## Split to leave 10% for validation 

## Now get the metadata info 
# Create a new dataframe 'actual_class_probs' from metadata_df
metadata_df <- readRDS("Metadata_df_all_with_corrected_cancer_subtype_Jan2025.rds")
actual_class_probs <- metadata_df %>%
  select(sample_id, endmotif_name, nucleosome_peak_sampleid, 
         CN_classifier, cancer_type_corrected_updated, cancer_subtype_corrected)

saveRDS(actual_class_probs, file = "PCA_results/Run_PCA_on_train_fit_on_validation/actual_class_probs.rds")


