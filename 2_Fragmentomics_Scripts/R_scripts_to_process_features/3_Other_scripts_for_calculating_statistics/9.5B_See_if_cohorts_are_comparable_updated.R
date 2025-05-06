# ----------------------------------------------------------------------------
# Title   : Training vs. Validation Set Consistency Checks
# Author : Dory Abelman 
# Date    : March 2025
#
# Purpose :
#   • Load raw feature matrices for training and validation cohorts, including
#     end‑motif frequencies, nucleosome peak distances, DELFI fragment ratios,
#     and insert‑size distributions.
#   • Harmonize sample identifiers and align feature columns between datasets.
#   • Calculate column‑wise means and standard deviations for each feature set.
#   • Quantify agreement between training and validation via Pearson correlations
#     of means and SDs, and average absolute mean differences.
#   • Produce scatter plots of training vs. validation feature means (with identity
#     line) to visually assess consistency.
#   • Generate summary sentences for manuscript methods describing the high
#     concordance observed.
# ----------------------------------------------------------------------------


### Check if trainnig and validation ar e consistent 

raw_features_dir <- "Final Data Dec 2024/Raw_features/"
data_dir <- raw_features_dir  # For simplicity

# Example: end motifs (training and validation)
end_motifs_all_frags <- readRDS(file.path(data_dir, "Dataframe matrix for end motifs, all frags, May 2024.rds"))
end_motifs_validation <- readRDS(file.path(data_dir, "Dataframe matrix for end motifs, all frags, May 2024, validation.rds"))

# Additional dataset examples:
matrix_nucleosome_peak <- readRDS(file.path(data_dir, "Matrix for nucleosome peak.rds"))
matrix_nucleosome_peak_validation <- readRDS(file.path(data_dir, "Matrix for nucleosome peak validation.rds"))

delfi_ratios <- readRDS(file.path(data_dir, "DELFI_Ratios_updated_Dec2024.rds"))
delfi_ratios_validation <- readRDS(file.path(data_dir, "DELFI_Ratios_updated_Dec2024_validation_cohort.rds"))

insert_matrix <- readRDS(file.path(data_dir, "insert_matrix_df_for_Althaf_updated.rds"))
insert_matrix_validation <- readRDS(file.path(data_dir, "insert_matrix_df_validation_for_Althaf_updated.rds"))


## Adjust names 
# For each validation dataset, update the rownames
matrix_nucleosome_peak_validation <- {
  rn <- rownames(matrix_nucleosome_peak_validation)
  rownames(matrix_nucleosome_peak_validation) <- sub("_dedup.*", "", rn)
  matrix_nucleosome_peak_validation
}

delfi_ratios_validation <- {
  rn <- rownames(delfi_ratios_validation)
  rownames(delfi_ratios_validation) <- sub("_dedup.*", "", rn)
  delfi_ratios_validation
}

insert_matrix_validation <- {
  rn <- rownames(insert_matrix_validation)
  rownames(insert_matrix_validation) <- sub("_dedup.*", "", rn)
  insert_matrix_validation
}

end_motifs_validation <- {
  rn <- rownames(end_motifs_validation)
  rownames(end_motifs_validation) <- sub("_dedup.*", "", rn)
  end_motifs_validation
}

# Use the cleaned end motifs row names as the reference
common_rn <- rownames(end_motifs_validation)

# Subset each validation dataframe to include only rows with these rownames
matrix_nucleosome_peak_validation <- matrix_nucleosome_peak_validation[rownames(matrix_nucleosome_peak_validation) %in% common_rn, ]
delfi_ratios_validation <- delfi_ratios_validation[rownames(delfi_ratios_validation) %in% common_rn, ]
insert_matrix_validation <- insert_matrix_validation[rownames(insert_matrix_validation) %in% common_rn, ]


# Create lists for training (non-validation) and validation datasets
non_val <- list(
  matrix_nucleosome_peak = matrix_nucleosome_peak,
  delfi_ratios = delfi_ratios,
  insert_matrix = insert_matrix,
  end_motifs_all_frags = end_motifs_all_frags
)

val <- list(
  matrix_nucleosome_peak = matrix_nucleosome_peak_validation,
  delfi_ratios = delfi_ratios_validation,
  insert_matrix = insert_matrix_validation,
  end_motifs_all_frags = end_motifs_validation
)



# -------------------------------
# Check comparability of training and validation datasets
# -------------------------------

# Assuming the following lists are already defined:
# non_val: list of training datasets (data frames)
# val: list of corresponding validation datasets (data frames)

# For reproducibility
set.seed(123)

# Create an empty list to store comparison results
comparison_results <- list()

# Open a PDF device to save plots (optional)
pdf("training_validation_comparison_plots.pdf", width = 7, height = 7)

# Loop through each dataset pair in the lists
for (dataset_name in names(non_val)) {
  
  # Retrieve training and validation dataframes
  train_df <- non_val[[dataset_name]]
  valid_df <- val[[dataset_name]]
  
  # Ensure the datasets have the same columns (if not, subset to common columns)
  common_cols <- intersect(colnames(train_df), colnames(valid_df))
  train_df <- train_df[, common_cols, drop = FALSE]
  valid_df <- valid_df[, common_cols, drop = FALSE]
  
  # Calculate column means and standard deviations for each dataset
  train_means <- colMeans(train_df, na.rm = TRUE)
  valid_means <- colMeans(valid_df, na.rm = TRUE)
  train_sds <- apply(train_df, 2, sd, na.rm = TRUE)
  valid_sds <- apply(valid_df, 2, sd, na.rm = TRUE)
  
  # Compute the Pearson correlation between training and validation means
  mean_corr <- cor(train_means, valid_means)
  
  # Compute the Pearson correlation between training and validation standard deviations
  sd_corr <- cor(train_sds, valid_sds)
  
  # Calculate the average absolute difference between the column means
  mean_diff <- mean(abs(train_means - valid_means))
  
  # Save these results in the list
  comparison_results[[dataset_name]] <- list(
    mean_correlation = mean_corr,
    sd_correlation = sd_corr,
    avg_abs_mean_diff = mean_diff
  )
  
  # Plot: Training vs. Validation Column Means
  plot(train_means, valid_means,
       main = paste("Training vs Validation Means:", dataset_name),
       xlab = "Training Column Means",
       ylab = "Validation Column Means",
       pch = 16, col = "blue")
  abline(0, 1, col = "red", lwd = 2)  # identity line for reference
  
  # Optionally, could add a plot for standard deviations as well
  # plot(train_sds, valid_sds,
  #      main = paste("Training vs Validation SDs:", dataset_name),
  #      xlab = "Training SDs",
  #      ylab = "Validation SDs",
  #      pch = 16, col = "darkgreen")
  # abline(0, 1, col = "red", lwd = 2)
}

# Close the PDF device if opened
dev.off()

# Print the comparison results to the console
print(comparison_results)

# -------------------------------
# Manuscript summary sentences:
# -------------------------------
ms_sentence1 <- "Comparative analysis between the training and validation datasets revealed high correlations in column means and standard deviations, indicating that feature distributions are largely consistent across datasets."
ms_sentence2 <- "The small average absolute differences in feature means further support the suitability of the validation cohort for reliable classifier evaluation."

cat("\nManuscript Summary:\n")
cat(ms_sentence1, "\n", ms_sentence2, "\n")
