## Running classifier on Yong's methylation signature 

## Load RDS file for paired end samples
Paired_end_methylation_sig <- readRDS("/Users/dabelman/Documents/Thesis_work/R/M4/Projects/Fragmentation_CfMeDIP_Yong/Methylation data from Yong/Common_pan_cancer_hyper_bins_adjusted_cnt_in_PE.RDS")

# Transpose the matrix
Paired_end_methylation_sig <- t(Paired_end_methylation_sig)

# Add '_dedup' to the end of every row name
rownames(Paired_end_methylation_sig) <- paste0(rownames(Paired_end_methylation_sig), "_dedup")
Paired_end_methylation_sig <- as.matrix(Paired_end_methylation_sig)

## Now work on dimensionality reduction

### This is the original way originally run where we take the columns that provide the top 2.5% of features
# Run PCA
pca_result <- prcomp(Paired_end_methylation_sig, center = FALSE, scale. = FALSE)

# Extract loadings
loadings <- abs(pca_result$rotation)

# Calculate the contribution of each region to the first few principal components
# Here, we sum across the first few PCs as an example
region_contributions <- rowSums(loadings[, 1:2])  # Adjust based on how many PCs you want to consider

# Identify the top 2.5% of regions
top_2_5_percent_index <- sort(region_contributions, decreasing = TRUE)[1:(length(region_contributions) * 0.025)]

# Get the names of the top 5% contributing regions
top_regions <- names(top_2_5_percent_index)

# Subset the matrix to include only columns in top_regions
Paired_end_methylation_sig_subset <- Paired_end_methylation_sig[, top_regions]

# Export it
saveRDS(Paired_end_methylation_sig, file = "Paired end methylation signatures, Jan 2024.rds")
saveRDS(Paired_end_methylation_sig_subset, file = "Paired end methylation signatures subset top 2.5% of regions, Jan 2024.rds")



## Now this time export just the PCAs
# Perform PCA
pca_result <- prcomp(Paired_end_methylation_sig, center = FALSE, scale. = FALSE)
pca_result_methylation <- pca_result
# Determine the number of components to retain
# Visualize the variance explained by each principal component
fviz_eig(pca_result, addlabels = TRUE, ylim = c(0, 70))

# Extract the proportion of variance explained by each component
explained_variance <- summary(pca_result_methylation)$importance[2,]
cum_variance <- cumsum(explained_variance)

# Find the number of components that explain at least 97% of the variance
num_components <- which(cum_variance >= 0.90)[1]

# Extract the principal components
pca_df_cfmedip <- as.data.frame(pca_result$x[, 1:num_components])
rownames(pca_df_cfmedip) <- rownames(Paired_end_methylation_sig)

# Create the variance explained plot with customized theme
variance_plot <- fviz_eig(pca_result, addlabels = TRUE, ylim = c(0, 70)) +
  theme_classic() +
  ggtitle("Scree plot for CfMeDIP methylation signature") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

# Save the plot to a file
ggsave(filename = file.path(outdir, "PCA_Variance_Explained_methylation_sig.png"), plot = variance_plot, width = 10, height = 6, units = "in", dpi = 300)

# Save the PCA results to a file 
write.csv(pca_df_cfmedip, file.path(outdir, "PCA_results_methylation.csv"), row.names = TRUE)
saveRDS(pca_df_cfmedip, file = "Methylation_sig_updated_June2024_PCAs_90_percent_variance.rds")
saveRDS(pca_result_methylation, file = "PCA results for methylation.rds")

## Now make an integrated classifier that takes any PCA explainnig more than 5% of the variance?

