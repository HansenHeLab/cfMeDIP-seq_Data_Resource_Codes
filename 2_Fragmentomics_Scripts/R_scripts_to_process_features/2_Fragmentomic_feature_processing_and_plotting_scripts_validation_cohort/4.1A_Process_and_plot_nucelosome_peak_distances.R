# Process and plot nucleosome peak distances

# authors: Dory Abelman and Derek Wong
# adapted from scripts by Derek Wong
# Original date: Aug 2023
# Updated for new samples Oct-Nov 2023
# Updated May 2024 for final set of adjusted prostate cancer samples

# Publication adapted from:
# Link: https://github.com/pughlab/TGL48_Uveal_Melanoma/blob/main/R_scripts_figures/Figure%202/Figure%202%20-%20Fragment%20Ratio.R
# Title: Integrated, Longitudinal Analysis of Cell-free DNA in Uveal Melanoma (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9973415/)

library(tidyverse)
library(dplyr)
library(matrixStats)
library(reshape2)
library(ggh4x)
library(ggpubr)
library(RColorBrewer)
library(factoextra) 


### Set paths
path <- "Nucleosome Peak/distance/"
outdir <- "Final Data Dec 2024/Nucelosome_peak/"
healthy_path <- "Nucleosome Peak/healthy_distance/"

### Find paths
data_files <- list.files(path, "\\.txt$", full.names = TRUE)
data_normal_files <- list.files(healthy_path, "\\.txt$", full.names = TRUE)

### Import data 
data_list <- lapply(data_files, read.delim)
data_normal_list <- lapply(data_normal_files, read.delim)

### Extract filenames without extension
data_filenames <- sapply(data_files, function(x) tools::file_path_sans_ext(basename(x)))
data_normal_filenames <- sapply(data_normal_files, function(x) tools::file_path_sans_ext(basename(x)))

### Add new column 'sample' to each dataframe in the list
data_list <- lapply(1:length(data_list), function(i) {
  df <- data_list[[i]]
  df$sample <- data_filenames[i]
  return(df)
})

data_normal_list <- lapply(1:length(data_normal_list), function(i) {
  df <- data_normal_list[[i]]
  df$sample <- data_normal_filenames[i]
  return(df)
})

### Convert list of data.frames to a single data.frame
data <- do.call(rbind, data_list)
data_normal <- do.call(rbind, data_normal_list)

### Format data
# Aggregating the values across all chromosomes, excluding 'distance' and 'sample'
data$aggregate_values <- rowSums(data[, 2:(ncol(data) - 1)])

# Reshaping the data to have a single row for each sample
# Select columns of interest and form shape
subset_data <- data %>% dplyr::select(sample, aggregate_values, distance)
reshaped_data_subset <- subset_data %>% 
  pivot_wider(names_from = distance, values_from = aggregate_values)

data_sum <- reshaped_data_subset

# Now for normal / healthy files
# Aggregating the values across all chromosomes, excluding 'distance' and 'sample'
data_normal$aggregate_values <- rowSums(data_normal[, 2:(ncol(data_normal) - 1)])

# Use a subset of the data to reshape - splitting by normal and cancer values

# Select a subset of the data to focus on samples, aggregate values, and distances
subset_data <- data_normal %>% dplyr::select(sample, aggregate_values, distance)

# Reshape the data using pivot_wider to spread aggregate values across columns based on distances
reshaped_data_subset_normal <- subset_data %>% 
  pivot_wider(names_from = distance, values_from = aggregate_values)

# Assign the reshaped data to a new data frame
data_normal_sum <- reshaped_data_subset_normal

# Extract unique sample names from the data
names <- unique(data$sample)

# Rename columns of data_sum to include 'sample' and a range of distance values from -1000 to 1000
colnames(data_sum) <- c("sample", -1000:1000)

# Filter columns of data_sum to only include 'sample' and distances from -300 to 300
data_sum <- data_sum[, colnames(data_sum) %in% c("sample", -300:300)] # take the regions of interest

# Normalize the aggregate values across each row by converting them into percentages of the row sum
data_sum[,2:ncol(data_sum)] <- data_sum[,2:ncol(data_sum)]/rowSums(data_sum[,2:ncol(data_sum)])*100

names <- unique(data_normal$sample) # Repeat extraction of unique sample names for the normal data subset
colnames(data_normal_sum) <- c("sample", -1000:1000) # Repeat column renaming for the data_normal_sum to match the new structure
data_normal_sum <- data_normal_sum[, colnames(data_normal_sum) %in% c("sample", -300:300)]  # Filter data_normal_sum to only include columns within the range of -300 to 300

# Normalize the values in data_normal_sum similarly by converting them into percentages of the row sum
data_normal_sum[,2:ncol(data_normal_sum)] <- data_normal_sum[,2:ncol(data_normal_sum)]/rowSums(data_normal_sum[,2:ncol(data_normal_sum)])*100

## Filter to only include patients of interest 
## Load and format metadata df 
# Now join to metadata df 
metadata_df <- read.csv("metadata_updated_Dec2024.csv")

## Annotate the metadata df
metadata_df$sample_id <- paste(metadata_df$sequencing_id, "_dedup", sep = "")
metadata_df$cancer_type_corrected_updated <- tolower(gsub(" ", "_", metadata_df$cancer_type))
metadata_df$cancer_type_corrected_updated <- as.character(metadata_df$cancer_type_corrected_updated)
metadata_df$cancer_type_corrected_updated <- gsub('blood_cancer', 'Aml', metadata_df$cancer_type_corrected_updated)
metadata_df$cancer_type_corrected_updated <- gsub('normal', 'healthy', metadata_df$cancer_type_corrected_updated)
metadata_df$cancer_type_title_case <- tools::toTitleCase(metadata_df$cancer_type_corrected_updated)
metadata_df$cancer_type_title_case  <- gsub("Lfs", "LFS", metadata_df$cancer_type_title_case)
metadata_df$cancer_type_title_case  <- gsub("Aml", "AML", metadata_df$cancer_type_title_case)

## Add nucleosome peak sample ID
metadata_df$nucleosome_peak_sampleid <- paste0(metadata_df$sequencing_id, "_dedup_peak_distance", sep = "")

# Filter to include only those in the metadata df 
data_normal_sum <- data_normal_sum %>%
  dplyr::filter(sample %in% metadata_df$nucleosome_peak_sampleid)

data_sum <- data_sum %>%
  dplyr::filter(sample %in% metadata_df$nucleosome_peak_sampleid)

## Include normal in data_sum to see how variance is affected there
data_sum <- bind_rows(data_sum, data_normal_sum)

## Only keep validation and resplit 
data_sum <- data_sum %>%
  dplyr::filter(
    sample %in% metadata_df$nucleosome_peak_sampleid &
      metadata_df$validation[match(sample, metadata_df$nucleosome_peak_sampleid)] == 1 &
      metadata_df$type[match(sample, metadata_df$nucleosome_peak_sampleid)] == "PE" &
      metadata_df$cancer_type_title_case[match(sample, metadata_df$nucleosome_peak_sampleid)] != "Liver_cancer"
  )


## Resplit normal seperately 
data_normal_sum <- data_sum %>%
  dplyr::filter(
    metadata_df$cancer_type[match(sample, metadata_df$nucleosome_peak_sampleid)] == "Normal")
      

### Calculate the healthy median
normal_median <- colMedians(as.matrix(data_normal_sum[,2:ncol(data_normal_sum)]))
normal_sd <- colSds(as.matrix(data_normal_sum[,2:ncol(data_normal_sum)]))


### Calculate the cancer median
## Initialize empty list to hold data frames for each cancer type
list_of_data_frames <- list()

cancer_types <- unique(metadata_df_validation$cancer_type_corrected_updated)

## Make edit 
metadata_df$cancer_type_corrected_updated <- gsub("aml", "AML", metadata_df$cancer_type_corrected_updated)

# Get zscore for each type
for (cancer_type in cancer_types) {
  
  # Calculate median and standard deviation
  type_median <- colMedians(as.matrix(
    data_sum[data_sum$sample %in%
               metadata_df$nucleosome_peak_sampleid[metadata_df$cancer_type_corrected_updated == cancer_type],
             2:ncol(data_sum)]))
  type_sd <- colSds(as.matrix(
    data_sum[data_sum$sample %in%
               metadata_df$nucleosome_peak_sampleid[metadata_df$cancer_type_corrected_updated == cancer_type],
             2:ncol(data_sum)]))
  
  # Calculate p-values
  p_type <- ks.test(normal_median[210:390], type_median[210:390])$p.value
  
  # Calculate z-scores
  z_normal <- (normal_median - normal_median) / normal_sd
  z_type <- (type_median - normal_median) / normal_sd
  
  type_dif <- type_median - normal_median
  
  # Make a data frame for this cancer type and add it to the list
  temp_df <- data.frame(
    length = c(-300:300),
    score = z_type,
    cancer_status = cancer_type
  )
  list_of_data_frames[[cancer_type]] <- temp_df
}

# Combine all data frames into one
data_plot <- do.call(rbind, list_of_data_frames)


## Now melt
### Make plotting table (individuals)
data_melt <- reshape2::melt(data_sum, id = "sample")
data_melt$variable <- as.numeric(as.character(data_melt$variable))
data_melt$diag <- "Cancer"
data_melt <- merge(data_melt, metadata_df, by.x = "sample", by.y = "nucleosome_peak_sampleid")

### Make plotting table (difference to healthy)
list_of_data_frames <- list()

# Initialize list to hold p-values
p_values_list <- list()

# Get diff and p-value for each type
for (cancer_type in cancer_types) {
  
  # Calculate median and standard deviation
  type_median <- colMedians(as.matrix(
    data_sum[data_sum$sample %in%
               metadata_df$nucleosome_peak_sampleid[metadata_df$cancer_type_corrected_updated == cancer_type],
             2:ncol(data_sum)]))
  
  type_sd <- colSds(as.matrix(
    data_sum[data_sum$sample %in%
               metadata_df$nucleosome_peak_sampleid[metadata_df$cancer_type_corrected_updated == cancer_type],
             2:ncol(data_sum)]))
  
  # Calculate p-values
  p_type <- ks.test(normal_median[210:390], type_median[210:390])$p.value
  
  # Store p-value in the list
  p_values_list[[cancer_type]] <- round(p_type, 4)
  
  # Calculate z-scores
  z_type <- (type_median - normal_median) / normal_sd
  
  # Calculate the difference from normal_median
  type_dif <- type_median - normal_median
  
  # Make a data frame for this cancer type and add it to the list
  temp_df <- data.frame(
    length = c(-300:300),
    score = type_dif,
    cancer_status = cancer_type
  )
  list_of_data_frames[[cancer_type]] <- temp_df
}

data_diff <- do.call(rbind, list_of_data_frames)

### Merge plots
data_melt <- data_melt[, c("sample", "variable", "value", "cancer_type_corrected_updated")]
colnames(data_melt) <- c("sample", "length", "score", "cancer_status")
data_melt$analysis <- "Frequency"

data_plot$sample <- data_plot$cancer_status
data_plot$analysis <- "Z-score"

data_diff$sample <- data_diff$cancer_status
data_diff$analysis <- "Difference"

data_plot <- bind_rows(data_plot, data_melt, data_diff)
data_plot$analysis <- factor(data_plot$analysis, levels = c("Frequency", "Z-score", "Difference"))
data_plot$cancer_status <- factor(data_plot$cancer_status, levels = c("head_and_neck_cancer", "breast_cancer", "ovarian_cancer", "melanoma",    "mixed_cancer", "liver_cancer", "healthy"))

### Set hlines
data_lines <- data.frame(cancer_status = rep(c("head_and_neck_cancer", "breast_cancer", "ovarian_cancer", "melanoma",    "mixed_cancer", "liver_cancer", "healthy"), 2),
                         analysis = c(rep("Difference", 10), rep("Z-score", 10)),
                         line = c(rep(0, 20)))
data_lines$analysis <- factor(data_lines$analysis, levels = c("Frequency", "Z-score", "Difference"))
data_lines$cancer_status <- factor(data_lines$cancer_status, levels = c("head_and_neck_cancer", "breast_cancer", "ovarian_cancer", "melanoma",    "mixed_cancer", "liver_cancer", "healthy"))
data_median <- data.frame(length = c(-300:300),
                          median = normal_median,
                          analysis = "Frequency")
data_median$analysis <- factor(data_median$analysis, levels = c("Frequency", "Z-score", "Difference"))

### Set p-value
# Create a data frame for p-values
data_pvalue <- data.frame(
  cancer_status = names(p_values_list),
  analysis = rep("Difference", length(p_values_list)),
  value = unlist(p_values_list)
)
data_pvalue$cancer_status <- factor(data_pvalue$cancer_status, levels = c("brain_cancer", "head_and_neck_cancer", "AML", "prostate_cancer", "lung_cancer", "lfs_survivor", "lfs_previvor", "lfs_positive",  "eye_cancer", "healthy"))
data_pvalue$analysis <- factor(data_pvalue$analysis, levels = c("Frequency", "Z-score", "Difference"))

### Set theme
theme <- theme(plot.title = element_text(hjust = 0.5, size = 16), 
               axis.line = element_line(colour = "black"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_rect(fill = NA),
               panel.background = element_blank(),
               legend.position = "bottom",
               legend.key = element_rect(fill = "white"),
               legend.title = element_text(size = 12),
               legend.text = element_text(size = 12),
               strip.background = element_blank(),
               strip.text = element_text(size = 13),
               axis.text = element_text(size = 13),
               axis.title = element_text(size = 13))


### Plot curves

## Get colors 

# With underscores 
col_cancer_type_split <- c('Healthy' = '#33a02c',
                     'Brain_cancer' = '#1f78b4',
                     'Lung_cancer' = '#b2df8a',
                     'Prostate_cancer' = '#a6cee3',
                     'AML' = '#fb9a99',
                     'Pancreatic_cancer' = '#e31a1c',
                     'Eye_cancer' = '#fdbf6f',
                     'Head_and\nneck_cancer' = '#ff7f00',
                     'Breast_cancer' = '#cab2d6',
                     'Colorectal_cancer' = '#6a3d9a',
                     'Bladder_cancer' = '#fb6a4a',
                     'Renal_cancer' = '#b15928',
                     'LFS_survivor' = '#bdbdbd',
                     'LFS_previvor' = '#969696',
                     'LFS_positive' = '#737373',
                     'Liver_cancer' = '#ffeda0',
                     'Melanoma' = '#41b6c4',
                     'Mixed_cancer' = '#e31a1c', # Same as pancreatic (mixed nature)
                     'Ovarian_cancer' = '#ffffb3')


## Plot curve 
data_plot_filtered <- data_plot %>% dplyr::filter(!is.na(cancer_status))


## Order by the order established in the insert size prop short
## First make title case
data_plot_filtered_reordered <- data_plot_filtered
data_plot_filtered_reordered$cancer_type_title_case <-as.character(data_plot_filtered_reordered$cancer_status)
data_plot_filtered_reordered$cancer_type_title_case <- sub("^(.)", "\\U\\1", data_plot_filtered_reordered$cancer_type_title_case, perl=TRUE)
data_plot_filtered_reordered$cancer_type_title_case  <- gsub("Lfs", "LFS", data_plot_filtered_reordered$cancer_type_title_case)
data_plot_filtered_reordered$cancer_type_title_case  <- gsub("Uveal_melanoma", "Eye_cancer", data_plot_filtered_reordered$cancer_type_title_case)

## Plot without difference 
data_plot_filtered_reordered <- data_plot_filtered_reordered %>% dplyr::filter(analysis != "Difference")
data_lines_no_diff <- data_lines %>% dplyr::filter(analysis != "Difference")

## Order based on prop short frags

# Replace the specific label in your data frame that doesn't fit
data_plot_filtered_reordered$cancer_type_title_case <- gsub("Head_and_neck_cancer", "Head_and\nneck_cancer", data_plot_filtered_reordered$cancer_type_title_case)

# Assuming ordered_cancer_types is defined as follows
ordered_cancer_types_peak <- c("Healthy", "Mixed_cancer", "Breast_cancer", "Ovarian_cancer", "Head_and\nneck_cancer", "Melanoma")

# Ensure that the cancer_type_title_case column is a factor and set the levels according to ordered_cancer_types
data_plot_filtered_reordered$cancer_type_title_case <- factor(data_plot_filtered_reordered$cancer_type_title_case, levels = ordered_cancer_types_peak)

fig_nuc_peaks <- ggplot(data_plot_filtered_reordered) +
  geom_line(aes(length, score, color = cancer_type_title_case, group = sample), linewidth = 0.5, alpha = 0.5) +
  geom_line(data = data_median, aes(length, median), linewidth = 0.5) +
  # geom_text(data = data_pvalue, aes(-200, -0.005, label = paste0("p-value\n", value)), size = 4, vjust = 1) +
  geom_hline(data = data_lines_no_diff, aes(yintercept = line), linetype = "dashed", size = 0.5) +
  geom_vline(xintercept = c(-83, 83), linetype = "dashed", size = 0.5) +
  xlab("Distance") + 
  ylab("") +
  labs(color = "") +
  facet_grid(analysis~cancer_type_title_case, scales = "free_y") +
  scale_color_manual(values = col_cancer_type_split, name = "Cancer type") + # Add name for legend title
  ggtitle(paste0("Distance to Nearest Peak")) + 
  theme # Enlarge the title
fig_nuc_peaks

fig <- fig_nuc_peaks
ggsave(file.path(outdir, paste0("nucleosome_peaks_167bp cfMeDIP with pvalue, updated normals, Mar 2025, updated no pvalue, 3, validation.pdf")), fig, width = 15, height = 7, units = "in", dpi = 500)
saveRDS(data_plot_filtered_reordered, file = "Dataframe for nucleosome peaks figure Dec 2024 validation.rds")



### Compare between patients
## Get batch data 

get_comparison_data <- function(type1, type2, data_plot) {
  data1 <- data_plot[data_plot$analysis == "Z-score" & data_plot$cancer_status == type1,]
  data2 <- data_plot[data_plot$analysis == "Z-score" & data_plot$cancer_status == type2,]
  
  data_z <- data.frame(length = data1$length,  # Assuming length is the same for all
                       score_diff = data1$score - data2$score,
                       comp = paste(type1, "vs", type2))
  
  p_value <- ks.test(data1$score, data2$score)$p.value
  
  if (p_value > 0.01) {
    p_value <- round(p_value, 2)
  } else {
    p_value <- sprintf("%.2e", p_value)
  }
  
  data_p <- data.frame(value = p_value, comp = paste(type1, "vs", type2))
  
  return(list(data_z = data_z, data_p = data_p))
}

# Assuming data_plot is your merged data frame containing all types
list_data_z <- list()
list_data_p <- list()



# Loop through all unique pairs of types
for (i in 1:(length(cancer_types) - 1)) {
  for (j in (i + 1):length(cancer_types)) {
    type1 <- cancer_types[i]
    type2 <- cancer_types[j]
    
    result <- get_comparison_data(type1, type2, data_plot)
    
    list_data_z[[paste(type1, "vs", type2)]] <- result$data_z
    list_data_p[[paste(type1, "vs", type2)]] <- result$data_p
  }
}

# Combine all the individual data frames into single data frames
data_z_all <- do.call(rbind, list_data_z)
data_p_all <- do.call(rbind, list_data_p)


fig <- ggplot(data_z_all) +
  geom_density(aes(score_diff), fill = "grey75", color = NA) +
  geom_text(data = data_p_all, aes(x = 0, y = Inf, label = paste0("KS-test = ", value)), vjust = 2) +
  xlab("Delta (Z-score)") + 
  ylab("Frequency") +
  ggtitle(paste0("Z-score Comparison 167 basepair fragments")) +
  facet_grid(.~comp) +
  theme +
  theme(legend.position = "none")
fig
ggsave(file.path(outdir, paste0("nucleosome_peaks_zscore_comparison, updated normals, May 2024, updated.pdf")), fig, width = 49, height = 8)

## Make the above for each comparison set 
# Loop through each unique cancer type
for (cancer_type in cancer_types) {
  
  # Subset the data to only include comparisons involving this cancer type
  data_z_subset <- data_z_all[grep(cancer_type, data_z_all$comp), ]
  data_p_subset <- data_p_all[grep(cancer_type, data_p_all$comp), ]
  
  # Check if the subset is empty, continue to next iteration if it is
  if (nrow(data_z_subset) == 0) next
  
  # Create a ggplot object
  fig <- ggplot(data_z_subset) +
    geom_density(aes(score_diff), fill = "grey75", color = NA) +
    geom_text(data = data_p_subset, aes(x = 0, y = Inf, label = paste0("KS-test = ", value)), vjust = 2) +
    xlab("Delta (Z-score)") + 
    ylab("Frequency") +
    ggtitle(paste0("Z-score Comparison 167 basepair fragments: ", cancer_type)) +
    facet_grid(.~comp) +
    theme +
    theme(legend.position = "none")
  
  # Show the plot (optional)
  print(fig)
  
  # Save the plot
  ggsave(file.path(outdir, paste0("nucleosome_peaks_zscore_comparison_", cancer_type, "updated May 2024.pdf")), fig, width = 25, height = 8)
  
}



## Now do the same, but for the differences 
get_comparison_data_diff <- function(type1, type2, data_plot) {
  data1 <- data_plot[data_plot$analysis == "Difference" & data_plot$cancer_status == type1,]
  data2 <- data_plot[data_plot$analysis == "Difference" & data_plot$cancer_status == type2,]
  
  data_diff <- data.frame(length = data1$length,  # Assuming length is the same for all
                          score_diff = data1$score - data2$score,
                          comp = paste(type1, "vs", type2))
  
  p_value <- ks.test(data1$score, data2$score)$p.value
  
  if (p_value > 0.01) {
    p_value <- round(p_value, 2)
  } else {
    p_value <- sprintf("%.2e", p_value)
  }
  
  data_p_diff <- data.frame(value = p_value, comp = paste(type1, "vs", type2))
  
  return(list(data_diff = data_diff, data_p_diff = data_p_diff))
}

# Assuming data_plot is your merged data frame containing all types
list_data_diff <- list()
list_data_p_diff <- list()


# Loop through all unique pairs of types
for (i in 1:(length(cancer_types) - 1)) {
  for (j in (i + 1):length(cancer_types)) {
    type1 <- cancer_types[i]
    type2 <- cancer_types[j]
    
    result <- get_comparison_data_diff(type1, type2, data_plot)
    
    list_data_diff[[paste(type1, "vs", type2)]] <- result$data_diff
    list_data_p_diff[[paste(type1, "vs", type2)]] <- result$data_p_diff
  }
}

# Combine all the individual data frames into single data frames
data_diff_all <- do.call(rbind, list_data_diff)
data_p_diff_all <- do.call(rbind, list_data_p_diff)

fig <- ggplot(data_diff_all) +
  geom_density(aes(score_diff), fill = "grey75", color = NA) +
  geom_text(data = data_p_diff_all, aes(x = -0, y = Inf, label = paste0("KS-test = ", value)), vjust = 2) +
  xlab("Delta (Difference)") + 
  ylab("Frequency") +
  ggtitle(paste0("Absolute Difference Comparison in 167bp fragments")) +
  facet_grid(.~comp) +
  theme +
  theme(legend.position = "none")
fig
ggsave(file.path(outdir, paste0("nucleosome_peaks_difference_absolute, May 2024 updated.pdf")), fig, width = 49, height = 8)


# Loop through each 
for (cancer_type in cancer_types) {
  
  # Subset the data to only include comparisons involving this cancer type
  data_diff_subset <- data_diff_all[grep(cancer_type, data_diff_all$comp), ]
  data_p_subset <- data_p_diff_all[grep(cancer_type, data_p_diff_all$comp), ]
  
  # Check if the subset is empty, continue to next iteration if it is
  if (nrow(data_diff_subset) == 0) next
  
  # Create a ggplot object
  fig <- ggplot(data_diff_subset) +
    geom_density(aes(score_diff), fill = "grey75", color = NA) +
    geom_text(data = data_p_subset, aes(x = 0, y = Inf, label = paste0("KS-test = ", value)), vjust = 2) +
    xlab("Delta (Difference)") + 
    ylab("Frequency") +
    ggtitle(paste0("Absolute Difference Comparison in 167 basepair fragments: ", cancer_type)) +
    facet_grid(.~comp) +
    theme +
    theme(legend.position = "none")
  
  # Show the plot (optional)
  print(fig)
  
  # Save the plot
  ggsave(file.path(outdir, paste0("nucleosome_peaks_absolute_diff_comparison_", cancer_type, ", May 2024 updated.pdf")), fig, width = 25, height = 8)
  
}




### Export the data 
Nucleosome_peak_matrix <- data_sum # previously added with data_normal

# Identify missing samples if any 
missing_samples <- setdiff(metadata_df$nucleosome_peak_sampleid, Nucleosome_peak_matrix$sample)
print(missing_samples)

# remove sample column but keep as rowname 
# Store row names
row_names_to_keep <- Nucleosome_peak_matrix$sample

# Remove the 'sample' column
Nucleosome_peak_matrix <- Nucleosome_peak_matrix %>% dplyr::select(-sample)

# Reassign row names
rownames(Nucleosome_peak_matrix) <- row_names_to_keep

## Export the dataframe 
saveRDS(Nucleosome_peak_matrix, file = "Matrix for nucleosome peak.rds")



## Export the PCAs on the dataframe 

# Perform PCA
pca_result_nucleosome_peak <- prcomp(Nucleosome_peak_matrix, center = TRUE, scale. = TRUE)

# Determine the number of components to retain
# Visualize the variance explained by each principal component
fviz_eig(pca_result_nucleosome_peak, addlabels = TRUE, ylim = c(0, 25))

# Extract the proportion of variance explained by each component
explained_variance <- summary(pca_result_nucleosome_peak)$importance[2,]

# Calculate cumulative variance
cum_variance <- cumsum(explained_variance)

# Plot cumulative explained variance
plot(cum_variance, xlab = "Number of Principal Components", ylab = "Cumulative Explained Variance", type = "b")
abline(h = 0.90, col = "red", lty = 2)
abline(h = 0.95, col = "purple", lty = 2)
abline(h = 0.99, col = "blue", lty = 2)

## Plot nicer 

# Create a data frame for ggplot
cum_variance_df <- tibble(Principal_Component = 1:length(cum_variance),
                          Cumulative_Variance = cum_variance)

# Create the plot
p <- ggplot(cum_variance_df, aes(x = Principal_Component, y = Cumulative_Variance)) +
  geom_line(color = "black") +
  geom_point(color = "black") +
  geom_hline(yintercept = 0.90, color = "red", linetype = "dashed") +
  geom_hline(yintercept = 0.95, color = "purple", linetype = "dashed") +
  geom_hline(yintercept = 0.99, color = "blue", linetype = "dashed") +
  labs(title = "Cumulative Explained Variance by Principal Components",
       subtitle = "Nucleosome peak. Red: 90%, Purple: 95%, Blue: 99%",
       x = "Number of Principal Components",
       y = "Cumulative Explained Variance") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10))

# Print the plot
print(p)
ggsave(filename = file.path(outdir, "Cumulative Explained Variance by Principal Components nucleosome peak.png"), plot = p, width = 7, height = 5, units = "in", dpi = 300)


# Find the number of components that explain at least 90% of the variance
num_components <- which(cum_variance >= 0.90)[1]

# Extract the principal components
pca_df_nucleosome_peak <- as.data.frame(pca_result_nucleosome_peak$x[, 1:num_components])
rownames(pca_df_nucleosome_peak) <- rownames(Nucleosome_peak_matrix)

# Create the variance explained plot with customized theme
variance_plot <- fviz_eig(pca_result_nucleosome_peak, addlabels = TRUE, ylim = c(0, 30)) +
  theme_classic() +
  ggtitle("Scree plot for nucleosome peak") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

# Save the plot to a file
ggsave(filename = file.path(outdir, "PCA_Variance_Explained_Nucleosome_peak.png"), plot = variance_plot, width = 10, height = 6, units = "in", dpi = 300)

# Save the PCA results to a file 
write.csv(pca_df_nucleosome_peak, file.path(outdir, "PCA_results_Nucleosome_peak.csv"), row.names = TRUE)
saveRDS(pca_df_nucleosome_peak, file = "Nucleosome_peaks_updated_June2024_PCAs_90_percent_variance.rds")
saveRDS(pca_result_nucleosome_peak, file = "Nucleosome_peaks_PCA_results.rds")

## Now calculate the z-score dataframe 

# Filter metadata_df for 'Normal' subtype
normal_ids <- metadata_df %>% 
  dplyr::filter(cancer_type_corrected_updated == "healthy") %>% 
  dplyr::select(nucleosome_peak_sampleid)

# Keep only rows in Nucleosome_peak_matrix that match the normal_ids
normal_nucleosome_peak <- Nucleosome_peak_matrix %>% 
  dplyr::filter(row.names(.) %in% normal_ids$nucleosome_peak_sampleid)

# Transpose matrix 
normal_nucleosome_peak <- t(normal_nucleosome_peak)

### Make healthy median
Normal_median <- rowMedians(as.matrix(normal_nucleosome_peak), na.rm = T)
Normal_sd <- rowSds(as.matrix(normal_nucleosome_peak), na.rm = T)


## Now get cancer median across all cancers
# Filter metadata_df for 'Normal' subtype
cancer_ids <- metadata_df %>% 
  dplyr::filter(cancer_type_corrected_updated != "healthy") %>% 
  dplyr::select(nucleosome_peak_sampleid)

# Keep only rows in Nucleosome_peak_matrix that match the normal_ids
cancer_nucleosome_peak <- Nucleosome_peak_matrix %>% 
  dplyr::filter(row.names(.) %in% cancer_ids$nucleosome_peak_sampleid)

# Transpose matrix 
cancer_nucleosome_peak <- t(cancer_nucleosome_peak)

### Make cancer median
cancer_median <- rowMedians(as.matrix(cancer_nucleosome_peak))
cancer_sd <- rowSds(as.matrix(cancer_nucleosome_peak))

## Get distance from normal median 
Cancer_nucleosome_peak_ratio <- (cancer_nucleosome_peak - Normal_median)/Normal_sd
Normal_nucleosome_peak_ratio <- (normal_nucleosome_peak - Normal_median)/Normal_sd

### Calculate Normal Z-scores (Normal relative to Normal)
Normal_sum <- colSums(abs(Normal_nucleosome_peak_ratio), na.rm = TRUE)
Normal_sum_median <- median(Normal_sum, na.rm = TRUE)
Normal_MAD <- mad(Normal_sum, na.rm = TRUE) #Median Absolute Deviation

Normal_zscores <- abs((Normal_sum - Normal_sum_median))/Normal_MAD
Normal_mean <- mean(Normal_zscores)
z_limit <- quantile(Normal_zscores, 0.9)

### Calculate Cancer Z-scores (Cancer relative to Normal)
Cancer_sum <- colSums(abs(Cancer_nucleosome_peak_ratio), na.rm = TRUE)
Cancer_zscore <- (Cancer_sum - Normal_sum_median)/Normal_MAD

### Make Z-score matrix
all_zscores <- c(Normal_zscores,Cancer_zscore)
zscore_df <- as.data.frame(all_zscores)
zscore_df$limit <- z_limit
zscores_mat_nucleosome_peak <- as.matrix(zscore_df) # save

### Set Normal median 
lower <- Normal_median - Normal_sd
upper <- Normal_median + Normal_sd
Normal_median <- as.matrix(cbind(Normal_median, lower, upper))

## Save the zscores to a seperate dataframe 
Nucleosome_peak_zscores <- zscore_df 
Nucleosome_peak_zscores <- rownames_to_column(Nucleosome_peak_zscores, "Sample")

# Join with metadata df
Nucleosome_peak_zscores <- left_join(Nucleosome_peak_zscores, metadata_df, by = c("Sample" = "nucleosome_peak_sampleid"))

### Save objects
saveRDS(Normal_ratio, file=file.path(outdir, paste0("Normal_ratio_Nucelosome_peak.rds")))
saveRDS(Cancer_ratio, file=file.path(outdir, paste0("Cancer_ratio_Nucelosome_peak.rds")))
saveRDS(Normal_median, file=file.path(outdir, paste0("Normal_median_Nucelosome_peak.rds")))
saveRDS(zscores_mat_nucleosome_peak, file = file.path(outdir, paste0("Zscore_mat_Nucelosome_peak.rds")))
saveRDS(Nucleosome_peak_zscores, file = "Zscores df for Nucelosome peak.rds")




## Export metadata df with endomotif_name also included to have everything together
metadata_df$endmotif_name <- paste(metadata_df$sequencing_id, "_dedup_motifs", sep = "")

saveRDS(metadata_df, file = "Metadata df Dec 2024 all samples.rds")





## Export the final dataframe 
write.table(Nucleosome_peak_matrix, file = "Nucleosome_peak_matrix_dataframe.txt", sep = "\t", row.names = TRUE, quote = F)
