# Delfi medians - 5Mb - heatmaps

# author: Dory Abelman
# adapted from Derek Wong
# Original date: Aug 2023
# Updated for new samples Oct-Nov 2023
# Updated May 2024 for final set of adjusted prostate cancer samples

# Publication adapted from:
# Link: https://github.com/pughlab/TGL48_Uveal_Melanoma/blob/main/R_scripts_figures/Figure%202/Figure%202%20-%20Fragment%20Ratio.R
# Title: Integrated, Longitudinal Analysis of Cell-free DNA in Uveal Melanoma (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9973415/)

# Load packages
library(tidyverse)
library(GenomicRanges)
library(ComplexHeatmap)
library(circlize)
library(matrixStats)
library(reshape2)
library(RColorBrewer)
library(factoextra)  # For visualizing PCA results


# Read in R objects
filedir <- '~/Documents/Thesis_work/R/M4/Projects/Fragmentation_CfMeDIP_Yong/Final Data May 2024/DELFI/'

basedir <- '~/Documents/Thesis_work/R/M4/Projects/Fragmentation_CfMeDIP_Yong/Final Data May 2024/DELFI/'
datasets <- "~/Documents/Thesis_work/R/M4/Projects/Fragmentation_CfMeDIP_Yong/Final Data May 2024/DELFI/"
outdir <- '~/Documents/Thesis_work/R/M4/Projects/Fragmentation_CfMeDIP_Yong/Final Data May 2024/DELFI/output'

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

## These are loaded from the 'Sample medians updated May 2024.R' file, which should be run prior
load(file=file.path(datasets, paste0("tcge.normal.median_fr_df.hg38.rda")))
load(file=file.path(datasets, paste0("tcge.cancer.median_fr_df.hg38.rda")))

# Add the sample metadata to the file 
# Adding a new column sample_id with _deduped at the end of sample_name
metadata_df <- read.csv("Updated metadata file/TCGE_cfDM_samples_only_fixed_after_QC_PE.csv")
metadata_df$sample_id <- paste(metadata_df$sequencing_id, "_dedup", sep = "")


### Reformat data to calculate zscores
Cancer_ratio <- dcast(df.fr.tcge.cancer, seqnames + arm + start + end + bin ~ sample_id, value.var="ratio_centered")
data_chr <- Cancer_ratio[, c("seqnames", "arm", "start", "end", "bin")] # Make chromosome data frame
row.names(Cancer_ratio) <- with(Cancer_ratio, bin)
row.names(data_chr) <- with(data_chr, bin)
Cancer_ratio <- as.matrix(Cancer_ratio[, -c(1:5)])

## Get ratios and calculate genome-wide zscores
Normal_ratio <- dcast(df.fr.tcge.normal, seqnames + arm + start + end + bin ~ sample_id, value.var="ratio_centered") 
row.names(Normal_ratio) <- with(Normal_ratio, bin)
Normal_ratio <- as.matrix(Normal_ratio[, -c(1:5)])

### Calculate the Normal median
Normal_median <- rowMedians(Normal_ratio, na.rm = TRUE)
Normal_sd <- rowSds(Normal_ratio, na.rm = TRUE)

### Calculate the distance from the Normal median
Cancer_ratio <- (Cancer_ratio - Normal_median)/Normal_sd
Normal_ratio <- (Normal_ratio - Normal_median)/Normal_sd

### Calculate Normal Z-scores (Normal relative to Normal)
Normal_sum <- colSums(abs(Normal_ratio), na.rm = TRUE)
Normal_sum_median <- median(Normal_sum, na.rm = TRUE)
Normal_MAD <- mad(Normal_sum, na.rm = TRUE) #Median Absolute Deviation

Normal_zscores <- abs((Normal_sum - Normal_sum_median))/Normal_MAD
Normal_mean <- mean(Normal_zscores)
z_limit <- quantile(Normal_zscores, 0.9)

### Calculate Cancer Z-scores (Cancer relative to Normal)
Cancer_sum <- colSums(abs(Cancer_ratio), na.rm = TRUE)
Cancer_zscore <- (Cancer_sum - Normal_sum_median)/Normal_MAD

### Make Z-score matrix
all_zscores <- c(Normal_zscores,Cancer_zscore)
zscore_df <- as.data.frame(all_zscores)
zscore_df$limit <- z_limit
zscores_mat <- as.matrix(zscore_df)
zscores_mat_delfi <- as.matrix(zscore_df)

### Set Normal median 
lower <- Normal_median - Normal_sd
upper <- Normal_median + Normal_sd
Normal_median <- as.matrix(cbind(Normal_median, lower, upper))
row.names(Normal_median) <- row.names(data_chr)

### Set alteration frequencies
data_freq <- rowMeans2(Cancer_ratio, na.rm = TRUE)

## Save the zscores to a seperate dataframe 
Delfi_df_zscore <- zscore_df 
Delfi_df_zscore <- rownames_to_column(Delfi_df_zscore, "Sample")

# Join with metadata df
Delfi_df_zscore <- left_join(Delfi_df_zscore, metadata_df, by = c("Sample" = "sample_id"))

### Save objects
saveRDS(Normal_ratio, file=file.path(outdir, paste0("Normal_ratio_DELFI.rds")))
saveRDS(Cancer_ratio, file=file.path(outdir, paste0("Cancer_ratio_DELFI.rds")))
saveRDS(Normal_median, file=file.path(outdir, paste0("Normal_median_DELFI.rds")))
saveRDS(zscores_mat_delfi, file = file.path(outdir, paste0("Zscore_mat_DELFI.rds")))
saveRDS(Delfi_df_zscore, file = file.path(outdir, paste0("Zscore_df_DELFI.rds")))

### Prepare for plotting
plot_df <- cbind(Normal_ratio,Cancer_ratio)
plot_df <- t(plot_df)


## Prepare to make heatmap
data_samples <- data.frame(sample_id = rownames(plot_df))

# Get order 
# Put in correct order
## Arrange by z-score and keep within group 
tmp <- tibble(sample_id = rownames(zscore_df), all_zscores = zscore_df$all_zscores)

## Create new column for ease of naming
# Replace spaces with underscores and convert to lowercase, ensure is charachter
metadata_df$cancer_type_corrected_updated <- tolower(gsub(" ", "_", metadata_df$cancer_type))
metadata_df$cancer_type_corrected_updated <- as.character(metadata_df$cancer_type_corrected_updated)

## Edit naming to correspond with Yong's update 
# Use gsub to replace values in the cancer_type_corrected_updated column
metadata_df$cancer_type_corrected_updated <- gsub('blood_cancer', 'AML', metadata_df$cancer_type_corrected_updated)
metadata_df$cancer_type_corrected_updated <- gsub('normal', 'healthy', metadata_df$cancer_type_corrected_updated)

# Join with metadata to get cancer types
ordering_df <- left_join(tmp, metadata_df, by = "sample_id")

# Arrange by cancer type and then by descending z-scores within each cancer type group
ordering_df <- ordering_df %>%
  dplyr::arrange(cancer_type_corrected_updated, desc(all_zscores))

# Now extract the sample order
sample_order <- ordering_df$sample_id

# Apply this order to your metadata_df and any other data frames that will be used for the heatmap
metadata_df <- metadata_df %>%
  dplyr::filter(sample_id %in% sample_order) %>%
  dplyr::mutate(sample_id = factor(sample_id, levels = sample_order)) %>%
  dplyr::arrange(sample_id)

# Apply the same ordering to your plot_df and any annotation data frames
plot_df <- plot_df[sample_order, ]  

# Set cancer types plotting
cancer_types <- c("healthy", "prostate_cancer", "head_and_neck_cancer", "brain_cancer", 
                  "AML", "eye_cancer", "lfs_survivor", "lfs_previvor", 
                  "lfs_positive", "lung_cancer")


## Get title case version 
metadata_df$cancer_type_title_case <- tools::toTitleCase(metadata_df$cancer_type_corrected_updated)
metadata_df$cancer_type_title_case  <- gsub("Lfs", "LFS", metadata_df$cancer_type_title_case)

### Set parameters for heatmap

# getting IDs
data_id <- as.matrix(metadata_df$sequencing_id)
data_cancer_type <- as.matrix(metadata_df$cancer_type_title_case)
data_age <- as.matrix(metadata_df$age)
data_cohort <- as.matrix(metadata_df$project_id)
#data_sex <- as.matrix(metadata_df$sex)
#data_ichorCNA <- as.matrix(metadata_df$sex_ichorCNA)


row.names(data_id) <- sample_order
#row.names(data_sex) <- sample_order
#row.names(data_ichorCNA) <- sample_order
row.names(data_cancer_type) <- sample_order
row.names(data_age) <- sample_order
row.names(data_cohort) <- sample_order



#data_sex <- factor(data_sex, levels = c("Male", "Female", ""),
 #                  labels = c("Male", "Female", "Not available"))
data_cancer_type <- factor(data_cancer_type, levels = c("Healthy", "Prostate_cancer", "Head_and_neck_cancer", "Brain_cancer", 
                                                        "AML", "Eye_cancer", "LFS_survivor", "LFS_previvor", 
                                                        "LFS_positive", "Lung_cancer"),
                           labels = c("Healthy", "Prostate_cancer", "Head_and_neck_cancer", "Brain_cancer", 
                                      "AML", "Eye_cancer", "LFS_survivor", "LFS_previvor", 
                                      "LFS_positive", "Lung_cancer"))
data_cohort <- factor(data_cohort, levels = c("TCGE-CFMe-HBC", "TCGE-CFMe-MCA", "TCGE-CFMe-BCA", "TCGE-CFMe-PRAD", "TCGE-CFMe-HNSC", "TCGE-CFMe-SCLC", "TCGE-CFMe-UM", "TCGE-CFMe-AML", "TCGE-CFMe-LFS"),
                      labels = c(levels = c("TCGE-CFMe-HBC", "TCGE-CFMe-MCA", "TCGE-CFMe-BCA", "TCGE-CFMe-PRAD", "TCGE-CFMe-HNSC", "TCGE-CFMe-SCLC", "TCGE-CFMe-UM", "TCGE-CFMe-AML", "TCGE-CFMe-LFS")))

## Set chromosome order
armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
               "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
               "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
               "19p", "19q","20p","20q","21q","22q")
data_chr$arm <- factor(data_chr$arm, levels=armlevels)
data_chr$arm <- factor(data_chr$arm, levels=armlevels,
                       labels = c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
                                  "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
                                  "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
                                  "19p", "19q","20p","20q","21q","22q"))


### Set annotation and heatmap colours
col_heat <- colorRamp2(c(-10, 0, 10), 
                       c("#1f78b4", "white", "#e31a1c"))
col_fun <- colorRamp2(c(-2, 0, 2), 
                      c("#1f78b4", "white", "#e31a1c"))
#col_sex <- c("Male" = "#1b9e77",
#             "Female" = "#d95f02",
#             "Not available" = "lightgrey")


# Assign specific colors to cancer types directly
cancer_type_col <- c('Healthy' = '#33a02c',
                     'Brain_cancer' = '#1f78b4',
                     'Lung_cancer' = '#b2df8a',
                     'Prostate_cancer' = '#a6cee3',
                     'AML' = '#fb9a99',
                     'Pancreatic_cancer' = '#e31a1c',
                     'Eye_cancer' = '#fdbf6f',
                     'Head_and_neck_cancer' = '#ff7f00',
                     'Breast_cancer' = '#cab2d6',
                     'Colorectal_cancer' = '#6a3d9a',
                     'Bladder_cancer' = '#fb6a4a',
                     'Renal_cancer' = '#b15928',
                     'LFS_survivor' = '#bdbdbd',
                     'LFS_previvor' = '#969696',
                     'LFS_positive' = '#737373')

col_cohort <- c('TCGE-CFMe-HBC' = '#4daf4a',
                 'TCGE-CFMe-MCA' = '#7fcdbb',
                 'TCGE-CFMe-BCA' = '#984ea3',
                 'TCGE-CFMe-PRAD' = '#377eb8',
                 'TCGE-CFMe-HNSC' = '#ff7f00',
                 'TCGE-CFMe-SCLC' = '#807dba',
                 'TCGE-CFMe-UM' = '#a65628',
                 'TCGE-CFMe-AML' = '#f781bf',
                 'TCGE-CFMe-LFS' = '#999999')

## Set column/row order and column/row splits
#sample_order <- rownames(plot_df)
chr_order <- rownames(data_chr)
column_split <- data_chr$arm
#row_split <- factor(data_id, levels = c("Normal", unique(data_patients$STUDY_ID))) #split by subject
#row_split <- factor(data_time, levels = c("Normal", "Baseline","3 months")) #split by time point

## Debugging 
# Ensure data_cancer_type and data_cohort are correctly factorized with corresponding labels
data_cancer_type <- factor(data_cancer_type, levels = names(cancer_type_col), labels = names(cancer_type_col))
data_cohort <- factor(data_cohort, levels = names(col_cohort), labels = names(col_cohort))


missing_colors <- setdiff(levels(data_cancer_type), names(cancer_type_col))
if (length(missing_colors) > 0) {
  print("Missing color mappings for the following cancer types:")
  print(missing_colors)
} else {
  print("All cancer types have corresponding color mappings.")
}

na_cancer_types <- sum(is.na(data_cancer_type))
if (na_cancer_types > 0) {
  print(paste("There are", na_cancer_types, "NA values in data_cancer_type"))
} else {
  print("There are no NA values in data_cancer_type")
}


### Generate annotations
left_annotation <- rowAnnotation(# "Sex" = data_sex,
                                 "Type" = data_cancer_type,
                                 "Cohort" = data_cohort, # added cohort info
                                 show_annotation_name = TRUE,
                                 #border = TRUE,
                                 col = list("Type" = cancer_type_col, "Cohort" = col_cohort), #"Sex" = col_sex, <- can add if including sex
                                 annotation_name_gp = gpar(fontsize = 12),
                                 annotation_name_side = "top",
                                 annotation_name_rot = 90,
                                 annotation_legend_param = list(title_gp = gpar(fontsize = 15), labels_gp = gpar(fontsize = 14)), # Font size for the legend
                                 #annotation_legend_param = list(border = TRUE),
                                 simple_anno_size = unit(0.35, "cm"),
                                 gap = unit(c(1, 1), "mm")) #  gap = unit(c(1, 1, 1), "mm")) <- add 3 if including sex


## Edit zscore order 
# Order the zscores_mat according to sample_order
zscores_mat_delfi <- zscores_mat_delfi[sample_order, ]

zscores_ordered <- zscores_mat_delfi[match(sample_order, rownames(zscores_mat_delfi)), ]

# Make sure that the row names match the sample order
rownames(zscores_ordered) <- sample_order

right_annotation <- rowAnnotation("Genome-Wide\nZ-score" = anno_lines(zscores_ordered,
                                                                      add_points = TRUE,
                                                                      pch = c(16, NA),
                                                                      ylim = c(-1,14),
                                                                      gp = gpar(col = c("black", "red"),
                                                                                lty = c("solid", "dashed")),
                                                                      axis_param = list(side = "bottom",
                                                                                        labels_rot = 0,
                                                                                        gp = gpar(fontsize = 10))),
                                  #height = unit(5, "in"),
                                  width = unit(2,"cm"),
                                  show_annotation_name = TRUE,
                                  border = TRUE,
                                  annotation_name_gp = gpar(fontsize = 12),
                                  annotation_name_side = "bottom",
                                  annotation_name_rot = 0,
                                  annotation_legend_param = list(title_gp = gpar(fontsize = 15), labels_gp = gpar(fontsize = 14)))  # Font size for the legend
                                  #annotation_legend_param = list(border = TRUE))

top_annotation <- HeatmapAnnotation("Healthy\nMedian" = anno_lines(Normal_median,
                                                                   ylim = c(-0.08, 0.08),
                                                                   axis_param = list(at = c(-0.05, 0, 0.05),
                                                                                     side = "right",
                                                                                     labels_rot = 0,
                                                                                     gp = gpar(fontsize = 10)),
                                                                   border = FALSE),
                                    height = unit(1.2, "cm"),
                                    show_annotation_name = FALSE,
                                    annotation_name_rot = 0,
                                    annotation_name_gp = gpar(fontsize = 13))

## Generate oncoprint
row_split <- data_cancer_type

Heatmap_delfi <- ComplexHeatmap::Heatmap(plot_df,
                   col = col_heat,
                   show_row_names = FALSE,
                   show_column_names = FALSE,
                   #heatmap_legend_param = heatmap_legend_param,
                   right_annotation = right_annotation,
                   top_annotation = top_annotation,
                   left_annotation = left_annotation,
                   column_order = chr_order,
                   row_order = sample_order,
                   row_labels = NULL,
                   row_split = row_split, # For split 
                   column_split = column_split,
                   column_title_gp = gpar(fontsize = 12),
                   row_title_gp = gpar(fontsize = 14),
                   cluster_row_slices = FALSE,
                   row_title_rot = 0,
                   border = FALSE,
                   heatmap_legend_param = list(title = "Ratio zscore", 
                                               title_gp = gpar(fontsize = 15),  # Font size for the legend title
                                               labels_gp = gpar(fontsize = 14)))  # Font size for the legend labels
    #               heatmap_legend_param = list(title = "Ratio zscore"))  # Set the title of the legend to "Proportion"


ht_opt$COLUMN_ANNO_PADDING = unit(5, "mm")
#draw(Heatmap, heatmap_legend_side = "right", merge_legend = TRUE)

png(filename = file.path(outdir, "Heatmap_Z-score_DELFI_updated_in_cfMedip_5Mb_updated_legend_updated_colors_May_2024_new_colors_ordered_split_updated_3.png"), 
    width = 30, height = 25, units = "in", res = 350)
draw(Heatmap, heatmap_legend_side = "right", merge_legend = TRUE)
dev.off()

png(filename = file.path(outdir, "Heatmap_Z-score_DELFI_updated_in_cfMedip_5Mb_updated_legend_updated_colors_May_2024_new_colors_ordered_split_updated_3_narrow2.png"), 
    width = 24, height = 26, units = "in", res = 350)
draw(Heatmap_delfi, heatmap_legend_side = "right", merge_legend = TRUE)
dev.off()

## Now get the coefficient of variance 

# Prepare data
data_df <- as.data.frame(plot_df)
data_df$Sample_ID <- rownames(data_df)
data_df <- left_join(data_df, metadata_df %>% dplyr::select(sample_id, cancer_type_corrected_updated, cancer_type_title_case), by = c("Sample_ID" = "sample_id"))
measurement_columns <- colnames(data_df)[1:(ncol(data_df)-3)] # Set measurement columns

# Calculate the mean and standard deviation for each measurement column within each cancer type
group_stats <- data_df %>%
  group_by(cancer_type_corrected_updated, cancer_type_title_case) %>%
  summarise(across(all_of(measurement_columns),
                   list(mean = ~ mean(.x, na.rm = TRUE), sd = ~ sd(.x, na.rm = TRUE)),
                   .names = "{.col}_{.fn}"),
            .groups = 'drop')

# Calculate the Coefficient of Variation (CV) for each group and each column
# CV is the standard deviation divided by the mean
group_stats_cv <- group_stats %>%
  rowwise() %>%
  dplyr::mutate(across(ends_with("_mean"), 
                       ~ get(str_replace(cur_column(), "_mean", "_sd")) / .x,
                       .names = "{.col}_cv"))

# Since the CV calculation produces NaN when the mean is 0, you might want to handle these cases
# Replace NaN CV values with 0 (assuming a mean of 0 implies no variation and hence a CV of 0)
group_stats_cv <- dplyr::mutate(group_stats_cv, across(ends_with("_cv"), ~ifelse(is.nan(.x), 0, .x)))

# Now, group_stats_cv contains the CV for each measurement column within each cancer type
# You might want to drop the mean and sd columns to keep only CV values
cv_columns <- grep("_cv$", colnames(group_stats_cv), value = TRUE)
cv_data <- dplyr::select(group_stats_cv, cancer_type_corrected_updated, cancer_type_title_case, all_of(cv_columns))

# View the CV data
head(cv_data)

## Correct naming 
cv_data$cancer_type_title_case  <- gsub("Head_and_neck_cancer", "Head_and\nneck_cancer", cv_data$cancer_type_title_case)

## Correct order 
cv_data$cancer_type_title_case <- factor(cv_data$cancer_type_title_case, c("Healthy", "AML", "Lung_cancer", 
                                                                           "Eye_cancer", "Prostate_cancer", "Head_and\nneck_cancer",
                                                                           "Brain_cancer", "LFS_positive", "LFS_survivor", 
                                                                           "LFS_previvor"))

## Pivot 
long_df_cv <- cv_data %>%
  pivot_longer(
    cols = -c(cancer_type_corrected_updated, cancer_type_title_case), # Exclude the non-numeric columns
    names_to = c("bin"), # Create 'bin' from part of the column names
    names_pattern = "^(\\d+)_mean_cv$", # Correct pattern to capture the bin number
    values_to = "cv" # Define where the values should go
  ) %>%
  # Extract only the number part of the column name to use as the bin
  # The number part is in the column named 'bin' after pivot_longer
  dplyr::mutate(bin = as.numeric(bin))

## Export 
write.table(long_df_cv, file = "Coefficient of variance long df, June 2024.txt", sep = "\t", quote = F, row.names = F)

## Add color 
col_cancer_type <- c('Healthy' = '#33a02c',
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
                     'LFS_positive' = '#737373')

# Now make the boxplot
# Cap the CV values at 100 and -100
long_df_cv$CV_capped <- ifelse(long_df_cv$cv > 100, 100, ifelse(long_df_cv$cv < -100, -100, long_df_cv$cv))

theme <- theme(plot.title = element_text(hjust = 0.5, size = 13), 
               axis.line = element_line(colour = "black"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               legend.position = "bottom",
               panel.border = element_rect(colour = "black", fill=NA, size=1), # Add a border
               legend.key = element_rect(fill = "white"),
               legend.title = element_text(size = 12),
               legend.text = element_text(size = 12),
               strip.background = element_blank(),
               strip.text = element_text(size = 13),
               axis.text = element_text(size = 13),
               axis.title = element_text(size = 13))

fig_delfi <- ggplot(long_df_cv, aes(x = cancer_type_title_case, y = CV_capped, color = cancer_type_title_case)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(pch = 16, alpha = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", size = 0.5) +
  xlab("Samples") + 
  ylab("Coefficient of Variation") +
  labs(color = "Cancer type") +
  scale_color_manual(values = col_cancer_type) + # Define color scale for cancer types
  ggtitle("Coefficient of Variation of sample DELFI ratios") + 
  theme +
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(
    limits = c(-100, 100), 
    breaks = c(-100, -50, 0, 50, 100),
    labels = c("<-100", "-50", "0", "50", ">100")
  )

ggsave(file.path(outdir, paste("Delfi fig, May 2024 with border.pdf", sep = "")), fig_delfi, width = 10, height = 3, dpi = 500, units = "in")





### Export the dataframes
saveRDS(plot_df, file = "DELFI_Ratios_updated_June2024.rds")
write.table(plot_df, file = "DELFI_Ratios_updated_June2024.txt", sep = "\t", quote = F, row.names = T)


## Export the PCAs on the dataframe 

# Perform PCA
pca_result <- prcomp(plot_df, center = TRUE, scale. = TRUE)

# Determine the number of components to retain
# Visualize the variance explained by each principal component
fviz_eig(pca_result, addlabels = TRUE, ylim = c(0, 15))

# Extract the proportion of variance explained by each component
explained_variance <- summary(pca_result)$importance[2,]
cum_variance <- cumsum(explained_variance)

# Find the number of components that explain at least 90% of the variance
num_components <- which(cum_variance >= 0.90)[1]

# Extract the principal components
pca_df <- as.data.frame(pca_result$x[, 1:num_components])
rownames(pca_df) <- rownames(plot_df)

# Create the variance explained plot with customized theme
variance_plot <- fviz_eig(pca_result, addlabels = TRUE, ylim = c(0, 15)) +
  theme_classic() +
  ggtitle("Scree plot for fragment proportions") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

# Save the plot to a file
ggsave(filename = file.path(outdir, "PCA_Variance_Explained_Delfi_ratios.png"), plot = variance_plot, width = 10, height = 6, units = "in", dpi = 300)

# Save the PCA results to a file 
write.csv(pca_df, file.path(outdir, "PCA_results_Delfi_ratios.csv"), row.names = TRUE)
saveRDS(pca_df, file = "DELFI_Ratios_updated_June2024_PCAs_90_percent_variance.rds")
saveRDS(pca_result, file = "PCA_result_Delfi_June2024.rds")

## Export the metadata df 
write.table(metadata_df, file = "Metadata df June 2024 annotated.txt", sep = "\t", quote = F, row.names = F)

## Export the table 
write.table(plot_df, file = "Fragment_proportions_matrix_dataframe.txt", sep = "\t", row.names = TRUE, quote = F)

