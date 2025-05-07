# Process and plot insert sizes

# author: Dory Abelman
# adapted from scripts by Derek Wong
# Original date: Aug 2023
# Updated for new samples Oct-Nov 2023
# Updated May 2024 for final set of adjusted prostate cancer samples

library(dplyr)
library(matrixStats)
library(ComplexHeatmap)
library(circlize)
library(cowplot)
library(ggpubr)
library(gridExtra)
library(RColorBrewer)
library(circlize)
library(reshape2)
library(factoextra)  # For visualizing PCA results



### Set paths
path <- "Insert_size/Data/"
outdir <- "Final Data Dec 2024/Insert_size/"

insert_path <- list.files(path = path, pattern = "insert_size.txt", full.names = TRUE)

names_path <- list.files(path = path, pattern = "insert_size.txt", full.names = FALSE)

# Edit the names
names_path <- sub(".insert_size.txt", "", names_path) # remove ending

# Make samples list from names
samples <- as.data.frame(names_path)
samples <- samples %>% dplyr::rename("Sample_id" = "names_path")

datalist <- lapply(insert_path, function(x){read.delim(file = x, skip = 13, colClasses = c(rep("integer", 2), rep("NULL", 2)), 
                                                       header = TRUE, col.names = c("length", "freq", "X", "Y"))})
insert_path <- Reduce(function(x,y) {merge(x,y, by = "length", all = TRUE, sort = TRUE)}, datalist)

### Format fragment data and save - need to redo this for longer fragment frequencies
row.names(insert_path) <- insert_path$length
insert_path <- insert_path[,-1]
colnames(insert_path) <- names_path
lengths <- c(10:600) # getting the single, di and tri peaks
insert_path <- insert_path[row.names(insert_path) %in% lengths, ]

## Now change from counts to frequencies
insert_path_long <- as.data.frame(t(insert_path))

# Make index a column value (sample type)
insert_path_long <- cbind(Sample_ID = rownames(insert_path_long), insert_path_long) 
rownames(insert_path_long) <- 1:nrow(insert_path_long)

# Then change from long to tall
insert_path_tall <- melt(insert_path_long, id.vars="Sample_ID") # From the most recent batch

# Make numeric
insert_path_tall$variable <- as.numeric(as.character(insert_path_tall$variable)) # change if needed
insert_path_tall$value <- as.numeric(as.character(insert_path_tall$value)) # change if needed

# Set NA to 0
insert_path_tall[is.na(insert_path_tall)] <- 0

# Next get the proportion instaed of the value
insert_path_tall <- insert_path_tall %>%
  dplyr::group_by(Sample_ID) %>%
  dplyr::mutate(proportion = value / sum(value))

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

## Join to insert size df
insert_size_df <- left_join(insert_path_tall, metadata_df, by = c("Sample_ID" = "sample_id"))

insert_size_df_validation <- insert_size_df_validation %>%
  filter(validation == "1", cancer_type != "Liver Cancer")

#insert_size_df <- insert_size_df %>% filter(validation == "0")
healthy_insert <- insert_size_df_validation %>% dplyr::filter(cancer_type_corrected_updated == "healthy") %>% 
  dplyr::select(variable, proportion, Sample_ID)

# Now form a matrix with each column as sample ID and each row as length for the normal samples
healthy_insert_matrix <- pivot_wider(healthy_insert, names_from = "variable", values_from = "proportion")

row_names_to_keep <- healthy_insert_matrix$Sample_ID

# Remove the 'sample' column
healthy_insert_matrix <- subset(healthy_insert_matrix, select = -1)

# Reassign row names
rownames(healthy_insert_matrix) <- row_names_to_keep

## Now reformat to flipped matrix
healthy_insert_matrix2 <- pivot_wider(healthy_insert, names_from = "Sample_ID", values_from = "proportion")
row_names_to_keep <- healthy_insert_matrix2$variable
healthy_insert_matrix2 <- subset(healthy_insert_matrix2, select = -1)
rownames(healthy_insert_matrix2) <- row_names_to_keep

### Make healthy median
Normal_median <- rowMedians(as.matrix(healthy_insert_matrix2), na.rm = T)
Normal_sd <- rowSds(as.matrix(healthy_insert_matrix2), na.rm = T)

## Now get cancer median across all cancers
cancer_insert <- insert_size_df_validation %>% dplyr::filter(cancer_type_corrected_updated != "healthy") %>% 
  dplyr::filter(!is.na(cancer_type_corrected_updated)) %>% 
  dplyr::select(variable, proportion, Sample_ID)

# Now form a matrix with each column as sample ID and each row as length for the normal samples
cancer_insert_matrix <- pivot_wider(cancer_insert, names_from = "variable", values_from = "proportion")

row_names_to_keep <- cancer_insert_matrix$Sample_ID

# Remove the 'sample' column
cancer_insert_matrix <- subset(cancer_insert_matrix, select = -1)

# Reassign row names
rownames(cancer_insert_matrix) <- row_names_to_keep

## Now reformat to flipped matrix
cancer_insert_matrix2 <- pivot_wider(cancer_insert, names_from = "Sample_ID", values_from = "proportion")
row_names_to_keep <- cancer_insert_matrix2$variable
cancer_insert_matrix2 <- subset(cancer_insert_matrix2, select = -1)
rownames(cancer_insert_matrix2) <- row_names_to_keep


### Make cancer median
cancer_median <- rowMedians(as.matrix(cancer_insert_matrix2))
cancer_sd <- rowSds(as.matrix(cancer_insert_matrix2))

## Get distance from normal median 
Cancer_insert_ratio <- (cancer_insert_matrix2 - Normal_median)/Normal_sd
Normal_insert_ratio <- (healthy_insert_matrix2 - Normal_median)/Normal_sd

### Calculate Normal Z-scores (Normal relative to Normal)
Normal_sum <- colSums(abs(Normal_insert_ratio), na.rm = TRUE)
Normal_sum_median <- median(Normal_sum, na.rm = TRUE)
Normal_MAD <- mad(Normal_sum, na.rm = TRUE) #Median Absolute Deviation

sum(is.na(Normal_sd))

Normal_zscores <- abs((Normal_sum - Normal_sum_median))/Normal_MAD
Normal_mean <- mean(Normal_zscores)
z_limit <- quantile(Normal_zscores, 0.9)

### Calculate Cancer Z-scores (Cancer relative to Normal)
Cancer_sum <- colSums(abs(Cancer_insert_ratio), na.rm = TRUE)
Cancer_zscore <- (Cancer_sum - Normal_sum_median)/Normal_MAD

### Make Z-score matrix
all_zscores <- c(Normal_zscores,Cancer_zscore)
zscore_df <- as.data.frame(all_zscores)
zscore_df$limit <- z_limit
zscores_mat <- as.matrix(zscore_df)
zscores_mat_insertsize <- zscores_mat # save

### Set Normal median 
lower <- Normal_median - Normal_sd
upper <- Normal_median + Normal_sd
Normal_median <- as.matrix(cbind(Normal_median, lower, upper))

## Save the zscores to a seperate dataframe 
Insert_df_zscores_validation <- zscore_df 
Insert_df_zscores_validation <- rownames_to_column(Insert_df_zscores_validation, "Sample")

# Join with metadata df
Insert_df_zscores_validation <- left_join(Insert_df_zscores_validation, metadata_df, by = c("Sample" = "sample_id"))

### Save objects
saveRDS(Normal_ratio, file=file.path(outdir, paste0("Normal_ratio_Insert_size_validation.rds")))
saveRDS(Cancer_ratio, file=file.path(outdir, paste0("Cancer_ratio_Insert_size_validation.rds")))
saveRDS(Normal_median, file=file.path(outdir, paste0("Normal_median_Insert_size_validation.rds")))
saveRDS(zscores_mat_insertsize, file = file.path(outdir, paste0("Zscore_mat_Insert_size_validation.rds")))
saveRDS(Insert_df_zscores_validation, file = "Zscores df for insert size_validation.rds")


# Check for missing samples
# Identify values in metadata_df$sample_id that are not in Insert_df_zscores_validation$Sample
unique_sample_ids <- setdiff(metadata_df$sample_id, Insert_df_zscores_validation$Sample)

# Print the unique sample IDs
print(unique_sample_ids)

## Join matrices together 
insert_matrix_df_validation <- bind_rows(healthy_insert_matrix, cancer_insert_matrix)


## Make a heatmap of the insert size matrix 
insert_matrix_validation <-as.matrix(insert_matrix_df_validation)

# getting IDs
sample_order <- rownames(insert_matrix_validation)
metadata_df <- metadata_df[match(sample_order, metadata_df$sample_id),]

data_id <- as.matrix(metadata_df$sequencing_id)
#data_sex <- as.matrix(metadata_df$sex)
#data_ichorCNA <- as.matrix(metadata_df$Tumor_fraction)
data_cancer_type <- as.matrix(metadata_df$cancer_type_title_case)
data_age <- as.matrix(metadata_df$age)
data_cohort <- as.matrix(metadata_df$project_id)


row.names(data_id) <- sample_order
#row.names(data_sex) <- sample_order
#row.names(data_ichorCNA) <- sample_order
row.names(data_cancer_type) <- sample_order
row.names(data_age) <- sample_order
row.names(data_cohort) <- sample_order


#data_sex <- factor(data_sex, levels = c("Male", "Female", ""),
#                   labels = c("Male", "Female", "Not available"))
# Updated factor levels and labels for data_cancer_type
data_cancer_type <- factor(data_cancer_type, 
                           levels = c("Healthy", 
                                      "Prostate_cancer", 
                                      "Head_and_neck_cancer", 
                                      "Brain_cancer", 
                                      "AML", 
                                      "Eye_cancer", 
                                      "LFS_survivor", 
                                      "LFS_previvor", 
                                      "LFS_positive", 
                                      "Lung_cancer", 
                                      "Liver_cancer", 
                                      "Melanoma", 
                                      "Mixed_cancer", 
                                      "Breast_cancer", 
                                      "Ovarian_cancer"),
                           labels = c("Healthy", 
                                      "Prostate_cancer", 
                                      "Head_and_neck_cancer", 
                                      "Brain_cancer", 
                                      "AML", 
                                      "Eye_cancer", 
                                      "LFS_survivor", 
                                      "LFS_previvor", 
                                      "LFS_positive", 
                                      "Lung_cancer", 
                                      "Liver_cancer", 
                                      "Melanoma", 
                                      "Mixed_cancer", 
                                      "Breast_cancer", 
                                      "Ovarian_cancer"))

# Updated factor levels and labels for data_cohort
data_cohort <- factor(data_cohort, 
                      levels = c("TCGE-CFMe-HBC", 
                                 "TCGE-CFMe-MCA", 
                                 "TCGE-CFMe-BCA", 
                                 "TCGE-CFMe-PRAD", 
                                 "TCGE-CFMe-HNSC", 
                                 "TCGE-CFMe-SCLC", 
                                 "TCGE-CFMe-UM", 
                                 "TCGE-CFMe-AML", 
                                 "TCGE-CFMe-LFS", 
                                 "TCGE-CFMe-HCC", 
                                 "TCGE-CFMe-INSPIRE"),
                      labels = c("TCGE-CFMe-HBC", 
                                 "TCGE-CFMe-MCA", 
                                 "TCGE-CFMe-BCA", 
                                 "TCGE-CFMe-PRAD", 
                                 "TCGE-CFMe-HNSC", 
                                 "TCGE-CFMe-SCLC", 
                                 "TCGE-CFMe-UM", 
                                 "TCGE-CFMe-AML", 
                                 "TCGE-CFMe-LFS", 
                                 "TCGE-CFMe-HCC", 
                                 "TCGE-CFMe-INSPIRE"))


### Set annotation and heatmap colours
col_heat <- colorRamp2(c(0, 0.01, 0.025), 
                       c("#1f78b4", "white", "#e31a1c"))
col_fun <- colorRamp2(c(-2, 0, 2), 
                      c("#1f78b4", "white", "#e31a1c"))
#col_sex <- c("Male" = "#1b9e77",
#             "Female" = "#d95f02",
#             "Not available" = "lightgrey")


# Match color scheme
# Match color scheme for cancer types (adding Liver_cancer, Melanoma, Mixed_cancer, Ovarian_cancer)
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
                     'LFS_positive' = '#737373',
                     'Liver_cancer' = '#ffeda0',
                     'Melanoma' = '#41b6c4',
                     'Mixed_cancer' = '#e31a1c', # Same as pancreatic (mixed nature)
                     'Ovarian_cancer' = '#ffffb3')


# Match color scheme for cohorts (adding TCGE-CFMe-HCC and TCGE-CFMe-INSPIRE)
col_cohort <- c('TCGE-CFMe-HBC' = '#4daf4a',
                'TCGE-CFMe-MCA' = '#7fcdbb',
                'TCGE-CFMe-BCA' = '#984ea3',
                'TCGE-CFMe-PRAD' = '#377eb8',
                'TCGE-CFMe-HNSC' = '#ff7f00',
                'TCGE-CFMe-SCLC' = '#807dba',
                'TCGE-CFMe-UM' = '#a65628',
                'TCGE-CFMe-AML' = '#f781bf',
                'TCGE-CFMe-LFS' = '#999999',
                'TCGE-CFMe-HCC' = '#fc9272',  # New for Liver_cancer cohort
                'TCGE-CFMe-INSPIRE' = '#bc80bd')  # New for INSPIRE cohort


## Set column/row order and column/row splits
sample_order <- rownames(insert_matrix_validation)
column_order <- colnames(insert_matrix_validation)

## Add some labels for nucleotide peaks

# Initialize a vector to hold the labels
nucleotide_regions <- rep(NA, ncol(insert_matrix_validation))

# Convert column names to numeric
column_nums <- as.numeric(colnames(insert_matrix_validation))

# Iterate over the column names and assign labels
for (i in seq_along(column_nums)) {
  if (column_nums[i] < 152) {
    nucleotide_regions[i] <- "<151bp"
  } else if (column_nums[i] > 151 & column_nums[i] < 183) {
    nucleotide_regions[i] <- "167+-15bp"
  } else if (column_nums[i] > 182 & column_nums[i] < 319) {
    nucleotide_regions[i] <- "183-318bp"
  } else if (column_nums[i] > 318 & column_nums[i] < 350) {
    nucleotide_regions[i] <- "334+-15bp"
  } else if (column_nums[i] > 349 & column_nums[i] < 486) {
    nucleotide_regions[i] <- "350-485bp"
  } else if (column_nums[i] > 485 & column_nums[i] < 517) {
    nucleotide_regions[i] <- "501+-15bp"
  } else {
    nucleotide_regions[i] <- ">517bp"
  }
}

# Convert to a factor with specified levels
#levels_order <- c("<167bp", "167bp", "168-333bp", "334bp", "335-500bp", "501bp", ">502bp")
levels_order <- c("<151bp", "167+-15bp", "183-318bp", "334+-15bp", "350-485bp", "501+-15bp", ">517bp")
nucleotide_regions <- factor(nucleotide_regions, levels = levels_order)

### Generate annotations
left_annotation <- rowAnnotation(# "Sex" = data_sex,
                                 "Type" = data_cancer_type,
                                 "Cohort" = data_cohort, # added cohort info
                                 show_annotation_name = TRUE,
                                 #border = TRUE,
                                 col = list("Type" = cancer_type_col, "Cohort" = col_cohort),
                                 annotation_name_gp = gpar(fontsize = 12),
                                 annotation_name_side = "top",
                                 annotation_name_rot = 90,
                                 annotation_legend_param = list(title_gp = gpar(fontsize = 15), labels_gp = gpar(fontsize = 14)), # Font size for the legend
                                 #annotation_legend_param = list(border = TRUE),
                                 simple_anno_size = unit(0.35, "cm"),
                                 gap = unit(c(1, 1), "mm"))

right_annotation <- rowAnnotation("Genome-Wide\nZ-score" = anno_lines(zscores_mat_insertsize,
                                                                      add_points = TRUE,
                                                                      pch = c(16, NA),
                                                                      ylim = c(-1,40), 
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
                             #     annotation_legend_param = list(border = TRUE))

top_annotation <- HeatmapAnnotation("Healthy\nMedian" = anno_lines(Normal_median,
                                                                   ylim = c(-0.03, 0.03),
                                                                   axis_param = list(at = c(-0.03, 0, 0.03),
                                                                                     side = "right",
                                                                                     labels_rot = 0,
                                                                                     gp = gpar(fontsize = 10)),
                                                                   border = FALSE,
                                    gp = gpar(lty = c(rep("solid", length(Normal_median) - 1), "dotted"))),  # Set last line as dotted
                                    height = unit(1.2, "cm"),
                                    show_annotation_name = FALSE,
                                    annotation_name_rot = 0,
                                    annotation_name_gp = gpar(fontsize = 13))

# Ensure order of insert_matrix_validation is correct 
insert_matrix_validation <- insert_matrix_validation[, column_order]

## Generate heatmap
Heatmap_insert <- Heatmap(insert_matrix_validation,
                   col = col_heat,
                   show_row_names = FALSE,
                   show_column_names = FALSE,
                   #heatmap_legend_param = heatmap_legend_param,
                   right_annotation = right_annotation,
                   top_annotation = top_annotation,
                   left_annotation = left_annotation,
                   column_order = column_order,
                   row_order = sample_order,
                   row_labels = NULL,
                   column_split = nucleotide_regions,
                   column_title_gp = gpar(fontsize = 10),
                   row_title_gp = gpar(fontsize = 12),
                   #row_split = row_split,
                   cluster_row_slices = FALSE,
                   row_title_rot = 0,
                   border = FALSE,
                   heatmap_legend_param = list(title = "Proportion", 
                                               title_gp = gpar(fontsize = 15),  # Font size for the legend title
                                               labels_gp = gpar(fontsize = 14)))  # Font size for the legend labels
 #     heatmap_legend_param = list(title = "Proportion"))  # Set the title of the legend to "Proportion"


ht_opt$COLUMN_ANNO_PADDING = unit(5, "mm")
# draw(Heatmap, heatmap_legend_side = "right", merge_legend = TRUE)

png(filename = "Heatmap of insert size by cancer type, updated legends, updated colors, with breaks, March 2025, validation.png", width = 30, height = 25, units = "in", res = 350)
draw(Heatmap_insert, heatmap_legend_side = "right", merge_legend = TRUE)
dev.off()


### Now create seperate matrix, ordered by zscore and with the cohort label information added - Jan 2024

# Order the zscores_mat_insertsize according to sample_order
# Suppose your Z-scores are stored in zscores_mat_insertsize
# Make a "capped" version, replacing Inf (and anything above 45) with 45
zscores_mat_insertsize_capped <- pmin(zscores_mat_insertsize, 45)
zscores_mat_insertsize_capped[is.infinite(zscores_mat_insertsize_capped)] <- 45


zscores_ordered <- zscores_mat_insertsize_capped[match(sample_order, rownames(zscores_mat_insertsize_capped)), ]

# Make sure that the row names match the sample order
rownames(zscores_ordered) <- sample_order

# Then supply zscores_mat_insertsize_capped to the anno_lines call:
right_annotation <- rowAnnotation(
  "Genome-Wide\nZ-score" = anno_lines(
    zscores_mat_insertsize_capped,
    add_points = TRUE,
    pch = c(16, NA),
    ylim = c(-5, 45),
    gp = gpar(col = c("black", "red"), lty = c("solid", "dashed")),
    axis_param = list(
      side = "bottom",
      # Pick whichever breaks you want to show
      at = c(-5, 0, 20, 40), 
      # Label 40 as ">40"
      labels = c("", "0", "20", ">40"), 
      labels_rot = 0,
      gp = gpar(fontsize = 10)
    )
  ),
  width = unit(2, "cm"),
  show_annotation_name = TRUE,
  border = TRUE,
  annotation_name_gp = gpar(fontsize = 12),
  annotation_name_side = "bottom",
  annotation_name_rot = 0,
  annotation_legend_param = list(title_gp = gpar(fontsize = 15), labels_gp = gpar(fontsize = 14))
)

## Generate heatmap
row_split <- data_cancer_type

Heatmap_insert2 <- Heatmap(insert_matrix_validation,
                   col = col_heat,
                   show_row_names = FALSE,
                   show_column_names = FALSE,
                   #heatmap_legend_param = heatmap_legend_param,
                   right_annotation = right_annotation,
                   top_annotation = top_annotation,
            #       left_annotation = left_annotation,
                   column_order = column_order,
                   row_order = sample_order,
                   row_labels = NULL,
                   column_split = nucleotide_regions,
                   column_title_gp = gpar(fontsize = 10),
                   row_title_gp = gpar(fontsize = 12),
                   row_split = row_split, # for split 
                   cluster_row_slices = FALSE,
                   row_title_rot = 0,
                   border = FALSE,
                   heatmap_legend_param = list(title = "Insert size proportion", 
                                               title_gp = gpar(fontsize = 15),  # Font size for the legend title
                                               labels_gp = gpar(fontsize = 14)))  # Font size for the legend labels


ht_opt$COLUMN_ANNO_PADDING = unit(5, "mm")
# draw(Heatmap, heatmap_legend_side = "right", merge_legend = TRUE)

png(filename = "Heatmap of insert size by cancer type, updated legends, updated colors, with breaks, split, March 2025 validation.png", width = 30, height = 25, units = "in", res = 350)
draw(Heatmap_insert2, heatmap_legend_side = "right", merge_legend = TRUE)
dev.off()

## Make heatmap but only up to 350bp. 
# Convert column names to numeric to ensure accurate selection
numeric_colnames <- as.numeric(colnames(insert_matrix_validation))

# Subset the matrix based on the numeric column names
insert_matrix_validation_subset <- insert_matrix_validation[, numeric_colnames >= 10 & numeric_colnames < 350]

# Initialize a vector to hold the labels
nucleotide_regions_150bp <- rep(NA, ncol(insert_matrix_validation_subset))

# Convert column names to numeric
column_nums <- as.numeric(colnames(insert_matrix_validation_subset))

# Iterate over the column names and assign labels
for (i in seq_along(column_nums)) {
  if (column_nums[i] < 152) {
    nucleotide_regions_150bp[i] <- "10-151bp"
  } else if (column_nums[i] > 151 & column_nums[i] < 183) {
    nucleotide_regions_150bp[i] <- "167+-15bp"
  } else if (column_nums[i] > 182 & column_nums[i] < 319) {
    nucleotide_regions_150bp[i] <- "183-318bp"
  } else if (column_nums[i] > 318 & column_nums[i] < 350) {
    nucleotide_regions_150bp[i] <- "334+-15bp"
  } else {
    nucleotide_regions_150bp[i] <- "Other"
  }
}

# Convert to a factor with specified levels
#levels_order <- c("<167bp", "167bp", "168-333bp", "334bp", "335-500bp", "501bp", ">502bp")
levels_order_short <- c("10-151bp", "167+-15bp", "183-318bp", "334+-15bp")
nucleotide_regions_150bp <- factor(nucleotide_regions_150bp, levels = levels_order_short)

# Generate heatmap 

# Subset Normal_median to keep only rows 0 to 311, corresponsing to insert sizes 10-320, or 341 for 350
Normal_median_subset <- Normal_median[0:340, ]

top_annotation <- HeatmapAnnotation("Healthy\nMedian" = anno_lines(Normal_median_subset,
                                                                   ylim = c(-0.01, 0.03),
                                                                   axis_param = list(at = c(0, 0.03),
                                                                                     side = "right",
                                                                                     labels_rot = 0,
                                                                                     gp = gpar(fontsize = 10)),
                                                                   border = FALSE,
                                                                   gp = gpar(lty = c(rep("solid", length(Normal_median) - 1), "dotted"))),  # Set last line as dotted
                                    height = unit(1.2, "cm"),
                                    show_annotation_name = FALSE,
                                    annotation_name_rot = 0,
                                    annotation_name_gp = gpar(fontsize = 13))

## Generate heatmap
column_order_narrow <- colnames(insert_matrix_validation_subset)

# Define a new sequential color scale
col_heat_seq <- colorRamp2(c(0, 0.007, 0.025), 
                           c("#dcedf6", "#92b4e8", "#1f78b4"))

Heatmap_insert_narrow <- Heatmap(insert_matrix_validation_subset,
                   col = col_heat_seq,
                   show_row_names = FALSE,
                   show_column_names = FALSE,
                   #heatmap_legend_param = heatmap_legend_param,
                   right_annotation = right_annotation,
                   top_annotation = top_annotation,
             #      left_annotation = left_annotation,
                   column_order = column_order_narrow,
                   row_order = sample_order,
                   row_labels = NULL,
                   column_split = nucleotide_regions_150bp,
                   column_title_gp = gpar(fontsize = 10),
                   row_title_gp = gpar(fontsize = 12),
                   row_split = row_split,
                   cluster_row_slices = FALSE,
                   row_title_rot = 0,
                   border = FALSE,
                   heatmap_legend_param = list(title = "Proportion", 
                                               title_gp = gpar(fontsize = 15),  # Font size for the legend title
                                               labels_gp = gpar(fontsize = 14)))  # Font size for the legend labels


ht_opt$COLUMN_ANNO_PADDING = unit(5, "mm")
#draw(Heatmap_insert_narrow, heatmap_legend_side = "right", merge_legend = TRUE)

png(filename = "Heatmap of insert size by cancer type, updated legends, updated colors, with breaks, March 2025, narrow, 2.png", width = 18, height = 25, units = "in", res = 350)
draw(Heatmap_insert_narrow, heatmap_legend_side = "right", merge_legend = TRUE)
dev.off()


png(filename = "Heatmap of insert size by cancer type, updated legends, updated colors, with breaks, March 2025, narrow, 3.png", width = 18, height = 23, units = "in", res = 350)
draw(Heatmap_insert_narrow, heatmap_legend_side = "right", merge_legend = TRUE)
dev.off()

png(filename = "Heatmap of insert size by cancer type, updated legends, updated colors, with breaks, March 2025, narrow, 4.png", width = 15, height = 20, units = "in", res = 350)
draw(Heatmap_insert_narrow, heatmap_legend_side = "right", merge_legend = TRUE)
dev.off()

png(filename = "Heatmap of insert size by cancer type, updated legends, updated colors, with breaks, March 2025, narrow, 5.png", width = 13, height = 17, units = "in", res = 500)
draw(Heatmap_insert_narrow, heatmap_legend_side = "right", merge_legend = TRUE)
dev.off()

### Now do boxplot showing the proportion of short reads
### Calculate fragment proportions

short <- rownames_to_column(insert_matrix_df_validation, "Sample")
# Select columns of interest
selected_columns <- c("Sample", paste(20:150))
short_subset <- short[, selected_columns]

# Calculate the column sums and convert to proportions
prop <- rowSums(as.matrix(short_subset[, -1])) / 100
prop <- as.data.frame(prop)
prop$Sample <- short_subset$Sample

## Merge with metadata df
prop_short_df_validation <- left_join(prop, metadata_df, by = c("Sample" = "sample_id"))

## Add column 
prop_short_df_validation$Twentyto150_pct <- prop_short_df_validation$prop*100

## Save
write.table(prop_short_df_validation, "Prop 20 to 150 bp in cfMeDip validation cohort.txt", row.names = FALSE, sep = "\t")
saveRDS(prop_short_df_validation, file = "Prop 20 to 150 bp in cfMeDip_validation.rds")



## Plot 
prop_short_df_validation$cancer_type_title_case <- gsub(
  pattern = "Head_and_neck_cancer", 
  replacement = "Head_and\nneck_cancer", 
  x = prop_short_df_validation$cancer_type_title_case
)

data_stats_cancer_timepoint_prop_insert <-  prop_short_df_validation %>%
  group_by(cancer_type_title_case)%>% 
  dplyr::summarise(Median=median(Twentyto150_pct, na.rm = TRUE),
                   Mean=mean(Twentyto150_pct, na.rm = TRUE),
                   SD=sd(Twentyto150_pct, na.rm = TRUE),
                   N=n()) %>% 
  dplyr::arrange(desc(Median)) # Arrange by decreasing median

# Get the ordered cancer types
ordered_cancer_types <- data_stats_cancer_timepoint_prop_insert$cancer_type_title_case

## Custom order reccomended by Yong
# Remove 'Normal' and 'LFS' from the vector, if they exist
ordered_cancer_types <- ordered_cancer_types[ordered_cancer_types != "Healthy" & ordered_cancer_types != "LFS_survivor" & ordered_cancer_types != "LFS_positive" & ordered_cancer_types != "LFS_previvor"]


# Prepend 'Normal' at the beginning and append 'LFS' at the end
ordered_cancer_types <- c("Healthy", ordered_cancer_types)

## Fixing head and neck naming 
ordered_cancer_types_peak <- ordered_cancer_types

# Create an ordered factor with levels in the order of decreasing median
#prop_short_df_validation$cancer_type_title_case <- gsub("Head_and_neck_cancer", "Head_and\nneck_cancer", prop_short_df_validation$cancer_type_title_case)
prop_short_df_validation$cancer_type_ordered <- factor(prop_short_df_validation$cancer_type_title_case, levels = ordered_cancer_types_peak)

data_stats_cancer_timepoint_prop_insert$cancer_type_title_case <- gsub("Head_and_neck_cancer", "Head_and\nneck_cancer", data_stats_cancer_timepoint_prop_insert$cancer_type_title_case)

## This is the main figure - updated with new samples
## Plot from here
png(filename = "Proportion of short reads in cfMos, March 2025, reordered, validation.png", width = 13, height = 5, units = "in", res = 500)
prop_short_df_validation %>%
  ggplot(aes(x = cancer_type_ordered, y = Twentyto150_pct)) + geom_boxplot(aes(fill=cancer_type_title_case), outlier.shape = NA) +
  geom_jitter(aes(x = cancer_type_ordered, y = Twentyto150_pct), width=0.3, size = 0.4, alpha = 0.5)+
  geom_hline(yintercept = 0.189, linetype = "dashed", color = "black") + # added intercept from table
  theme(axis.text.x = element_text(angle = 90)) + 
  geom_text(data = data_stats_cancer_timepoint_prop_insert, aes(x = cancer_type_title_case, y = 0.03, label = N), size = 4) +
  labs(title="Proportion of short fragments between 20-150bp across cancer types", subtitle = "CfMeDIP-seq samples") + 
  ylab("Proportion of fragments between 20-150bp\n") + xlab("Cancer type") +
  scale_fill_manual(values = cancer_type_col, name = "Cancer type") + # Add name for legend title
  theme(panel.grid.major = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )
dev.off()


## Make text bigger 

# Assuming your cancer_type_ordered is a factor, you can modify its levels
# to automatically include line breaks if needed
#levels(prop_short_df_validation$cancer_type_ordered) <- str_wrap(levels(prop_short_df_validation$cancer_type_ordered), width = 10) # Adjust 'width' as needed

png(filename = "Proportion of short reads in cfMos, Dec 2024, reordered, bigger text, 2, validation.png", width = 13, height = 4.5, units = "in", res = 500)
ggplot(prop_short_df_validation, aes(x = cancer_type_ordered, y = Twentyto150_pct)) +
  geom_boxplot(aes(fill=cancer_type_title_case), outlier.shape = NA) +
  geom_jitter(aes(x = cancer_type_ordered, y = Twentyto150_pct), width=0.3, size = 0.4, alpha = 0.5) +
  geom_hline(yintercept = 0.1487008, linetype = "dashed", color = "black") + 
  theme(axis.text.x = element_text(angle = 45, size = 12, vjust = 1, hjust=1)) + # Adjust size as needed
  geom_text(data = data_stats_cancer_timepoint_prop_insert, aes(x = cancer_type_title_case, y = 0.03, label = N), size = 4) +
  labs(title="Proportion of short reads between 20-150bp across cancer types", subtitle = "CfMeDIP-seq samples") + 
  ylab("Proportion of reads under 150bp\n") + xlab("Cancer type") +
  scale_fill_manual(values = cancer_type_col, name = "Cancer type") +
  theme(panel.grid.major = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()


## Updated theme 

theme_updated <- theme(plot.title = element_text(hjust = 0.5, size = 13), 
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

png(filename = "Proportion of short reads in cfMos, March 2025, reordered, bigger text, 5, validation.png", width = 13, height = 4.5, units = "in", res = 500)
ggplot(prop_short_df_validation, aes(x = cancer_type_ordered, y = Twentyto150_pct)) +
  geom_boxplot(aes(fill=cancer_type_title_case), outlier.shape = NA) +
  geom_jitter(aes(x = cancer_type_ordered, y = Twentyto150_pct), width=0.3, size = 0.4, alpha = 0.5) +
  geom_hline(yintercept = 0.189, linetype = "dashed", color = "black") + 
  theme(axis.text.x = element_text(angle = 45, size = 12, vjust = 1, hjust=1)) + # Adjust size as needed
  geom_text(data = data_stats_cancer_timepoint_prop_insert, aes(x = cancer_type_title_case, y = 0.03, label = N), size = 4) +
  labs(title="") + 
  ylab("Proportion of fragments under 150bp\n") + xlab("Cancer type") +
  scale_fill_manual(values = cancer_type_col_split, name = "Cancer type") +
  theme_updated +
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



## Now export the dataframes required for machine learning 
### Export the dataframes
saveRDS(insert_matrix_df_validation, file = "Prop_insert_sizes_updated_Dec2024_validation.rds")
write.table(insert_matrix_df_validation, file = "Prop_insert_sizes_updated_Dec2024_validation.txt", sep = "\t", quote = F, row.names = T)

## Export the PCAs on the dataframe 

# Perform PCA
pca_result <- prcomp(insert_matrix, center = TRUE, scale. = TRUE)
pca_result_insertsize <- pca_result
# Determine the number of components to retain
# Visualize the variance explained by each principal component
fviz_eig(pca_result, addlabels = TRUE, ylim = c(0, 40))

# Extract the proportion of variance explained by each component
explained_variance <- summary(pca_result_insertsize)$importance[2,]
cum_variance <- cumsum(explained_variance)

# Find the number of components that explain at least 97% of the variance
num_components <- which(cum_variance >= 0.97)[1]

# Extract the principal components
pca_df <- as.data.frame(pca_result$x[, 1:num_components])
rownames(pca_df) <- rownames(insert_matrix_df_validation)

# Create the variance explained plot with customized theme
variance_plot <- fviz_eig(pca_result, addlabels = TRUE, ylim = c(0, 40)) +
  theme_classic() +
  ggtitle("Scree plot for insert sizes") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

# Save the plot to a file
ggsave(filename = file.path(outdir, "PCA_Variance_Explained_Insert_sizes_validation.png"), plot = variance_plot, width = 10, height = 6, units = "in", dpi = 300)

# Save the PCA results to a file 
write.csv(pca_df, file.path(outdir, "PCA_results_Insert_sizes_validation.csv"), row.names = TRUE)
saveRDS(pca_df, file = "Insert_sizes_updated_Dec2024_PCAs_97_percent_variance_validation.rds")

## Export the metadata df 
write.table(metadata_df, file = "Metadata df June 2024 annotated.txt", sep = "\t", quote = F, row.names = F)
saveRDS(metadata_df, file = "Metadata df June 2024 annotated_validation.rds")
saveRDS(pca_result_insertsize, file = "PCA result for insert size_validation.rds")



## Export for Althaf 

# Export as .txt
write.table(insert_matrix_df_validation, file = "insert_matrix_df_validation_for_Althaf_updated_validation.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Export as _validation.rds
saveRDS(insert_matrix_df_validation, file = "insert_matrix_df_validation_for_Althaf_updated_validation_validation.rds")

## Export for the supplemental tables
write.table(insert_matrix_df_validation, file = "Insert_matrix_dataframe_validation.txt", sep = "\t", row.names = TRUE, quote = F)





## Do some stats for prop short frags

# Check if assumptions for normality met 
prop_short_df_validation %>% 
  group_by(cancer_type_corrected_updated) %>% 
  summarise(shapiro_test = shapiro.test(prop)$p.value)
## Not met

# Check variance 
library(car)

leveneTest(prop ~ cancer_type_corrected_updated, data = prop_short_df_validation)

#difference between groups 
anova_result <- aov(prop ~ cancer_type_corrected_updated, data = prop_short_df_validation)
summary(anova_result)

kruskal_result_prop <- kruskal.test(prop ~ cancer_type_corrected_updated, data = prop_short_df_validation)
print(kruskal_result_prop)

# what is different to normal 
# Conduct pairwise t-tests comparing each group to "normal", with p-value adjustment
pairwise_results <- pairwise.t.test(x = prop_short_df_validation$prop,
                                    g = prop_short_df_validation$cancer_type_corrected_updated,
                                    p.adjust.method = "holm",
                                    paired = FALSE,
                                    pool.sd = !is.factor(prop_short_df_validation$cancer_type_corrected_updated))

# Print the results
print(pairwise_results)

# Extract the p-value matrix
p_values_matrix <- pairwise_results$p.value


# Do non-parametric test
# Perform Dunn's test for post-hoc analysis
library(dunn.test)

dunn_result <- dunn.test(x = prop_short_df_validation$prop, 
                         g = prop_short_df_validation$cancer_type_corrected_updated, 
                         method = "bh") # Adjusting p-values using the Benjamini-Hochberg method


## Extract Dunn files 
# Create a data frame for the Dunn test results
dunn_df <- data.frame(
  Comparison = dunn_result$comparisons,
  Adjusted_P_Value = dunn_result$P.adjusted
)

# Now, filter for comparisons involving the "normal" group
normal_comparisons <- dunn_df[grep("healthy", dunn_df$Comparison), ]

# Sort by Adjusted_P_Value to highlight the most significant differences
normal_comparisons_sorted <- normal_comparisons[order(normal_comparisons$Adjusted_P_Value), ]

# If you want to export this as a CSV file
write.csv(normal_comparisons_sorted, "normal_comparisons_dunn_test_proportion_short_frags_Dec2024_validation.csv", row.names = FALSE)



## Get direction
# Calculate median or mean for each group
group_medians <- prop_short_df_validation %>%
  group_by(cancer_type_corrected_updated) %>%
  summarise(median_prop = median(prop))



## Do some additional stats on the proportion of short fragments

# Check if assumptions for normality met 
prop_short_df_validation %>% 
  group_by(cancer_type_corrected_updated) %>% 
  summarise(shapiro_test = shapiro.test(prop)$p.value)
## Not met

# Check variance 
library(car)

leveneTest(prop ~ cancer_type_corrected_updated, data = prop_short_df_validation)

#difference between groups 
anova_result <- aov(prop ~ cancer_type_corrected_updated, data = prop_short_df_validation)
summary(anova_result)

kruskal_result_prop <- kruskal.test(prop ~ cancer_type_corrected_updated, data = prop_short_df_validation)
print(kruskal_result_prop)

# what is different to normal 
# Conduct pairwise t-tests comparing each group to "normal", with p-value adjustment
pairwise_results <- pairwise.t.test(x = prop_short_df_validation$prop,
                                    g = prop_short_df_validation$cancer_type_corrected_updated,
                                    p.adjust.method = "holm",
                                    paired = FALSE,
                                    pool.sd = !is.factor(prop_short_df_validation$cancer_type_corrected_updated))

# Print the results
print(pairwise_results)

# Extract the p-value matrix
p_values_matrix <- pairwise_results$p.value


# Do non-parametric test
# Perform Dunn's test for post-hoc analysis
library(dunn.test)

dunn_result <- dunn.test(x = prop_short_df_validation$prop, 
                         g = prop_short_df_validation$cancer_type_corrected_updated, 
                         method = "bh") # Adjusting p-values using the Benjamini-Hochberg method


## Extract Dunn files 
# Create a data frame for the Dunn test results
dunn_df <- data.frame(
  Comparison = dunn_result$comparisons,
  Adjusted_P_Value = dunn_result$P.adjusted
)

# Now, filter for comparisons involving the "normal" group
normal_comparisons <- dunn_df[grep("healthy", dunn_df$Comparison), ]

# Sort by Adjusted_P_Value to highlight the most significant differences
normal_comparisons_sorted <- normal_comparisons[order(normal_comparisons$Adjusted_P_Value), ]

# If you want to export this as a CSV file
write.csv(normal_comparisons_sorted, "normal_comparisons_dunn_test_proportion_short_frags_Dec2024_validation_2.csv", row.names = FALSE)


# Export the ANOVA results 
write.table(p_values_matrix, "Pairwise results for the proportion of short fragments Dec 2024 validation.txt", quote = F, sep = "\t")



# Group by cancer_type_title_case and project_id, then count the number of rows
row_counts <- metadata_df %>%
  group_by(cancer_type_title_case, project_id) %>%
  summarize(row_count = n(), .groups = "drop")

# Print the result
print(row_counts)

# Optionally, save the results to a CSV file
write.csv(row_counts, "Row_counts_by_cancer_type_and_project_validation.csv", row.names = FALSE)
