## Creating the heatmap of all metrics together
## Author: Dory Abelman
## date: Jan 2024. Updated May 2024.

## Things to load in
# - Nucleosome peak zscore
# - Insert size zscore
# - End motif zscore 
# - Delfi zscore 
# - Prop short fragments
# Potentially TF 

## Combine everything together into one dataframe 

## Load libraries
library(tidyverse)
library(GGally)
library(dplyr)
library(tidyr)
library(circlize)
library(ComplexHeatmap)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(psych)
library(car)
library(multcomp)
library(FSA)

## Load in the input dataframes 
# Define the directory where the dataframes are stored
data_dir <- "Final Data May 2024/Zscores/"

# Load the dataframes
# These are generated from the previous analysis scripts
Endmotif_df_zscores <- readRDS(file.path(data_dir, "Zscores df for end motifs.rds"))
Insert_df_zscores <- readRDS(file.path(data_dir, "Zscores df for insert size.rds"))
DELFI_df_zscores <- readRDS(file.path(data_dir, "Zscore_df_DELFI.rds"))
Nucleosome_peak_df_zscores <- readRDS(file.path(data_dir, "Zscores df for Nucelosome peak.rds"))
prop_short_df <- readRDS(file.path(data_dir, "Prop 20 to 150 bp in cfMeDip.rds"))

# Give each new column a different name so they will all join together
# These are generated from the previous scripts - for clarity can also export and re-upload here
Endmotif_df_zscores$Endmotif_zscore <- Endmotif_df_zscores$all_zscores
Insert_df_zscores$Insertsize_zscore <- Insert_df_zscores$all_zscores
DELFI_df_zscore$Delfi_zscore <- DELFI_df_zscore$all_zscores
Nucleosome_peak_df_zscores$Nucleosome_peak_zscore <- Nucleosome_peak_df_zscores$all_zscores

## Now join together 
Merged_zscore_df <- Endmotif_df_zscores %>% dplyr::select(sequencing_id, Endmotif_zscore) %>%
  left_join(Insert_df_zscores %>% dplyr::select(sequencing_id, Insertsize_zscore)) %>%
  left_join(DELFI_df_zscore %>% dplyr::select(sequencing_id, Delfi_zscore)) %>%
  left_join(Nucleosome_peak_df_zscores %>% dplyr::select(sequencing_id, Nucleosome_peak_zscore)) %>%
  left_join(prop_short_df %>% dplyr::select(sequencing_id, Twentyto150_pct))

## Now add to metadata df
metadata_df <- readRDS("Final Data May 2024/Machine learning/Metadata df June 2024 all samples.rds")
Merged_zscore_df <- Merged_zscore_df %>% 
  left_join(metadata_df, by = "sequencing_id")
Merged_zscore_df$Twentyto150_pct <- Merged_zscore_df$Twentyto150_pct*100 # Convert to percent

## Add in Althaf's data 
Fragment_scores <- readRDS("Fragment scores from Althaf/All_samples_FS.rds")

## Join fo the zscore df 
Merged_zscore_df <- Merged_zscore_df %>% 
  left_join(Fragment_scores %>% dplyr::select(sequencing_id, FS), by = "sequencing_id")

# Save
saveRDS(Merged_zscore_df, file = "Merged zscore df May 2024.rds")

## Create matrix for heatmap
# Select the relevant columns and set sequencing_id as row names
matrix_heatmap_fig <- Merged_zscore_df %>%
  dplyr::select(sequencing_id, Endmotif_zscore, Insertsize_zscore, Delfi_zscore, Nucleosome_peak_zscore) %>%
  column_to_rownames(var = "sequencing_id")

# Convert to matrix
matrix_heatmap_fig <- as.matrix(matrix_heatmap_fig)

## Now make the heatmap
sample_order <- rownames(matrix_heatmap_fig)

# getting IDs
data_id <- as.matrix(Merged_zscore_df$sequencing_id)
#data_sex <- as.matrix(Merged_zscore_df$sex)
#data_ichorCNA <- as.matrix(Merged_zscore_df$Tumor_fraction)
data_cancer_type <- as.matrix(Merged_zscore_df$cancer_type_title_case)
data_age <- as.matrix(Merged_zscore_df$age)
data_cohort <- as.matrix(Merged_zscore_df$project_id)
data_propshortfrag <- as.matrix(Merged_zscore_df$Twentyto150_pct)
data_fragment_score <- as.matrix(Merged_zscore_df$FS)

row.names(data_id) <- sample_order
#row.names(data_sex) <- sample_order
#row.names(data_ichorCNA) <- sample_order
row.names(data_cancer_type) <- sample_order
row.names(data_age) <- sample_order
row.names(data_cohort) <- sample_order
row.names(data_propshortfrag) <- sample_order
row.names(data_fragment_score) <- sample_order



#data_sex <- factor(data_sex, levels = c("Male", "Female", ""),
#                   labels = c("Male", "Female", "Not available"))
data_cancer_type <- factor(data_cancer_type, levels = c("Healthy", "Prostate_cancer", "Head_and_neck_cancer", "Brain_cancer", 
                                                        "AML", "Eye_cancer", "LFS_survivor", "LFS_previvor", 
                                                        "LFS_positive", "Lung_cancer"),
                           labels = c("Healthy", "Prostate_cancer", "Head_and_neck_cancer", "Brain_cancer", 
                                      "AML", "Eye_cancer", "LFS_survivor", "LFS_previvor", 
                                      "LFS_positive", "Lung_cancer"))
data_cohort <- factor(data_cohort, levels = c("TCGE-CFMe-HBC", "TCGE-CFMe-MCA", "TCGE-CFMe-BCA", "TCGE-CFMe-PRAD", "TCGE-CFMe-HNSC", "TCGE-CFMe-SCLC", "TCGE-CFMe-UM", "TCGE-CFMe-AML", "TCGE-CFMe-LFS"),
                      labels = c(levels = c("TCGE-CFMe-HBC", "TCGE-CFMe-MCA", "TCGE-CFMe-BCA", "TCGE-CFMe-PRAD", "TCGE-CFMe-HNSC", "TCGE-CFMe-SCLC", "TCGE-CFMe-UM", "TCGE-CFMe-AML", "TCGE-CFMe-LFS")))


### Set annotation and heatmap colours
col_heat_together <- colorRamp2(c(-2, 0, 10), 
                       c("#1f78b4", "white", "#e31a1c"))
col_fun <- colorRamp2(c(-2, 0, 40), 
                      c("#1f78b4", "white", "#e31a1c"))
#col_sex <- c("Male" = "#1b9e77",
#             "Female" = "#d95f02",
#             "Not available" = "lightgrey")


#cancer_colors <- color_func(15)
# Match Yong's 
cancer_type_col <- c('Healthy' = '#33a02c',
                     'Brain_cancer' = '#1f78b4',
                     'Lung_cancer' = '#b2df8a',
                     'Prostate_cancer' = '#a6cee3',
                     'AML' = '#fb9a99',
                     'Pancreatic_cancer' = '#e31a1c',
                     'Eye_cancer' = '#fdbf6f',
                     'Uveal_melanoma' = '#fdbf6f',
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


### Generate annotations
## Keep only top annotation 

right_annotation <- rowAnnotation("Fragment score" = anno_lines(data_fragment_score,
                                                                      add_points = TRUE,
                                                                      pch = c(16, NA),
                                                                      ylim = c(-1.5,1.5), 
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
                        #          annotation_legend_param = list(border = TRUE))

left_annotation <- rowAnnotation(#"Sex" = data_sex,
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

data_technology <- colnames(matrix_heatmap_fig)
column_split <- data_technology

top_annotation <- HeatmapAnnotation("Technology" = data_technology)

# Remove unnecessary columns if still included, ie 'Twentyto150_pct'
#matrix_heatmap_fig <- subset(matrix_heatmap_fig, select = -Twentyto150_pct)
row_split <- data_cancer_type

Heatmap_altogether_narrow_May2024 <- Heatmap(matrix_heatmap_fig,
                              col = col_heat_together,
                              show_row_names = FALSE,
                              show_column_names = FALSE,
                              right_annotation = right_annotation,
                           #   left_annotation = left_annotation, # remove this annotation
                #              top_annotation = top_annotation,
                              #heatmap_legend_param = heatmap_legend_param,
                              #       row_order = sample_order,
                              row_labels = NULL,
                              column_title_gp = gpar(fontsize = 12),
                              row_title_gp = gpar(fontsize = 14),
                              row_split = row_split,
                              column_split = column_split,
                              cluster_row_slices = FALSE,
                              row_title_rot = 0,
                              border = FALSE,
                              heatmap_legend_param = list(title = "Zscore to normal", 
                                                          title_gp = gpar(fontsize = 15),  # Font size for the legend title
                                                          labels_gp = gpar(fontsize = 14),
                                                          at = seq(-2, 40, by = 20))) # can remove the at and seq for other breaks
ht_opt$COLUMN_ANNO_PADDING = unit(8, "mm")

# What had before 
#                    heatmap_legend_param = list(title = "Score"))  # Set the title of the legend to "Proportion"


png(filename = "Heatmap of frag metrics all together, May 2024, narrow 3.5.png", width = 13, height = 15, units = "in", res = 350)
draw(Heatmap_altogether_narrow_May2024, heatmap_legend_side = "right", merge_legend = TRUE)
dev.off()


## Now see the proportion of difference in z-score between each cancer type 
cutoff <- 2 # for p of ~0.05

# Calculate the proportion of significant z-scores for each cancer type and z-score column
significant_proportions <- Merged_zscore_df %>%
  gather(key = "zscore_type", value = "zscore", Endmotif_zscore, Insertsize_zscore, Delfi_zscore, Nucleosome_peak_zscore) %>%
  dplyr::mutate(is_significant = abs(zscore) > cutoff) %>%
  dplyr::group_by(cancer_type_corrected_updated, zscore_type) %>%
  dplyr::summarise(proportion_significant = mean(is_significant, na.rm = TRUE)) %>%
  ungroup()

# View the results
significant_proportions


### Do stats to see if correlated 
# Shapiro-Wilk test for normality
shapiro.test(Merged_zscore_df$Endmotif_zscore)
shapiro.test(Merged_zscore_df$Insertsize_zscore)
shapiro.test(Merged_zscore_df$Delfi_zscore)
shapiro.test(Merged_zscore_df$Nucleosome_peak_zscore)
shapiro.test(Merged_zscore_df$Twentyto150_pct)
shapiro.test(Merged_zscore_df$FS)

# See what is correlated 
# Not normally distributed, based on Shapiro-Wilk test
subset_df_expanded <- Merged_zscore_df[, c("Endmotif_zscore", "Insertsize_zscore", "Delfi_zscore", "Nucleosome_peak_zscore", "Twentyto150_pct", "FS", "cancer_type_corrected_updated", "sex")]
subset_df <- Merged_zscore_df[, c("Endmotif_zscore", "Insertsize_zscore", "Delfi_zscore", "Nucleosome_peak_zscore", "Twentyto150_pct", "FS")]

# Calculate correlation matrix using Spearman method
cor_matrix <- cor(subset_df, method = "spearman", use = "complete.obs")

# Calculate p-values and adjust them
p_matrix <- psych::corr.test(subset_df, method = "spearman", adjust = "none")
adj_p_values <- p.adjust(p_matrix$p, method = "BH")

# Create adjusted p-value matrix (same structure as p_matrix$p)
adj_p_matrix <- matrix(adj_p_values, nrow = nrow(p_matrix$p), ncol = ncol(p_matrix$p))
rownames(adj_p_matrix) <- colnames(subset_df)
colnames(adj_p_matrix) <- colnames(subset_df)

# Visualize the correlation matrix
ggcorr(subset_df, method = c("everything", "spearman"), label = TRUE, label_round = 2)
ggpairs(subset_df)

ggpairs(subset_df_expanded)

# Visualize with ggally and plot
corr_plot <- ggpairs(subset_df, 
        lower = list(continuous = wrap("points", alpha = 0.5)),
        upper = list(continuous = wrap("cor", method = "spearman")),
        diag = list(continuous = "barDiag"))

ggsave(filename = "Correlation_matrix_all_features_May2024.pdf", plot = corr_plot, width = 12, height = 12, units = "in", dpi = 500)

## With colors 

subset_df <- Merged_zscore_df[, c("Endmotif_zscore", "Insertsize_zscore", "Delfi_zscore", "Nucleosome_peak_zscore", "Twentyto150_pct", "FS", "cancer_type_title_case")]

corr_plot_colored <- ggpairs(subset_df[,c(1:7)],
                             aes(colour = cancer_type_title_case),
                             lower = list(continuous = wrap("points", alpha = 0.5)),
                             upper = list(continuous = wrap("cor", method = "spearman")),
                             diag = list(continuous = "barDiag")) +
  scale_colour_manual(values = cancer_type_col)

ggsave(filename = "Correlation_matrix_all_features_colored_May2024.pdf", plot = corr_plot_colored, width = 14, height = 12, units = "in", dpi = 500)

## Trying manually 

# set dplyr functions
select <- dplyr::select; rename <- dplyr::rename; mutate <- dplyr::mutate; 
summarize <- dplyr::summarize; arrange <- dplyr::arrange; slice <- dplyr::slice; filter <- dplyr::filter; recode<-dplyr::recode

## Based on https://stackoverflow.com/questions/55657386/how-to-display-coloured-group-correlations-with-scale-colour-manual-in-ggpairs
mycorrelations <- function(data,mapping,...){
  data2 = data
  data2$x = as.numeric(data[,as_label(mapping$x)])
  data2$y = as.numeric(data[,as_label(mapping$y)])
  data2$group = data[,as_label(mapping$colour)]
  
  correlation_df = data2 %>% 
    bind_rows(data2 %>% mutate(group="Overall Corr")) %>%
    group_by(group) %>% 
    filter(sum(!is.na(x),na.rm=T)>1) %>%
    filter(sum(!is.na(y),na.rm=T)>1) %>%
    summarize(estimate = round(as.numeric(cor.test(x,y,method="spearman")$estimate),2),
              pvalue = cor.test(x,y,method="spearman")$p.value,
              pvalue_star = as.character(symnum(pvalue, corr = FALSE, na = FALSE, 
                                                cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                                symbols = c("***", "**", "*", "'", " "))))%>%
    group_by() %>%
    mutate(group = factor(group, levels=c(as.character(unique(sort(data[,as_label(mapping$colour)]))), "Overall Corr")))
  
  ggplot(data = correlation_df, aes(x = 1, y = group, color = group)) +
    geom_text(aes(label = paste0(group, ": ", estimate, pvalue_star))) +
    theme_minimal() + # Apply the minimal theme
    theme(panel.grid.major = element_blank(), # Remove major gridlines
          panel.grid.minor = element_blank())
  }

## Works!! For diagonal color
mydiagonals <- function(data, mapping, ...){
  # Extract the variable name for x and the group (coloring) variable
  x_var <- as_label(mapping$x)
  group_var <- as_label(mapping$colour)
  
  # Check if 'group_var' is correctly specified and exists in 'data'
  if(!group_var %in% names(data)) {
    stop("The specified group variable does not exist in the data.")
  }
  
  # Create the ggplot for the diagonal with appropriate coloring
  p <- ggplot(data, aes(x = !!sym(x_var), fill = !!sym(group_var))) +
    geom_density(alpha = 0.5) +  # Use only fill for the densities
    scale_fill_manual(values = cancer_type_col) +  # Set fill colors for the densities
    theme_minimal() +
    theme(legend.position = "none",  # No legend needed
          panel.grid.major = element_blank(),  # Remove major gridlines
          panel.grid.minor = element_blank(),  # Remove minor gridlines
          axis.title.x = element_blank(),  # Remove x axis title
          axis.text.x = element_blank(),  # Remove x axis text
          axis.ticks.x = element_blank())  # Remove x axis ticks
  
  # Print the group levels for debugging
  print(levels(as.factor(data[[group_var]])))
  
  return(p)
}

## Without new theme 
mydiagonals <- function(data, mapping, ...){
  # Extract the variable name for x and the group (coloring) variable
  x_var <- as_label(mapping$x)
  group_var <- as_label(mapping$colour)
  
  # Check if 'group_var' is correctly specified and exists in 'data'
  if(!group_var %in% names(data)) {
    stop("The specified group variable does not exist in the data.")
  }
  
  # Create the ggplot for the diagonal with appropriate coloring
  p <- ggplot(data, aes(x = !!sym(x_var), fill = !!sym(group_var))) +
    geom_density(alpha = 0.5) +  # Use only fill for the densities
    scale_fill_manual(values = cancer_type_col) # Set fill colors for the densities
 
  # Print the group levels for debugging
  print(levels(as.factor(data[[group_var]])))
  
  return(p)
}

### Getting diagonal to be correct color 
corr_plot_colored_manual2 <- ggpairs(subset_df,columns=1:6,
                                    upper = list(continuous = mycorrelations), 
                                    diag = list(continuous = mydiagonals),
                                    mapping = ggplot2::aes(color=cancer_type_title_case)) +
  scale_colour_manual(values = cancer_type_col) 

ggsave(filename = "Correlation_matrix_all_features_colored_manual 3 may 2024.pdf", plot = corr_plot_colored_manual2, width = 14, height = 12, units = "in", dpi = 500)


# Create a new dataframe with the renamed columns
subset_df_renamed <- subset_df %>%
  rename(
    Prop_short_frags = Twentyto150_pct,
    Fragment_size_score = FS
  )


corr_plot_colored_manual3 <- ggpairs(subset_df_renamed,columns=1:6,
                                     upper = list(continuous = mycorrelations), 
                                     diag = list(continuous = mydiagonals),
                                     mapping = ggplot2::aes(color=cancer_type_title_case)) +
  scale_colour_manual(values = cancer_type_col) 

ggsave(filename = "Correlation_matrix_all_features_colored_manual 4 may 2024.pdf", plot = corr_plot_colored_manual3, width = 14, height = 12, units = "in", dpi = 500)

## Correct naming 
# Replace "Eye_cancer" with "Uveal_melanoma" in the subset_df_renamed dataframe
subset_df_renamed <- subset_df_renamed %>%
  dplyr::mutate(cancer_type_title_case = ifelse(cancer_type_title_case == "Eye_cancer", "Uveal_melanoma", cancer_type_title_case))


## Export the z-score table as a supplementary table 

Z_score_table_export <- Merged_zscore_df[, c("Endmotif_zscore", "Insertsize_zscore", "Delfi_zscore", "Nucleosome_peak_zscore", "Twentyto150_pct", "FS", "cancer_type_title_case", "sample_id")]
Z_score_table_export <- Z_score_table_export %>%
  dplyr::mutate(cancer_type_title_case = ifelse(cancer_type_title_case == "Eye_cancer", "Uveal_melanoma", cancer_type_title_case))
Z_score_table_export <- Z_score_table_export %>%
  rename(
    Prop_short_frags = Twentyto150_pct,
    Fragment_size_score = FS
  )

## Export 
write.table(Z_score_table_export, sep = "\t", quote = F, row.names = F, file = "Genome-wide Z-score table per feature.txt")

## Rerun 
corr_plot_colored_manual4 <- ggpairs(subset_df_renamed,columns=1:6,
                                     upper = list(continuous = mycorrelations), 
                                     diag = list(continuous = mydiagonals),
                                     mapping = ggplot2::aes(color=cancer_type_title_case)) +
  scale_colour_manual(values = cancer_type_col) 

ggsave(filename = "Correlation_matrix_all_features_colored_manual 5, naming corrected, may 2024.pdf", plot = corr_plot_colored_manual4, width = 14, height = 12, units = "in", dpi = 500)



## Get stats for fragment scores - which are significantly different? 

# Subset data for analysis
df <- Merged_zscore_df %>%
  select(cancer_type_title_case, FS) %>%
  filter(!is.na(FS))

# Check if data meets ANOVA assumptions
# 1. Homogeneity of variances
leveneTest(FS ~ cancer_type_title_case, data = df)

# 2. Normality of residuals
anova_model <- aov(FS ~ cancer_type_title_case, data = df)
qqnorm(residuals(anova_model))
qqline(residuals(anova_model)) # assumptions not met

# If assumptions are met, proceed with ANOVA; otherwise, use Kruskal-Wallis test
# Adjust this based on the outcome of the assumption tests

# Assuming ANOVA assumptions are met
anova_results <- aov(FS ~ cancer_type_title_case, data = df)
summary(anova_results)

# Post-hoc test to compare each cancer type against "Healthy"
post_hoc_results <- TukeyHSD(anova_results)
post_hoc_results

# If ANOVA assumptions are not met, use Kruskal-Wallis and Dunn's test
kruskal_test <- kruskal.test(FS ~ cancer_type_title_case, data = df)

# Dunn's test for post-hoc comparisons
dunn_test <- dunnTest(df$FS, df$cancer_type_title_case, method="bonferroni")

# Extract rows where 'Healthy' is part of the comparison
normal_comparisons <- dunn_test$res[grep("Healthy", dunn_test$res$Comparison), ]

# If you want to specifically compare against 'Normal', you might need to ensure that 'Normal' is in every comparison string
# This step is just to ensure that the grep function above should work as intended if 'Normal' is indeed a part of the comparisons

# Exporting the results to a CSV file for further inspection
write.csv(normal_comparisons, "normal_comparisons_dunn_test_FS_June2024.csv", row.names = FALSE)

# Print to check if any rows are caught
print(normal_comparisons)

