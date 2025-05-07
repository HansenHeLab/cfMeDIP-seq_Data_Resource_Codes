# End motifs - heatmaps, plots and stats

# author: Dory Abelman and Derek Wong
# adapted from scripts by Derek Wong
# Original date: Aug 2023
# Updated for new samples Oct-Nov 2023
# Updated May 2024 for final set of adjusted prostate cancer samples

# Publication adapted from:
# Link: https://github.com/pughlab/TGL48_Uveal_Melanoma/blob/main/R_scripts_figures/Figure%202/Figure%202%20-%20Fragment%20Ratio.R
# Title: Integrated, Longitudinal Analysis of Cell-free DNA in Uveal Melanoma (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9973415/)

# Load packages
library(tidyverse)
library(dplyr)
library(matrixStats)
library(cowplot)
library(ggh4x)
# library(zplyr)
library(patchwork)
library(ggpubr)
library(chisq.posthoc.test)
library(broom)
library(factoextra) 
library(ComplexHeatmap)
library(colorRamp2)



### Set paths
path <- "end_motifs/Data_20_600bp_frags/"
outdir <- "Final Data Dec 2024/End_motifs/"
healthy_path <- "end_motifs/Healthy_data_20_600bp_frags/"

### Find paths
data_files <- list.files(path, "5prime_motifs\\.txt$", full.names = TRUE)
data_normal_files <- list.files(healthy_path, "5prime_motifs\\.txt$", full.names = TRUE)

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

# Now join to metadata df 
metadata_df <- readRDS("Metadata_df_all_samples_with_CN_classifier_Jan_2025.rds")

## Annotate the metadata df
metadata_df$sample_id <- paste(metadata_df$sequencing_id, "_dedup", sep = "")
metadata_df$cancer_type_corrected_updated <- tolower(gsub(" ", "_", metadata_df$cancer_type))
metadata_df$cancer_type_corrected_updated <- as.character(metadata_df$cancer_type_corrected_updated)
metadata_df$cancer_type_corrected_updated <- gsub('blood_cancer', 'Aml', metadata_df$cancer_type_corrected_updated)
metadata_df$cancer_type_corrected_updated <- gsub('normal', 'healthy', metadata_df$cancer_type_corrected_updated)
metadata_df$cancer_type_title_case <- tools::toTitleCase(metadata_df$cancer_type_corrected_updated)
metadata_df$cancer_type_title_case  <- gsub("Lfs", "LFS", metadata_df$cancer_type_title_case)
metadata_df$cancer_type_title_case  <- gsub("Aml", "AML", metadata_df$cancer_type_title_case)

metadata_df$endmotif_name <- paste(metadata_df$sequencing_id, "_dedup_motifs", sep = "")



# Filter and summarize the data based on conditions
metadata_summary <- metadata_df %>%
  filter(type %in% c("PE", "SE") | validation %in% c(1, 0)) %>%
  group_by(cancer_type_title_case, type, validation) %>%
  summarise(row_count = n(), .groups = "drop")

## Bind rows
data_endmotif <- bind_rows(data, data_normal)

# Remove the 5prime from naming and add to metadata df
data_endmotif$sample <- gsub("5prime_", "", data_endmotif$sample)
data_endmotif <- left_join(data_endmotif, metadata_df, by = c("sample" = "endmotif_name"))

# Get frequency
data_endmotif$frequency <- data_endmotif$frequency*100

# Remove ambiguous or low quality motifs and recalculate frequencies
data_endmotif_filtered <- data_endmotif %>%
  dplyr::filter(!str_detect(motif, "[KYSWRBMV]")) %>% 
  dplyr::group_by(sample) %>%
  dplyr::mutate(total_freq = sum(frequency)) %>%
  ungroup() %>%
  dplyr::mutate(normalized_frequency = frequency / total_freq) %>%
  dplyr::select(-total_freq)

#test <- data_endmotif_filtered %>% filter(Sample)
# Keep as percent to stay consistent
data_endmotif_filtered$normalized_frequency <- data_endmotif_filtered$normalized_frequency*100

# Remove the samples not in metadata df 
data_endmotif_filtered <- data_endmotif_filtered %>% dplyr::filter(!is.na(cancer_type))

### Now remove the SE and validation for time being 
data_endmotif_filtered_SE <- data_endmotif_filtered %>% filter(type == "SE")
data_endmotif_filtered_validation <- data_endmotif_filtered %>% filter(validation == "1") %>% filter(type == "PE")
#data_endmotif_filtered <- data_endmotif_filtered %>% filter(validation == "0") %>% filter(type == "PE")


### Do the SE 
data_endmotif_filtered <- data_endmotif_filtered_SE


### Calculate stats
data_stats <- data_endmotif_filtered %>%
  dplyr::group_by(motif, cancer_type_corrected_updated) %>%
  dplyr::summarise(mean=mean(normalized_frequency),
                   sd=sd(normalized_frequency))


motifs <- unique(data_stats$motif)
patients <- unique(data_stats$cancer_type_corrected_updated)


# Initialize vectors for storing t-test results
x <- c()
y <- c()
motif_vec <- c()
patient_vec <- c()

# Do t-test if sufficient data on differences in end motifs
for (motif in motifs) {
  for (patient in patients) {
    n1 <- sum(data_endmotif_filtered$cancer_type_corrected_updated == "healthy" & data_endmotif_filtered$motif == motif)
    n2 <- sum(data_endmotif_filtered$cancer_type_corrected_updated == patient & data_endmotif_filtered$motif == motif)
    
    if(n1 >= 2 && n2 >= 2) {
      a <- t.test(data_endmotif_filtered$normalized_frequency[data_endmotif_filtered$cancer_type_corrected_updated == "healthy" & data_endmotif_filtered$motif == motif], 
                  data_endmotif_filtered$normalized_frequency[data_endmotif_filtered$cancer_type_corrected_updated == patient & data_endmotif_filtered$motif == motif])$p.value
      x <- c(x, a)
      
      b <- mean(data_endmotif_filtered$normalized_frequency[data_endmotif_filtered$cancer_type_corrected_updated == patient & data_endmotif_filtered$motif == motif]) / 
        mean(data_endmotif_filtered$normalized_frequency[data_endmotif_filtered$cancer_type_corrected_updated == "healthy" & data_endmotif_filtered$motif == motif])
      y <- c(y, b)
    } else {
      x <- c(x, NA)
      y <- c(y, NA)
    }
    motif_vec <- c(motif_vec, motif)
    patient_vec <- c(patient_vec, patient)
  }
}

# Create a new dataframe to store the t-test results
t_test_results <- data.frame(
  motif = motif_vec,
  cancer_type_corrected_updated = patient_vec,
  pvalue = p.adjust(x),
  foldchange = y
)

# Merge it back into data_stats
data_stats <- merge(data_stats, t_test_results, by = c("motif", "cancer_type_corrected_updated"))



### Calculate change of DNASE1L3 motifs
data_dnase <- data_endmotif_filtered[data_endmotif_filtered$motif %in% c("CCCA", "CCAG", "CCTG", "CCAA", "CCCT", "CCAT"), ]

data_stats_dnase <- data_dnase %>%
  group_by(motif, cancer_type_corrected_updated, cancer_type_title_case) %>%
  dplyr::summarise(mean=mean(normalized_frequency),
                   sd=sd(normalized_frequency))

motifs <- unique(data_stats_dnase$motif)
patients <- unique(data_stats_dnase$cancer_type_corrected_updated)

# Initialize vectors for storing t-test results
x <- c()
y <- c()
motif_vec <- c()
patient_vec <- c()

# Do t-test if sufficient data
for (motif in motifs) {
  for (patient in patients) {
    n1 <- sum(data_endmotif_filtered$cancer_type_corrected_updated == "healthy" & data_endmotif_filtered$motif == motif)
    n2 <- sum(data_endmotif_filtered$cancer_type_corrected_updated == patient & data_endmotif_filtered$motif == motif)
    
    if(n1 >= 2 && n2 >= 2) {
      a <- t.test(data_endmotif_filtered$normalized_frequency[data_endmotif_filtered$cancer_type_corrected_updated == "healthy" & data_endmotif_filtered$motif == motif], 
                  data_endmotif_filtered$normalized_frequency[data_endmotif_filtered$cancer_type_corrected_updated == patient & data_endmotif_filtered$motif == motif])$p.value
      x <- c(x, a)
      
      b <- mean(data_endmotif_filtered$normalized_frequency[data_endmotif_filtered$cancer_type_corrected_updated == patient & data_endmotif_filtered$motif == motif]) / 
        mean(data_endmotif_filtered$normalized_frequency[data_endmotif_filtered$cancer_type_corrected_updated == "healthy" & data_endmotif_filtered$motif == motif])
      y <- c(y, b)
    } else {
      x <- c(x, NA)
      y <- c(y, NA)
    }
    motif_vec <- c(motif_vec, motif)
    patient_vec <- c(patient_vec, patient)
  }
}

# Create a new dataframe to store the t-test results
t_test_results <- data.frame(
  motif = motif_vec,
  cancer_type_corrected_updated = patient_vec,
  pvalue = p.adjust(x),
  foldchange = y
)

# Merge it back into data_stats
data_stats_dnase <- merge(data_stats_dnase, t_test_results, by = c("motif", "cancer_type_corrected_updated"))

data_stats_dnase$annot <- ifelse(data_stats_dnase$pvalue < 0.05 & data_stats_dnase$pvalue > 0.01, "*",
                                 ifelse(data_stats_dnase$pvalue < 0.01 & data_stats_dnase$pvalue > 0.001, "**",
                                        ifelse(data_stats_dnase$pvalue < 0.001, "***", "")))





### Calculate median fold changes
# Calculate median frequency for each motif in the 'normal' data
normal_medians <- data_endmotif_filtered %>%
  dplyr::filter(cancer_type_corrected_updated == "healthy") %>%
  group_by(motif) %>%
  dplyr::summarise(median_normal = median(normalized_frequency))

# Initialize empty data frame to store fold changes for each cancer type
fold_change_all <- data.frame()

# Loop through each unique cancer type in data_endmotif_filtered
unique_cancer_types <- unique(data_endmotif_filtered$cancer_type_corrected_updated)
data_endmotif_filtered$cancer_type_corrected_updated <- as.factor(data_endmotif_filtered$cancer_type_corrected_updated)

for (cancer_type in unique_cancer_types) {
  
  # Subset data for current cancer type
  data_cancer <- data_endmotif_filtered[data_endmotif_filtered$cancer_type_corrected_updated == cancer_type, ]
  
  # Calculate median frequency for each motif in the cancer data
  cancer_medians <- data_cancer %>%
    group_by(motif) %>%
    summarise(median_cancer = median(normalized_frequency))
  
  # Merge the normal and cancer medians for each motif
  merged_data <- merge(normal_medians, cancer_medians, by = "motif", all = TRUE)
  
  # Calculate fold change and add cancer type as a new column
  merged_data$fold_change <- merged_data$median_cancer / merged_data$median_normal
  merged_data$cancer_type <- cancer_type
  
  # Get standard deviation
  std_dev_data <- merged_data %>%
    summarise(std_dev_fold_change = sd(fold_change, na.rm = TRUE))
  
  merged_data$sd <- std_dev_data$std_dev_fold_change
  
  # Append to the overall data frame
  fold_change_all <- rbind(fold_change_all, merged_data)
}


# Order by fold_change
fold_change_all <- fold_change_all[order(fold_change_all$fold_change, decreasing = TRUE),]
order <- unique(fold_change_all$motif)
fold_change_all$motif <- factor(fold_change_all$motif, levels = order)

data_freq <- fold_change_all

### Make motif table
data_motif <- as.data.frame(str_split_fixed(order, "", 4))
colnames(data_motif) <- c("5' Base 1", "5' Base 2", "5' Base 3", "5' Base 4")
data_motif$bin <- c(1:nrow(data_motif))
data_motif[data_motif == "A" | data_motif == "T"] <- "A/T"
data_motif[data_motif == "G" | data_motif == "C"] <- "G/C"

### Make frequency change comparison table
data_freq$change <- ifelse(data_freq$fold_change > 1, "Increase", "Decrease")
data_freq$change <- factor(data_freq$change, levels = c("Increase", "Decrease"))
data_freq <- data_freq[order(data_freq$change,
                             data_freq$cancer, decreasing = TRUE), ]
order <- unique(data_freq$motif)
data_freq$motif <- factor(data_freq$motif, levels = order)


### Compare frequencies based on expected frequencies (shannon)
expected_freq <- 100/length(unique(normal_medians$motif)) # 100 since multiplied by it at the beginning
normal_medians$expected <- ifelse(normal_medians$median_normal > expected_freq, "Above", "Below")

data_expected_freq <- merge(normal_medians, fold_change_all)
data_expected_freq$change <- ifelse(data_expected_freq$fold_change > 1, "Increase", "Decrease")



## Function  
# Create an empty list to store results for each cancer type
list_motif2 <- list()
list_reg <- list()
list_shannon <- list()
list_shannon2 <- list()
list_vectors <- list()

# Loop through each unique cancer type
for(cancer_type in unique_cancer_types) {
  
  # Subset the data to only include the current cancer type
  sub_data <- data_freq[data_freq$cancer_type == cancer_type, ]
  
  ### Make motif table 2
  data_motif2 <- as.data.frame(str_split_fixed(order, "", 4))
  data_motif2$sum <- sub_data$change
  colnames(data_motif2) <- c("5' Base 1", "5' Base 2", "5' Base 3", "5' Base 4", "sum")
  data_motif2$bin <- c(1:nrow(data_motif2))
  data_motif2[data_motif2 == "A" | data_motif2 == "T"] <- "A/T"
  data_motif2[data_motif2 == "G" | data_motif2 == "C"] <- "G/C"
  list_motif2[[cancer_type]] <- data_motif2
  
  ### Code for Shannon's index and Chi-squared test
  sub_data2 <- data_expected_freq[data_expected_freq$cancer_type == cancer_type, ]
  stat_shannon <- as.data.frame(table(sub_data2$expected, sub_data2$change))
  
  # Check if expected and change have at least 2 levels
  if(length(unique(sub_data2$expected)) > 1 & length(unique(sub_data2$change)) > 1) {
    chisq_shannon <- chisq.test(sub_data2$expected, sub_data2$change)$p.value
  } else {
    chisq_shannon <- NA
  }
  
  post_hoc <- as.data.frame(table(sub_data2$expected, sub_data2$change))
  
  if(length(unique(post_hoc$Var2)) > 1) {
    chisq_shannon_decrease <- chisq.test(post_hoc$Freq[post_hoc$Var2 == "Decrease"])$p.value
    chisq_shannon_increase <- chisq.test(post_hoc$Freq[post_hoc$Var2 == "Increase"])$p.value
  } else {
    chisq_shannon_decrease <- NA
    chisq_shannon_increase <- NA
  }
  
  my_chars <- strsplit(sub_data2$motif, "")
  char_frequency <- lapply(my_chars, function(x) {table(x)})
  char_frequency <- do.call(bind_rows, char_frequency)
  char_frequency <- as.data.frame(char_frequency)
  char_frequency[is.na(char_frequency)] <- 0
  sub_data2$AT <- char_frequency$A + char_frequency$T
  
  sub_data2$expected <- as.factor(sub_data2$expected)
  sub_data2$change <- as.factor(sub_data2$change)
  sub_data2$AT <- as.factor(sub_data2$AT)
  
  stat_shannon2 <- as.data.frame(table(sub_data2$expected, sub_data2$change, sub_data2$AT))
  
  input <- sub_data2[sub_data2$expected == "Above", ]
  input <- table(input$change, input$AT)
  if(length(unique(input)) > 1) {
    chisq_above <- chisq.test(input)$p.value
  } else {
    chisq_above <- NA
  }
  
  input <- sub_data2[sub_data2$expected == "Below", ]
  input <- table(input$change, input$AT)
  if(length(unique(input)) > 1) {
    chisq_below <- chisq.test(input)$p.value
  } else {
    chisq_below <- NA
  }
  
  list_vectors[[cancer_type]] <- list(chisq_shannon, chisq_shannon_decrease, chisq_shannon_increase, chisq_above, chisq_below)
  list_shannon[[cancer_type]] <- stat_shannon
  list_shannon2[[cancer_type]] <- stat_shannon2
}

# Combine the results into a single DataFrame
final_motif2 <- do.call(rbind, list_motif2)
final_shannon <- do.call(rbind, list_shannon)
final_shannon2 <- do.call(rbind, list_shannon2)
final_vectors <- do.call(rbind, list_vectors)

# Edit so rownames are in column
final_motif2$cancer_type <- rownames(final_motif2)
final_motif2$cancer_type <- sub("\\..*$", "", final_motif2$cancer_type)
final_shannon$cancer_type <- rownames(final_shannon)
final_shannon$cancer_type <- sub("\\..*$", "", final_shannon$cancer_type)
final_shannon2$cancer_type <- rownames(final_shannon2)
final_shannon2$cancer_type <- sub("\\..*$", "", final_shannon2$cancer_type)







### Make regression table

# Initialize an empty data frame to hold the regression data
data_reg <- data.frame()
# Loop through each unique cancer type
for (cancer1 in unique_cancer_types) {
  
  # Select data for the first cancer type
  data1 <- fold_change_all[fold_change_all$cancer_type == cancer1, ]
  
  # Loop through each other unique cancer type
  for (cancer2 in unique_cancer_types) {
    
    # Skip if it's the same as the first type
    if (cancer1 == cancer2) next
    
    # Select data for the second cancer type
    data2 <- fold_change_all[fold_change_all$cancer_type == cancer2, ]
    
    # Find common motifs between the two types
    common_motifs <- intersect(data1$motif, data2$motif)
    
    # Filter data by common motifs
    data1_common <- data1[data1$motif %in% common_motifs, ]
    data2_common <- data2[data2$motif %in% common_motifs, ]
    
    # Create regression data for these pairs
    reg_data <- data.frame(var1 = data1_common$fold_change,
                           var2 = data2_common$fold_change,
                           comp = paste0(cancer1, " vs\n", cancer2))
    
    # Append to the overall data frame
    data_reg <- rbind(data_reg, reg_data)
  }
}

# Convert 'comp' to a factor with the levels defined by the unique comparisons
data_reg$comp <- factor(data_reg$comp, levels = unique(data_reg$comp))


### Set theme
theme <- theme(plot.title = element_text(hjust = 0.5, size = 13), 
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

### Plot comparisons

## Set color scheme 
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


fig <- ggplot(data_dnase) +
  geom_boxplot(aes(cancer_type_title_case, normalized_frequency, fill = cancer_type_title_case), outlier.size = 1, alpha = 1) +
  geom_text(data = data_stats_dnase, aes(cancer_type_title_case, 2, label = annot), size = 5) +
  xlab("Frequency (%)") + 
  ylab("Frequency (%)") +
  labs(color = "", fill = "") +
  scale_fill_manual(values = cancer_type_col) +  # Apply custom color palette
  facet_wrap(.~motif, scales = "free_y", nrow = 1) +
  labs(title = "DNASE1L3 motifs") +
  theme +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
fig

ggsave(file.path(outdir, paste0("fragment_end_contexts_dnase1l3, Mar 2025, updated to all frags 5prime only, 3, SE.pdf")), fig, width = 15, height = 8, dpi = 500, units = "in")








## --- NEW CODE FOR DNASE1L3 FOLD CHANGE EXCLUDING HIGH Z-SCORES --- ##

# Identify healthy samples with acceptable z-scores (<= 3)
healthy_samples_to_include <- End_motif_df_zscores %>%
  filter(cancer_type_corrected_updated == "healthy", all_zscores <= 3) %>%
  pull(Sample)

# Filter DNASE1L3 motif data from the overall filtered data
data_dnase_filtered <- data_endmotif_filtered %>%
  filter(motif %in% c("CCCA", "CCAG", "CCTG", "CCAA", "CCCT", "CCAT"))

# Summarize the DNASE1L3 data (mean and sd) by motif and cancer type
data_stats_dnase_filtered <- data_dnase_filtered %>%
  group_by(motif, cancer_type_corrected_updated, cancer_type_title_case) %>%
  summarise(mean = mean(normalized_frequency),
            sd = sd(normalized_frequency),
            .groups = "drop")

motifs <- unique(data_stats_dnase_filtered$motif)
patients <- unique(data_stats_dnase_filtered$cancer_type_corrected_updated)

# Initialize vectors for storing t-test results and fold change values
pvals <- c()
foldchanges <- c()
motif_vec <- c()
patient_vec <- c()

# Loop over each motif and cancer type (patient group)
for (motif in motifs) {
  for (patient in patients) {
    # For the healthy group, filter out samples with high z-scores
    healthy_values <- data_dnase_filtered$normalized_frequency[
      data_dnase_filtered$cancer_type_corrected_updated == "healthy" &
        data_dnase_filtered$motif == motif &
        data_dnase_filtered$sample %in% healthy_samples_to_include
    ]
    
    # Values for the current patient group
    patient_values <- data_dnase_filtered$normalized_frequency[
      data_dnase_filtered$cancer_type_corrected_updated == patient &
        data_dnase_filtered$motif == motif
    ]
    
    n1 <- length(healthy_values)
    n2 <- length(patient_values)
    
    if (n1 >= 2 && n2 >= 2) {
      # Conduct t-test between healthy (filtered) and current patient group
      p_val <- t.test(healthy_values, patient_values)$p.value
      pvals <- c(pvals, p_val)
      
      # Calculate fold change: (patient mean) / (healthy mean)
      fc <- mean(patient_values) / mean(healthy_values)
      foldchanges <- c(foldchanges, fc)
    } else {
      pvals <- c(pvals, NA)
      foldchanges <- c(foldchanges, NA)
    }
    motif_vec <- c(motif_vec, motif)
    patient_vec <- c(patient_vec, patient)
  }
}


# Create a data frame to store the t-test results with multiple hypothesis adjustment
t_test_results_filtered <- data.frame(
  motif = motif_vec,
  cancer_type_corrected_updated = patient_vec,
  pvalue = p.adjust(pvals),
  foldchange = foldchanges
)

# Merge the t-test results back into the summarized DNASE1L3 statistics
data_stats_dnase_filtered_SE <- merge(data_stats_dnase_filtered, t_test_results_filtered, 
                                   by = c("motif", "cancer_type_corrected_updated"))


data_stats_dnase_filtered_SE$annot <- ifelse(data_stats_dnase_filtered_SE$pvalue < 0.05 & data_stats_dnase_filtered_SE$pvalue > 0.01, "*",
                                 ifelse(data_stats_dnase_filtered_SE$pvalue < 0.01 & data_stats_dnase_filtered_SE$pvalue > 0.001, "**",
                                        ifelse(data_stats_dnase_filtered_SE$pvalue < 0.001, "***", "")))


### replot DNASE-1L3
tmp <- data_dnase %>%
  filter(cancer_type_corrected_updated != "healthy" |
           (cancer_type_corrected_updated == "healthy" & sample %in% healthy_samples_to_include))

fig <- ggplot(tmp) +
  geom_boxplot(aes(cancer_type_title_case, normalized_frequency, fill = cancer_type_title_case), outlier.size = 1, alpha = 1) +
  geom_text(data = data_stats_dnase_filtered_SE, aes(cancer_type_title_case, 2, label = annot), size = 5) +
  xlab("Frequency (%)") + 
  ylab("Frequency (%)") +
  labs(color = "", fill = "") +
  scale_fill_manual(values = cancer_type_col) +  # Apply custom color palette
  facet_wrap(.~motif, scales = "free_y", nrow = 1) +
  labs(title = "DNASE1L3 motifs") +
  theme +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
fig

ggsave(file.path(outdir, paste0("fragment_end_contexts_dnase1l3, Mar 2025, updated to all frags 5prime only, 3, SE, but only good healthies.pdf")), fig, width = 15, height = 8, dpi = 500, units = "in")
















### Plot fold changes

# If head and neck cancer split 
cancer_type_col_split <- c('Healthy' = '#33a02c',
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


## Do the above ordered by variance 
# Calculate the variance for each cancer type

## First make capital 
fold_change_all$cancer_type_title_case <- fold_change_all$cancer_type
fold_change_all$cancer_type_title_case <- tools::toTitleCase(fold_change_all$cancer_type)
fold_change_all$cancer_type_title_case  <- gsub("Lfs", "LFS", fold_change_all$cancer_type_title_case)
fold_change_all$cancer_type_title_case  <- gsub("Uveal_melanoma", "Eye_cancer", fold_change_all$cancer_type_title_case)
fold_change_all$cancer_type_title_case  <- gsub("Head_and_neck_cancer", "Head_and\nneck_cancer", fold_change_all$cancer_type_title_case)
fold_change_all$cancer_type_title_case  <- gsub("Aml", "AML", fold_change_all$cancer_type_title_case)

variance_by_cancer <- fold_change_all %>%
  dplyr::group_by(cancer_type_title_case) %>%
  dplyr::summarise(variance = var(fold_change)) %>%
  dplyr::arrange(desc(variance))

# Reorder cancer_type in fold_change_all based on variance
fold_change_all$cancer_type_title_case <- factor(fold_change_all$cancer_type_title_case, levels = variance_by_cancer$cancer_type_title_case)

# Now create the plot
fig_fold_cancer <- ggplot(fold_change_all, aes(x = cancer_type_title_case, y = fold_change, col = cancer_type_title_case)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(pch = 16, alpha = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", size = 0.5) +
  xlab("Samples") + 
  ylab("Fold Change") +
  scale_color_manual(values = cancer_type_col_split) +
  ggtitle("Fold change of cancers vs normal samples, by motif") + 
  theme +
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

fig_fold_cancer
ggsave(file.path(outdir, paste("fragment_fold_change_all_cancers, Dec 2024, updated by Variance, all frags.pdf", sep = "")), fig_fold_cancer, width = 10, height = 4, dpi = 500, units = "in")


## Now reorder
# Get the ordered cancer types to match the insert size
## since don't have for SE change to without
ordered_cancer_types <- unique(fold_change_all$cancer_type_title_case)

## Custom order reccomended by Yong
# Remove 'Healthy' and 'LFS' from the vector, if they exist
ordered_cancer_types <- ordered_cancer_types[ordered_cancer_types != "Healthy"]


# Prepend 'Healthy' at the beginning and append 'LFS' at the end
ordered_cancer_types <- c("Healthy", ordered_cancer_types)

## Fixing head and neck naming 
ordered_cancer_types_peak <- ordered_cancer_types

# Create an ordered factor with levels in the order of decreasing median
fold_change_all$cancer_type_title_case <- gsub("Head_and_neck_cancer", "Head_and\nneck_cancer", fold_change_all$cancer_type_title_case)
fold_change_all$cancer_type_ordered <- factor(fold_change_all$cancer_type_title_case, levels = ordered_cancer_types_peak)

## Replot 
# Now create the plot
fig_fold_cancer_ordered <- ggplot(fold_change_all, aes(x = cancer_type_ordered, y = fold_change, col = cancer_type_ordered)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(pch = 16, alpha = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", size = 0.5) +
  xlab("Samples") + 
  ylab("Fold Change") +
  scale_color_manual(values = cancer_type_col_split) +
  ggtitle("Fold change of cancers vs normal samples, by motif") + 
  theme +
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

fig_fold_cancer_ordered
ggsave(file.path(outdir, paste("fragment_fold_change_all_cancers, Dec 2024, ordered by prop short insert size, all frags, narrower.pdf", sep = "")), fig_fold_cancer_ordered, width = 10, height = 3, dpi = 500, units = "in")






## Now plot the A/T and C/G ratio for each sample type
## As loop

## With A/T and C/G being the same color 
for(cancer_type in unique_cancer_types){
  # Filter fold_change_all to only the current cancer_type
  fold_change_filtered <- fold_change_all[fold_change_all$cancer_type == cancer_type, ]
  
  ## Set order
  order <- fold_change_filtered$motif
  fold_change_filtered$motif <- factor(fold_change_filtered$motif, levels = order)
  
  # Create fig_fold_original
  fig_fold_original <- ggplot(fold_change_filtered) +
    geom_point(aes(motif, fold_change), pch = 16, alpha = 0.5, color = "#1F78B4") +
    geom_hline(yintercept = 1, linetype = "dashed", size = 0.5) +
    xlab("") + 
    ylab("Fold Change") +
    labs(title = paste("Cancer vs normal:", cancer_type)) +
    theme +
    theme(legend.position = c(0.25, 0.4),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  ### Make motif table
  data_motif <- as.data.frame(str_split_fixed(order, "", 4))
  colnames(data_motif) <- c("5' Base 1", "5' Base 2", "5' Base 3", "5' Base 4")
  data_motif$bin <- c(1:nrow(data_motif))
  data_motif[data_motif == "A" | data_motif == "T"] <- "A/T"
  data_motif[data_motif == "G" | data_motif == "C"] <- "G/C"
  data_motif_melted <- reshape2::melt(data_motif, id = "bin")
  fig_motif <- ggplot(data_motif_melted) +
    geom_tile(aes(bin, variable, fill = value)) +
    scale_fill_manual(values = c("A/T" = "#FC766AFF", "G/C" = "#5B84B1FF")) +
    xlab("Tetranucleotide") +
    ylab("Position") +
    labs(fill = "") +
    theme +
    theme(legend.position = "right",  # moved the legend to the right
          legend.background = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(-1,0.5,0.5,0.5), "lines")) +
    scale_x_continuous(expand = c(0,0))
  
  # Combine the two figures
  figure <- fig_fold_original / fig_motif + plot_layout(heights = c(2, 1))
  
  # Save the combined plot
  ggsave(file.path(outdir, paste("fragment_end_contexts_change_", cancer_type, "May 2024.pdf", sep = "")), figure, width = 10, height = 6, dpi = 500, units = "in")
}

## With A/T and C/G being different colors: 
for(cancer_type in unique_cancer_types){
  # Filter fold_change_all to only the current cancer_type
  fold_change_filtered <- fold_change_all[fold_change_all$cancer_type == cancer_type, ]
  
  ## Set order
  order <- fold_change_filtered$motif
  fold_change_filtered$motif <- factor(fold_change_filtered$motif, levels = order)
  
  # Create fig_fold_original
  fig_fold_original <- ggplot(fold_change_filtered) +
    geom_point(aes(motif, fold_change), pch = 16, alpha = 0.5, color = "#1F78B4") +
    geom_hline(yintercept = 1, linetype = "dashed", size = 0.5) +
    xlab("") + 
    ylab("Fold Change") +
    labs(title = paste("Cancer vs normal:", cancer_type)) +
    theme +
    theme(legend.position = c(0.25, 0.4),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  ### Make motif table
  data_motif <- as.data.frame(str_split_fixed(order, "", 4))
  colnames(data_motif) <- c("5' Base 1", "5' Base 2", "5' Base 3", "5' Base 4")
  data_motif$bin <- c(1:nrow(data_motif))
  # data_motif[data_motif == "A" | data_motif == "T"] <- "A/T"
  # data_motif[data_motif == "G" | data_motif == "C"] <- "G/C"
  data_motif_melted <- reshape2::melt(data_motif, id = "bin")
  fig_motif <- ggplot(data_motif_melted) +
    geom_tile(aes(bin, variable, fill = value)) +
    scale_fill_manual(values = c("A" = "#5F9E6E", "T" = "#E58686", "C" = "#6DAEDB", "G" = "#F6C85F")) +
    xlab("Tetranucleotide") +
    ylab("Position") +
    labs(fill = "") +
    theme +
    theme(legend.position = "right",  # moved the legend to the right
          legend.background = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(-1,0.5,0.5,0.5), "lines")) +
    scale_x_continuous(expand = c(0,0))
  
  # Combine the two figures
  figure <- fig_fold_original / fig_motif + plot_layout(heights = c(2, 1))
  
  # Save the combined plot
  ggsave(file.path(outdir, paste("fragment_end_contexts_change_different_colors_", cancer_type, "May 2024.pdf", sep = "")), figure, width = 10, height = 6, dpi = 500, units = "in")
}

### Plot regressions
for(cancer_type in unique_cancer_types){
  # Filter to only the current cancer_type in comp
  data_reg_filtered <- data_reg[grepl(paste0("\\b", cancer_type, "\\b"), data_reg$comp), ]
  
  fig_reg <- ggplot(data_reg_filtered, aes(var1, var2)) +
    geom_point(pch = 16, alpha = 0.5) +
    stat_regline_equation(label.y = 1.5, aes(label = ..rr.label..)) +
    facet_grid(.~comp) +
    xlab("Fold Change (Var1)") + 
    ylab("Fold Change (Var2)") +
    labs(color = "", fill = "") +
    #  scale_color_manual(name = " ", values =c("#1F78B4", "#33A02C", "#E31A1C")) +
    ggtitle(paste0("Fold Change in ", cancer_type, " vs Other Samples")) + 
    theme +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          strip.text.x = element_text(size = 8)) # for smaller text labels
  
  ggsave(file.path(outdir, paste0("fragment_end_contexts_change_reg, ", cancer_type, "Spring 2024.pdf")), fig_reg, width = 16, height = 8)
}






### Shannon plots
## Edit for classifier
theme <- theme(plot.title = element_text(hjust = 0.5, size = 13), 
               axis.line = element_line(colour = "black"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank(),
               legend.position = "bottom",
               legend.key = element_rect(fill = "white"),
               legend.title = element_text(size = 12),
               legend.text = element_text(size = 12),
               strip.background = element_blank(),
               strip.text = element_text(size = 13),
               axis.text = element_text(size = 13),
               axis.title = element_text(size = 13))

### Plot regressions
for(cancer_type in unique_cancer_types){
  # Filter to only the current cancer_type in comp
  stat_shannon <- final_shannon[final_shannon$cancer_type == cancer_type, ]
  
  ## Edit it to be clearer 
  stat_shannon$Var1
  
  plot_shan <- ggplot(stat_shannon) +
    geom_bar(aes(Var2, Freq, fill = Var2), color = NA, stat = "identity") +
    facet_grid(.~Var1, scales = "free_x", space = "free_x") +
    xlab("Change vs normal") + 
    ylab("# of Motifs") +
    ggtitle(paste0("Expected Motif Frequecy: ", cancer_type)) + 
    scale_fill_manual(name = " ", values = c("#5B84B1FF", "#FC766AFF")) +
    theme +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  plot_shan
  ggsave(file.path(outdir, paste0("fragment_end_shannon_expected, ", cancer_type, ".pdf")), plot_shan, width = 5, height = 8)
  
  stat_shannon2 <- final_shannon2[final_shannon2$cancer_type == cancer_type, ]
  
  plot_shan2 <- ggplot(stat_shannon2) +
    geom_bar(aes(Var3, Freq, fill = Var2), color = NA, stat = "identity") +
    facet_grid(Var2~Var1, scales = "free_x", space = "free_x") +
    xlab("# of A/T Nucleotides") + 
    ylab("# of Motifs") +
    ggtitle(paste0("Expected Motif Frequecy: ", cancer_type)) + 
    scale_fill_manual(name = " ", values = c("#5B84B1FF", "#FC766AFF")) +
    theme +
    theme(legend.position = "none",
          legend.background = element_blank())
  plot_shan2
  ggsave(file.path(outdir, paste0("fragment_end_shannon_AT, ", cancer_type, "F2023.pdf")), plot_shan2, width = 7, height = 7)
  
}

## Reformat the data for the classifier 

# Reformat the dataframe to wide format
tmp <- data_endmotif_filtered %>% dplyr::select(sample, motif, normalized_frequency)
data_motif_wide <- spread(tmp, key = motif, value = normalized_frequency)

row_names_to_keep <- data_motif_wide$sample

# Remove the 'sample' column
data_motif_wide <- data_motif_wide %>% dplyr::select(-sample)

# Reassign row names
rownames(data_motif_wide) <- row_names_to_keep

# Export the dataframe
saveRDS(data_motif_wide, file = "Dataframe matrix for end motifs, all frags, May 2024, validation.rds")
saveRDS(metadata_df, file = "Metadata df for end motifs, May 2024 SE.rds")

# Save other dataframes
saveRDS(data_endmotif_filtered, file = file.path(outdir, "filtered_endmotif_data_normalized_frequencies.rds"))
saveRDS(data_stats, file = file.path(outdir, "motif_statistics_by_cancer_type_SE.rds"))
saveRDS(data_stats_dnase, file = file.path(outdir, "DNASE1L3_motif_statistics_by_cancer_type_SE.rds"))
saveRDS(fold_change_all, file = file.path(outdir, "fold_change_all_cancer_vs_normal_SE.rds"))
saveRDS(data_reg, file = file.path(outdir, "regression_data_cancer_comparisons_SE.rds"))
saveRDS(final_shannon, file = file.path(outdir, "shannon_index_statistics_by_cancer_type_SE.rds"))
saveRDS(final_shannon2, file = file.path(outdir, "detailed_shannon_statistics_AT_content_by_cancer_type_SE.rds"))






### Now make heatmap and get z-scores

## Get cancer median across all cancers
normal_endmotif <- data_endmotif_filtered %>% dplyr::filter(cancer_type_corrected_updated == "healthy") %>% 
  dplyr::filter(!is.na(cancer_type_corrected_updated)) %>% 
  dplyr::select(motif, normalized_frequency, sample)

normal_endmotif_matrix <- pivot_wider(normal_endmotif, names_from = "sample", values_from = "normalized_frequency")

row_names_to_keep <- normal_endmotif_matrix$motif

# Remove the 'sample' column
normal_endmotif_matrix <- subset(normal_endmotif_matrix, select = -1)

# Reassign row names
rownames(normal_endmotif_matrix) <- row_names_to_keep

### Make healthy median
Normal_median <- rowMedians(as.matrix(normal_endmotif_matrix), na.rm = T)
Normal_sd <- rowSds(as.matrix(normal_endmotif_matrix), na.rm = T)


## Now get cancer median across all cancers
cancer_endmotif <- data_endmotif_filtered %>% dplyr::filter(cancer_type_corrected_updated != "healthy") %>% 
  dplyr::filter(!is.na(cancer_type_corrected_updated)) %>% 
  dplyr::select(motif, normalized_frequency, sample)

cancer_endmotif_matrix <- pivot_wider(cancer_endmotif, names_from = "sample", values_from = "normalized_frequency")

row_names_to_keep <- cancer_endmotif_matrix$motif

# Remove the 'sample' column
cancer_endmotif_matrix <- subset(cancer_endmotif_matrix, select = -1)

# Reassign row names
rownames(cancer_endmotif_matrix) <- row_names_to_keep

### Make cancer median
cancer_median <- rowMedians(as.matrix(cancer_endmotif_matrix))
cancer_sd <- rowSds(as.matrix(cancer_endmotif_matrix))

## Get distance from normal median 
Cancer_endmotif_ratio <- (cancer_endmotif_matrix - Normal_median)/Normal_sd
Normal_endmotif_ratio <- (normal_endmotif_matrix - Normal_median)/Normal_sd

### Calculate Normal Z-scores (Normal relative to Normal)
Normal_sum <- colSums(abs(Normal_endmotif_ratio), na.rm = TRUE)
Normal_sum_median <- median(Normal_sum, na.rm = TRUE)
Normal_MAD <- mad(Normal_sum, na.rm = TRUE) #Median Absolute Deviation

Normal_zscores <- abs((Normal_sum - Normal_sum_median))/Normal_MAD
Normal_mean <- mean(Normal_zscores)
z_limit <- quantile(Normal_zscores, 0.9)

### Calculate Cancer Z-scores (Cancer relative to Normal)
Cancer_sum <- colSums(abs(Cancer_endmotif_ratio), na.rm = TRUE)
Cancer_zscore <- (Cancer_sum - Normal_sum_median)/Normal_MAD

### Make Z-score matrix
all_zscores <- c(Normal_zscores,Cancer_zscore)
zscore_df <- as.data.frame(all_zscores)
zscore_df$limit <- z_limit
zscores_mat <- as.matrix(zscore_df)
zscores_mat_endmotif <- zscores_mat # save

### Set Normal median 
lower <- Normal_median - Normal_sd
upper <- Normal_median + Normal_sd
Normal_median <- as.matrix(cbind(Normal_median, lower, upper))

## Save the zscores to a seperate dataframe 
End_motif_df_zscores <- zscore_df 
End_motif_df_zscores <- rownames_to_column(End_motif_df_zscores, "Sample")

# Join with metadata df
End_motif_df_zscores <- left_join(End_motif_df_zscores, metadata_df, by = c("Sample" = "endmotif_name"))

### Save objects
saveRDS(Normal_ratio, file=file.path(outdir, paste0("Normal_ratio_Endmotif_SE.rds")))
saveRDS(Cancer_ratio, file=file.path(outdir, paste0("Cancer_ratio_Endmotif_SE.rds")))
saveRDS(Normal_median, file=file.path(outdir, paste0("Normal_median_Endmotif_SE.rds")))
saveRDS(zscores_mat_endmotif, file = file.path(outdir, paste0("Zscore_mat_Endmotif_SE.rds")))
saveRDS(End_motif_df_zscores, file = "Zscores df for end motifs_SE.rds")

# Now use this for the heatmap - although consider showing the end motif deconvolution signatures rather than heatmap


### Get significant proportions 
# Calculate the proportion of significant z-scores for each cancer type and z-score column
cutoff <- 2
significant_proportions_SE <- End_motif_df_zscores %>%
  dplyr::mutate(is_significant = abs(all_zscores) > cutoff) %>%
  dplyr::group_by(cancer_type_corrected_updated) %>%
  dplyr::summarise(proportion_significant = mean(is_significant, na.rm = TRUE)) %>%
  ungroup()

print(significant_proportions_SE)

# View the results
significant_proportions


## Make heatmap 

## Make a heatmap of the insert size matrix 

endmotif_matrix <-as.matrix(data_motif_wide)


# getting IDs
sample_order <- rownames(endmotif_matrix)
metadata_df <- metadata_df[match(sample_order, metadata_df$endmotif_name),]

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
                                      "Ovarian_cancer", 
                                      "Pancreatic_cancer", 
                                      "Colorectal_cancer", 
                                      "Bladder_cancer", 
                                      "Renal_cancer"),
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
                                      "Ovarian_cancer", 
                                      "Pancreatic_cancer", 
                                      "Colorectal_cancer", 
                                      "Bladder_cancer", 
                                      "Renal_cancer"))

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
col_heat <- colorRamp2(c(0, 0.3, 2), 
                       c("#1f78b4", "white", "#e31a1c"))
col_fun <- colorRamp2(c(-2, 0, 2), 
                      c("#1f78b4", "white", "#e31a1c"))
#col_sex <- c("Male" = "#1b9e77",
#             "Female" = "#d95f02",
#             "Not available" = "lightgrey")


color_func <- colorRampPalette(brewer.pal(12, "Set3"))
#cancer_colors <- color_func(15)
# Match Yong's 
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
sample_order <- rownames(endmotif_matrix)
column_order <- colnames(endmotif_matrix)

## Add some lables for endmotif regions
# Initialize a vector to hold the labels
endmotif_regions <- rep(NA, ncol(endmotif_matrix))

# Convert column names to numeric
column_nums <- colnames(endmotif_matrix)

# Iterate over the column names and assign labels
labels <- c("A***", "C***", "G***", "T***")
for (i in seq_along(column_nums)) {
  label_index <- ((i - 1) %/% 64) + 1
  endmotif_regions[i] <- labels[label_index]
}

# Convert to a factor with specified levels
#levels_order <- c("<167bp", "167bp", "168-333bp", "334bp", "335-500bp", "501bp", ">502bp")
levels_order_endmotif <- c("A***", "C***", "G***", "T***")
endmotif_regions <- factor(endmotif_regions, levels = levels_order_endmotif)

### Generate annotations
left_annotation <- rowAnnotation(# "Sex" = data_sex,
  "Type" = data_cancer_type,
  "Cohort" = data_cohort, # added cohort info
  show_annotation_name = TRUE,
  #border = TRUE,
  col = list("Sex" = col_sex, "Type" = cancer_type_col, "Cohort" = col_cohort),
  annotation_name_gp = gpar(fontsize = 12),
  annotation_name_side = "top",
  annotation_name_rot = 90,
  annotation_legend_param = list(title_gp = gpar(fontsize = 15), labels_gp = gpar(fontsize = 14)), # Font size for the legend
  #annotation_legend_param = list(border = TRUE),
  simple_anno_size = unit(0.35, "cm"),
  gap = unit(c(1, 1), "mm"))

# Reorder the rows of zscores_mat_endmotif based on sample_order
zscores_mat_endmotif <- zscores_mat_endmotif[sample_order, ]

right_annotation <- rowAnnotation("Genome-Wide\nZ-score" = anno_lines(zscores_mat_endmotif,
                                                                      add_points = TRUE,
                                                                      pch = c(16, NA),
                                                                      ylim = c(-2,20), #c(-1,7) for Butler, 
                                                                      #c(-1,40) for Cancer
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
#    annotation_legend_param = list(border = TRUE))

top_annotation <- HeatmapAnnotation("Healthy\nMedian" = anno_lines(Normal_median,
                                                                   ylim = c(0.02, 3),
                                                                   axis_param = list(at = c(0, 1.5, 3),
                                                                                     side = "right",
                                                                                     labels_rot = 0,
                                                                                     gp = gpar(fontsize = 10)),
                                                                   border = FALSE,
                                                                   gp = gpar(lty = c(rep("solid", length(Normal_median) - 1), "dotted"))),  # Set last line as dotted
                                    height = unit(1.2, "cm"),
                                    show_annotation_name = FALSE,
                                    annotation_name_rot = 0,
                                    annotation_name_gp = gpar(fontsize = 13))

# Ensure order of endmotif_matrix is correct 
endmotif_matrix <- endmotif_matrix[, column_order]

## Generate heatmap

Heatmap_endmotif_freq <- Heatmap(endmotif_matrix,
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
                                 column_split = endmotif_regions,
                                 column_title_gp = gpar(fontsize = 10),
                                 row_title_gp = gpar(fontsize = 12),
                                 column_title = "End motif proportions",
                                 row_split = row_split,
                                 cluster_row_slices = FALSE,
                                 row_title_rot = 0,
                                 border = FALSE,
                                 heatmap_legend_param = list(title = "Motif proportions", 
                                                             title_gp = gpar(fontsize = 15),  # Font size for the legend title
                                                             labels_gp = gpar(fontsize = 14)))  # Font size for the legend labels


ht_opt$COLUMN_ANNO_PADDING = unit(7, "mm")
#draw(Heatmap, heatmap_legend_side = "right", merge_legend = TRUE)

png(filename = "Heatmap of endmotif by cancer type, updated legends, updated colors, with breaks, May 2024, all frags 5' ends only, 20 to 600bp, SE.png", width = 20, height = 25, units = "in", res = 350)
draw(Heatmap_endmotif_freq, heatmap_legend_side = "right", merge_legend = TRUE)
dev.off()


# Order the zscores_mat according to sample_order
zscores_ordered_endmotif <- zscores_mat_endmotif[match(sample_order, rownames(zscores_mat_endmotif)), ]

# Make sure that the row names match the sample order
rownames(zscores_ordered_endmotif) <- sample_order

right_annotation <- rowAnnotation("Genome-Wide\nZ-score" = anno_lines(zscores_ordered_endmotif,
                                                                      add_points = TRUE,
                                                                      pch = c(16, NA),
                                                                      ylim = c(-2,20), 
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
#        annotation_legend_param = list(border = TRUE))
## Generate heatmap
row_split <- data_cancer_type

Heatmap_endmotif_freq <- Heatmap(endmotif_matrix,
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
                                 column_split = endmotif_regions,
                                 column_title_gp = gpar(fontsize = 10),
                                 row_title_gp = gpar(fontsize = 12),
                                 row_split = row_split, # for split 
                                 cluster_row_slices = FALSE,
                                 row_title_rot = 0,
                                 border = FALSE,
                                 heatmap_legend_param = list(title = "Motif proportions", 
                                                             title_gp = gpar(fontsize = 15),  # Font size for the legend title
                                                             labels_gp = gpar(fontsize = 14)))  # Font size for the legend labels


ht_opt$COLUMN_ANNO_PADDING = unit(6, "mm")
#draw(Heatmap, heatmap_legend_side = "right", merge_legend = TRUE)

png(filename = "Heatmap of end motif proportions by cancer type, updated legends, updated colors, with breaks, split, May 2024, updated 5' only all frags, 20-600bp.png", width = 18, height = 25, units = "in", res = 350)
draw(Heatmap_endmotif_freq, heatmap_legend_side = "right", merge_legend = TRUE)
dev.off()


## Zscore 

### Prepare for plotting
zscores_mat_endmotif_plot <- cbind(Normal_endmotif_ratio,Cancer_endmotif_ratio)
zscores_mat_endmotif_plot <- t(zscores_mat_endmotif_plot)

zscores_mat_endmotif_plot_ordered <- zscores_mat_endmotif_plot[match(sample_order, rownames(zscores_mat_endmotif_plot)), ]

# Change col
col_heat_zscore <- colorRamp2(c(-5, 0, 10), 
                              c("#1f78b4", "white", "#e31a1c"))

Heatmap_endmotif_freq_zscore <- Heatmap(zscores_mat_endmotif_plot_ordered,
                                        col = col_heat_zscore,
                                        show_row_names = FALSE,
                                        show_column_names = FALSE,
                                        #heatmap_legend_param = heatmap_legend_param,
                                        right_annotation = right_annotation,
                                        top_annotation = top_annotation,
                                        left_annotation = left_annotation,
                                        column_order = column_order,
                                        row_order = sample_order,
                                        row_labels = NULL,
                                        column_split = endmotif_regions,
                                        column_title_gp = gpar(fontsize = 10),
                                        row_title_gp = gpar(fontsize = 12),
                                        row_split = row_split, # for split 
                                        cluster_row_slices = FALSE,
                                        row_title_rot = 0,
                                        border = FALSE,
                                        heatmap_legend_param = list(title = "Zscore to normal", 
                                                                    title_gp = gpar(fontsize = 15),  # Font size for the legend title
                                                                    labels_gp = gpar(fontsize = 14)))  # Font size for the legend labels


ht_opt$COLUMN_ANNO_PADDING = unit(6, "mm")

png(filename = "Heatmap of end motif proportions by cancer type, updated legends, updated colors, with breaks, split, May 2024, updated 5' only all frags, 20-600bp, zscores, updated color.png", width = 18, height = 25, units = "in", res = 350)
draw(Heatmap_endmotif_freq_zscore, heatmap_legend_side = "right", merge_legend = TRUE)
dev.off()

png(filename = "Heatmap of end motif proportions by cancer type, updated legends, updated colors, with breaks, split, May 2024, updated 5' only all frags, 20-600bp, zscores, updated color, narrow 2.png", width = 18, height = 22, units = "in", res = 350)
draw(Heatmap_endmotif_freq_zscore, heatmap_legend_side = "right", merge_legend = TRUE)
dev.off()




## Export files for machine learning 

saveRDS(endmotif_matrix, file=file.path(outdir, paste0("End_motif_matrix_June2024_SE.rds")))
saveRDS(metadata_df, file=file.path(outdir, paste0("Metadata_df_endmotif_June2024_SE.rds")))


## Export the PCAs on the dataframe 

# Perform PCA
pca_result_endmotif <- prcomp(endmotif_matrix, center = TRUE, scale. = TRUE)

# Determine the number of components to retain
# Visualize the variance explained by each principal component
fviz_eig(pca_result_endmotif, addlabels = TRUE, ylim = c(0, 30))

# Extract the proportion of variance explained by each component
explained_variance <- summary(pca_result_endmotif)$importance[2,]

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
       subtitle = "End motifs. Red: 90%, Purple: 95%, Blue: 99%",
       x = "Number of Principal Components",
       y = "Cumulative Explained Variance") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10))

# Print the plot
print(p)
ggsave(filename = file.path(outdir, "Cumulative Explained Variance by Principal Components end motifs SE.png"), plot = p, width = 7, height = 5, units = "in", dpi = 300)


# Find the number of components that explain at least 97% of the variance
num_components <- which(cum_variance >= 0.97)[1]

# Extract the principal components
pca_df_endmotif <- as.data.frame(pca_result_endmotif$x[, 1:num_components])
rownames(pca_df_endmotif) <- rownames(endmotif_matrix)

# Create the variance explained plot with customized theme
variance_plot <- fviz_eig(pca_result_endmotif, addlabels = TRUE, ylim = c(0, 45)) +
  theme_classic() +
  ggtitle("Scree plot for end motifs") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

# Save the plot to a file
ggsave(filename = file.path(outdir, "PCA_Variance_Explained_Endmotif_SE.png"), plot = variance_plot, width = 10, height = 6, units = "in", dpi = 300)

# Save the PCA results to a file 
write.csv(pca_df_endmotif, file.path(outdir, "PCA_results_End_motif SE.csv"), row.names = TRUE)
saveRDS(pca_df_endmotif, file = "End_motifs_updated_June2024_PCAs_97_percent_variance_SE.rds")
saveRDS(pca_result_endmotif, file = "PCA results for end motifs_SE.rds")
























## Below here is testing - did not use for manuscript 

## Now see if can identify the most important motifs which are differentially present 
# Can comment on this in the stats part of the paper

# Perform Kruskal-Wallis test for each motif
results_motif_stats <- data_endmotif_filtered %>%
  group_by(motif) %>%
  do(tidy(kruskal.test(normalized_frequency ~ cancer_type_corrected_updated, data = .)))

# Get significant motifs
significant_motifs <- results_motif_stats %>% filter(p.value < 0.05) %>% pull(motif)

# Do p-value adjustment for multiple hypothesis correction
pairwise_adjusted_results <- list()
for (motif in significant_motifs) {
  data_subset <- data_endmotif_filtered %>% filter(motif == !!motif)
  
  # Perform pairwise Wilcoxon tests (as a substitute for post-hoc Kruskal-Wallis analysis)
  pairwise_results <- pairwise.wilcox.test(data_subset$normalized_frequency,
                                           data_subset$cancer_type_corrected_updated,
                                           p.adjust.method = "BH")
  
  # Store adjusted p-values and motif identifier for later use
  pairwise_adjusted_results[[motif]] <- pairwise_results$p.value
}

# Extract the output 

# Initialize an empty data frame to store the tidied results
pairwise_tidied <- data.frame(motif = character(),
                              group1 = character(),
                              group2 = character(),
                              p_value = numeric(),
                              p_adjusted = numeric(),
                              stringsAsFactors = FALSE)

# Loop through each motif in the list of pairwise adjusted results
for (motif in names(pairwise_adjusted_results)) {
  # Extract the matrix of adjusted p-values for the current motif
  p_matrix <- pairwise_adjusted_results[[motif]]
  
  # Convert the matrix to a tidy data frame
  tidy_df <- expand.grid(group1 = rownames(p_matrix),
                         group2 = colnames(p_matrix)) %>%
    mutate(p_value = as.vector(p_matrix),
           motif = motif)
  
  # Adjust p-values (if not already adjusted, depending on your prior steps)
  tidy_df$p_adjusted <- p.adjust(tidy_df$p_value, method = "BH")
  
  # Bind this motif's tidied and adjusted results to the collective data frame
  pairwise_tidied <- rbind(pairwise_tidied, tidy_df)
}

# Filter to retain only significant comparisons
significant_pairs <- pairwise_tidied %>% 
  filter(p_adjusted < 0.05) %>% select(group1, group2) %>% unique()

# Get only significant motifs - seems nearly all are - can try clustering everything together
significant_motifs <- pairwise_tidied %>% 
  filter(p_adjusted < 0.05) %>% select(motif) %>% unique()


## Heatmap with clustering - needs some work

# Convert categorical metadata into factors with specified levels and labels
data_sex <- factor(metadata_df$sex, levels = c("Male", "Female", ""), labels = c("Male", "Female", "Not available"))
data_cancer_type <- factor(metadata_df$cancer_type_title_case, 
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
                                      "Ovarian_cancer", 
                                      "Pancreatic_cancer", 
                                      "Colorectal_cancer", 
                                      "Bladder_cancer", 
                                      "Renal_cancer"),
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
                                      "Ovarian_cancer", 
                                      "Pancreatic_cancer", 
                                      "Colorectal_cancer", 
                                      "Bladder_cancer", 
                                      "Renal_cancer"))
# Prepare annotations
# Note: It's important to create these annotations before drawing the heatmap to ensure they align with the original sample order
left_annotation <- rowAnnotation(df = data.frame("Sex" = data_sex, "Cancer Type" = data_cancer_type),
                                 col = list("Sex" = col_sex, "Cancer Type" = cancer_type_col),
                                 width = unit(2, "cm"))

# Generate the heatmap with enabled clustering on columns (samples)
Heatmap_endmotif_freq <- Heatmap(endmotif_matrix,
                                 col = col_heat,
                                 name = "End Motif Frequency",
                                 #top_annotation = top_annotation,
                                 bottom_annotation = NULL,
                                 left_annotation = left_annotation,
                                 show_row_names = FALSE,
                                 show_column_names = FALSE,
                                 cluster_rows = FALSE, # Disable clustering of rows
                                 clustering_distance_rows = "euclidean",
                                 clustering_distance_columns = "euclidean",
                                 clustering_method_rows = "complete",
                                 clustering_method_columns = "complete",
                                 row_split = row_split, # for split 
                                 #        column_split = endmotif_regions,
                                 heatmap_legend_param = list(title = "Motif proportions", legend_direction = "horizontal"))

# Draw the heatmap
draw(Heatmap_endmotif_freq, heatmap_legend_side = "right", merge_legend = TRUE)

# Save to file
png(filename = "Clustered_Heatmap_Endmotif_Frequencies_with_split.png", width = 20, height = 25, units = "in", res = 300)
draw(Heatmap_endmotif_freq, heatmap_legend_side = "right", merge_legend = TRUE)
dev.off()

## Export table of the end motifs in a matrix 
write.table(endmotif_matrix, file = "End_motif_matrix_dataframe.txt", sep = "\t", row.names = TRUE, quote = F)







#### Export 
# Export the significant proportions for SE
write.csv(significant_proportions_SE, file = file.path(outdir, "significant_proportions_SE_Endmotif.csv"), row.names = FALSE)

# Export the new DNASE1L3 table (with fold change, p-values, and annotations)
write.csv(data_stats_dnase_filtered_SE, file = file.path(outdir, "DNASE1L3_statistics_SE_Endmotif.csv"), row.names = FALSE)

write.csv(data_dnase, file = file.path(outdir, "DNASE1L3_all_fold_changes_SE.csv"), row.names = FALSE)
write.csv(data_dnase_filtered, file = file.path(outdir, "DNASE1L3_all_fold_changes_SE_filtered_without_healthy_controls_that_are_problematic.csv"), row.names = FALSE)
