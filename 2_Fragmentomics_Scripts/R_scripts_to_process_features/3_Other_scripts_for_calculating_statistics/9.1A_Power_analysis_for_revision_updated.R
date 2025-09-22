# ----------------------------------------------------------------------------
# Title   : Combined Power Analysis for DMRs and Fragmentation Features
# Authors : Dory Abelman 
# Date    : March 2025
#
# Purpose :
#   • Compute Cohen’s d for differentially methylated regions (DMRs) in PE and SE
#     cfDNA samples and for six fragmentomic metrics (end‑motifs, insert‑size,
#     DELFI ratios, nucleosome‑peak distances, % short fragments, fragmentation score)
#     contrasting Cancer vs Normal.
#   • Perform two‑sample t‑test power analyses using the observed effect sizes
#     to estimate sample sizes needed for 80% power.
#   • Conduct one‑way ANOVA power analyses to assess ability to detect differences
#     among multiple cancer subtypes for each feature, including aggregate DMR metrics.
#   • Generate sensitivity curves (power vs. sample size) for all features and DMRs.
#   • Create and save: effect‑size distributions, sensitivity plots, ANOVA power tables,
#     pairwise t‑test power tables, and combined power summaries for manuscript figures
#     and supplementary materials.
# ----------------------------------------------------------------------------


# -------------------------------
# Load Required Libraries
# -------------------------------

# ============================
# LOAD LIBRARIES & DATASETS
# ============================
library(tidyverse)
library(pwr)
library(effsize)  # for cohen.d
library(pwr)
library(vegan)      # for adonis2
library(mvnormtest)



# --- Load Metadata ---
metadata_df_tmp <- readRDS("Metadata_df_all_with_corrected_cancer_subtype_Jan2025.rds")
metadata_df_tmp <- metadata_df_tmp %>% 
  filter(!grepl("lfs", cancer_type_corrected, ignore.case = TRUE)) %>%
  mutate(group = ifelse(tolower(CN_classifier) == "healthy", "Normal", "Cancer"))
metadata_df_tmp$group <- factor(metadata_df_tmp$group, levels = c("Normal", "Cancer"))

# --- Load Raw Features for Fragmentation ---
data_dir <- "Final Data Dec 2024/Raw_features/"

matrix_nucleosome_peak <- readRDS(file.path(data_dir, "Matrix for nucleosome peak.rds"))
matrix_nucleosome_peak_validation <- readRDS(file.path(data_dir, "Matrix for nucleosome peak validation.rds"))

delfi_ratios <- readRDS(file.path(data_dir, "DELFI_Ratios_updated_Dec2024.rds"))
delfi_ratios_validation <- readRDS(file.path(data_dir, "DELFI_Ratios_updated_Dec2024_validation_cohort.rds"))

insert_matrix <- readRDS(file.path(data_dir, "insert_matrix_df_for_Althaf_updated.rds"))
insert_matrix_validation <- readRDS(file.path(data_dir, "insert_matrix_df_validation_for_Althaf_updated.rds"))

end_motifs_all_frags <- readRDS(file.path(data_dir, "Dataframe matrix for end motifs, all frags, May 2024.rds"))
end_motifs_validation <- readRDS(file.path(data_dir, "Dataframe matrix for end motifs, all frags, May 2024, validation.rds"))

# --- Load DMR Data ---
# PE DMR data
norm_counts <- readRDS("Final Data Dec 2024/Methylation_data/Updated_Feb3/PanCancer_hyper_DMRs_combat_deseq2_normalized_count_with_blood_cell_age_sex_filtering_for_PE_samples_after_QC.RDS")
common_samples <- intersect(colnames(norm_counts), metadata_df_tmp$sequencing_id)
norm_counts <- norm_counts[, common_samples]
dmr_metadata <- metadata_df_tmp %>% filter(sequencing_id %in% common_samples)

# SE DMR data
norm_counts_se <- readRDS("Final Data Dec 2024/Methylation_data/Updated_Feb3/PanCancer_hyper_DMRs_combat_deseq2_normalized_count_with_blood_cell_age_sex_filtering_for_SE_samples_after_QC.RDS")
common_samples_se <- intersect(colnames(norm_counts_se), metadata_df_tmp$sequencing_id)
norm_counts_se <- norm_counts_se[, common_samples_se]
dmr_metadata_se <- metadata_df_tmp %>% filter(sequencing_id %in% common_samples_se)


# ============================
# FUNCTIONS
# ============================
# Function to compute Cohen's d for each feature (column) in a matrix/dataframe
compute_effects <- function(feature_matrix, sample_ids, metadata, group_col = "group") {
  # Subset feature matrix to only those samples present in metadata
  common_ids <- intersect(colnames(feature_matrix), metadata$sequencing_id)
  submat <- feature_matrix[, common_ids, drop = FALSE]
  meta_sub <- metadata %>% filter(sequencing_id %in% common_ids)
  
  effects <- apply(submat, 1, function(x) {
    tryCatch({
      cohen.d(x, meta_sub[[group_col]], na.rm = TRUE)$estimate
    }, error = function(e) NA)
  })
  effects <- effects[!is.na(effects)]
  return(effects)
}

# Function to simulate power for a given effect size over a range of sample sizes
simulate_power <- function(effect_size, sample_sizes, sig.level = 0.05) {
  sapply(sample_sizes, function(n) {
    pwr.t.test(n = n, d = effect_size, sig.level = sig.level,
               type = "two.sample", alternative = "two.sided")$power
  })
}


# ============================
# SECTION 1: DMR (Methylation) ANALYSIS
# ============================

# --- 1.2 Calculate Effect Sizes for DMRs (PE) ---
dmr_effects <- apply(norm_counts, 1, function(x) {
  tryCatch({
    cohen.d(x, dmr_metadata$group, na.rm = TRUE)$estimate
  }, error = function(e) NA)
})
dmr_effects <- dmr_effects[!is.na(dmr_effects)]
median_dmr_d <- median(abs(dmr_effects))
mean_dmr_d <- mean(abs(dmr_effects))
cat("PE DMRs - Median Cohen's d:", median_dmr_d, "\n")
cat("PE DMRs - Mean Cohen's d:", mean_dmr_d, "\n")

# --- 1.3 Sensitivity Analysis for PE DMRs ---
n_cancer <- sum(df_power$group == "Cancer")  # df_power should be defined later from fragmentation features metadata
# If df_power is not defined yet, we can use metadata_df_tmp (assuming it covers all samples)
n_cancer <- sum(metadata_df_tmp$group == "Cancer")
max_sample_size <- n_cancer  # use available number per group
dmr_sample_sizes <- seq(10, max_sample_size, by = 10)
dmr_power_values <- simulate_power(mean_dmr_d, dmr_sample_sizes)

sensitivity_dmr <- data.frame(Data_Type = "DMR",
                              Sample_Size = dmr_sample_sizes,
                              Power = dmr_power_values,
                              Feature = "PE DMR (mean)",
                              Analysis = "DMR",
                              stringsAsFactors = FALSE)

# --- 1.3b Sensitivity Analysis for PE DMRs using median effect size ---
dmr_power_values_median <- simulate_power(median_dmr_d, dmr_sample_sizes)

sensitivity_dmr_median <- data.frame(Data_Type = "DMR",
                                     Sample_Size = dmr_sample_sizes,
                                     Power = dmr_power_values_median,
                                     Feature = "PE DMR (median)",
                                     Analysis = "DMR",
                                     stringsAsFactors = FALSE)

# --- 1.4 Plot Distribution of DMR Effect Sizes (PE) ---
df_dmr_effects <- data.frame(Effect_Size = abs(dmr_effects))
p_dmr_effects <- ggplot(df_dmr_effects, aes(x = Effect_Size)) +
  geom_histogram(bins = 50, fill = "skyblue", color = "black") +
  labs(title = "Distribution of Absolute Cohen's d for PE DMRs",
       x = "Absolute Cohen's d", y = "Count") +
  theme_classic(base_size = 14)
ggsave("DMR_Effect_Size_Distribution_PE.pdf", p_dmr_effects, width = 8, height = 5)
print(p_dmr_effects)


# --- 1B. DMR Analysis for SE Data ---
dmr_effects_se <- apply(norm_counts_se, 1, function(x) {
  tryCatch({
    cohen.d(x, dmr_metadata_se$group, na.rm = TRUE)$estimate
  }, error = function(e) NA)
})
dmr_effects_se <- dmr_effects_se[!is.na(dmr_effects_se)]
median_dmr_d_se <- median(abs(dmr_effects_se))
mean_dmr_d_se <- mean(abs(dmr_effects_se))
cat("SE DMRs - Median Cohen's d:", median_dmr_d_se, "\n")
cat("SE DMRs - Mean Cohen's d:", mean_dmr_d_se, "\n")

dmr_sample_sizes_se <- seq(10, max_sample_size, by = 10)
dmr_power_values_se <- simulate_power(median_dmr_d_se, dmr_sample_sizes_se)
sensitivity_dmr_se <- data.frame(Data_Type = "SE DMR",
                                 Sample_Size = dmr_sample_sizes_se,
                                 Power = dmr_power_values_se,
                                 Feature = "SE DMR (median)",
                                 Analysis = "SE DMR",
                                 stringsAsFactors = FALSE)
dmr_power_values_se_mean <- simulate_power(mean_dmr_d_se, dmr_sample_sizes_se)
sensitivity_dmr_se_mean <- data.frame(Data_Type = "SE DMR",
                                      Sample_Size = dmr_sample_sizes_se,
                                      Power = dmr_power_values_se_mean,
                                      Feature = "SE DMR (mean)",
                                      Analysis = "SE DMR",
                                      stringsAsFactors = FALSE)


# ============================
# SECTION 2: FRAGMENTATION FEATURES ANALYSIS (RAW MATRICES)
# ============================

# Define a list of fragmentation feature matrices and labels.
frag_list <- list(
  "Nucleosome_peak" = matrix_nucleosome_peak,
  "DELFI_Ratios" = delfi_ratios,
  "Insert_Size" = insert_matrix,
  "End_Motifs" = end_motifs_all_frags
)

# For each fragmentation dataset, compute Cohen's d per feature.
frag_effects_list <- list()
frag_sensitivity_list <- list()

for(feature_name in names(frag_list)) {
  mat <- frag_list[[feature_name]]
  
  # Clean rownames: remove "_dedup" and following characters.
  rownames(mat) <- sub("_dedup.*", "", rownames(mat))
  
  # Rows are samples; get common sample IDs.
  common_ids <- intersect(rownames(mat), metadata_df_tmp$sequencing_id)
  submat <- mat[common_ids, , drop = FALSE]
  meta_sub <- metadata_df_tmp %>% filter(sequencing_id %in% common_ids)
  
  # Compute Cohen's d for each column (feature) in the matrix.
  effect_sizes <- apply(submat, 2, function(x) {
    tryCatch({
      cohen.d(x, meta_sub$group, na.rm = TRUE)$estimate
    }, error = function(e) NA)
  })
  effect_sizes <- effect_sizes[!is.na(effect_sizes)]
  
  # Store individual effect sizes.
  df_effect <- data.frame(Feature = paste(feature_name, names(effect_sizes), sep = "_"),
                          cohen_d = effect_sizes,
                          Data_Type = feature_name,
                          stringsAsFactors = FALSE)
  frag_effects_list[[feature_name]] <- df_effect
  
  # Simulate power for each feature individually, then average.
  sample_sizes_frag <- seq(5, max(nrow(meta_sub[meta_sub$group=="Cancer", ]), na.rm = TRUE), by = 5)
  power_mat <- sapply(effect_sizes, function(d_val) {
    simulate_power(abs(d_val), sample_sizes_frag)
  })
  avg_power <- rowMeans(power_mat, na.rm = TRUE)
  
  df_sens <- data.frame(Data_Type = feature_name,
                        Sample_Size = sample_sizes_frag,
                        Power = avg_power,
                        Feature = paste(feature_name, "(avg)", sep = " "),
                        Analysis = "Fragmentation",
                        stringsAsFactors = FALSE)
  frag_sensitivity_list[[feature_name]] <- df_sens
}

# --- Aggregate Fragmentation Features ---
# For each fragmentation data type, compute the median effect size (aggregated metric).
frag_effects_agg <- bind_rows(frag_effects_list) %>%
  group_by(Data_Type) %>%
  summarize(agg_cohen_d = median(cohen_d, na.rm = TRUE), .groups = "drop") %>%
  mutate(Feature = paste(Data_Type, "(median)"),
         Analysis = "Fragmentation")

# Simulate power curves for aggregated fragmentation features.
frag_sample_sizes <- seq(5, max(nrow(metadata_df_tmp[metadata_df_tmp$group=="Cancer", ]), na.rm = TRUE), by = 5)
frag_sensitivity_agg <- frag_effects_agg %>%
  rowwise() %>%
  mutate(Power_list = list(simulate_power(abs(agg_cohen_d), frag_sample_sizes)),
         Sample_Size_list = list(frag_sample_sizes)) %>%
  ungroup() %>%
  unnest(cols = c(Power_list, Sample_Size_list)) %>%
  rename(Power = Power_list, Sample_Size = Sample_Size_list) %>%
  select(Data_Type, Sample_Size, Power, Feature, Analysis)


###############################
# SECTION 3: COMBINE SENSITIVITY RESULTS
###############################
# Use median values for all features.
sensitivity_dmr_all <- bind_rows(sensitivity_dmr_median, sensitivity_dmr_se)
sensitivity_all <- bind_rows(sensitivity_dmr_all, frag_sensitivity_agg)

# Clean up feature names
sensitivity_all <- sensitivity_all %>%
  mutate(Feature = recode(Feature,
                          "DELFI_Ratios (median)" = "Fragment Ratios (median)",
                          "End_Motifs (median)"   = "End Motifs (median)",
                          "Nucleosome_peak (median)" = "Nucleosome Peak (median)",
                          "Insert_Size (median)" = "Insert Size (median)"
                          
  ))

# Define your custom color palette again (for safety)
color_palette <- c(
  "End Motifs (median)"                              = "#E69F00",
  "Insert Size (median)"                             = "#56B4E9",
  "PE DMR (median)"                             = "#CC79A7",
  "SE DMR (median)" = "#9932CC",
  "Nucleosome Peak (median)"                         = "#0072B2",
  "Fragment Ratios (median)"                          = "#009E73"
)

# Plot combined sensitivity curves.
p_combined_sens <- ggplot(sensitivity_all, aes(x = Sample_Size, y = Power, color = Feature)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(title = "Sensitivity Analysis",
       subtitle = "DMRs (PE & SE) and Aggregated Fragmentation Features",
       x = "Sample Size per Group", y = "Statistical Power",
       color = "Feature", linetype = "Data Category") +
  scale_color_manual(values = color_palette) +
  theme_classic(base_size = 13)
ggsave("Combined_Sensitivity_Power_All_Median_updated.pdf", p_combined_sens, width = 7, height = 5)
print(p_combined_sens)

###############################
# SECTION 4: COMBINED EFFECT SIZE BAR PLOT
###############################
# Create combined effect size plot using aggregated (median) metrics.
dmr_effects_PE_df <- data.frame(Feature = "PE DMR (median)", cohen_d = median_dmr_d, Data_Type = "Methylation", stringsAsFactors = FALSE)
dmr_effects_SE_df <- data.frame(Feature = "SE DMR (median)", cohen_d = median_dmr_d_se, Data_Type = "Methylation", stringsAsFactors = FALSE)
combined_effects_all <- bind_rows(frag_effects_agg %>% select(Feature, cohen_d = agg_cohen_d, Data_Type),
                                  dmr_effects_PE_df, dmr_effects_SE_df)

# Rename features for plotting
# Clean up feature names to match the sensitivity plot
combined_effects_all <- combined_effects_all %>%
  mutate(Feature = recode(Feature,
                          "DELFI_Ratios (median)" = "Fragment Ratios (median)",
                          "End_Motifs (median)"   = "End Motifs (median)",
                          "Nucleosome_peak (median)" = "Nucleosome Peak (median)",
                          "Insert_Size (median)"  = "Insert Size (median)"
  ))

p_combined_effects <- ggplot(combined_effects_all, aes(x = reorder(Feature, cohen_d), y = cohen_d, fill = Feature)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = color_palette) +
  labs(title = "Median Effect Sizes",
       x = "Feature", y = "Cohen's d") +
  theme_classic() +
  theme(legend.position = "none")

ggsave("Combined_Effect_Sizes_median.pdf", p_combined_effects, width = 5, height = 5)
print(p_combined_effects)


# ============================
# OPTIONAL: Save Summary Tables
# ============================
# Save power simulation summaries for DMR and fragmentation features if needed
write.csv(sensitivity_all, "Combined_Sensitivity_Power_All.csv", row.names = FALSE)
write.csv(combined_effects_all, "Combined_Effect_Sizes_All.csv", row.names = FALSE)

# ============================
# Manuscript Reporting Summary (Console Output)
# ============================
cat("\n--- Manuscript Reporting Summary ---\n")
cat("For PE DMRs, the mean absolute Cohen's d is", round(mean_dmr_d, 3), "and for SE DMRs it is", round(mean_dmr_d_se, 3), ".\n")
# For power, we can report required sample size for 80% power for a chosen feature group:
required_PE <- dmr_sample_sizes[which(dmr_power_values >= 0.8)[1]]
required_SE <- dmr_sample_sizes_se[which(dmr_power_values_se >= 0.8)[1]]
cat("Approximately", required_PE, "samples per group are required for 80% power for PE DMR detection,\n")
cat("while for SE DMRs (using median effect size) it is", required_SE, "samples per group.\n")


#### Now redo for the kept principal components as aggregate 

# ================================
# SET PATHS AND FILE NAMES
# ================================
pca_dir <- "/Users/dabelman/Documents/Thesis_work/R/M4/Projects/Fragmentation_CfMeDIP_Yong/PCA_results/Run_PCA_on_train_fit_on_validation_March2025/"
file_names <- c(
  "Combined_EndMotif_Methylation_1pct_PCs_training.rds",
  "Combined_ratios_Methylation_1pct_PCs_training.rds",
  "delfi_ratios_Training_PCA_80percent.rds",
  "end_motifs_Training_PCA_97percent.rds",
  "insert_matrix_Training_PCA_97percent.rds",
  "methylation_Training_PCA_90percent.rds",
  "nucleosome_peak_Training_PCA_80percent.rds",
  "Training_combined_1pct_PCs.rds",
  "Training_combined_5pct_PCs.rds",
  "Combined_EndMotif_Methylation_FragmentRatios_1pct_PCs_training.rds",
  "Combined_EndMotif_Methylation_FragmentRatios_1pct_PCs_training.rds"
)

# ================================
# FUNCTION: PERMANOVA POWER SIMULATION
# ================================
simulate_permanova_power <- function(pc_data, pc_columns = NULL, sample_sizes = seq(10, 500, by = 10),
                                     n_reps = 100, dist_method = "euclidean", sig.level = 0.05) {
  # If PC columns not provided, assume all columns except 'group' and 'sample' are PCs.
  if (is.null(pc_columns)) {
    pc_columns <- setdiff(colnames(pc_data), c("group", "sample"))
  }
  
  # Data matrix and groups.
  data_matrix <- as.matrix(pc_data[, pc_columns])
  groups <- pc_data$group
  
  # Determine the minimum number of samples available in any group.
  min_group_size <- min(table(groups))
  allowed_sample_sizes <- sample_sizes[sample_sizes <= min_group_size]
  
  if(length(allowed_sample_sizes) == 0) {
    stop("No sample size in the provided range is less than or equal to the smallest group size.")
  }
  
  power_results <- data.frame(Sample_Size = allowed_sample_sizes, Power = NA)
  
  for (n in allowed_sample_sizes) {
    rep_results <- replicate(n_reps, {
      # For each group, sample n indices (without replacement)
      group_levels <- unique(groups)
      idxs <- unlist(lapply(group_levels, function(g) {
        sample(which(groups == g), n, replace = FALSE)
      }))
      
      # Subset the data matrix and groups.
      data_sub <- data_matrix[idxs, ]
      group_sub <- groups[idxs]
      
      # Run PERMANOVA (adonis2) on the subset.
      adonis_res <- adonis2(data_sub ~ group_sub, permutations = 999, method = dist_method)
      p_val <- adonis_res$`Pr(>F)`[1]
      
      return(p_val < sig.level)
    })
    
    power_results$Power[power_results$Sample_Size == n] <- mean(rep_results)
  }
  
  return(power_results)
}

# ================================
# FUNCTION: CHECK MULTIVARIATE NORMALITY
# ================================
check_mv_normality <- function(pc_data, pc_columns) {
  groups <- unique(pc_data$group)
  normality_results <- list()
  for(g in groups) {
    group_data <- pc_data %>% filter(group == g)
    test_mat <- t(as.matrix(group_data[, pc_columns]))
    test_result <- tryCatch({
      mvnormtest::mshapiro.test(test_mat)
    }, error = function(e) {
      warning("Error in mshapiro.test for group ", g, ": ", e$message)
      return(NULL)
    })
    normality_results[[g]] <- test_result
    cat("Multivariate Normality for group:", g, "\n")
    print(test_result)
  }
  return(normality_results)
}

# ================================
# INITIALIZE LISTS FOR RESULTS
# ================================
permanova_results_list <- list()
power_simulation_list <- list()

# ================================
# LOOP THROUGH EACH PCA FILE
# ================================
for(file in file_names) {
  file_path <- file.path(pca_dir, file)
  pc_data <- readRDS(file_path)
  
  # If "group" column doesn't exist, create it by cleaning rownames.
  if(!"group" %in% colnames(pc_data)) {
    pc_data <- pc_data %>%
      mutate(sample = sub("_dedup.*", "", rownames(pc_data)))
    
    # Join with metadata to get group assignment.
    pc_data <- left_join(pc_data, metadata_df_tmp %>% select(sequencing_id, group),
                         by = c("sample" = "sequencing_id"))
    
    if(any(is.na(pc_data$group))) {
      warning("Some samples in file ", file, " could not be matched to metadata for group assignment.")
    }
  }
  
  # Assume every column except "group" and "sample" is a PC.
  pc_cols <- setdiff(colnames(pc_data), c("group", "sample"))
  
  # Check multivariate normality (for diagnostic purposes).
  cat("Checking multivariate normality for file:", file, "\n")
  normality_results <- check_mv_normality(pc_data, pc_cols)
  
  # Run PERMANOVA-based power simulation.
  power_results <- simulate_permanova_power(pc_data, pc_columns = pc_cols,
                                            sample_sizes = seq(10, 500, by = 10),
                                            n_reps = 100, dist_method = "euclidean", sig.level = 0.05)
  
  permanova_results_list[[file]] <- data.frame(
    File = file,
    # Here you could store summary statistics from adonis2 if desired.
    # For example, the R^2 value:
    R2 = NA,  # placeholder; you could compute an overall R2 from the full dataset if needed.
    stringsAsFactors = FALSE
  )
  
  power_simulation_list[[file]] <- data.frame(
    File = file,
    Sample_Size = power_results$Sample_Size,
    Power = power_results$Power,
    stringsAsFactors = FALSE
  )
}

# ================================
# AGGREGATE AND PLOT RESULTS
# ================================
permanova_summary_df <- bind_rows(permanova_results_list)
power_simulation_df <- bind_rows(power_simulation_list)

p_power <- ggplot(power_simulation_df, aes(x = Sample_Size, y = Power, color = File)) +
  geom_line(size = 1) +
  geom_point() +
  labs(title = "PERMANOVA Power Curves for PCA-Transformed Data",
       subtitle = "Using permutation-based MANOVA (adonis2)",
       x = "Sample Size per Group",
       y = "Power") +
  theme_classic(base_size = 14)

ggsave("PERMANOVA_PCA_Power_Curves_All.pdf", p_power, width = 10, height = 6)
print(p_power)

# Print PERMANOVA summary table (if any additional summary statistics are desired)
print(permanova_summary_df)


## Helper function to map raw metric names to pretty labels
## Define the mapping from raw file names to pretty labels
pretty_mapping <- c(
  "Combined_EndMotif_Methylation_1pct_PCs_training.rds" = "Combined Motif + Methylation (1%)",
  "Combined_ratios_Methylation_1pct_PCs_training.rds" = "Combined Ratios + Methylation (1%)",
  "delfi_ratios_Training_PCA_80percent.rds" = "Fragment Ratio",
  "end_motifs_Training_PCA_97percent.rds" = "End Motifs",
  "insert_matrix_Training_PCA_97percent.rds" = "Insert Size",
  "methylation_Training_PCA_90percent.rds" = "Methylation",
  "nucleosome_peak_Training_PCA_80percent.rds" = "Nucleosome Peak",
  "Training_combined_1pct_PCs.rds" = "Combined All Features (1%)",
  "Training_combined_5pct_PCs.rds" = "Combined All Features (5%)",
  "Combined_EndMotif_Methylation_FragmentRatios_1pct_PCs_training.rds" = "Combined Motif + Methylation + Ratios (1%)"
)

## Define your custom color palette (names must match pretty labels exactly)
my_color_palette <- c(
  "Combined All Features (1%)"              = "#6A3D9A",
  "Combined All Features (5%)"              = "#8B4513",
  "Combined Motif + Methylation (1%)"         = "#32CD32",
  "Combined Motif + Methylation + Ratios (1%)" = "#FFD700",
  "Combined Ratios + Methylation (1%)"        = "#FF4500",
  "End Motifs"                              = "#E69F00",
  "Insert Size"                             = "#56B4E9",
  "Methylation"                             = "#CC79A7",
  "Nucleosome Peak"                         = "#0072B2",
  "Fragment Ratio"                          = "#009E73"
)

## Map the raw file names (from the File column) to pretty labels.
## Here we use basename() to ensure we match the file name.
power_simulation_df <- power_simulation_df %>%
  mutate(Pretty_Label = pretty_mapping[basename(File)])

## Now create the power curve plot using the pretty labels and the custom color palette.
p_power <- ggplot(power_simulation_df, aes(x = Sample_Size, y = Power, color = Pretty_Label)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "gray30") +
  labs(title = "PERMANOVA Power Curves for PCA-Transformed Data",
       subtitle = "Using permutation-based MANOVA",
       x = "Sample Size per Group",
       y = "Power",
       color = "Metric") +
  scale_color_manual(values = my_color_palette) +
  theme_classic(base_size = 14)

ggsave("PERMANOVA_PCA_Power_Curves_All_Pretty.pdf", p_power, width = 10, height = 6)
print(p_power)

# Export PERMANOVA summary results and power simulation data
write.csv(permanova_summary_df, "PERMANOVA_PCA_Manova_Summary.csv", row.names = FALSE)
write.csv(power_simulation_df, "PERMANOVA_PCA_Power_Simulation.csv", row.names = FALSE)

# Export session information for reproducibility
writeLines(capture.output(sessionInfo()), "session_info.txt")




























#------------------------------
# Extended Sensitivity Analysis for DMR
#------------------------------
# (Assuming mean_dmr_d has been computed as the mean absolute Cohen's d for DMRs)
# We'll simulate power estimates up to 1,000 samples per group.
dmr_sample_sizes_extended <- seq(5, 1000, by = 10)

dmr_power_values_extended <- sapply(dmr_sample_sizes_extended, function(n) {
  pwr.t.test(n = n, d = mean_dmr_d, sig.level = 0.05, 
             type = "two.sample", alternative = "two.sided")$power
})

sensitivity_dmr_extended <- data.frame(
  Data_Type = "DMR",
  Sample_Size = dmr_sample_sizes_extended,
  Power = dmr_power_values_extended,
  Feature = "DMR (mean)",
  Analysis = "DMR",
  stringsAsFactors = FALSE
)

#------------------------------
# Extended Sensitivity Analysis for Fragmentation Features
#------------------------------
# (Assuming frag_effects is a data frame with columns "feature" and "cohen_d")
features_frag <- c("Endmotif_zscore", "Insertsize_zscore", "Delfi_zscore", 
                   "Nucleosome_peak_zscore", "Twentyto150_pct", "FS")
frag_sample_sizes_extended <- seq(5, 1000, by = 5)

sensitivity_frag_extended <- expand.grid(Data_Type = features_frag, Sample_Size = frag_sample_sizes_extended, stringsAsFactors = FALSE)
sensitivity_frag_extended$Power <- NA

for(feat in features_frag) {
  d_val <- frag_effects$cohen_d[frag_effects$feature == feat]
  for(n in frag_sample_sizes_extended) {
    sensitivity_frag_extended$Power[sensitivity_frag_extended$Data_Type == feat & sensitivity_frag_extended$Sample_Size == n] <-
      pwr.t.test(n = n, d = d_val, sig.level = 0.05, type = "two.sample", alternative = "two.sided")$power
  }
}

sensitivity_frag_extended <- sensitivity_frag_extended %>% 
  mutate(Feature = paste("Fragment:", Data_Type),
         Analysis = "Fragment")

#------------------------------
# Combine Extended Sensitivity Data
#------------------------------
sensitivity_dmr_extended <- sensitivity_dmr_extended %>% select(Feature, Sample_Size, Power, Analysis)
sensitivity_frag_extended <- sensitivity_frag_extended %>% select(Feature, Sample_Size, Power, Analysis)
combined_sensitivity_extended <- bind_rows(sensitivity_dmr_extended, sensitivity_frag_extended)

# Define the simulation maximum sample size (we simulated up to 1000)
max_simulated_sample <- 365

#------------------------------
# Summarize Required Sample Size for 80% Power per Feature
#------------------------------
power_summary <- combined_sensitivity_extended %>%
  group_by(Feature) %>%
  summarize(
    min_n_for_80 = if (any(Power >= 0.8)) min(Sample_Size[Power >= 0.8]) else NA_real_,
    power_at_max = Power[Sample_Size == max_simulated_sample][1]
  ) %>%
  mutate(
    comment = ifelse(is.na(min_n_for_80),
                     paste("80% power not reached even at", max_simulated_sample, "samples (max power =", round(power_at_max, 3), ")"),
                     ifelse(min_n_for_80 > 356,
                            paste("80% power reached at", min_n_for_80, "samples, which is above the current maximum of 356"),
                            paste("80% power reached at", min_n_for_80, "samples")
                     )
    )
  )

print(power_summary)



































#### Below here is old testing 

# -------------------------------
# Data Preparation
# -------------------------------
# We assume that cancer_type_corrected_updated == "healthy" means normal.
df_power <- Merged_zscore_df %>%
  mutate(group = ifelse(cancer_type_corrected_updated == "healthy", "Normal", "Cancer"))

## Remove the LFS from this since used in validation cohort 
df_power <- df_power %>% filter(cancer_type_corrected_updated != "lfs_previvor")
df_power <- df_power %>% filter(cancer_type_corrected_updated != "lfs_survivor")
df_power <- df_power %>% filter(cancer_type_corrected_updated != "lfs_positive")

# Calculate sample sizes for overall Cancer vs. Normal comparison:
n_normal <- sum(df_power$group == "Normal")
n_cancer <- sum(df_power$group == "Cancer")
cat("Sample sizes: Normal =", n_normal, "Cancer =", n_cancer, "\n")

## Remove the LFS from the model 


# -------------------------------
# 1. Calculate Effect Sizes (Cohen's d) and T-test Power Analysis
# -------------------------------
features <- c("Endmotif_zscore", "Insertsize_zscore", "Delfi_zscore", 
              "Nucleosome_peak_zscore", "Twentyto150_pct", "FS")

effect_sizes <- data.frame(feature = character(), cohen_d = numeric(), stringsAsFactors = FALSE)
power_results_ttest <- data.frame(feature = character(), power = numeric(), 
                                  required_n_equal = numeric(), stringsAsFactors = FALSE)

for (feat in features) {
  # Compute Cohen's d comparing the two groups (Normal vs Cancer)
  d_result <- cohen.d(df_power[[feat]], df_power$group, na.rm = TRUE)
  effect_sizes <- rbind(effect_sizes, data.frame(feature = feat, cohen_d = d_result$estimate))
  
  # Compute achieved power using current (unequal) sample sizes
  power_test <- pwr.t2n.test(n1 = n_normal, n2 = n_cancer, d = d_result$estimate, sig.level = 0.05)
  
  # Approximate required sample size per group (using equal groups)
  req <- pwr.t.test(d = d_result$estimate, power = 0.8, sig.level = 0.05, type = "two.sample", alternative = "two.sided")
  
  power_results_ttest <- rbind(power_results_ttest, data.frame(feature = feat,
                                                               power = power_test$power,
                                                               required_n_equal = ceiling(req$n)))
}

cat("Two-sample t-test power analysis results:\n")
print(power_results_ttest)

# Export the tables as CSV files
write.csv(effect_sizes, "Effect_Sizes_Cancer_vs_Normal.csv", row.names = FALSE)
write.csv(power_results_ttest, "TTest_Power_Analysis_Cancer_vs_Normal.csv", row.names = FALSE)

# Plot the Cohen's d effect sizes using theme_classic
p_effect_sizes <- ggplot(effect_sizes, aes(x = reorder(feature, cohen_d), y = cohen_d)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +
  labs(title = "Effect Sizes (Cohen's d) for Each Feature\nCancer vs. Normal", 
       x = "Feature", y = "Cohen's d") +
  theme_classic(base_size = 14)
ggsave("Effect_Sizes_Cancer_vs_Normal.pdf", p_effect_sizes, width = 8, height = 5)

# -------------------------------
# 2. ANOVA Power Analysis for Differences Among Cancer Types
# -------------------------------
anova_power_results <- data.frame(feature = character(), f_effect = numeric(), 
                                  eta_sq = numeric(), power = numeric(), 
                                  avg_n = numeric(), k = numeric(), stringsAsFactors = FALSE)

for (feat in features) {
  # Run one-way ANOVA: feature ~ cancer_type_corrected_updated
  aov_model <- aov(as.formula(paste(feat, "~ cancer_type_corrected_updated")), data = df_complete)
  anova_summary <- summary(aov_model)
  
  # Compute eta-squared: eta_sq = SSeffect / SStotal
  ss_effect <- anova_summary[[1]]["cancer_type_corrected_updated", "Sum Sq"]
  ss_total <- sum(anova_summary[[1]][, "Sum Sq"])
  eta_sq <- ss_effect / ss_total
  
  # Convert eta-squared to Cohen's f: f = sqrt(eta_sq / (1 - eta_sq))
  f_effect <- sqrt(eta_sq / (1 - eta_sq))
  
  # Number of groups and average sample size per group for the factor
  group_counts <- table(df_complete$cancer_type_corrected_updated)
  k <- length(group_counts)
  avg_n <- mean(group_counts)
  
  # Power analysis using pwr.anova.test
  power_anova <- pwr.anova.test(k = k, n = avg_n, f = f_effect, sig.level = 0.05)
  
  anova_power_results <- rbind(anova_power_results, data.frame(feature = feat, 
                                                               f_effect = f_effect, 
                                                               eta_sq = eta_sq, 
                                                               power = power_anova$power, 
                                                               avg_n = avg_n, 
                                                               k = k))
}
cat("ANOVA power analysis results for differences among cancer types:\n")
print(anova_power_results)

write.csv(anova_power_results, "ANOVA_Power_Analysis_CancerTypes.csv", row.names = FALSE)

# -------------------------------
# 3. Sensitivity Analysis: Power vs. Sample Size for All Features
# -------------------------------
# We'll vary the sample size per group from a minimum (e.g., 20) up to the max observed sample size in the Cancer vs. Normal comparison.
sample_sizes <- seq(20, 1000, by = 5)

# Create an empty data frame to store sensitivity results for all features
sensitivity_df <- expand.grid(feature = features, n = sample_sizes)
sensitivity_df$power <- NA

# Fill in the power values using a two-sample t-test for each feature (assuming equal group sizes)
for (feat in features) {
  # Get Cohen's d for the feature
  d_val <- effect_sizes$cohen_d[effect_sizes$feature == feat]
  for (i in seq_along(sample_sizes)) {
    n_temp <- sample_sizes[i]
    power_temp <- pwr.t.test(n = n_temp, d = d_val, sig.level = 0.05, 
                             type = "two.sample", alternative = "two.sided")$power
    sensitivity_df$power[sensitivity_df$feature == feat & sensitivity_df$n == n_temp] <- power_temp
  }
}

# Plot sensitivity curves for all features on the same plot
p_sensitivity <- ggplot(sensitivity_df, aes(x = n, y = power, color = feature)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(title = "Sensitivity Analysis: Power vs. Sample Size\n(Cancer vs. Normal Comparison)",
       x = "Sample Size per Group",
       y = "Statistical Power",
       color = "Feature") +
  theme_classic(base_size = 14)
ggsave("Sensitivity_Power_All_Features.pdf", p_sensitivity, width = 10, height = 6)

# -------------------------------
# 4. Manuscript Paragraphs for Reporting
# -------------------------------
cat("\n--- Manuscript Reporting ---\n")
cat("We conducted a comprehensive power analysis to assess the sensitivity of our study design for detecting differences between cancer and normal groups as well as among multiple cancer types. \n\n")

cat("For the cancer versus normal comparison, effect sizes were computed using Cohen's d for key features (e.g., Endmotif_zscore, Insertsize_zscore, Delfi_zscore, Nucleosome_peak_zscore, Twentyto150_pct, and FS). The observed effect sizes ranged from ", 
    round(min(effect_sizes$cohen_d), 2), " to ", round(max(effect_sizes$cohen_d), 2), 
    ". Based on the current sample sizes (Normal: ", n_normal, ", Cancer: ", n_cancer, 
    "), the achieved power for these features varied (see Table 1), and the estimated required sample sizes to achieve 80% power ranged from ", 
    min(power_results_ttest$required_n_equal), " to ", max(power_results_ttest$required_n_equal), " per group. \n\n")

cat("In addition, one-way ANOVA was performed to assess differences among individual cancer types. Effect sizes were quantified using eta-squared (η²) and converted to Cohen's f, yielding values that suggest a moderate effect (see Table 2). The power analysis based on these parameters indicates that, given the average group sizes and the number of cancer type categories, the study is moderately powered to detect differences among cancer types. \n\n")

cat("A sensitivity analysis was conducted to explore how the power of a two-sample t-test varies with sample size for each feature (Figure 1). The combined sensitivity curves demonstrate that as the sample size per group increases, the power to detect significant differences improves across all features. This analysis also provides insight into the potential impact of reducing the sample size, where power would decrease and the risk of type II error would increase. \n\n")

cat("Overall, our analysis suggests that the current sample sizes are sufficient for the primary cancer versus normal comparison, while also being moderately powered to detect differences among the individual cancer types. However, for features with lower effect sizes, increasing the sample size would further enhance our ability to detect significant differences.\n")




### Get power summary 
# Suppose you tested up to 355 samples in combined_sensitivity
max_sample_tested <- 355

power_summary <- sensitivity_df %>%
  group_by(Feature) %>%
  summarize(
    # If there's any point where Power >= 0.8, get the smallest such Sample_Size
    min_n_for_80 = if (any(Power >= 0.8)) {
      min(Sample_Size[Power >= 0.8])
    } else {
      NA_real_
    },
    # Check power at the max sample size tested (e.g., 355)
    power_at_max = if (any(Sample_Size == max_sample_tested)) {
      Power[Sample_Size == max_sample_tested]
    } else {
      NA_real_
    }
  ) %>%
  # Optionally, add a note if 80% power wasn't reached
  mutate(
    comment = ifelse(is.na(min_n_for_80),
                     paste("80% power not reached; max power =", round(power_at_max, 3)),
                     paste("80% power reached at n =", min_n_for_80))
  )

print(power_summary)





#### Between cancer types 
# -------------------------------
# Load Required Libraries
# -------------------------------

df_complete <- Merged_zscore_df %>%
  mutate(cancer_type_corrected_updated = as.factor(cancer_type_corrected_updated))

# Count sample sizes for each cancer type
cancer_type_counts <- df_complete %>%
  group_by(cancer_type_corrected_updated) %>%
  summarise(n = n(), .groups = "drop")

print(cancer_type_counts)  # Print table of sample sizes per group

# -------------------------------
# 2. One-Way ANOVA Power Analysis for All Cancer Types
# -------------------------------
anova_power_results <- data.frame(feature = character(), eta_sq = numeric(), 
                                  f_effect = numeric(), power = numeric(), 
                                  required_n = numeric(), k = numeric(), stringsAsFactors = FALSE)

features <- c("Endmotif_zscore", "Insertsize_zscore", "Delfi_zscore", 
              "Nucleosome_peak_zscore", "Twentyto150_pct", "FS")

for (feat in features) {
  # Run one-way ANOVA
  aov_model <- aov(as.formula(paste(feat, "~ cancer_type_corrected_updated")), data = df_complete)
  anova_summary <- summary(aov_model)
  
  # Compute eta-squared (η² = SS_effect / SS_total)
  ss_effect <- anova_summary[[1]]["cancer_type_corrected_updated", "Sum Sq"]
  ss_total <- sum(anova_summary[[1]][, "Sum Sq"])
  eta_sq <- ss_effect / ss_total
  
  # Convert eta-squared to Cohen's f: f = sqrt(η² / (1 - η²))
  f_effect <- sqrt(eta_sq / (1 - eta_sq))
  
  # Get the number of groups (cancer types)
  k <- length(unique(df_complete$cancer_type_corrected_updated))
  
  # Compute power with current average sample size per group
  avg_n <- mean(cancer_type_counts$n)
  power_anova <- pwr.anova.test(k = k, n = avg_n, f = f_effect, sig.level = 0.05)
  
  # Estimate required sample size per group to achieve 80% power
  required_n <- pwr.anova.test(k = k, f = f_effect, power = 0.8, sig.level = 0.05)$n
  
  anova_power_results <- rbind(anova_power_results, data.frame(feature = feat, 
                                                               eta_sq = eta_sq, 
                                                               f_effect = f_effect, 
                                                               power = power_anova$power, 
                                                               required_n = ceiling(required_n), 
                                                               k = k))
}

print(anova_power_results)  # Print the power analysis table
write.csv(anova_power_results, "ANOVA_Power_Analysis_CancerTypes.csv", row.names = FALSE)

# -------------------------------
# 3. Pairwise t-test Power Analysis (Cancer Type Comparisons)
# -------------------------------
pairwise_power_results <- data.frame(feature = character(), group1 = character(),
                                     group2 = character(), cohen_d = numeric(),
                                     power = numeric(), required_n = numeric(),
                                     stringsAsFactors = FALSE)

for (feat in features) {
  # Get all cancer type combinations
  cancer_types <- unique(df_complete$cancer_type_corrected_updated)
  for (i in 1:(length(cancer_types) - 1)) {
    for (j in (i + 1):length(cancer_types)) {
      g1 <- cancer_types[i]
      g2 <- cancer_types[j]
      
      # Subset data for the two groups
      subset_data <- df_complete %>%
        filter(cancer_type_corrected_updated %in% c(g1, g2)) %>%
        mutate(group = as.factor(cancer_type_corrected_updated))
      
      # Compute Cohen's d
      d_result <- cohen.d(subset_data[[feat]], subset_data$group, na.rm = TRUE)
      
      # Get sample sizes
      n1 <- sum(subset_data$group == g1)
      n2 <- sum(subset_data$group == g2)
      
      # Compute achieved power
      power_test <- pwr.t2n.test(n1 = n1, n2 = n2, d = d_result$estimate, sig.level = 0.05)
      
      # Estimate required sample size per group for 80% power
      required_n <- pwr.t.test(d = d_result$estimate, power = 0.8, sig.level = 0.05, 
                               type = "two.sample", alternative = "two.sided")$n
      
      pairwise_power_results <- rbind(pairwise_power_results, 
                                      data.frame(feature = feat, 
                                                 group1 = g1, 
                                                 group2 = g2, 
                                                 cohen_d = d_result$estimate, 
                                                 power = power_test$power, 
                                                 required_n = ceiling(required_n)))
    }
  }
}

print(pairwise_power_results)  # Print the power analysis results
write.csv(pairwise_power_results, "Pairwise_Power_Analysis_CancerTypes.csv", row.names = FALSE)

# -------------------------------
# 4. Plot Power vs. Sample Size for All Features (ANOVA)
# -------------------------------
sample_sizes <- seq(20, max(cancer_type_counts$n), by = 5)
sensitivity_df <- expand.grid(feature = features, n = sample_sizes)
sensitivity_df$power <- NA

for (feat in features) {
  f_effect <- anova_power_results$f_effect[anova_power_results$feature == feat]
  k <- anova_power_results$k[anova_power_results$feature == feat]
  for (i in seq_along(sample_sizes)) {
    n_temp <- sample_sizes[i]
    power_temp <- pwr.anova.test(k = k, n = n_temp, f = f_effect, sig.level = 0.05)$power
    sensitivity_df$power[sensitivity_df$feature == feat & sensitivity_df$n == n_temp] <- power_temp
  }
}

p_sensitivity <- ggplot(sensitivity_df, aes(x = n, y = power, color = feature)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(title = "Sensitivity Analysis: Power vs. Sample Size (ANOVA)",
       x = "Sample Size per Group",
       y = "Statistical Power",
       color = "Feature") +
  theme_classic(base_size = 14)
ggsave("ANOVA_Sensitivity_Power.pdf", p_sensitivity, width = 8, height = 6)

# -------------------------------
# Final Interpretation
# -------------------------------
cat("\nInterpretation:\n")
cat("1. The one-way ANOVA power analysis estimates the power to detect differences among multiple cancer types based on eta-squared effect sizes.\n")
cat("2. Pairwise t-test power analysis assesses whether we are sufficiently powered to detect differences between individual cancer types.\n")
cat("3. The sensitivity analysis plots show how power improves with increasing sample size, guiding how many samples are needed for stronger conclusions.\n")

