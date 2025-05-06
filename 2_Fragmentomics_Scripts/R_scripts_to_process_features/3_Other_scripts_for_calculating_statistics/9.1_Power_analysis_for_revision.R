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
library(dplyr)
library(effsize)    # for computing Cohen's d
library(ggplot2)
library(pwr)        # for power analysis
library(tidyr)      # for data manipulation
library(gridExtra)  # for arranging plots


##############################################
# SECTION 0: Load and Prepare Metadata
##############################################
# Load the metadata containing sequencing_id, CN_classifier, cancer_type_corrected, etc.
metadata_df_tmp <- readRDS("Metadata_df_all_with_corrected_cancer_subtype_Jan2025.rds")

# Exclude samples with "lfs" (case-insensitive) in cancer_type_corrected
metadata_df_tmp <- metadata_df_tmp %>% 
  filter(!grepl("lfs", cancer_type_corrected, ignore.case = TRUE))

# Define group based on CN_classifier:
# Assume that CN_classifier == "healthy" means Normal and anything else is Cancer.
metadata_df_tmp <- metadata_df_tmp %>% 
  mutate(group = ifelse(tolower(CN_classifier) == "healthy", "Normal", "Cancer"))

# For consistency, convert group to a factor with levels Normal then Cancer.
metadata_df_tmp$group <- factor(metadata_df_tmp$group, levels = c("Normal", "Cancer"))


##############################################
# Combined & Extended Power Analysis for 
#   (1) DMRs (Methylation), 
#   (2) Fragmentation Features USING Z-SCORES, and 
#   (3) ANOVA Power Analysis for Differences Among Cancer Types
##############################################

##############################################
# SECTION 1: DMR (Methylation) Analysis
##############################################

# 1.1 Load Data for DMR Analysis
# norm_counts is assumed to be an RDS file with rows = DMRs and columns = samples.
norm_counts <- readRDS("Final Data Dec 2024/Methylation_data/Updated_Feb3/PanCancer_hyper_DMRs_combat_deseq2_normalized_count_with_blood_cell_age_sex_filtering_for_PE_samples_after_QC.RDS")

# Match samples between norm_counts and metadata using sequencing_id
common_samples <- intersect(colnames(norm_counts), metadata_df_tmp$sequencing_id)
norm_counts <- norm_counts[, common_samples]
dmr_metadata <- metadata_df_tmp %>% filter(sequencing_id %in% common_samples)

# 1.2 Calculate Effect Sizes for Each DMR (Cohen's d)
dmr_effects <- apply(norm_counts, 1, function(x) {
  tryCatch({
    cohen.d(x, dmr_metadata$group, na.rm = TRUE)$estimate
  }, error = function(e) NA)
})
dmr_effects <- dmr_effects[!is.na(dmr_effects)]

# Use the median of the absolute effect sizes as the representative value
median_dmr_d <- median(abs(dmr_effects))
cat("Median absolute Cohen's d for DMRs:", median_dmr_d, "\n")

# 1.3 Sensitivity Analysis for DMRs: Power vs. Sample Size (using median effect size)
dmr_sample_sizes <- seq(10, 1000, by = 10)
dmr_power_values <- sapply(dmr_sample_sizes, function(n) {
  pwr.t.test(n = n, d = median_dmr_d, sig.level = 0.05,
             type = "two.sample", alternative = "two.sided")$power
})
sensitivity_dmr <- data.frame(Data_Type = "DMR", 
                              Sample_Size = dmr_sample_sizes, 
                              Power = dmr_power_values)
sensitivity_dmr$Feature <- "PE DMR (median)"
sensitivity_dmr$Analysis <- "DMR"



# 1.4 Plot Distribution of Absolute DMR Effect Sizes
df_dmr_effects <- data.frame(Effect_Size = abs(dmr_effects))
p_dmr_effects <- ggplot(df_dmr_effects, aes(x = Effect_Size)) +
  geom_histogram(bins = 50, fill = "skyblue", color = "black") +
  labs(title = "Distribution of Absolute Cohen's d for DMRs in PE samples",
       x = "Absolute Cohen's d", y = "Count") +
  theme_classic(base_size = 14)
ggsave("DMR_Effect_Size_Distribution_PE.pdf", p_dmr_effects, width = 8, height = 5)
print(p_dmr_effects)

##############################################
# SECTION 1B: DMR (Methylation) Analysis from SE Data
##############################################

# 1B.1 Load Data for SE DMR Analysis
norm_counts_se <- readRDS("Final Data Dec 2024/Methylation_data/Updated_Feb3/PanCancer_hyper_DMRs_combat_deseq2_normalized_count_with_blood_cell_age_sex_filtering_for_SE_samples_after_QC.RDS")

# Match samples between norm_counts_se and metadata using sequencing_id
common_samples_se <- intersect(colnames(norm_counts_se), metadata_df_tmp$sequencing_id)
norm_counts_se <- norm_counts_se[, common_samples_se]
dmr_metadata_se <- metadata_df_tmp %>% filter(sequencing_id %in% common_samples_se)

# 1B.2 Calculate Effect Sizes for Each DMR in SE Data (Cohen's d)
dmr_effects_se <- apply(norm_counts_se, 1, function(x) {
  tryCatch({
    cohen.d(x, dmr_metadata_se$group, na.rm = TRUE)$estimate
  }, error = function(e) NA)
})
dmr_effects_se <- dmr_effects_se[!is.na(dmr_effects_se)]
# Use the median of the absolute effect sizes as a representative metric for SE data
median_dmr_d_se <- median(abs(dmr_effects_se))
cat("Median absolute Cohen's d for SE DMRs:", median_dmr_d_se, "\n")

mean_dmr_d_se <- mean(abs(dmr_effects_se))
cat("Mean absolute Cohen's d for SE DMRs:", mean_dmr_d_se, "\n")

# 1B.3 Sensitivity Analysis for SE DMRs: Power vs. Sample Size
# Determine maximum available sample size per group (using df_power from fragmentation features analysis)
# Determine maximum available sample size per group (using df_power from fragmentation features analysis)
n_cancer_se <- sum(df_power$group == "Cancer")
max_sample_size_se <- min(n_cancer_se)
#dmr_sample_sizes_se <- seq(10, max_sample_size_se, by = 10)
dmr_sample_sizes_se <- seq(10, 1000, by = 10)

# Calculate power using the median effect size for SE DMRs
dmr_power_values_se <- sapply(dmr_sample_sizes_se, function(n) {
  pwr.t.test(n = n, d = median_dmr_d_se, sig.level = 0.05, 
             type = "two.sample", alternative = "two.sided")$power
})
sensitivity_dmr_se <- data.frame(Data_Type = "SE DMR", 
                                 Sample_Size = dmr_sample_sizes_se, 
                                 Power = dmr_power_values_se)
sensitivity_dmr_se$Feature <- "SE DMR (median)"
sensitivity_dmr_se$Analysis <- "SE DMR"

# Calculate power using the mean effect size for SE DMRs
dmr_power_values_se_mean <- sapply(dmr_sample_sizes_se, function(n) {
  pwr.t.test(n = n, d = mean_dmr_d_se, sig.level = 0.05, 
             type = "two.sample", alternative = "two.sided")$power
})
sensitivity_dmr_se_mean <- data.frame(Data_Type = "SE DMR", 
                                      Sample_Size = dmr_sample_sizes_se, 
                                      Power = dmr_power_values_se_mean)
sensitivity_dmr_se_mean$Feature <- "SE DMR (mean)"
sensitivity_dmr_se_mean$Analysis <- "SE DMR"


# Optionally, print a summary table for the SE DMR analysis
se_dmr_summary <- data.frame(
  Median_Cohen_d = median_dmr_d_se,
  Required_Sample_Size_for_80_Power = dmr_sample_sizes_se[which(dmr_power_values_se >= 0.8)[1]]
)
cat("SE DMR summary:\n")
print(se_dmr_summary)

# Create a data frame for the SE DMR effect sizes
df_dmr_effects_se <- data.frame(Effect_Size = abs(dmr_effects_se))

# Plot histogram of the absolute Cohen's d for SE DMRs
p_dmr_effects_se <- ggplot(df_dmr_effects_se, aes(x = Effect_Size)) +
  geom_histogram(bins = 50, fill = "skyblue", color = "black") +
  labs(title = "Distribution of Absolute Cohen's d for DMRs in SE samples",
       x = "Absolute Cohen's d", y = "Count") +
  theme_classic(base_size = 14)
ggsave("SE_DMR_Effect_Size_Distribution_PE.pdf", p_dmr_effects_se, width = 8, height = 5)
print(p_dmr_effects_se)


##############################################
# SECTION 2: Fragmentation Features Analysis
##############################################

# 2.1 Load and Prepare Data for Fragment Features
# Assume Merged_zscore_df is already available; otherwise, load it (e.g., readRDS("Merged_zscore_df.RDS"))
# Merge with metadata using sequencing_id.
df_power <- Merged_zscore_df %>%
  mutate(group = ifelse(cancer_type_corrected_updated == "healthy", "Normal", "Cancer"))

# 2.2 Exclude any remaining LFS samples (if not already removed)
df_power <- df_power %>% filter(!grepl("lfs", cancer_type_corrected_updated, ignore.case = TRUE))

# 2.3 Calculate sample sizes for Cancer vs. Normal in fragmentation data
n_normal <- sum(df_power$group == "Normal")
n_cancer <- sum(df_power$group == "Cancer")
cat("Fragment features sample sizes: Normal =", n_normal, "Cancer =", n_cancer, "\n")

# 2.4 Calculate Effect Sizes (Cohen's d) for Fragmentation Features
features_frag <- c("Endmotif_zscore", "Insertsize_zscore", "Delfi_zscore", 
                   "Nucleosome_peak_zscore", "Twentyto150_pct", "FS")
frag_effects <- data.frame(feature = character(), cohen_d = numeric(), stringsAsFactors = FALSE)
for(feat in features_frag){
  d_res <- cohen.d(df_power[[feat]], df_power$group, na.rm = TRUE)
  frag_effects <- rbind(frag_effects, data.frame(feature = feat, cohen_d = d_res$estimate))
}
print(frag_effects)

# 2.5 Sensitivity Analysis for Fragment Features: Power vs. Sample Size
#frag_sample_sizes <- seq(5, max(n_normal, n_cancer), by = 5)
frag_sample_sizes <- seq(10, 1000, by = 10)
sensitivity_frag <- expand.grid(Data_Type = features_frag, Sample_Size = frag_sample_sizes)
sensitivity_frag$Power <- NA
for(feat in features_frag){
  d_val <- frag_effects$cohen_d[frag_effects$feature == feat]
  for(n in frag_sample_sizes){
    sensitivity_frag$Power[sensitivity_frag$Data_Type == feat & sensitivity_frag$Sample_Size == n] <-
      pwr.t.test(n = n, d = d_val, sig.level = 0.05, type = "two.sample", alternative = "two.sided")$power
  }
}
# Modify Data_Type to indicate Fragment features
sensitivity_frag <- sensitivity_frag %>% 
  mutate(Feature = paste("Fragment:", Data_Type),
         Analysis = "Fragment")

# 2.6 Plot Cohen's d Effect Sizes for Fragment Features
p_frag_effect <- ggplot(frag_effects, aes(x = reorder(feature, cohen_d), y = cohen_d)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +
  labs(title = "Effect Sizes (Cohen's d) for Fragment Features\nCancer vs. Normal", 
       x = "Feature", y = "Cohen's d") +
  theme_classic(base_size = 14)
ggsave("Fragment_Effect_Sizes.pdf", p_frag_effect, width = 8, height = 5)
print(p_frag_effect)

# Compute mean absolute DMR effect size (from previously computed dmr_effects)
mean_dmr_d <- mean(abs(dmr_effects))
cat("Mean absolute Cohen's d for DMRs:", mean_dmr_d, "\n")

# Create a combined data frame with fragmentation features and DMR
combined_effects <- frag_effects %>%
  mutate(feature = as.character(feature))  # ensure it's character for binding

# Add a new row for DMR with its mean effect size
combined_effects <- rbind(combined_effects,
                          data.frame(feature = "PE DMR (median)", cohen_d = median_dmr_d))

# Add a new row for SE DMR with its mean effect size
combined_effects <- rbind(combined_effects,
                          data.frame(feature = "SE DMR (median)", cohen_d = median_dmr_d_se))

# Step 1: Recode feature names with cleaned labels (no "zscore", no "median")
combined_effects$feature <- dplyr::recode(combined_effects$feature,
                                          "Endmotif_zscore"        = "End Motifs",
                                          "Insertsize_zscore"      = "Insert Size",
                                          "Delfi_zscore"           = "Fragment Ratios",
                                          "Nucleosome_peak_zscore" = "Nucleosome Peak",
                                          "FS"                     = "Fragment Size Score",
                                          "Twentyto150_pct"        = "Prop_short_frags",  # Temp for filtering
                                          .default = combined_effects$feature)

# Step 2: Remove 'Prop_short_frags' (Twentyto150_pct) from the dataframe
combined_effects <- combined_effects %>%
  filter(feature != "Prop_short_frags")

# Step 3: Define the updated color palette (no "median" in labels)
custom_color_palette <- c(
  "End Motifs"           = "#E69F00",
  "Insert Size"          = "#56B4E9",
  "Fragment Ratios"      = "#009E73",
  "Nucleosome Peak"      = "#0072B2",
  "Fragment Size Score"  = "#F0E442",
  "PE DMR (median)"      = "#CC79A7",  # Leave this one as-is if it's a median
  "SE DMR (median)"      = "#9932CC"   # Leave as-is unless you want to rename
)

# Step 4: Plot the cleaned feature names with consistent color scheme
p_combined_effects_cleaned <- ggplot(combined_effects, aes(x = reorder(feature, cohen_d), y = cohen_d, fill = feature)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = custom_color_palette) +
  labs(title = "Effect Sizes (Cohen's d)",
       x = "Feature", y = "Cohen's d") +
  theme_classic() +
  theme(legend.position = "none")

ggsave("Combined_Effect_Sizes_Cleaned_zscore_narrow.pdf", p_combined_effects_cleaned, width = 5, height = 5)
print(p_combined_effects_cleaned)

write.csv(combined_effects, "combined_effect_sizes.csv", row.names = FALSE)




# Plot all effect sizes in one bar plot
### old version of fig
p_combined_effects <- ggplot(combined_effects, aes(x = reorder(feature, cohen_d), y = cohen_d)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +
  labs(title = "Effect Sizes (Cohen's d) for Features", 
       x = "Feature", y = "Cohen's d") +
  theme_classic(base_size = 14)

ggsave("Combined_Effect_Sizes_DMR_with_SE.pdf", p_combined_effects, width = 8, height = 6)
print(p_combined_effects)


##############################################
# SECTION 3: ANOVA Power Analysis for Differences Among Cancer Types
##############################################

# Use the complete merged Z-score data as df_complete.
df_complete <- Merged_zscore_df %>%
  mutate(cancer_type_corrected_updated = as.factor(cancer_type_corrected_updated)) %>%
  filter(!grepl("lfs", cancer_type_corrected_updated, ignore.case = TRUE))
df_complete <- droplevels(df_complete)

# Count sample sizes for each cancer type.
cancer_type_counts <- df_complete %>%
  group_by(cancer_type_corrected_updated) %>%
  summarise(n = n(), .groups = "drop")
print(cancer_type_counts)

# 3.1 One-Way ANOVA Power Analysis for Each Fragmentation Feature Among Cancer Types
anova_power_results <- data.frame(feature = character(), f_effect = numeric(), 
                                  eta_sq = numeric(), power = numeric(), 
                                  avg_n = numeric(), k = numeric(), stringsAsFactors = FALSE)
for (feat in features_frag) {
  aov_model <- aov(as.formula(paste(feat, "~ cancer_type_corrected_updated")), data = df_complete)
  anova_summary <- summary(aov_model)
  ss_effect <- anova_summary[[1]]["cancer_type_corrected_updated", "Sum Sq"]
  ss_total <- sum(anova_summary[[1]][, "Sum Sq"])
  eta_sq <- ss_effect / ss_total
  f_effect <- sqrt(eta_sq / (1 - eta_sq))
  group_counts <- table(df_complete$cancer_type_corrected_updated)
  k <- length(group_counts)
  avg_n <- mean(group_counts)
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

### Add DMR 


# Create an aggregate DMR metric for each sample (e.g., median methylation across all DMRs)
# Note: norm_counts has rows = DMRs and columns = samples.
dmr_metric <- apply(norm_counts, 2, median, na.rm = TRUE)

# Build a data frame linking sample IDs, cancer types, and the aggregate DMR metric.
dmr_df <- dmr_metadata %>%
  dplyr::select(sequencing_id, cancer_type_corrected_updated) %>%
  mutate(DMR_metric = dmr_metric[match(sequencing_id, names(dmr_metric))])
dmr_df$cancer_type_corrected_updated <- as.factor(dmr_df$cancer_type_corrected_updated)

# Perform one-way ANOVA for the aggregate DMR metric across cancer types.
aov_model_dmr <- aov(DMR_metric ~ cancer_type_corrected_updated, data = dmr_df)
anova_summary_dmr <- summary(aov_model_dmr)

# Extract sums of squares and compute eta-squared.
ss_effect_dmr <- anova_summary_dmr[[1]]["cancer_type_corrected_updated", "Sum Sq"]
ss_total_dmr <- sum(anova_summary_dmr[[1]][, "Sum Sq"])
eta_sq_dmr <- ss_effect_dmr / ss_total_dmr

# Convert eta-squared to Cohen's f.
f_effect_dmr <- sqrt(eta_sq_dmr / (1 - eta_sq_dmr))

# Get group information and compute the average sample size.
group_counts_dmr <- table(dmr_df$cancer_type_corrected_updated)
k_dmr <- length(group_counts_dmr)
avg_n_dmr <- mean(group_counts_dmr)

# Estimate power using pwr.anova.test.
power_anova_dmr <- pwr.anova.test(k = k_dmr, n = avg_n_dmr, f = f_effect_dmr, sig.level = 0.05)

# Create a summary row for the DMR analysis.
dmr_anova_results <- data.frame(feature = "DMR (aggregate metric)",
                                f_effect = f_effect_dmr,
                                eta_sq = eta_sq_dmr,
                                power = power_anova_dmr$power,
                                avg_n = avg_n_dmr,
                                k = k_dmr,
                                stringsAsFactors = FALSE)

##############################################
# Combine the DMR ANOVA Results with the Fragmentation Features ANOVA Results
##############################################
# (Assuming your existing fragmentation features ANOVA results are stored in 'anova_power_results')
combined_power_results <- rbind(dmr_anova_results, anova_power_results)
cat("Combined ANOVA Power Analysis Results:\n")
print(combined_power_results)
write.csv(combined_power_results, "Combined_ANOVA_Power_Analysis.csv", row.names = FALSE)


### Add SE DMR ANOVA ###

# Create an aggregate SE DMR metric per sample (e.g., median methylation across all SE DMRs)
dmr_metric_se <- apply(norm_counts_se, 2, median, na.rm = TRUE)

# Build a data frame linking SE sample IDs, cancer types, and the aggregate DMR metric.
dmr_df_se <- dmr_metadata_se %>%
  dplyr::select(sequencing_id, cancer_type_corrected_updated) %>%
  mutate(DMR_metric = dmr_metric_se[match(sequencing_id, names(dmr_metric_se))])
dmr_df_se$cancer_type_corrected_updated <- as.factor(dmr_df_se$cancer_type_corrected_updated)

# Perform one-way ANOVA on SE DMR aggregate metric across cancer types
aov_model_dmr_se <- aov(DMR_metric ~ cancer_type_corrected_updated, data = dmr_df_se)
anova_summary_dmr_se <- summary(aov_model_dmr_se)

# Extract sums of squares and compute eta-squared for SE DMRs
ss_effect_dmr_se <- anova_summary_dmr_se[[1]]["cancer_type_corrected_updated", "Sum Sq"]
ss_total_dmr_se <- sum(anova_summary_dmr_se[[1]][, "Sum Sq"])
eta_sq_dmr_se <- ss_effect_dmr_se / ss_total_dmr_se

# Convert eta-squared to Cohen's f
f_effect_dmr_se <- sqrt(eta_sq_dmr_se / (1 - eta_sq_dmr_se))

# Get group information and average sample size
group_counts_dmr_se <- table(dmr_df_se$cancer_type_corrected_updated)
k_dmr_se <- length(group_counts_dmr_se)
avg_n_dmr_se <- mean(group_counts_dmr_se)

# Estimate power for SE DMRs using pwr.anova.test
power_anova_dmr_se <- pwr.anova.test(k = k_dmr_se, n = avg_n_dmr_se, f = f_effect_dmr_se, sig.level = 0.05)

# Create summary row for SE DMR ANOVA result
dmr_anova_results_se <- data.frame(feature = "SE DMR (median)",
                                   f_effect = f_effect_dmr_se,
                                   eta_sq = eta_sq_dmr_se,
                                   power = power_anova_dmr_se$power,
                                   avg_n = avg_n_dmr_se,
                                   k = k_dmr_se,
                                   stringsAsFactors = FALSE)

# Export
combined_power_results <- rbind(combined_power_results, dmr_anova_results_se)
cat("Combined ANOVA Power Analysis Results:\n")
print(combined_power_results)
write.csv(combined_power_results, "Combined_ANOVA_Power_Analysis.csv", row.names = FALSE)


##############################################
# SECTION 4: Combine and Plot Sensitivity Curves Together
##############################################

sensitivity_dmr <- sensitivity_dmr %>% select(Feature, Sample_Size, Power, Analysis)
sensitivity_dmr_se <- sensitivity_dmr_se %>% select(Feature, Sample_Size, Power, Analysis)
sensitivity_frag <- sensitivity_frag %>% select(Feature, Sample_Size, Power, Analysis)
combined_sensitivity <- bind_rows(sensitivity_dmr, sensitivity_frag, sensitivity_dmr_se)

# Recode feature labels in combined_sensitivity
combined_sensitivity <- combined_sensitivity %>%
  mutate(Feature = dplyr::recode(Feature,
                                 "Fragment: Endmotif_zscore"        = "End Motifs",
                                 "Fragment: Insertsize_zscore"      = "Insert Size",
                                 "Fragment: Delfi_zscore"           = "Fragment Ratios",
                                 "Fragment: Nucleosome_peak_zscore" = "Nucleosome Peak",
                                 "Fragment: FS"                     = "Fragment Size Score",
                                 "Fragment: Twentyto150_pct"        = "Prop_short_frags",  # will be removed
                                 "DMR (median)"                       = "PE DMR (median)",
                                 .default = Feature
  )) %>%
  filter(Feature != "Prop_short_frags")  # remove if needed

# Define cleaned color palette
custom_color_palette <- c(
  "End Motifs"           = "#E69F00",
  "Insert Size"          = "#56B4E9",
  "Fragment Ratios"      = "#009E73",
  "Nucleosome Peak"      = "#0072B2",
  "Fragment Size Score"  = "#F0E442",
  "PE DMR (median)"        = "#CC79A7",
  "SE DMR (median)"      = "#9932CC"
)

# Plot with cleaned labels and custom colors
p_combined <- ggplot(combined_sensitivity, aes(x = Sample_Size, y = Power, color = Feature)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "Sensitivity Analysis: Power vs. Sample Size",
    subtitle = "Comparison of DMR (Methylation) and Fragment Features",
    x = "Sample Size per Group",
    y = "Statistical Power",
    color = "Feature"
  ) +
  scale_color_manual(values = custom_color_palette) +
  theme_classic(base_size = 14)

ggsave("Combined_Sensitivity_Power_Cleaned.pdf", p_combined, width = 7, height = 5)
print(p_combined)


write.csv(combined_sensitivity, "Combined_Sensitivity_Table_Cleaned.csv", row.names = FALSE)


##############################################
# SECTION 5: Manuscript Reporting & Tables Summary
##############################################

cat("\n--- Manuscript Reporting ---\n")
cat("For the methylation (DMR) analysis, normalized counts (samples as columns) and metadata (based on sequencing_id and CN_classifier) were used. ",
    "After excluding LFS samples, Cohen's d was computed for each DMR, yielding a median absolute effect size of", round(median_dmr_d, 2), ". ",
    "The sensitivity analysis (Figure 1) indicates that approximately", dmr_sample_sizes[which(dmr_power_values >= 0.8)[1]], 
    "samples per group are needed to achieve 80% power for DMR detection. Additionally, the distribution of absolute DMR effect sizes is shown in 'DMR_Effect_Size_Distribution.pdf'.\n\n")

cat("For the fragmentation features, six key metrics (", paste(features_frag, collapse = ", "), ") were analyzed using merged Z-scores. ",
    "The effect sizes (see 'Fragment_Effect_Sizes.pdf') and sensitivity curves (Figure 2) show that power increases with sample size. ",
    "The ANOVA power analysis among cancer types (Table saved as 'ANOVA_Power_Analysis_CancerTypes.csv') indicates that differences among cancer types are moderately powered given the current group sizes. ",
    "Overall, while the current sample sizes may be adequate for some features, increasing the sample number could further enhance detection power for features with smaller effect sizes.\n")








## Testing 
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
