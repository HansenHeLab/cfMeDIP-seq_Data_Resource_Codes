# ----------------------------------------------------------------------------
# Title   : Correlation of DNASE1L3 Motif Changes with Fragmentation Score
# Authors : Dory Abelman
# Date    : February 2025 (Revision)
#
# Purpose :
#   • Compare cancer‑type–level averages of the weighted Fragmentation Score (FS)
#     with DNASE1L3‑associated motif fold‑changes.
#   • Generate a scatter plot and fit a linear model across cancer types,
#     testing Pearson and Spearman correlations.
#   • At the sample level, merge individual FS values with per‑sample DNASE1L3
#     median fold‑change scores, compute Pearson correlation, and visualize
#     colored by cancer subtype.
#   • Save correlation statistics and publication‑ready plots for revision figures.
# ----------------------------------------------------------------------------


##### Added for revision Feb 2025
## Look at the correlation with DNASE1L3 and the fragment score 
## Using merged zscore df from other part in later stages
# --- Aggregate FS by cancer type from Merged_zscore_df ---
fs_summary <- Merged_zscore_df %>%
  group_by(cancer_type_corrected_updated) %>%
  summarise(mean_FS = mean(FS, na.rm = TRUE),
            sd_FS = sd(FS, na.rm = TRUE),
            n = n()) %>%
  mutate(cancer_type_corrected_updated = tolower(cancer_type_corrected_updated))

# --- Aggregate DNASE1L3 fold change by cancer type from data_stats_dnase ---
tmp <- readRDS("Final Data Dec 2024/End_motifs/Without validation/DNASE1L3_motif_statistics_by_cancer_type.rds")
dnase_summary <- tmp %>%
  group_by(cancer_type_corrected_updated) %>%
  summarise(mean_dnaseFC = mean(foldchange, na.rm = TRUE),
            sd_dnaseFC = sd(foldchange, na.rm = TRUE),
            n = n()) %>%
  mutate(cancer_type_corrected_updated = tolower(as.character(cancer_type_corrected_updated)))


# Merge the summaries by cancer type
combined_df <- merge(fs_summary, dnase_summary, by = "cancer_type_corrected_updated")

# Create scatterplot
p <- ggplot(combined_df, aes(x = mean_FS, y = mean_dnaseFC, color = cancer_type_corrected_updated)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(x = "Average Weighted Fragment Score (FS)",
       y = "Average DNASE1L3 Fold Change",
       title = "Correlation between FS and DNASE1L3 Variation by Cancer Type",
       color = "Cancer Type") +
  theme_classic(base_size = 14)

print(p)
ggsave("FS_vs_DNASE1L3_Correlation_colored.png", plot = p, width = 8, height = 5, dpi = 300)

# Perform a correlation test
cor_test <- cor.test(combined_df$mean_FS, combined_df$mean_dnaseFC)
print(cor_test)

# Perform a Spearman's rank correlation test
cor_test <- cor.test(combined_df$mean_FS, combined_df$mean_dnaseFC, method = "spearman")
print(cor_test)

# Export the correlation test result to a text file
sink("Spearman_Correlation_Test_DNASE1L3.txt")
print(cor_test)
sink()




#### Now do per sample 
# --- Aggregate DNASE1L3 fold change by cancer type from data_stats_dnase ---
dnase_summary <- read_csv("Final Data Dec 2024/End_motifs/DNASE1L3_score_PE.csv")
dnase_summary <- merge(metadata_df, dnase1l3_score, by.x = "endmotif_name", by.y = "sample", all.x = FALSE)
dnase_summary <- dnase_summary %>%
  mutate(cancer_type_corrected_updated = tolower(as.character(cancer_type_corrected_updated)))


# Merge the summaries by cancer type
combined_df <- merge(all_FS, dnase_summary)

# Calculate Pearson correlation and p-value
cor_test <- cor.test(combined_df$FS, combined_df$DNASE1L3_median_score, method = "pearson")

# Extract correlation and p-value
cor_value <- round(cor_test$estimate, 3)
p_value <- signif(cor_test$p.value, 3)
cor_text <- paste0("Pearson r = ", cor_value, "\nP-value = ", p_value)

# Update label 
combined_df <- combined_df %>%
  mutate(cancer_type_title_case = ifelse(cancer_type_title_case == "Eye_cancer", "Uveal_melanoma", cancer_type_title_case))

cancer_type_col_split <- c(
  "Healthy" = "#33a02c",
  "Brain_cancer" = "#1f78b4",
  "Lung_cancer" = "#b2df8a",
  "Prostate_cancer" = "#a6cee3",
  "AML" = "#fb9a99",
  "Pancreatic_cancer" = "#e31a1c",
  "Eye_cancer" = "#fdbf6f",
  "Uveal_melanoma" = "#fdbf6f",  # same as Eye_cancer
  "Head_and\nneck_cancer" = "#ff7f00",
  "Breast_cancer" = "#cab2d6",
  "Colorectal_cancer" = "#6a3d9a",
  "Bladder_cancer" = "#fb6a4a",
  "Renal_cancer" = "#b15928",
  "LFS_survivor" = "#bdbdbd",
  "LFS_previvor" = "#969696",
  "LFS_positive" = "#737373",
  "Liver_cancer" = "#ffeda0",
  "Melanoma" = "#41b6c4",
  "Mixed_cancer" = "#e31a1c",
  "Ovarian_cancer" = "#ffffb3"
)

# Scatter plot with regression line and correlation text
Corr_DNASE <- ggplot(combined_df, aes(x = FS, y = DNASE1L3_mean_score, col = cancer_type_title_case)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  annotate("text", 
           x = max(combined_df$FS, na.rm = TRUE), 
           y = max(combined_df$DNASE1L3_mean_score, na.rm = TRUE), 
           label = cor_text, 
           hjust = 1, 
           vjust = 1, 
           size = 3) +
  xlab("Fragment Score") + 
  ylab("Mean DNASE1L3 Fold Change") +
  scale_color_manual(name = "Cancer Type", values = cancer_type_col_split) +
  ggtitle("Mean DNASE1L3 Fold Change vs Fragment Score") + 
  theme_classic() +
  theme(
    legend.position = "right")

Corr_DNASE

ggsave("FS_vs_DNASE1L3_Correlation_colored_per_sample_mean.png", plot = Corr_DNASE, width = 6, height = 3.5, dpi = 500)



