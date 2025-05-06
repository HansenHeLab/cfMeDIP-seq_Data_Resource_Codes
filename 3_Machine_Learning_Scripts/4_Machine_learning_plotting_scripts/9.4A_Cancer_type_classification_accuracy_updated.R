###############################################################################
# Script:    Classification Accuracy and Probability by Cancer Type
# Author:    Dory Abelman
# Date:      2025‑03‑27
# Purpose:   
#   1) Load per‑sample classifier outputs (Original or Validation cohort)
#   2) Merge predictions with sample metadata (to get cancer_type)
#   3) Calculate overall accuracy by cancer type
#   4) Plot and export barplots of accuracy and boxplots of predicted probabilities
#   5) (Optional) Perform ANOVA + Tukey HSD on predicted probabilities across cancer types
#
# Usage:     
#   - Edit `cohort` ("Original" or "Validation") and input file paths as needed
#   - Run the entire script to produce CSV summaries and PDF/PNG figures
###############################################################################


# Load required libraries
library(dplyr)
library(ggplot2)
library(readr)
library(cowplot)


### Updated Mar 27, 2024

# 1. Load the per-sample classification table.
path <- "Final Data Dec 2024/ML_data/PE/Cancer_vs_normal/"
class_df <- read.csv(paste0(path, "classification_details_cancer_vs_healthy_Original.csv"))

# 2. Filter for the specific file.
## Using the one that did the best for validation 
filtered_class_df <- class_df %>%
  filter(Analysis == "Methylation")

# 3. Load metadata.
metadata_df <- readRDS("Metadata_df_all_with_corrected_cancer_subtype_Jan2025.rds")

# 4. Merge classification data with metadata by matching sample (classification) to sequencing_id (metadata).
#    (Assuming the 'sample' column in the classification data corresponds to the 'sequencing_id' in metadata.)
merged_df <- filtered_class_df %>%
  left_join(metadata_df, by = c("sample" = "sequencing_id"))

# 5. Create a new column indicating whether the classification was correct.
merged_df <- merged_df %>%
  mutate(Correct = ifelse(ActualClass == PredictedClass, TRUE, FALSE))

merged_df_original <- merged_df

# 6. Compute classification accuracy by cancer type.
#    Here, we assume that the metadata has a column named 'cancer_type' that specifies the cancer type.
accuracy_by_cancer <- merged_df %>%
  group_by(cancer_type) %>%
  summarise(
    Accuracy = mean(Correct, na.rm = TRUE),
    Count = n()
  )

# 7. Create a bar plot of accuracy by cancer type.
accuracy_plot <- ggplot(accuracy_by_cancer, aes(x = reorder(cancer_type, Accuracy), y = Accuracy)) +
  geom_bar(stat = "identity", fill = "#0072B2") +
  coord_flip() +
  xlab("Cancer Type") +
  ylab("Classification Accuracy") +
  ggtitle("Classification Accuracy by Cancer Type") +
  theme_minimal()

print(accuracy_plot)
ggsave("Accuracy_by_Cancer_Type.pdf", accuracy_plot, width = 8, height = 6)

# 8. (Optional) Create a boxplot of predicted cancer probability by cancer type, colored by correctness.
boxplot_prob <- ggplot(merged_df, aes(x = cancer_type, y = prob_cancer, fill = Correct)) +
  geom_boxplot() +
  coord_flip() +
  xlab("Cancer Type") +
  ylab("Predicted Probability of Cancer") +
  ggtitle("Distribution of Predicted Cancer Probability by Cancer Type") +
  theme_minimal()

print(boxplot_prob)
ggsave("Prob_Cancer_by_Cancer_Type.pdf", boxplot_prob, width = 8, height = 6)



# 2. Filter out normal/healthy samples. 
# Adjust the filtering condition based on how the 'cancer_type' (or similar) is coded in your metadata.
non_normal_df <- merged_df %>%
  filter(cancer_type != "healthy" & cancer_type != "Normal")

# 3. Create a boxplot of the predicted probability of cancer (prob_cancer) by cancer type with jittered points.
boxplot_prob <- ggplot(non_normal_df, aes(x = cancer_type, y = prob_cancer)) +
  geom_boxplot(fill = "lightblue", outlier.color = "NA") +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.7, color = "darkblue") +
  coord_flip() +
  xlab("Cancer Type") +
  ylab("Predicted Probability of Cancer") +
  ggtitle("Predicted Probability of Cancer by Cancer Type (Excluding Normal)") +
  theme_classic()


print(boxplot_prob)
ggsave("Boxplot_Predicted_Probability_Excluding_Normal.pdf", boxplot_prob, width = 8, height = 5)
ggsave("Boxplot_Predicted_Probability_Excluding_Normal.png", boxplot_prob, width = 8, height = 4, dpi = 500)

# 4. Statistical Testing: Perform one-way ANOVA to assess differences among cancer types.
anova_res <- aov(prob_cancer ~ cancer_type, data = non_normal_df)
summary(anova_res)

# 5. If the ANOVA is significant, perform Tukey's HSD post hoc test.
tukey_res <- TukeyHSD(anova_res)
print(tukey_res)

# Optionally, save the Tukey results to a CSV file.
tukey_df <- as.data.frame(tukey_res$cancer_type)
write.csv(tukey_df, "TukeyHSD_Results_Cancer_Type.csv", row.names = TRUE)


# --- Summarize ANOVA Results ---
# Extract the ANOVA table from your anova_res object.
# Convert the summary ANOVA table (a matrix) to a data frame while preserving all columns
anova_table <- as.data.frame.matrix(summary(anova_res)[[1]])
# Add the row names (terms) as a new column
anova_table$Term <- rownames(anova_table)
# Reorder the columns to bring Term to the front
anova_table <- anova_table[, c("Term", "Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)")]
# Print the table to check
print(anova_table)

# Export the formatted ANOVA table to CSV
write.csv(anova_table, "ANOVA_Summary.csv", row.names = FALSE)


# --- Summarize Tukey HSD Results ---
# Convert the TukeyHSD results for the 'cancer_type' factor into a data frame.
tukey_df <- as.data.frame(tukey_res$cancer_type)
tukey_df$Comparison <- rownames(tukey_df)
rownames(tukey_df) <- NULL
# Reorder columns
tukey_df <- tukey_df[, c("Comparison", "diff", "lwr", "upr", "p adj")]
# Print the table.
kable(tukey_df, caption = "Tukey HSD Post Hoc Test Results for Differences Among Cancer Types (prob_cancer)")
# Export the Tukey results as a CSV.
write.csv(tukey_df, "TukeyHSD_Results_Cancer_Type.csv", row.names = FALSE)


# --- Create a Forest Plot for Tukey HSD Comparisons ---
# Order comparisons by the estimated difference (diff)
tukey_df$Comparison <- factor(tukey_df$Comparison, levels = tukey_df$Comparison[order(tukey_df$diff)])
tukey_plot <- ggplot(tukey_df, aes(x = diff, y = Comparison)) +
  geom_point(color = "darkblue", size = 2) +
  geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2, color = "darkblue") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  xlab("Difference in Mean (diff)") +
  ylab("Comparison") +
  ggtitle("Forest Plot: Tukey HSD Comparisons for Cancer Type Differences") +
  theme_minimal(base_size = 12)

print(tukey_plot)
ggsave("TukeyHSD_ForestPlot.pdf", tukey_plot, width = 10, height = 8)



### Other code 
# Compute the mean, standard deviation, and count of predicted probabilities for each cancer type
mean_prob_by_cancer <- non_normal_df %>%
  group_by(cancer_type) %>%
  summarise(
    mean_prob = mean(prob_cancer, na.rm = TRUE),
    sd_prob = sd(prob_cancer, na.rm = TRUE),
    n = n()
  )
print(mean_prob_by_cancer)

# --- Identify the Three Cancer Types with the Lowest Predicted Probability ---
lowest_three <- mean_prob_by_cancer %>% 
  arrange(mean_prob) %>% 
  head(3)
print(lowest_three)

# Extract the cancer type names as a vector
lowest_groups <- as.character(lowest_three$cancer_type)
cat("The three cancer types with the lowest predicted probability are:", paste(lowest_groups, collapse = ", "), "\n")

# --- Extract Significant Tukey HSD Comparisons Involving Any of the Three Lowest Groups ---
# Load the Tukey HSD results (if not already loaded)
# tukey_df <- read.csv("TukeyHSD_Results_Cancer_Type.csv", stringsAsFactors = FALSE)

# Filter for significant comparisons (p adj < 0.05)
sig_tukey <- tukey_df %>% filter(`p adj` < 0.05)

# Create a regex pattern matching any of the three lowest groups
pattern <- paste(lowest_groups, collapse = "|")

# Filter Tukey comparisons where the Comparison column contains any of the lowest groups
sig_lowest_three <- sig_tukey[grepl(pattern, sig_tukey$Comparison), ]
print(sig_lowest_three)

# Export the significant comparisons to a CSV file
write.csv(sig_lowest_three, "Significant_Tukey_Comparisons_Involving_3_Lowest_Groups.csv", row.names = FALSE)












### For validaiton now 

# 1. Load the per-sample classification table.
path <- "Final Data Dec 2024/ML_data/PE/Cancer_vs_normal/"
class_df <- read.csv(paste0(path, "classification_details_cancer_vs_healthy_Validation.csv"))

### see which class is correct 
most_correct <- class_df %>%
  mutate(correct = (ActualClass == PredictedClass)) %>%
  group_by(Analysis) %>%
  summarise(
    total_predictions = n(),
    correct_predictions = sum(correct),
    accuracy = mean(correct)
  ) %>%
  arrange(desc(accuracy))

conf_df <- class_df %>% filter(ActualClass == "cancer") %>%
  group_by(Analysis) %>%
  summarise(
    total = n(),
    TP = sum(ActualClass == "cancer" & PredictedClass == "cancer"),
    TN = sum(ActualClass == "healthy" & PredictedClass == "healthy"),
    FP = sum(ActualClass == "healthy" & PredictedClass == "cancer"),
    FN = sum(ActualClass == "cancer" & PredictedClass == "healthy"),
    accuracy = (TP + TN) / total
  ) %>%
  arrange(desc(accuracy))


# 2. Filter for the specific file.
filtered_class_df <- class_df %>%
  filter(File == "Combined_ratios_methylation_PCs_1pct_validation_outcomes_healthy_vs_cancer.rds")

# 3. Load metadata.
metadata_df <- readRDS("Metadata_df_all_with_corrected_cancer_subtype_Jan2025.rds")

# 4. Merge classification data with metadata by matching sample (classification) to sequencing_id (metadata).
#    (Assuming the 'sample' column in the classification data corresponds to the 'sequencing_id' in metadata.)
merged_df <- filtered_class_df %>%
  left_join(metadata_df, by = c("sample" = "sequencing_id"))

# 5. Create a new column indicating whether the classification was correct.
merged_df <- merged_df %>%
  mutate(Correct = ifelse(ActualClass == PredictedClass, TRUE, FALSE))

merged_df_validation <- merged_df

# 6. Compute classification accuracy by cancer type.
#    Here, we assume that the metadata has a column named 'cancer_type' that specifies the cancer type.
accuracy_by_cancer_validation <- merged_df %>%
  group_by(cancer_type) %>%
  summarise(
    Accuracy = mean(Correct, na.rm = TRUE),
    Count = n()
  )

# 7. Create a bar plot of accuracy by cancer type.
accuracy_plot <- ggplot(accuracy_by_cancer_validation, aes(x = reorder(cancer_type, Accuracy), y = Accuracy)) +
  geom_bar(stat = "identity", fill = "#0072B2") +
  coord_flip() +
  xlab("Cancer Type") +
  ylab("Classification Accuracy") +
  ggtitle("Classification Accuracy by Cancer Type") +
  theme_minimal()

print(accuracy_plot)
ggsave("Accuracy_by_Cancer_Type_validation_new_motif_methylation.pdf", accuracy_plot, width = 8, height = 6)

# 8. (Optional) Create a boxplot of predicted cancer probability by cancer type, colored by correctness.
boxplot_prob <- ggplot(merged_df, aes(x = cancer_type, y = prob_cancer, fill = Correct)) +
  geom_boxplot() +
  coord_flip() +
  xlab("Cancer Type") +
  ylab("Predicted Probability of Cancer") +
  ggtitle("Distribution of Predicted Cancer Probability by Cancer Type") +
  theme_minimal()

print(boxplot_prob)
ggsave("Prob_Cancer_by_Cancer_Type_validation_new_motif_methylation.pdf", boxplot_prob, width = 8, height = 6)



# 2. Filter out normal/healthy samples. 
# Adjust the filtering condition based on how the 'cancer_type' (or similar) is coded in your metadata.
non_normal_df <- merged_df %>%
  filter(cancer_type != "healthy" & cancer_type != "Normal")

# 3. Create a boxplot of the predicted probability of cancer (prob_cancer) by cancer type with jittered points.
boxplot_prob <- ggplot(merged_df, aes(x = cancer_type, y = prob_cancer)) +
  geom_boxplot(fill = "lightblue", outlier.color = "NA") +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.7, color = "darkblue") +
  coord_flip() +
  xlab("Cancer Type") +
  ylab("Predicted Probability of Cancer") +
  ggtitle("Predicted Probability of Cancer by Cancer Type (Excluding Normal)") +
  theme_classic()


print(boxplot_prob)
ggsave("Boxplot_Predicted_Probability_with_Normal_validation2_motif_methylation.pdf", boxplot_prob, width = 8, height = 5)
ggsave("Boxplot_Predicted_Probability_with_Normal_validation2_motif_methylation.png", boxplot_prob, width = 8, height = 4, dpi = 500)

# 4. Statistical Testing: Perform one-way ANOVA to assess differences among cancer types.
anova_res <- aov(prob_cancer ~ cancer_type, data = non_normal_df)
summary(anova_res)

# 5. If the ANOVA is significant, perform Tukey's HSD post hoc test.
tukey_res <- TukeyHSD(anova_res)
print(tukey_res)

# Optionally, save the Tukey results to a CSV file.
tukey_df <- as.data.frame(tukey_res$cancer_type)
write.csv(tukey_df, "TukeyHSD_Results_Cancer_Type.csv", row.names = TRUE)


# --- Summarize ANOVA Results ---
# Extract the ANOVA table from your anova_res object.
# Convert the summary ANOVA table (a matrix) to a data frame while preserving all columns
anova_table <- as.data.frame.matrix(summary(anova_res)[[1]])
# Add the row names (terms) as a new column
anova_table$Term <- rownames(anova_table)
# Reorder the columns to bring Term to the front
anova_table <- anova_table[, c("Term", "Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)")]
# Print the table to check
print(anova_table)

# Export the formatted ANOVA table to CSV
write.csv(anova_table, "ANOVA_Summary.csv", row.names = FALSE)


# --- Summarize Tukey HSD Results ---
# Convert the TukeyHSD results for the 'cancer_type' factor into a data frame.
tukey_df <- as.data.frame(tukey_res$cancer_type)
tukey_df$Comparison <- rownames(tukey_df)
rownames(tukey_df) <- NULL
# Reorder columns
tukey_df <- tukey_df[, c("Comparison", "diff", "lwr", "upr", "p adj")]
# Print the table.
kable(tukey_df, caption = "Tukey HSD Post Hoc Test Results for Differences Among Cancer Types (prob_cancer)")
# Export the Tukey results as a CSV.
write.csv(tukey_df, "TukeyHSD_Results_Cancer_Type.csv", row.names = FALSE)


# --- Create a Forest Plot for Tukey HSD Comparisons ---
# Order comparisons by the estimated difference (diff)
tukey_df$Comparison <- factor(tukey_df$Comparison, levels = tukey_df$Comparison[order(tukey_df$diff)])
tukey_plot <- ggplot(tukey_df, aes(x = diff, y = Comparison)) +
  geom_point(color = "darkblue", size = 2) +
  geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2, color = "darkblue") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  xlab("Difference in Mean (diff)") +
  ylab("Comparison") +
  ggtitle("Forest Plot: Tukey HSD Comparisons for Cancer Type Differences") +
  theme_minimal(base_size = 12)

print(tukey_plot)
ggsave("TukeyHSD_ForestPlot.pdf", tukey_plot, width = 10, height = 8)



### Other code 
# Compute the mean, standard deviation, and count of predicted probabilities for each cancer type
mean_prob_by_cancer <- non_normal_df %>%
  group_by(cancer_type) %>%
  summarise(
    mean_prob = mean(prob_cancer, na.rm = TRUE),
    sd_prob = sd(prob_cancer, na.rm = TRUE),
    n = n()
  )
print(mean_prob_by_cancer)

# --- Identify the Three Cancer Types with the Lowest Predicted Probability ---
lowest_three <- mean_prob_by_cancer %>% 
  arrange(mean_prob) %>% 
  head(3)
print(lowest_three)

# Extract the cancer type names as a vector
lowest_groups <- as.character(lowest_three$cancer_type)
cat("The three cancer types with the lowest predicted probability are:", paste(lowest_groups, collapse = ", "), "\n")

# --- Extract Significant Tukey HSD Comparisons Involving Any of the Three Lowest Groups ---
# Load the Tukey HSD results (if not already loaded)
# tukey_df <- read.csv("TukeyHSD_Results_Cancer_Type.csv", stringsAsFactors = FALSE)

# Filter for significant comparisons (p adj < 0.05)
sig_tukey <- tukey_df %>% filter(`p adj` < 0.05)

# Create a regex pattern matching any of the three lowest groups
pattern <- paste(lowest_groups, collapse = "|")

# Filter Tukey comparisons where the Comparison column contains any of the lowest groups
sig_lowest_three <- sig_tukey[grepl(pattern, sig_tukey$Comparison), ]
print(sig_lowest_three)

# Export the significant comparisons to a CSV file
write.csv(sig_lowest_three, "Significant_Tukey_Comparisons_Involving_3_Lowest_Groups.csv", row.names = FALSE)






#### Update for the MS 
### Start here


##############################################################################
# 1. Load Libraries and Define Color Palette
##############################################################################
library(dplyr)
library(ggplot2)
library(cowplot)

# Custom color palette (keys must match your actual factor levels)
cancer_type_col <- c(
  'Healthy'           = '#33a02c',
  'Brain_cancer'      = '#1f78b4',
  'Lung_cancer'       = '#b2df8a',
  'Prostate_cancer'   = '#a6cee3',
  'AML'               = '#fb9a99',
  'Pancreatic_cancer' = '#e31a1c',
  'Eye_cancer'        = '#fdbf6f',
  'Head_and_neck_cancer' = '#ff7f00',
  'Breast_cancer'     = '#cab2d6',
  'Colorectal_cancer' = '#6a3d9a',
  'Bladder_cancer'    = '#fb6a4a',
  'Renal_cancer'      = '#b15928',
  'LFS_survivor'      = '#bdbdbd',
  'LFS_previvor'      = '#969696',
  'LFS_positive'      = '#737373',
  'Liver_cancer'      = '#ffeda0',
  'Melanoma'          = '#41b6c4',
  'Mixed_cancer'      = '#e31a1c', # same as Pancreatic if truly "Mixed"
  'Ovarian_cancer'    = '#ffffb3'
)

##############################################################################
# 2. ORIGINAL DATASET: Load, Merge, Compute Accuracy, and Plot
##############################################################################

### 2.1 Load the per-sample classification table for "Original" ###
path_original <- "Final Data Dec 2024/ML_data/PE/Cancer_vs_normal/"
class_df_original <- read.csv(paste0(path_original, "classification_details_cancer_vs_healthy_Original.csv"))

# Merge with metadata; assume 'sample' matches 'sequencing_id'
merged_df_original <- class_df_original %>%
  left_join(metadata_df, by = c("sample" = "sequencing_id")) %>%
  mutate(
    Correct = (ActualClass == PredictedClass),
    # Example of renaming "Eye Cancer" -> "Uveal Melanoma" if needed
    cancer_type = ifelse(cancer_type == "Eye Cancer", "Uveal Melanoma", cancer_type)
  )

### 2.2 Compute Accuracy by Cancer Type (Original) ###
accuracy_by_cancer_original <- merged_df_original %>%
  group_by(Analysis, cancer_type, cancer_type_title_case) %>%
  summarise(Accuracy = mean(Correct, na.rm = TRUE), .groups = "drop")

### 2.3 Create Subsets & Extract Accuracy-Based Ordering ###

#### 2.3.a Original Methylation ####
accuracy_order_original_meth <- accuracy_by_cancer_original %>%
  filter(Analysis == "Methylation") %>%
  arrange(Accuracy) %>%    # ascending order; use desc(Accuracy) for descending
  pull(cancer_type)

merged_df_original_meth <- merged_df_original %>%
  filter(Analysis == "Methylation") %>%
  mutate(cancer_type = factor(cancer_type, levels = accuracy_order_original_meth))

#### 2.3.b Original Combined Ratios + Methylation ####
accuracy_order_original_ratios <- accuracy_by_cancer_original %>%
  filter(Analysis == "Combined_ratios_methylation_PCs_1pct") %>%
  arrange(Accuracy) %>%
  pull(cancer_type)

merged_df_original_ratios <- merged_df_original %>%
  filter(Analysis == "Combined_ratios_methylation_PCs_1pct") %>%
  mutate(cancer_type = factor(cancer_type, levels = accuracy_order_original_ratios))

### 2.4 Boxplots (Original) ###

# --- Methylation ---
boxplot_prob_original_meth <- ggplot(
  merged_df_original_meth,
  aes(x = cancer_type, y = prob_cancer, fill = cancer_type_title_case)
) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(width = 0.2, size = 1.2, alpha = 0.8, color = "#1a1a1a") +
  coord_flip() +
  labs(
    x = "Cancer Type",
    y = "Predicted Probability of Cancer",
    title = "Predicted Probability (Primary, Methylation)"
  ) +
  theme_classic(base_size = 12) +
  scale_fill_manual(values = cancer_type_col, guide = "none")

ggsave("Boxplot_Predicted_Probability_Original_Methylation_Ordered.pdf",
       boxplot_prob_original_meth, width = 7, height = 5, dpi = 600)

# --- Combined Ratios + Methylation ---
boxplot_prob_original_ratios <- ggplot(
  merged_df_original_ratios,
  aes(x = cancer_type, y = prob_cancer, fill = cancer_type_title_case)
) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(width = 0.2, size = 1.2, alpha = 0.8, color = "#1a1a1a") +
  coord_flip() +
  labs(
    x = "Cancer Type",
    y = "Predicted Probability of Cancer",
    title = "Predicted Probability (Primary, Ratios + Methylation)"
  ) +
  theme_classic(base_size = 12) +
  scale_fill_manual(values = cancer_type_col, guide = "none")

ggsave("Boxplot_Predicted_Probability_Original_Ratios_Methylation_Ordered.pdf",
       boxplot_prob_original_ratios, width = 7, height = 5, dpi = 600)

### 2.5 Accuracy Barplots (Original) ###

# --- Methylation ---
accuracy_plot_original_meth <- ggplot(
  accuracy_by_cancer_original %>% filter(Analysis == "Methylation"),
  aes(x = factor(cancer_type, levels = accuracy_order_original_meth),
      y = Accuracy,
      fill = cancer_type_title_case)
) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    x = "Cancer Type",
    y = "Classification Accuracy",
    title = "Classification Accuracy (Primary, Methylation)"
  ) +
  theme_classic(base_size = 12) +
  scale_fill_manual(values = cancer_type_col, guide = "none")

ggsave("Accuracy_by_Cancer_Type_Original_Methylation_Ordered.pdf",
       accuracy_plot_original_meth, width = 7, height = 5, dpi = 600)

# --- Combined Ratios + Methylation ---
accuracy_plot_original_ratios <- ggplot(
  accuracy_by_cancer_original %>% filter(Analysis == "Combined_ratios_methylation_PCs_1pct"),
  aes(x = factor(cancer_type, levels = accuracy_order_original_ratios),
      y = Accuracy,
      fill = cancer_type_title_case)
) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    x = "Cancer Type",
    y = "Classification Accuracy",
    title = "Classification Accuracy (Primary, Ratios + Methylation)"
  ) +
  theme_classic(base_size = 12) +
  scale_fill_manual(values = cancer_type_col, guide = "none")

ggsave("Accuracy_by_Cancer_Type_Original_Ratios_Methylation_Ordered.pdf",
       accuracy_plot_original_ratios, width = 7, height = 5, dpi = 600)

### 2.6 Example: Combine Methylation Boxplot & Accuracy (Side by Side) ###
accuracy_plot_original_meth_no_axis <- accuracy_plot_original_meth +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  )

combined_plot_original_meth <- plot_grid(
  boxplot_prob_original_meth,
  accuracy_plot_original_meth_no_axis,
  ncol  = 2,
  align = "v",
  axis  = "l"
)

ggsave("Combined_Plot_Original_Methylation_Ordered.png",
       combined_plot_original_meth, width = 10, height = 5, dpi = 600)


##############################################################################
# 3. VALIDATION DATASET: Load, Merge, Compute Accuracy, and Plot
##############################################################################

### 3.1 Load the classification table for "Validation" ###
path_val <- "Final Data Dec 2024/ML_data/PE/Cancer_vs_normal/"
class_df_validation <- read.csv(paste0(path_val, "classification_details_cancer_vs_healthy_Validation.csv"))

merged_df_validation <- class_df_validation %>%
  left_join(metadata_df, by = c("sample" = "sequencing_id")) %>%
  mutate(Correct = (ActualClass == PredictedClass))

### 3.2 Compute Accuracy by Cancer Type (Validation) ###
accuracy_by_cancer_validation <- merged_df_validation %>%
  group_by(Analysis, cancer_type, cancer_type_title_case) %>%
  summarise(Accuracy = mean(Correct, na.rm = TRUE), .groups = "drop")

### 3.3 Create Subsets & Extract Accuracy-Based Ordering ###

#### 3.3.a Validation Methylation ####
accuracy_order_validation_meth <- accuracy_by_cancer_validation %>%
  filter(Analysis == "Methylation") %>%
  arrange(Accuracy) %>%
  pull(cancer_type)

merged_df_validation_meth <- merged_df_validation %>%
  filter(Analysis == "Methylation") %>%
  mutate(cancer_type = factor(cancer_type, levels = accuracy_order_validation_meth))

#### 3.3.b Validation Combined Ratios + Methylation ####
accuracy_order_validation_ratios <- accuracy_by_cancer_validation %>%
  filter(Analysis == "Combined_ratios_methylation_PCs_1pct") %>%
  arrange(Accuracy) %>%
  pull(cancer_type)

merged_df_validation_ratios <- merged_df_validation %>%
  filter(Analysis == "Combined_ratios_methylation_PCs_1pct") %>%
  mutate(cancer_type = factor(cancer_type, levels = accuracy_order_validation_ratios))

### 3.4 Boxplots (Validation) ###

# --- Validation: Methylation ---
boxplot_prob_validation_meth <- ggplot(
  merged_df_validation_meth,
  aes(x = cancer_type, y = prob_cancer, fill = cancer_type_title_case)
) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(width = 0.2, size = 1.2, alpha = 0.8, color = "#1a1a1a") +
  coord_flip() +
  labs(
    x = "Cancer Type", 
    y = "Predicted Probability of Cancer",
    title = "Predicted Probability (Validation, Methylation)"
  ) +
  theme_classic(base_size = 12) +
  scale_fill_manual(values = cancer_type_col, guide = "none")

ggsave("Boxplot_Predicted_Probability_Validation_Methylation_Ordered.pdf",
       boxplot_prob_validation_meth, width = 7, height = 5, dpi = 600)

# --- Validation: Ratios + Methylation ---
boxplot_prob_validation_ratios <- ggplot(
  merged_df_validation_ratios,
  aes(x = cancer_type, y = prob_cancer, fill = cancer_type_title_case)
) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(width = 0.2, size = 1.2, alpha = 0.8, color = "#1a1a1a") +
  coord_flip() +
  labs(
    x = "Cancer Type", 
    y = "Predicted Probability of Cancer",
    title = "Predicted Probability (Validation, Ratios + Methylation)"
  ) +
  theme_classic(base_size = 12) +
  scale_fill_manual(values = cancer_type_col, guide = "none")

ggsave("Boxplot_Predicted_Probability_Validation_Ratios_Methylation_Ordered.pdf",
       boxplot_prob_validation_ratios, width = 7, height = 5, dpi = 600)

### 3.5 Accuracy Barplots (Validation) ###

# --- Methylation ---
accuracy_plot_validation_meth <- ggplot(
  accuracy_by_cancer_validation %>% filter(Analysis == "Methylation"),
  aes(x = factor(cancer_type, levels = accuracy_order_validation_meth),
      y = Accuracy,
      fill = cancer_type_title_case)
) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    x = "Cancer Type",
    y = "Classification Accuracy",
    title = "Classification Accuracy (Validation, Methylation)"
  ) +
  theme_classic(base_size = 12) +
  scale_fill_manual(values = cancer_type_col, guide = "none")

ggsave("Accuracy_by_Cancer_Type_Validation_Methylation_Ordered.pdf",
       accuracy_plot_validation_meth, width = 7, height = 5, dpi = 600)

# --- Ratios + Methylation ---
accuracy_plot_validation_ratios <- ggplot(
  accuracy_by_cancer_validation %>% filter(Analysis == "Combined_ratios_methylation_PCs_1pct"),
  aes(x = factor(cancer_type, levels = accuracy_order_validation_ratios),
      y = Accuracy,
      fill = cancer_type_title_case)
) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    x = "Cancer Type",
    y = "Classification Accuracy",
    title = "Classification Accuracy (Validation, Ratios + Methylation)"
  ) +
  theme_classic(base_size = 12) +
  scale_fill_manual(values = cancer_type_col, guide = "none")

ggsave("Accuracy_by_Cancer_Type_Validation_Ratios_Methylation_Ordered.pdf",
       accuracy_plot_validation_ratios, width = 7, height = 5, dpi = 600)


##############################################################################
# Additional Code for End_motifs Analysis
##############################################################################

### Original Dataset: End Motifs ###

# Create ordering vector for End_motifs (Original)
accuracy_order_original_endmotifs <- accuracy_by_cancer_original %>%
  filter(Analysis == "End_motifs") %>%
  arrange(Accuracy) %>%
  pull(cancer_type)

# Accuracy Barplot: Original End Motifs
accuracy_plot_original_endmotifs <- ggplot(
  accuracy_by_cancer_original %>% filter(Analysis == "End_motifs"),
  aes(x = factor(cancer_type, levels = accuracy_order_original_endmotifs),
      y = Accuracy,
      fill = cancer_type_title_case)
) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    x = "Cancer Type",
    y = "Classification Accuracy",
    title = "Classification Accuracy (Primary, End Motifs)"
  ) +
  theme_classic(base_size = 12) +
  scale_fill_manual(values = cancer_type_col, guide = "none")

ggsave("Accuracy_by_Cancer_Type_Original_End_Motifs_Ordered.pdf",
       accuracy_plot_original_endmotifs, width = 7, height = 5, dpi = 600)


### Validation Dataset: End Motifs ###

# Create ordering vector for End_motifs (Validation)
accuracy_order_validation_endmotifs <- accuracy_by_cancer_validation %>%
  filter(Analysis == "End_motifs") %>%
  arrange(Accuracy) %>%
  pull(cancer_type)

# Subset merged validation data for End_motifs
merged_df_validation_endmotifs <- merged_df_validation %>%
  filter(Analysis == "End_motifs") %>%
  mutate(cancer_type = factor(cancer_type, levels = accuracy_order_validation_endmotifs))

# Accuracy Barplot: Validation End Motifs
accuracy_plot_validation_endmotifs <- ggplot(
  accuracy_by_cancer_validation %>% filter(Analysis == "End_motifs"),
  aes(x = factor(cancer_type, levels = accuracy_order_validation_endmotifs),
      y = Accuracy,
      fill = cancer_type_title_case)
) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    x = "Cancer Type",
    y = "Classification Accuracy",
    title = "Classification Accuracy (Validation, End Motifs)"
  ) +
  theme_classic(base_size = 12) +
  scale_fill_manual(values = cancer_type_col, guide = "none")

ggsave("Accuracy_by_Cancer_Type_Validation_End_Motifs_Ordered.pdf",
       accuracy_plot_validation_endmotifs, width = 7, height = 5, dpi = 600)

# Boxplot: Validation Predicted Probability for End Motifs
boxplot_prob_validation_endmotifs <- ggplot(
  merged_df_validation_endmotifs,
  aes(x = cancer_type, y = prob_cancer, fill = cancer_type_title_case)
) +
  geom_boxplot(outlier.color = NA) +
  scale_y_continuous(limits = c(0, 1)) +
  geom_jitter(width = 0.2, size = 1.2, alpha = 0.8, color = "#1a1a1a") +
  coord_flip() +
  labs(
    x = "Cancer Type", 
    y = "Predicted Probability of Cancer",
    title = "Predicted Probability (Validation, End Motifs)"
  ) +
  theme_classic(base_size = 12) +
  scale_fill_manual(values = cancer_type_col, guide = "none")

ggsave("Boxplot_Predicted_Probability_Validation_End_Motifs_Ordered_updated.pdf",
       boxplot_prob_validation_endmotifs, width = 7, height = 5, dpi = 600)

# --- Primary ---
### ORIGINAL DATASET: End Motifs ###

# Create ordering vector for End_motifs (Original)
accuracy_order_original_endmotifs <- accuracy_by_cancer_original %>%
  filter(Analysis == "End_motifs") %>%
  arrange(Accuracy) %>%
  pull(cancer_type)

# Accuracy Barplot: Original End Motifs
accuracy_plot_original_endmotifs <- ggplot(
  accuracy_by_cancer_original %>% filter(Analysis == "End_motifs"),
  aes(x = factor(cancer_type, levels = accuracy_order_original_endmotifs),
      y = Accuracy,
      fill = cancer_type_title_case)
) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    x = "Cancer Type",
    y = "Classification Accuracy",
    title = "Classification Accuracy (Original, End Motifs)"
  ) +
  theme_classic(base_size = 12) +
  scale_fill_manual(values = cancer_type_col, guide = "none")

ggsave("Accuracy_by_Cancer_Type_Original_End_Motifs_Ordered.pdf",
       accuracy_plot_original_endmotifs, width = 7, height = 5, dpi = 600)

# Subset original data for End_motifs boxplot
merged_df_original_endmotifs <- merged_df_original %>% 
  filter(Analysis == "End_motifs") %>% 
  mutate(cancer_type = factor(cancer_type, levels = accuracy_order_original_endmotifs))

# Boxplot: Original Predicted Probability for End Motifs
boxplot_prob_original_endmotifs <- ggplot(
  merged_df_original_endmotifs,
  aes(x = cancer_type, y = prob_cancer, fill = cancer_type_title_case)
) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(width = 0.2, size = 1.2, alpha = 0.8, color = "#1a1a1a") +
  coord_flip() +
  labs(
    x = "Cancer Type", 
    y = "Predicted Probability of Cancer",
    title = "Predicted Probability (Original, End Motifs)"
  ) +
  theme_classic(base_size = 12) +
  scale_fill_manual(values = cancer_type_col, guide = "none")

ggsave("Boxplot_Predicted_Probability_Original_End_Motifs_Ordered.pdf",
       boxplot_prob_original_endmotifs, width = 7, height = 5, dpi = 600)




##############################################################################
# 4. (Optional) Export Accuracy Tables
##############################################################################
accuracy_by_cancer_original_export <- accuracy_by_cancer_original %>%
  mutate(PercentCorrect = Accuracy * 100)
write.csv(accuracy_by_cancer_original_export, 
          "Accuracy_By_Cancer_Type_Original.csv", 
          row.names = FALSE)

accuracy_by_cancer_validation_export <- accuracy_by_cancer_validation %>%
  mutate(PercentCorrect = Accuracy * 100)
write.csv(accuracy_by_cancer_validation_export, 
          "Accuracy_By_Cancer_Type_Validation.csv", 
          row.names = FALSE)



#### As loop
# Define a function that creates the plots and exports the accuracy table
plot_feature <- function(feature_name) {
  # ----- Primary (Original) Dataset -----
  # Subset the merged data and compute ordering based on accuracy for the feature
  df_orig <- merged_df_original %>% filter(Analysis == feature_name)
  order_orig <- accuracy_by_cancer_original %>% 
    filter(Analysis == feature_name) %>% 
    arrange(Accuracy) %>% 
    pull(cancer_type)
  df_orig <- df_orig %>% mutate(cancer_type = factor(cancer_type, levels = order_orig))
  
  # Boxplot for Primary Dataset
  bp_orig <- ggplot(df_orig, aes(x = cancer_type, y = prob_cancer, fill = cancer_type_title_case)) +
    geom_boxplot(outlier.color = NA) +
    geom_jitter(width = 0.2, size = 1.2, alpha = 0.8, color = "#1a1a1a") +
    coord_flip() +
    labs(x = "Cancer Type", 
         y = "Predicted Probability of Cancer", 
         title = paste("Predicted Probability (Primary,", feature_name, ")")) +
    theme_classic(base_size = 12) +
    scale_fill_manual(values = cancer_type_col, guide = "none")
  
  ggsave(paste0("Boxplot_Predicted_Probability_Original_", feature_name, "_Ordered.pdf"),
         bp_orig, width = 7, height = 5, dpi = 600)
  
  # Accuracy Barplot for Primary Dataset
  acc_plot_orig <- ggplot(accuracy_by_cancer_original %>% filter(Analysis == feature_name), 
                          aes(x = factor(cancer_type, levels = order_orig),
                              y = Accuracy,
                              fill = cancer_type_title_case)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(x = "Cancer Type", 
         y = "Classification Accuracy", 
         title = paste("Classification Accuracy (Primary,", feature_name, ")")) +
    theme_classic(base_size = 12) +
    scale_fill_manual(values = cancer_type_col, guide = "none")
  
  ggsave(paste0("Accuracy_by_Cancer_Type_Original_", feature_name, "_Ordered.pdf"),
         acc_plot_orig, width = 7, height = 5, dpi = 600)
  
  # Export accuracy table for Primary Dataset
  accuracy_table_orig <- accuracy_by_cancer_original %>% 
    filter(Analysis == feature_name) %>% 
    mutate(PercentCorrect = Accuracy * 100)
  
  write.csv(accuracy_table_orig, 
            paste0("Accuracy_Table_Original_", feature_name, ".csv"), 
            row.names = FALSE)
  
  # ----- Validation Dataset -----
  # Subset the merged data and compute ordering based on accuracy for the feature
  df_val <- merged_df_validation %>% filter(Analysis == feature_name)
  order_val <- accuracy_by_cancer_validation %>% 
    filter(Analysis == feature_name) %>% 
    arrange(Accuracy) %>% 
    pull(cancer_type)
  df_val <- df_val %>% mutate(cancer_type = factor(cancer_type, levels = order_val))
  
  # Boxplot for Validation Dataset
  bp_val <- ggplot(df_val, aes(x = cancer_type, y = prob_cancer, fill = cancer_type_title_case)) +
    geom_boxplot(outlier.color = NA) +
    geom_jitter(width = 0.2, size = 1.2, alpha = 0.8, color = "#1a1a1a") +
    coord_flip() +
    labs(x = "Cancer Type", 
         y = "Predicted Probability of Cancer", 
         title = paste("Predicted Probability (Validation,", feature_name, ")")) +
    theme_classic(base_size = 12) +
    scale_fill_manual(values = cancer_type_col, guide = "none")
  
  ggsave(paste0("Boxplot_Predicted_Probability_Validation_", feature_name, "_Ordered.pdf"),
         bp_val, width = 7, height = 5, dpi = 600)
  
  # Accuracy Barplot for Validation Dataset
  acc_plot_val <- ggplot(accuracy_by_cancer_validation %>% filter(Analysis == feature_name), 
                         aes(x = factor(cancer_type, levels = order_val),
                             y = Accuracy,
                             fill = cancer_type_title_case)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(x = "Cancer Type", 
         y = "Classification Accuracy", 
         title = paste("Classification Accuracy (Validation,", feature_name, ")")) +
    theme_classic(base_size = 12) +
    scale_fill_manual(values = cancer_type_col, guide = "none")
  
  ggsave(paste0("Accuracy_by_Cancer_Type_Validation_", feature_name, "_Ordered.pdf"),
         acc_plot_val, width = 7, height = 5, dpi = 600)
  
  # Export accuracy table for Validation Dataset
  accuracy_table_val <- accuracy_by_cancer_validation %>% 
    filter(Analysis == feature_name) %>% 
    mutate(PercentCorrect = Accuracy * 100)
  
  write.csv(accuracy_table_val, 
            paste0("Accuracy_Table_Validation_", feature_name, ".csv"), 
            row.names = FALSE)
}

# Now, if you want to run the function for a set of features, you can loop over them:
features_to_plot <- unique(merged_df_original$Analysis)

for(feature in features_to_plot){
  plot_feature(feature)
}

# Save as RDS
saveRDS(merged_df_validation, "merged_df_validation.rds")
saveRDS(merged_df_original, "merged_df_original.rds")

# Save as CSV
write.csv(merged_df_validation, "merged_df_validation.csv", row.names = FALSE)
write.csv(merged_df_original, "merged_df_original.csv", row.names = FALSE)

