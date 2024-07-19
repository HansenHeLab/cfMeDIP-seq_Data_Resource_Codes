# Additional stats on merged zscore df - June 2024 
library(tidyverse)
# Select relevant columns
selected_data <- Merged_zscore_df %>%
  select(cancer_type_corrected_updated, Endmotif_zscore, Insertsize_zscore, Delfi_zscore, Nucleosome_peak_zscore)



# Define a function to calculate summary statistics
calculate_stats <- function(df, zscore_col) {
  df %>%
    summarize(
      mean = mean(.data[[zscore_col]], na.rm = TRUE),
      min = min(.data[[zscore_col]], na.rm = TRUE),
      max = max(.data[[zscore_col]], na.rm = TRUE),
      median = median(.data[[zscore_col]], na.rm = TRUE),
      sd = sd(.data[[zscore_col]], na.rm = TRUE)
    )
}

# Calculate summary statistics for each z-score and combine into one table
summary_table <- selected_data %>%
  group_by(cancer_type_corrected_updated) %>%
  summarize(
    Endmotif_mean = calculate_stats(cur_data(), "Endmotif_zscore")$mean,
    Endmotif_min = calculate_stats(cur_data(), "Endmotif_zscore")$min,
    Endmotif_max = calculate_stats(cur_data(), "Endmotif_zscore")$max,
    Endmotif_median = calculate_stats(cur_data(), "Endmotif_zscore")$median,
    Endmotif_sd = calculate_stats(cur_data(), "Endmotif_zscore")$sd,
    Insertsize_mean = calculate_stats(cur_data(), "Insertsize_zscore")$mean,
    Insertsize_min = calculate_stats(cur_data(), "Insertsize_zscore")$min,
    Insertsize_max = calculate_stats(cur_data(), "Insertsize_zscore")$max,
    Insertsize_median = calculate_stats(cur_data(), "Insertsize_zscore")$median,
    Insertsize_sd = calculate_stats(cur_data(), "Insertsize_zscore")$sd,
    Delfi_mean = calculate_stats(cur_data(), "Delfi_zscore")$mean,
    Delfi_min = calculate_stats(cur_data(), "Delfi_zscore")$min,
    Delfi_max = calculate_stats(cur_data(), "Delfi_zscore")$max,
    Delfi_median = calculate_stats(cur_data(), "Delfi_zscore")$median,
    Delfi_sd = calculate_stats(cur_data(), "Delfi_zscore")$sd,
    Nucleosome_peak_mean = calculate_stats(cur_data(), "Nucleosome_peak_zscore")$mean,
    Nucleosome_peak_min = calculate_stats(cur_data(), "Nucleosome_peak_zscore")$min,
    Nucleosome_peak_max = calculate_stats(cur_data(), "Nucleosome_peak_zscore")$max,
    Nucleosome_peak_median = calculate_stats(cur_data(), "Nucleosome_peak_zscore")$median,
    Nucleosome_peak_sd = calculate_stats(cur_data(), "Nucleosome_peak_zscore")$sd
  )

## Get overall 
# Convert to long format
long_data <- selected_data %>%
  pivot_longer(cols = starts_with("Endmotif_zscore") | starts_with("Insertsize_zscore") | starts_with("Delfi_zscore") | starts_with("Nucleosome_peak_zscore"),
               names_to = "feature",
               values_to = "zscore")

# Calculate summary statistics for each z-score and overall
summary_table_2 <- long_data %>%
  group_by(cancer_type_corrected_updated) %>%
  summarize(
    Endmotif_mean = mean(zscore[feature == "Endmotif_zscore"], na.rm = TRUE),
    Endmotif_min = min(zscore[feature == "Endmotif_zscore"], na.rm = TRUE),
    Endmotif_max = max(zscore[feature == "Endmotif_zscore"], na.rm = TRUE),
    Endmotif_median = median(zscore[feature == "Endmotif_zscore"], na.rm = TRUE),
    Endmotif_sd = sd(zscore[feature == "Endmotif_zscore"], na.rm = TRUE),
    Insertsize_mean = mean(zscore[feature == "Insertsize_zscore"], na.rm = TRUE),
    Insertsize_min = min(zscore[feature == "Insertsize_zscore"], na.rm = TRUE),
    Insertsize_max = max(zscore[feature == "Insertsize_zscore"], na.rm = TRUE),
    Insertsize_median = median(zscore[feature == "Insertsize_zscore"], na.rm = TRUE),
    Insertsize_sd = sd(zscore[feature == "Insertsize_zscore"], na.rm = TRUE),
    Delfi_mean = mean(zscore[feature == "Delfi_zscore"], na.rm = TRUE),
    Delfi_min = min(zscore[feature == "Delfi_zscore"], na.rm = TRUE),
    Delfi_max = max(zscore[feature == "Delfi_zscore"], na.rm = TRUE),
    Delfi_median = median(zscore[feature == "Delfi_zscore"], na.rm = TRUE),
    Delfi_sd = sd(zscore[feature == "Delfi_zscore"], na.rm = TRUE),
    Nucleosome_peak_mean = mean(zscore[feature == "Nucleosome_peak_zscore"], na.rm = TRUE),
    Nucleosome_peak_min = min(zscore[feature == "Nucleosome_peak_zscore"], na.rm = TRUE),
    Nucleosome_peak_max = max(zscore[feature == "Nucleosome_peak_zscore"], na.rm = TRUE),
    Nucleosome_peak_median = median(zscore[feature == "Nucleosome_peak_zscore"], na.rm = TRUE),
    Nucleosome_peak_sd = sd(zscore[feature == "Nucleosome_peak_zscore"], na.rm = TRUE),
    overall_mean = mean(zscore, na.rm = TRUE),
    overall_min = min(zscore, na.rm = TRUE),
    overall_max = max(zscore, na.rm = TRUE),
    overall_median = median(zscore, na.rm = TRUE),
    overall_sd = sd(zscore, na.rm = TRUE)
  )

# Print the summary table
print(summary_table)



# Get feature summary 
# Calculate summary statistics for each feature and cancer type
feature_summary <- long_data %>%
  group_by(cancer_type_corrected_updated, feature) %>%
  summarise(
    mean = mean(zscore, na.rm = TRUE),
    min = min(zscore, na.rm = TRUE),
    max = max(zscore, na.rm = TRUE),
    median = median(zscore, na.rm = TRUE),
    sd = sd(zscore, na.rm = TRUE)
  ) %>%
  pivot_wider(names_from = feature, values_from = c(mean, min, max, median, sd))

overall_summary <- long_data %>%
  group_by(cancer_type_corrected_updated) %>%
  summarise(overall_mean = mean(zscore, na.rm = TRUE),
            overall_min = min(zscore, na.rm = TRUE),
            overall_max = max(zscore, na.rm = TRUE),
            overall_median = median(zscore, na.rm = TRUE),
            overall_sd = sd(zscore, na.rm = TRUE))

## Get based on feature 
mean_zscores <- Merged_zscore_df %>% filter(cancer_type_corrected_updated != "healthy") %>%
  summarise(mean_Endmotif_zscore = mean(Endmotif_zscore, na.rm = TRUE),
            mean_Insertsize_zscore = mean(Insertsize_zscore, na.rm = TRUE),
            mean_Delfi_zscore = mean(Delfi_zscore, na.rm = TRUE),
            mean_Nucleosome_peak_zscore = mean(Nucleosome_peak_zscore, na.rm = TRUE))

## Now see if sd is significant 
# Group by cancer type and calculate standard deviations for each feature and overall
sd_table <- long_data %>%
  group_by(cancer_type_corrected_updated) %>%
  summarize(
    Endmotif_sd = sd(zscore[feature == "Endmotif_zscore"], na.rm = TRUE),
    Insertsize_sd = sd(zscore[feature == "Insertsize_zscore"], na.rm = TRUE),
    Delfi_sd = sd(zscore[feature == "Delfi_zscore"], na.rm = TRUE),
    Nucleosome_peak_sd = sd(zscore[feature == "Nucleosome_peak_zscore"], na.rm = TRUE),
    overall_sd = sd(zscore, na.rm = TRUE)
  )


healthy_samples <- long_data %>%
  filter(cancer_type_corrected_updated == "healthy")

# Perform Levene's test
library(car)
levene_test <- leveneTest(zscore ~ cancer_type_corrected_updated, data = long_data)

# Summarize the results
summary_report <- sd_table %>%
  summarize(
    healthy_sd = overall_sd[cancer_type_corrected_updated == "healthy"],
    mean_cancer_sd = mean(overall_sd[cancer_type_corrected_updated != "healthy"], na.rm = TRUE),
    min_cancer_sd = min(overall_sd[cancer_type_corrected_updated != "healthy"], na.rm = TRUE),
    max_cancer_sd = max(overall_sd[cancer_type_corrected_updated != "healthy"], na.rm = TRUE)
  )

# Print the summary report
print(summary_report)

# Create a narrative based on the results
narrative <- paste(
  "Additionally, when considering these four fragmentomic features together, we observed that healthy samples exhibited greater consistency compared to other cancer types.",
  "This was quantified by analyzing the variation in z-scores across the cohort, where healthy samples showed lower standard deviations in z-scores (SD =",
  round(summary_report$healthy_sd, 3), ") than those observed in the cancer samples (Mean SD =",
  round(summary_report$mean_cancer_sd, 3), ", Min SD =",
  round(summary_report$min_cancer_sd, 3), ", Max SD =",
  round(summary_report$max_cancer_sd, 3), ").",
  "Levene's test confirmed that the differences in variance between healthy samples and other cancer types were statistically significant (p-value =",
  round(levene_test$`Pr(>F)`[1], 3), ")."
)

# Print the narrative
cat(narrative)




## See which mean higgest for which feature 
# Find the highest mean for each cancer type
mean_columns <- c("Endmotif_mean", "Insertsize_mean", "Delfi_mean", "Nucleosome_peak_mean", "overall_mean")
summary_table_mean <- summary_table_2 %>%
  rowwise() %>%
  mutate(Highest_mean_feature = mean_columns[which.max(c(Endmotif_mean, Insertsize_mean, Delfi_mean, Nucleosome_peak_mean, overall_mean))],
         Highest_mean_value = max(c(Endmotif_mean, Insertsize_mean, Delfi_mean, Nucleosome_peak_mean, overall_mean)))

# Select relevant columns for the result
result <- summary_table_mean %>%
  select(cancer_type_corrected_updated, Highest_mean_feature, Highest_mean_value)



### Testing
# Gather the data to long format for easier summarization
long_data <- selected_data %>%
  gather(key = "Method", value = "zscore", -cancer_type_corrected_updated)

# Summarize mean, range, median, and standard deviation of z-scores for each cancer type and method
summary_stats <- long_data %>%
  group_by(cancer_type_corrected_updated, Method) %>%
  summarize(
    mean_zscore = mean(zscore, na.rm = TRUE),
    range_zscore = range(zscore, na.rm = TRUE),
    median_zscore = median(zscore, na.rm = TRUE),
    sd_zscore = sd(zscore, na.rm = TRUE)
  )


# Identify the method that was best for the most samples based on the highest mean z-score
best_method <- summary_stats %>%
  group_by(Method, cancer_type_corrected_updated) %>%
  summarize(mean_zscore = mean(mean_zscore, na.rm = TRUE)) %>%
  arrange(desc(mean_zscore)) %>%
  slice(1)

