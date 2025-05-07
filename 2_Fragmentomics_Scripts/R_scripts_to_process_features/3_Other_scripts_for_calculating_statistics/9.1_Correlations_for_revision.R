# ----------------------------------------------------------------------------
# Title   : Correlation and Regression Analyses of Fragmentomic Metrics
# Authors : Dory Abelman
# Date    : March 2025
#
# Purpose :
#   1. Assess pairwise Pearson correlations between each fragmentation metric
#      (end‑motif, insert‑size, DELFI, nucleosome‑peak Z‑scores, percent short
#      fragments, fragmentation score) and patient age.
#   2. Visualize these relationships via scatterplots with fitted regression lines
#      and annotated correlation statistics.
#   3. Compare each metric by sex using boxplots and two‑sample t‑tests.
#   4. Extend comparisons across cancer subtypes:
#        • Boxplots for each metric stratified by cancer type with ANOVA.
#        • Multiple linear regression for each metric including age, sex, and
#          cancer subtype; extract coefficient estimates, confidence intervals,
#          and diagnostic plots.
#   5. Summarize and save all model outputs, diagnostic figures, and coefficient
#      tables for downstream reporting.
#
# ----------------------------------------------------------------------------

# Load necessary libraries
library(tidyverse)
library(ggpubr)
library(gridExtra)  # for arranging multiple plots

# (Optional) Ensure that 'sex' is a factor
Merged_zscore_df$sex <- as.factor(Merged_zscore_df$sex)

# Define the set of variables of interest:
# All columns with "zscore" in their name, plus "Twentyto150_pct" and "FS"
zscore_vars <- c("Endmotif_zscore", "Insertsize_zscore", "Delfi_zscore", 
                 "Nucleosome_peak_zscore", "Twentyto150_pct", "FS")

# 1. Correlation analysis with age -------------------------------------------

# Create an empty data frame to store correlation results
cor_results <- data.frame(variable = zscore_vars,
                          correlation = NA,
                          p_value = NA,
                          stringsAsFactors = FALSE)

# Loop over variables and perform Pearson correlation with age
for (i in seq_along(zscore_vars)) {
  var <- zscore_vars[i]
  test <- cor.test(Merged_zscore_df$age, Merged_zscore_df[[var]], method = "pearson")
  cor_results$correlation[i] <- test$estimate
  cor_results$p_value[i] <- test$p.value
}
# Print correlation results
print(cor_results)

# 2. Scatter plots for age correlations ---------------------------------------

# Create a list to store individual plots
age_plots <- list()

for (var in zscore_vars) {
  p <- ggplot(Merged_zscore_df, aes_string(x = "age", y = var)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "blue") +
    # Annotate with correlation coefficient and p-value (from ggpubr's stat_cor)
    stat_cor(method = "pearson",
             label.x = min(Merged_zscore_df$age, na.rm = TRUE) + 1,
             label.y = max(Merged_zscore_df[[var]], na.rm = TRUE) - 1,
             size = 5) +
    labs(title = paste("Correlation between", var, "and Age"),
         x = "Age",
         y = var) +
    theme_classic(base_size = 14)
  
  age_plots[[var]] <- p
}

# Arrange and display all age correlation plots in a grid
age_plot_grid <- grid.arrange(grobs = age_plots, ncol = 2)
# Optionally, save the figure:
ggsave("Age_Correlation_Plots.pdf", age_plot_grid, width = 14, height = 10)

# 3. Boxplots for sex comparisons ---------------------------------------------

# Create a temporary dataframe that excludes rows with NA in sex
df_sex <- Merged_zscore_df %>% filter(!is.na(sex))

# (Optional) Ensure that 'sex' is a factor in the temporary dataframe
df_sex$sex <- as.factor(df_sex$sex)

# Define the set of variables of interest (zscore columns, Twentyto150_pct, and FS)
zscore_vars <- c("Endmotif_zscore", "Insertsize_zscore", "Delfi_zscore", 
                 "Nucleosome_peak_zscore", "Twentyto150_pct", "FS")

# Create a list to store boxplots by sex using the filtered dataframe
sex_plots <- list()

for (var in zscore_vars) {
  p <- ggplot(df_sex, aes_string(x = "sex", y = var, fill = "sex")) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
    stat_compare_means(method = "t.test", label.y = max(df_sex[[var]], na.rm = TRUE) * 1.05) +
    labs(title = paste("Comparison of", var, "by Sex"),
         x = "Sex",
         y = var) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none")
  
  sex_plots[[var]] <- p
}

# Arrange and display all sex comparison plots in a grid
sex_plot_grid <- grid.arrange(grobs = sex_plots, ncol = 2)
# Optionally, save the figure:
ggsave("Sex_Comparison_Plots.pdf", sex_plot_grid, width = 14, height = 10)



#### Now consider cancer type 

# ---------------------------
# Comprehensive Analysis Script
# ---------------------------
# Load necessary libraries
library(tidyverse)
library(ggpubr)
library(gridExtra)  # for arranging plots
library(broom)      # for tidying model outputs
library(car)        # for VIF (multicollinearity check)

# ---------------------------
# 1. Data Preparation
# ---------------------------
# Convert key variables to factors
tmp <- Merged_zscore_df %>%
  mutate(
    sex = as.factor(sex),
    cancer_type_corrected_updated = as.factor(cancer_type_corrected_updated)
  )

# Create a filtered dataframe excluding rows with NA in key predictors (age, sex, and cancer type)
df_complete <- tmp %>% 
  filter(!is.na(age) & !is.na(sex) & !is.na(cancer_type_corrected_updated))

# Define the response variables of interest
response_vars <- c("Endmotif_zscore", "Insertsize_zscore", "Delfi_zscore", 
                   "Nucleosome_peak_zscore", "Twentyto150_pct", "FS")

# ---------------------------
# 2. Descriptive Analysis & Visualizations
# ---------------------------
# A. Boxplots by cancer type
cancer_plots <- list()
for (var in response_vars) {
  p <- ggplot(df_complete, aes_string(x = "cancer_type_corrected_updated", y = var, fill = "cancer_type_corrected_updated")) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
    stat_compare_means(method = "anova", label.y = max(df_complete[[var]], na.rm = TRUE)*1.05) +
    labs(title = paste("Comparison of", var, "by Cancer Type"),
         x = "Cancer Type", y = var) +
    theme_classic(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.position = "none")
  
  cancer_plots[[var]] <- p
}
# Arrange and display cancer type plots in a grid
cancer_plot_grid <- grid.arrange(grobs = cancer_plots, ncol = 2)
# Optionally save the plot grid:
ggsave("CancerType_Comparison_Plots.pdf", cancer_plot_grid, width = 14, height = 13)

# B. Scatter plots vs Age (as before)
age_plots <- list()
for (var in response_vars) {
  p <- ggplot(df_complete, aes_string(x = "age", y = var)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "blue") +
    stat_cor(method = "pearson",
             label.x = min(df_complete$age, na.rm = TRUE) + 1,
             label.y = max(df_complete[[var]], na.rm = TRUE) - 1,
             size = 5) +
    labs(title = paste("Correlation between", var, "and Age"),
         x = "Age", y = var) +
    theme_classic(base_size = 14)
  age_plots[[var]] <- p
}
age_plot_grid <- grid.arrange(grobs = age_plots, ncol = 2)
ggsave("Age_Correlation_Plots2.pdf", age_plot_grid, width = 14, height = 10)

# C. Boxplots by Sex (using filtered data)
sex_plots <- list()
for (var in response_vars) {
  p <- ggplot(df_complete, aes_string(x = "sex", y = var, fill = "sex")) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
    stat_compare_means(method = "t.test", label.y = max(df_complete[[var]], na.rm = TRUE) * 1.05) +
    labs(title = paste("Comparison of", var, "by Sex"),
         x = "Sex", y = var) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none")
  sex_plots[[var]] <- p
}
sex_plot_grid <- grid.arrange(grobs = sex_plots, ncol = 2)
ggsave("Sex_Comparison_Plot2s.pdf", sex_plot_grid, width = 14, height = 12)

# ---------------------------
# 3. Multiple Regression Analysis
# ---------------------------
# We now fit a linear model for each response variable with age, sex, and cancer type as predictors.
# This will help us determine whether the effect of cancer type remains significant
# after accounting for age and sex.
model_results <- list()  # to store model summaries
coeff_plots <- list()    # to store coefficient plots

for (var in response_vars) {
  # Construct the regression formula (response ~ age + sex + cancer_type)
  formula_str <- as.formula(paste(var, "~ age + sex + cancer_type_corrected_updated"))
  
  # Fit the linear model
  model <- lm(formula_str, data = df_complete)
  
  # Print and store the model summary
  cat("\n\n==============================\n")
  cat("Results for", var, "\n")
  model_summary <- summary(model)
  print(model_summary)
  model_results[[var]] <- model_summary
  
  # Check multicollinearity with Variance Inflation Factor (VIF)
  cat("\nVariance Inflation Factors:\n")
  print(vif(model))
  
  # Save diagnostic plots (residuals, Q-Q, etc.) directly to a PDF file
  pdf(paste0("Diagnostic_Plots_", var, ".pdf"), width = 8, height = 8)
  par(mfrow = c(2,2), mar = c(4, 4, 2, 1))  # set margins appropriately
  plot(model)
  dev.off()
  
  # Create a coefficient plot using broom::tidy() for a clean summary
  tidy_mod <- tidy(model)
  
  # Exclude the intercept for plotting
  cp <- tidy_mod %>% 
    filter(term != "(Intercept)") %>%
    ggplot(aes(x = term, y = estimate)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), 
                  width = 0.2, color = "darkblue") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(title = paste("Coefficient Estimates for", var),
         x = "Predictor", y = "Estimate") +
    theme_bw(base_size = 14) +
    coord_flip()  # flip for easier reading
  
  print(cp)
  coeff_plots[[var]] <- cp
  
  # Save each coefficient plot to file
  ggsave(paste0("Coeff_Plot_", var, ".pdf"), cp, width = 12, height = 4)
}


### Updated 
# Create a named vector to map variable names to desired labels
var_labels <- c(
  "Endmotif_zscore"       = "End Motifs",
  "Insertsize_zscore"     = "Insert Size",
  "Delfi_zscore"          = "Fragment Ratios",
  "Nucleosome_peak_zscore"= "Nucleosome Peak",
  "Twentyto150_pct"       = "Prop Short Frags",
  "FS"                    = "Fragment Size Score"
)

# If you haven’t already done so, define your vector of response variables
zscore_vars <- c("Endmotif_zscore", "Insertsize_zscore", "Delfi_zscore",
                 "Nucleosome_peak_zscore", "Twentyto150_pct", "FS")

model_results <- list()  # to store model summaries
coeff_plots   <- list()  # to store coefficient plots

for (var in zscore_vars) {
  # Construct the regression formula
  formula_str <- as.formula(paste(var, "~ age + sex + cancer_type_corrected_updated"))
  
  # Fit the linear model
  model <- lm(formula_str, data = df_complete)
  
  # Print and store the model summary
  cat("\n\n==============================\n")
  cat("Results for", var, "\n")
  model_summary <- summary(model)
  print(model_summary)
  model_results[[var]] <- model_summary
  
  # Check multicollinearity
  cat("\nVariance Inflation Factors:\n")
  print(vif(model))
  
  # Save diagnostic plots (residuals, Q-Q, etc.) directly to a PDF
  pdf(paste0("Diagnostic_Plots_", var, ".pdf"), width = 8, height = 8)
  par(mfrow = c(2,2), mar = c(4, 4, 2, 1))
  plot(model)
  dev.off()
  
  # Tidy the model output
  tidy_mod <- tidy(model)
  
  # Save tidy model results to CSV
  write.csv(tidy_mod, paste0("Model_Table_", var, ".csv"), row.names = FALSE)
  
  # Create coefficient plot
  cp <- tidy_mod %>% 
    filter(term != "(Intercept)") %>%
    # Remove the "cancer_type_corrected_updated" prefix from term labels
    mutate(term = str_remove(term, "^cancer_type_corrected_updated")) %>%
    ggplot(aes(x = term, y = estimate)) +
    geom_point(size = 3) +
    geom_errorbar(
      aes(ymin = estimate - std.error, ymax = estimate + std.error), 
      width = 0.2, color = "darkblue"
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    # Use your desired label for the plot title
    labs(
      title = paste("Coefficient Estimates for", var_labels[[var]]),
      x = "Predictor", 
      y = "Estimate"
    ) +
    theme_bw(base_size = 14) +
    coord_flip()  # flip for easier reading
  
  print(cp)
  coeff_plots[[var]] <- cp
  
  # Save each coefficient plot to file
  ggsave(paste0("Coeff_Plot_", var, "updated.pdf"), cp, width = 12, height = 4)
}


# ---------------------------
# 4. Interpretation Considerations
# ---------------------------
# The outputs above allow you to examine:
# 1. Whether cancer_type_corrected_updated is a significant predictor
#    of the response variable even after controlling for age and sex.
# 2. The magnitude and direction of the effect from cancer type versus age and sex.
# 3. Model diagnostics (residual plots and VIF) to verify that assumptions are met.
#
# In a publication, one would discuss:
# - The significance levels (p-values) for cancer_type_corrected_updated versus age and sex.
# - The model's adjusted R-squared, which shows how much variation is explained by these predictors.
# - Diagnostic plots ensuring that residuals are approximately normal and that there is no undue influence of multicollinearity.
#
# If cancer_type_corrected_updated remains statistically significant in multiple models,
# you have evidence that its effect is independent of age and sex.
#
# Conversely, if its significance disappears after accounting for age/sex,
# then it may be that the observed raw correlations are confounded by these factors.
#
# Additionally, consider exploring potential interactions:
# For example, you might test whether the effect of cancer type varies by sex:
#
#   lm(Endmotif_zscore ~ age + sex * cancer_type_corrected_updated, data = df_complete)
#
# which you could similarly analyze.
#
# ---------------------------
# End of Analysis
# ---------------------------



#### Now see the distribution of age and sex 

# -------------------------------
# 1. Bar Plot: Proportion of Sex by Cancer Type
# -------------------------------
sex_prop_df <- df_complete %>%
  group_by(cancer_type_corrected_updated, sex) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(cancer_type_corrected_updated) %>%
  mutate(prop = n / sum(n))

# Stacked bar plot showing proportion of each sex within each cancer type
p1 <- ggplot(sex_prop_df, aes(x = cancer_type_corrected_updated, y = prop, fill = sex)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = percent_format()) +
  labs(title = "Proportion of Sex by Cancer Type",
       x = "Cancer Type",
       y = "Proportion") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank())

# If you prefer side-by-side bars (instead of stacked), use: position = "dodge"

# -------------------------------
# 2. Bar Plot: Mean Age by Cancer Type
# -------------------------------
age_stats_df <- df_complete %>%
  group_by(cancer_type_corrected_updated) %>%
  summarise(mean_age = mean(age, na.rm = TRUE),
            sd_age = sd(age, na.rm = TRUE),
            n = n(),
            .groups = "drop")

p2 <- ggplot(age_stats_df, aes(x = cancer_type_corrected_updated, y = mean_age, fill = cancer_type_corrected_updated)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = mean_age - sd_age, ymax = mean_age + sd_age), width = 0.2) +
  labs(title = "Mean Age by Cancer Type",
       x = "Cancer Type",
       y = "Mean Age (+/- SD)") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# -------------------------------
# 3. Combine Plots Side by Side
# -------------------------------
combined_plot <- grid.arrange(p1, p2, ncol = 2)

# Optionally, save the combined figure to a file
ggsave("ProportionSex_and_MeanAge_byCancerType.pdf", combined_plot, width = 14, height = 6)






##### Now do power analysis 