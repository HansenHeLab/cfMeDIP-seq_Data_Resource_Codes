# ----------------------------------------------------------------------------
# Title   : Compare End‑Motif Profiles Between PE and SE Samples via MANOVA
# Authors : Dory Abelman
# Date    : March 2025
#
# Purpose :
#   • Load PE and SE end‑motif frequency matrices.
#   • Harmonize sample identifiers and annotate each row by sample type (PE vs. SE).
#   • Check matrix rank to account for compositional redundancy.
#   • Perform a multivariate analysis of variance (MANOVA) using Pillai’s trace
#     to test for global differences in end‑motif patterns between PE and SE.
#   • Extract and report Pillai’s trace, approximate F, degrees of freedom, and
#     p‑value for manuscript narrative.
# ----------------------------------------------------------------------------


### Check if end motifs are consistent between PE and SE 

endmotif_se <- readRDS("path/file.RDS")

end_motif_pe <- readRDS("Final Data Dec 2024/End_motifs/Without_validation_updated/End_motif_matrix_Mar2025.rds")


# Load necessary libraries
library(tidyverse)
library(car)

# --- Prepare Data -------------------------------------------------------------
# For SE data (assumed to be a tibble)
df_se <- end_motifs_se %>% 
  mutate(sample_type = "SE")

# For PE data (convert matrix to tibble and add sample_type)
df_pe <- as_tibble(end_motif_pe) %>% 
  mutate(sample_type = "PE")

# Combine both datasets into one data frame
df_all <- bind_rows(df_se, df_pe)

# --- Examine the Response Matrix ---------------------------------------------
# Extract the response matrix (columns 1:256 represent the 256 motif features)
response_matrix <- as.matrix(df_all[, 1:256])
orig_rank <- qr(response_matrix)$rank
cat("Original rank of response matrix:", orig_rank, "out of", ncol(response_matrix), "\n")

# Because the data are compositional (the motif frequencies sum to a constant per sample),
# one feature is redundant even if QR indicates full rank. To be conservative, we remove one column.
if (orig_rank == ncol(response_matrix)) {
  cat("The response matrix appears full rank, but one motif is redundant due to compositional constraints.\n")
  # Remove the last column (or any one column) to remove redundancy.
  response_matrix_adj <- response_matrix[, -ncol(response_matrix)]
  cat("Adjusted response matrix dimensions:", dim(response_matrix_adj), "\n")
} else {
  response_matrix_adj <- response_matrix
}

# --- Perform MANOVA on the Adjusted Response Matrix -------------------------
manova_fit <- manova(response_matrix_adj ~ sample_type, data = df_all)
manova_res <- summary(manova_fit, test = "Pillai")
print(manova_res)

# --- Extract MANOVA Statistics for Narrative Reporting -----------------------
manova_stats <- manova_res$stats
pillai <- manova_stats[1, "Pillai"]
approx_F <- manova_stats[1, "approx F"]
num_df <- manova_stats[1, "num Df"]
den_df <- manova_stats[1, "den Df"]
pvalue_manova <- manova_stats[1, "Pr(>F)"]

cat("\n--- MANOVA Results ---\n")
cat("Pillai's trace =", round(pillai, 3), "\n")
cat("Approximate F(", num_df, ",", den_df, ") =", round(approx_F, 3), "\n")
cat("p-value =", format(pvalue_manova, scientific = TRUE, digits = 3), "\n")
