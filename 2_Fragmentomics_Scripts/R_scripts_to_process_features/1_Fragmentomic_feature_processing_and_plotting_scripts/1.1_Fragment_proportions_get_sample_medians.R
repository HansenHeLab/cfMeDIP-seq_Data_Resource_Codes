# Delfi medians - 5Mb 

# Adapted from:
# file: M4_median.R
# author: Dory Abelman
# adapted from Derek Wong
# Original date: August 16th, 2023
# Updated for new samples Oct-Nov 2023
# Updated May 2024 for final set of adjusted prostate cancer samples
# Updated from Dec 2024 - April 2025 for final set of revisions
# Purpose:  Load DELFI ratio files, merge with metadata, compute medians, save RDS

library(dplyr)
library(GenomicRanges)
library(ggplot2)

# First load in the sample ID file 

#metadata_df <- read_csv("TCGE-CFMe-Only-Samples-light.csv")

# Read in files and combine into data frame
filedir <- '~/Documents/Thesis_work/R/M4/Projects/Fragmentation_CfMeDIP_Yong/DELFI_ratios/All_5Mb'
filedir2 <- '~/Documents/Thesis_work/R/M4/Projects/Fragmentation_CfMeDIP_Yong/Revisions Fall 2024/New_data/Paired_end/Output_Delfi/'
outdir <- '~/Documents/Thesis_work/R/M4/Projects/Fragmentation_CfMeDIP_Yong/Final Data Dec 2024/DELFI'
files <- list.files(filedir, pattern = "_5Mb_bins.txt", recursive = TRUE, full.names = TRUE)

# Combine file lists from both directories
files <- c(
  list.files(filedir, pattern = "_5Mb_bins.txt", recursive = TRUE, full.names = TRUE),
  list.files(filedir2, pattern = "_5Mb_bins.txt", recursive = TRUE, full.names = TRUE)
)

dir.create(outdir, showWarnings = FALSE)

bins.list <- lapply(files, read.delim)
tib.list <- lapply(bins.list, as_tibble)

# Get everything together
df.fr.tcge <- tib.list %>%
  bind_rows() %>% dplyr::select(everything())

# Add the sample metadata to the file 
# Adding a new column sample_id with _deduped at the end of sample_name
# Load metadata files into the environment
metadata_main <- read.csv("Revisions Fall 2024/Metadata tables/After QC/TCGE_cfDM_samples_only_fixed_after_QC_vialA_only.csv")
metadata_validation <- read.csv("Revisions Fall 2024/Metadata tables/After QC/Validation_TCGE_cfDM_samples_only_fixed_after_QC_vialA_only.csv")

# Combine metadata_main and metadata_validation
metadata_main <- full_join(metadata_main, metadata_validation)

# Load additional data files
file_PE <- read.csv("Revisions Fall 2024/Metadata tables/TCGE_cfDM_samples_only_fixed_after_QC_PE_vialA_only.csv")
file_SE <- read.csv("Revisions Fall 2024/Metadata tables/TCGE_cfDM_samples_only_fixed_after_QC_SE_vialA_only.csv")

# Add a new column 'type' based on sequencing_id
# Update metadata_main with 'type' and 'validation' columns
metadata_main$type <- ifelse(metadata_main$sequencing_id %in% metadata_validation$sequencing_id, "PE", 
                             ifelse(metadata_main$sequencing_id %in% file_PE$sequencing_id, "PE", 
                                    ifelse(metadata_main$sequencing_id %in% file_SE$sequencing_id, "SE", NA)))

metadata_main$validation <- ifelse(metadata_main$sequencing_id %in% metadata_validation$sequencing_id, 1, 0)

metadata_df <- metadata_main
metadata_df$cancer_type[metadata_df$cancer_type == "Healthy"] <- "Normal" # Edit change
metadata_df$sample_id <- paste(metadata_df$sequencing_id, "_dedup", sep = "")

metadata_df <- metadata_df %>% 
  mutate(cancer_type = ifelse(cancer_type == "Uveal Melanoma", "Eye Cancer", cancer_type))

metadata_df$cancer_type_corrected <- tolower(gsub(" ", "_", metadata_df$cancer_type)) # rename for ease of use 

#Export
write.csv(metadata_df, "metadata_updated_Dec2024.csv", row.names = FALSE)
saveRDS(metadata_df, "metadata_updated_Dec2024.rds")

# Keep only values in the metadata df
df.fr.tcge <- df.fr.tcge %>%
  inner_join(metadata_df, by = "sample_id")

## Check if any ids are missing 
# Extract unique sample_id from metadata_df
unique_metadata_sample_ids <- unique(metadata_df$sample_id)

# Extract unique sample_id from df.fr.tcge after merging
unique_merged_sample_ids <- unique(df.fr.tcge$sample_id)

# Find sample_ids present in metadata_df but not in the merged df.fr.tcge
missing_sample_ids <- setdiff(unique_metadata_sample_ids, unique_merged_sample_ids)

# Filter metadata_combined for the missing sample_ids and get their types
missing_samples_with_type <- metadata_df %>%
  filter(sample_id %in% missing_sample_ids) %>%
  select(sample_id, type)

# Print the results
print(missing_samples_with_type)

# Export 
#write.table(df.fr.tcge, file.path(outdir, paste0("tcge.median_fr_df.hg38.txt")), sep = "\t", row.names = FALSE)
save(df.fr.tcge, file=file.path(outdir, paste0("tcge.median_fr_df.hg38.rda")))

# Export cancer only 
df.fr.tcge.cancer <- df.fr.tcge %>% filter(cancer_type != "Normal") %>% filter(validation == "0")
df.fr.tcge.normal <- df.fr.tcge %>% filter(cancer_type == "Normal") %>% filter(validation == "0")

# Do this for the validation then plot seperately for those
df.fr.tcge.cancer.validation <- df.fr.tcge %>% filter(cancer_type != "Normal") %>% filter(validation == "1")
df.fr.tcge.normal.validation <- df.fr.tcge %>% filter(cancer_type == "Normal") %>% filter(validation == "1")

# Export 
#write.table(df.fr.tcge.cancer, file.path(outdir, paste0("tcge.cancer.median_fr_df.hg38.txt")), sep = "\t", row.names = FALSE)
save(df.fr.tcge.cancer, file=file.path(outdir, paste0("tcge.cancer.median_fr_df.hg38.rda")))
save(df.fr.tcge.normal, file=file.path(outdir, paste0("tcge.normal.median_fr_df.hg38.rda")))
save(df.fr.tcge.cancer.validation, file=file.path(outdir, paste0("tcge.cancer.median_fr_df.hg38.validation.rda")))
save(df.fr.tcge.normal.validation, file=file.path(outdir, paste0("tcge.normal.median_fr_df.hg38.validation.rda")))


## For ease of plotting 
df.fr.tcge <- df.fr.tcge %>% filter(validation == "0")

# Get coverage 
coverage <- df.fr.tcge %>% ungroup() %>% group_by(cancer_subtype) %>%
  dplyr::summarize(hqbases_analyzed = 100*sum(nfrags)*2,
                   depth = hqbases_analyzed/(504*5e6)
  )

# Generate tcge median
df.fr.tcge$cancer_type_corrected <- tolower(gsub(" ", "_", df.fr.tcge$cancer_type)) # rename for ease of use 

tcge_median <- df.fr.tcge %>%
  group_by(bin, cancer_type_corrected) %>% 
  dplyr::summarize(seqnames = unique(seqnames),
                   arm = unique(arm),
                   start = unique(start),
                   end = unique(end),
                   gc = unique(gc),
                   median_frag_gc = median(frag_gc, na.rm=TRUE),
                   median_ratio=median(short/long, na.rm=TRUE),
                   sd_ratio=sd(short/long, na.rm=TRUE),
                   median_ratio_corrected=median(short_corrected/long_corrected, na.rm=TRUE),
                   sd_ratio_corrected=sd(short_corrected/long_corrected, na.rm=TRUE),
                   median_ratio_centered=median(ratio_centered, na.rm=TRUE),
                   sd_ratio_centered=sd(ratio_centered, na.rm=TRUE),
                   median_coverage=median(coverage, na.rm=TRUE),
                   sd_coverage=sd(coverage, na.rm=TRUE),
                   median_coverage_corrected=median(coverage_corrected, na.rm=TRUE),
                   sd_coverage_corrected=sd(coverage_corrected, na.rm=TRUE),
                   median_coverage_centered=median(coverage_centered, na.rm=TRUE),
                   sd_coverage_centered=sd(coverage_centered, na.rm=TRUE),
                   median_combined=median(combined, na.rm=TRUE),
                   sd_combined=sd(combined, na.rm=TRUE),
                   median_combined_centered=median(combined_centered, na.rm=TRUE),
                   sd_combined_centered=sd(combined_centered, na.rm=TRUE),
                   median_mode_size=median(mode_size, na.rm=TRUE),
                   median_mean_size=median(mean_size, na.rm=TRUE),
                   sd_mean_size=sd(mean_size, na.rm=TRUE),
                   median_median_size=median(median_size, na.rm=TRUE),
                   sd_median_size=sd(median_size, na.rm=TRUE),
                   median_q25_size=median(q25_size, na.rm=TRUE),
                   median_q75_size=median(q75_size, na.rm=TRUE))

#write.table(tcge_median, file.path(outdir, paste0("tcge.median.hg38.txt")), sep = "\t", row.names = FALSE)
save(tcge_median, file=file.path(outdir, paste0("tcge.median.hg38.rda")))


# Plot profiles
mytheme <- theme_classic(base_size=12) + theme(
  axis.text.x = element_blank(),
  axis.ticks.x=element_blank(),
  strip.text.x = element_text(size=11),
  strip.text.y = element_text(size=12),
  axis.title.x = element_text(face="bold", size=17),
  axis.title.y = element_text(size=15),
  axis.text.y = element_text(size=15),
  plot.title = element_text(size=15),
  legend.position = "none",
  legend.title = element_text(size=10),
  legend.text = element_text(size=10),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background=element_rect(fill="white", color="white"),
  panel.spacing.x=unit(0.1, "lines"))

armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
               "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
               "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
               "19p", "19q","20p","20q","21q","22q")
df.fr.tcge$arm <- factor(df.fr.tcge$arm, levels=armlevels)
tcge_median$arm <- factor(tcge_median$arm, levels=armlevels)

arm <- df.fr.tcge %>% group_by(arm) %>%
  dplyr::summarize(n=n()) %>%
  dplyr::mutate(arm = as.character(arm))
small.arms <- setNames(c("", "10q", "", "12q", "", "16q",
                         "", "17q", "", "18q",
                         "", "", "", "",
                         "", ""),
                       c("10p", "10q", "12p", "12q", "16p", "16q",
                         "17p", "17q", "18p", "18q",
                         "19p", "19q", "20p", "20q",
                         "21q", "22q"))
arm.labels <- setNames(arm$arm, arm$arm)
arm.labels[names(small.arms)] <- small.arms

### Analyze by variable of interest

## Cancer Subtype 
## As function 
run_script_for_cancer_subtype <- function(cancer_subtype) {
  # Filter for cancer_subtype
  df.fr.tcge <- df.fr.tcge[grepl(paste0(cancer_subtype), df.fr.tcge$cancer_subtype), ]
  write.table(df.fr.tcge, file.path(outdir, paste0("tcge.median_fr_df.hg38", cancer_subtype, ".txt")), sep = "\t", row.names = FALSE)
  #  save(df.fr.tcge, file=file.path(outdir, paste0("tcge.median_fr_df.hg38", cancer_subtype, ".rda")))
  
  coverage <- df.fr.tcge %>% ungroup() %>% group_by(sample_id) %>%
    dplyr::summarize(hqbases_analyzed = 100*sum(nfrags)*2,
                     depth = hqbases_analyzed/(504*5e6)
    )
  
  # Generate tcge median
  tcge_median <- df.fr.tcge %>%
    group_by(bin) %>% 
    dplyr::summarize(sample_id = "median_tcge",
                     seqnames = unique(seqnames),
                     arm = unique(arm),
                     start = unique(start),
                     end = unique(end),
                     gc = unique(gc),
                     median_frag_gc = median(frag_gc, na.rm=TRUE),
                     median_ratio=median(short/long, na.rm=TRUE),
                     sd_ratio=sd(short/long, na.rm=TRUE),
                     median_ratio_corrected=median(short_corrected/long_corrected, na.rm=TRUE),
                     sd_ratio_corrected=sd(short_corrected/long_corrected, na.rm=TRUE),
                     median_ratio_centered=median(ratio_centered, na.rm=TRUE),
                     sd_ratio_centered=sd(ratio_centered, na.rm=TRUE),
                     median_coverage=median(coverage, na.rm=TRUE),
                     sd_coverage=sd(coverage, na.rm=TRUE),
                     median_coverage_corrected=median(coverage_corrected, na.rm=TRUE),
                     sd_coverage_corrected=sd(coverage_corrected, na.rm=TRUE),
                     median_coverage_centered=median(coverage_centered, na.rm=TRUE),
                     sd_coverage_centered=sd(coverage_centered, na.rm=TRUE),
                     median_combined=median(combined, na.rm=TRUE),
                     sd_combined=sd(combined, na.rm=TRUE),
                     median_combined_centered=median(combined_centered, na.rm=TRUE),
                     sd_combined_centered=sd(combined_centered, na.rm=TRUE),
                     median_mode_size=median(mode_size, na.rm=TRUE),
                     median_mean_size=median(mean_size, na.rm=TRUE),
                     sd_mean_size=sd(mean_size, na.rm=TRUE),
                     median_median_size=median(median_size, na.rm=TRUE),
                     sd_median_size=sd(median_size, na.rm=TRUE),
                     median_q25_size=median(q25_size, na.rm=TRUE),
                     median_q75_size=median(q75_size, na.rm=TRUE))
  
  write.table(tcge_median, file.path(outdir, paste0("tcge.median.hg38", cancer_subtype, ".txt")), sep = "\t", row.names = FALSE)
  #  save(tcge_median, file=file.path(outdir, paste0("tcge.median.hg38", cancer_subtype, ".rda")))
  
  # Generate correlations
  summary_df <- df.fr.tcge %>% ungroup() %>% group_by(sample_id) %>%
    dplyr::summarize(
      ratio_cor=cor(ratio, tcge_median$median_ratio, method="pearson", use="complete.obs"),
      ratio_corrected_cor=cor(ratio_corrected, tcge_median$median_ratio_corrected, method="pearson", use="complete.obs"),
      ratio_centered_cor=cor(ratio_centered, tcge_median$median_ratio_centered, method="pearson", use="complete.obs"),
      coverage_cor=cor(coverage, tcge_median$median_coverage, method="pearson", use="complete.obs"),
      coverage_corrected_cor=cor(coverage_corrected, tcge_median$median_coverage_corrected, method="pearson", use="complete.obs"),
      coverage_centered_cor=cor(coverage_centered, tcge_median$median_coverage_centered, method="pearson", use="complete.obs"),
      combined_cor=cor(combined, tcge_median$median_combined, method="pearson", use="complete.obs"),
      combined_centered_cor=cor(combined_centered, tcge_median$median_combined_centered, method="pearson", use="complete.obs"),
      nfrags = sum(nfrags),
      mode_size=unique(mode_size),
      mean_size=unique(mean_size),
      median_size=unique(median_size),
      q25_size=unique(q25_size),
      q75_size=unique(q75_size),
      hqbases_analyzed = 100*sum(nfrags)*2,
      coverage = hqbases_analyzed/(504*5e6)
    )
  
  write.table(summary_df, file.path(outdir, paste0("summary_correlations", cancer_subtype, ".txt")), sep = "\t", row.names = FALSE)
  
  # Plot profiles
  mytheme <- theme_classic(base_size=12) + theme(
    axis.text.x = element_blank(),
    axis.ticks.x=element_blank(),
    strip.text.x = element_text(size=11),
    strip.text.y = element_text(size=12),
    axis.title.x = element_text(face="bold", size=17),
    axis.title.y = element_text(size=15),
    axis.text.y = element_text(size=15),
    plot.title = element_text(size=15),
    legend.position = "none",
    legend.title = element_text(size=10),
    legend.text = element_text(size=10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background=element_rect(fill="white", color="white"),
    panel.spacing.x=unit(0.1, "lines"))
  
  armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
                 "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
                 "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
                 "19p", "19q","20p","20q","21q","22q")
  df.fr.tcge$arm <- factor(df.fr.tcge$arm, levels=armlevels)
  tcge_median$arm <- factor(tcge_median$arm, levels=armlevels)
  
  arm <- df.fr.tcge %>% group_by(arm) %>%
    dplyr::summarize(n=n()) %>%
    dplyr::mutate(arm = as.character(arm))
  small.arms <- setNames(c("", "10q", "", "12q", "", "16q",
                           "", "17q", "", "18q",
                           "", "", "", "",
                           "", ""),
                         c("10p", "10q", "12p", "12q", "16p", "16q",
                           "17p", "17q", "18p", "18q",
                           "19p", "19q", "20p", "20q",
                           "21q", "22q"))
  arm.labels <- setNames(arm$arm, arm$arm)
  arm.labels[names(small.arms)] <- small.arms
  
  # Generate Fragmentation and Coverage plots
  g1 <- ggplot(df.fr.tcge, aes(x=bin, y=ratio_centered, group=sample_id, color="red")) + 
    geom_line(linewidth=0.5, alpha=0.5)
  #geom_line(data=subset(df.fr.tcge, sample_id=="", size=0.75, alpha=1))
  g1 <- g1 + geom_line(data=tcge_median, aes(x=bin, y=median_ratio_centered), size=0.75, alpha=0.5, color="black")
  g1 <- g1 + labs(x="", y="Fragmentation profile\n", color="")
  g1 <- g1 + facet_grid(~arm, switch="x",space="free_x", scales="free_x", labeller=labeller(arm=arm.labels))
  g1 <- g1 + coord_cartesian(xlim = NULL, ylim=c(-.10,.12), expand = TRUE)
  g1 <- g1 + mytheme
  g1
  ggsave(file.path(outdir, paste0("tcge_fragment", cancer_subtype, ".pdf")), g1, width=15, height=3, units="in")
  
  c1 <- ggplot(df.fr.tcge, aes(x=bin, y=coverage_centered, group=sample_id, color="red")) + 
    geom_line(linewidth=0.5, alpha=0.5)
  #geom_line(data=subset(df.fr.tcge, sample_id==""), size=0.75, alpha=1)
  c1 <- c1 + geom_line(data=tcge_median, aes(x=bin, y=median_coverage_centered), size=0.75, alpha=0.5, color="black")
  c1 <- c1 + labs(x="", y="Coverage profile\n", color="")
  c1 <- c1 + facet_grid(~arm, switch="x",space="free_x", scales="free_x", labeller=labeller(arm=arm.labels))
  c1 <- c1 + coord_cartesian(xlim = NULL, ylim=c(-.10,.12), expand = TRUE)
  c1 <- c1 + mytheme
  c1
  ggsave(file.path(outdir, paste0("tcge_coverage", cancer_subtype, ".pdf")), c1, width=15, height=3, units="in")
  
  b1 <- ggplot(df.fr.tcge, aes(x=bin, y=combined_centered, group=sample_id, color="red")) + 
    geom_line(linewidth=0.5, alpha=0.5)
  #geom_line(data=subset(df.fr.tcge, sample_id==""), size=0.75, alpha=1)
  b1 <- b1 + geom_line(data=tcge_median, aes(x=bin, y=median_combined_centered), size=0.75, alpha=0.5, color="black")
  b1 <- b1 + labs(x="", y="Combined profile\n", color="")
  b1 <- b1 + facet_grid(~arm, switch="x",space="free_x", scales="free_x", labeller=labeller(arm=arm.labels))
  b1 <- b1 + coord_cartesian(xlim = NULL, ylim=c(-.10,.12), expand = TRUE)
  b1 <- b1 + mytheme
  b1
  ggsave(file.path(outdir, paste0("tcge_combined", cancer_subtype, ".pdf")), b1, width=15, height=3, units="in")
  
  
  return(paste0("Script run for subtype ", cancer_subtype))
}

# Now  can run the script for each subtype like this:
## Get subtypes 
# Filter metadata_df for only 'PE' type
metadata_df_pe <- metadata_df[metadata_df$type == "PE", ]

# Get unique cancer types for PE samples
cancer_types_delfi <- unique(metadata_df_pe$cancer_type)

# Create a table of cancer type counts for PE samples
type_counts <- as.data.frame(table(metadata_df_pe$cancer_type))

## Get subtypes for PE samples
cancer_subtypes <- unique(metadata_df_pe$cancer_subtype)

# First, create a table of counts for each cancer subtype for PE samples
subtype_counts <- table(metadata_df_pe$cancer_subtype)

# Filter for cancer subtypes with at least 3 samples for PE samples
valid_subtypes <- names(subtype_counts)[subtype_counts >= 3]

# Loop over the subtypes and run the script if there are 3 or more occurrences
for (cancer_subtype in valid_subtypes) {
  if (!is.na(cancer_subtype)) {
    print(run_script_for_cancer_subtype(cancer_subtype))
  }
}


## Now run the script to plot the ratios for the cancer types
run_script_for_cancer_type <- function(cancer_type) {
  # Filter for cancer_type
  df.fr.tcge <- df.fr.tcge[grepl(paste0(cancer_type), df.fr.tcge$cancer_type), ]
  write.table(df.fr.tcge, file.path(outdir, paste0("tcge.median_fr_df.hg38", cancer_type, ".txt")), sep = "\t", row.names = FALSE)
  #  save(df.fr.tcge, file=file.path(outdir, paste0("tcge.median_fr_df.hg38", cancer_type, ".rda")))
  
  coverage <- df.fr.tcge %>% ungroup() %>% group_by(sample_id) %>%
    dplyr::summarize(hqbases_analyzed = 100*sum(nfrags)*2,
                     depth = hqbases_analyzed/(504*5e6)
    )
  
  # Generate tcge median
  tcge_median <- df.fr.tcge %>%
    group_by(bin) %>% 
    dplyr::summarize(sample_id = "median_tcge",
                     seqnames = unique(seqnames),
                     arm = unique(arm),
                     start = unique(start),
                     end = unique(end),
                     gc = unique(gc),
                     median_frag_gc = median(frag_gc, na.rm=TRUE),
                     median_ratio=median(short/long, na.rm=TRUE),
                     sd_ratio=sd(short/long, na.rm=TRUE),
                     median_ratio_corrected=median(short_corrected/long_corrected, na.rm=TRUE),
                     sd_ratio_corrected=sd(short_corrected/long_corrected, na.rm=TRUE),
                     median_ratio_centered=median(ratio_centered, na.rm=TRUE),
                     sd_ratio_centered=sd(ratio_centered, na.rm=TRUE),
                     median_coverage=median(coverage, na.rm=TRUE),
                     sd_coverage=sd(coverage, na.rm=TRUE),
                     median_coverage_corrected=median(coverage_corrected, na.rm=TRUE),
                     sd_coverage_corrected=sd(coverage_corrected, na.rm=TRUE),
                     median_coverage_centered=median(coverage_centered, na.rm=TRUE),
                     sd_coverage_centered=sd(coverage_centered, na.rm=TRUE),
                     median_combined=median(combined, na.rm=TRUE),
                     sd_combined=sd(combined, na.rm=TRUE),
                     median_combined_centered=median(combined_centered, na.rm=TRUE),
                     sd_combined_centered=sd(combined_centered, na.rm=TRUE),
                     median_mode_size=median(mode_size, na.rm=TRUE),
                     median_mean_size=median(mean_size, na.rm=TRUE),
                     sd_mean_size=sd(mean_size, na.rm=TRUE),
                     median_median_size=median(median_size, na.rm=TRUE),
                     sd_median_size=sd(median_size, na.rm=TRUE),
                     median_q25_size=median(q25_size, na.rm=TRUE),
                     median_q75_size=median(q75_size, na.rm=TRUE))
  
  #write.table(tcge_median, file.path(outdir, paste0("tcge.median.hg38", cancer_type, ".txt")), sep = "\t", row.names = FALSE)
  save(tcge_median, file=file.path(outdir, paste0("tcge.median.hg38", cancer_type, ".rda")))
  
  safe_cor <- function(x, y, method = "pearson") {
    if (sum(!is.na(x) & !is.na(y)) > 1) {
      return(cor(x, y, method = method, use = "complete.obs"))
    } else {
      return(NA)
    }
  }
  
  # Generate correlations
  summary_df <- df.fr.tcge %>% ungroup() %>% group_by(sample_id) %>%
    dplyr::summarize(
      ratio_cor = safe_cor(ratio, tcge_median$median_ratio),
      ratio_corrected_cor = safe_cor(ratio_corrected, tcge_median$median_ratio_corrected),
      ratio_centered_cor = safe_cor(ratio_centered, tcge_median$median_ratio_centered),
      coverage_cor = safe_cor(coverage, tcge_median$median_coverage),
      coverage_corrected_cor = safe_cor(coverage_corrected, tcge_median$median_coverage_corrected),
      coverage_centered_cor = safe_cor(coverage_centered, tcge_median$median_coverage_centered),
      combined_cor = safe_cor(combined, tcge_median$median_combined),
      combined_centered_cor = safe_cor(combined_centered, tcge_median$median_combined_centered),
      nfrags = sum(nfrags),
      mode_size=unique(mode_size),
      mean_size=unique(mean_size),
      median_size=unique(median_size),
      q25_size=unique(q25_size),
      q75_size=unique(q75_size),
      hqbases_analyzed = 100*sum(nfrags)*2,
      coverage = hqbases_analyzed/(504*5e6)
    )
  
  write.table(summary_df, file.path(outdir, paste0("summary_correlations", cancer_type, ".txt")), sep = "\t", row.names = FALSE)
  
  # Plot profiles
  mytheme <- theme_classic(base_size=12) + theme(
    axis.text.x = element_blank(),
    axis.ticks.x=element_blank(),
    strip.text.x = element_text(size=11),
    strip.text.y = element_text(size=12),
    axis.title.x = element_text(face="bold", size=17),
    axis.title.y = element_text(size=15),
    axis.text.y = element_text(size=15),
    plot.title = element_text(size=15),
    legend.position = "none",
    legend.title = element_text(size=10),
    legend.text = element_text(size=10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background=element_rect(fill="white", color="white"),
    panel.spacing.x=unit(0.1, "lines"))
  
  armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
                 "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
                 "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
                 "19p", "19q","20p","20q","21q","22q")
  df.fr.tcge$arm <- factor(df.fr.tcge$arm, levels=armlevels)
  tcge_median$arm <- factor(tcge_median$arm, levels=armlevels)
  
  arm <- df.fr.tcge %>% group_by(arm) %>%
    dplyr::summarize(n=n()) %>%
    dplyr::mutate(arm = as.character(arm))
  small.arms <- setNames(c("", "10q", "", "12q", "", "16q",
                           "", "17q", "", "18q",
                           "", "", "", "",
                           "", ""),
                         c("10p", "10q", "12p", "12q", "16p", "16q",
                           "17p", "17q", "18p", "18q",
                           "19p", "19q", "20p", "20q",
                           "21q", "22q"))
  arm.labels <- setNames(arm$arm, arm$arm)
  arm.labels[names(small.arms)] <- small.arms
  
  # Generate Fragmentation and Coverage plots
  g1 <- ggplot(df.fr.tcge, aes(x=bin, y=ratio_centered, group=sample_id, color="red")) + 
    geom_line(linewidth=0.5, alpha=0.5)
  #geom_line(data=subset(df.fr.tcge, sample_id=="", size=0.75, alpha=1))
  g1 <- g1 + geom_line(data=tcge_median, aes(x=bin, y=median_ratio_centered), size=0.75, alpha=0.5, color="black")
  g1 <- g1 + labs(x="", y="Fragmentation profile\n", color="")
  g1 <- g1 + facet_grid(~arm, switch="x",space="free_x", scales="free_x", labeller=labeller(arm=arm.labels))
  g1 <- g1 + coord_cartesian(xlim = NULL, ylim=c(-.10,.12), expand = TRUE)
  g1 <- g1 + mytheme
  g1
  ggsave(file.path(outdir, paste0("tcge_fragment", cancer_type, ".pdf")), g1, width=15, height=3, units="in")
  
  c1 <- ggplot(df.fr.tcge, aes(x=bin, y=coverage_centered, group=sample_id, color="red")) + 
    geom_line(linewidth=0.5, alpha=0.5)
  #geom_line(data=subset(df.fr.tcge, sample_id==""), size=0.75, alpha=1)
  c1 <- c1 + geom_line(data=tcge_median, aes(x=bin, y=median_coverage_centered), size=0.75, alpha=0.5, color="black")
  c1 <- c1 + labs(x="", y="Coverage profile\n", color="")
  c1 <- c1 + facet_grid(~arm, switch="x",space="free_x", scales="free_x", labeller=labeller(arm=arm.labels))
  c1 <- c1 + coord_cartesian(xlim = NULL, ylim=c(-.10,.12), expand = TRUE)
  c1 <- c1 + mytheme
  c1
  ggsave(file.path(outdir, paste0("tcge_coverage", cancer_type, ".pdf")), c1, width=15, height=3, units="in")
  
  b1 <- ggplot(df.fr.tcge, aes(x=bin, y=combined_centered, group=sample_id, color="red")) + 
    geom_line(linewidth=0.5, alpha=0.5)
  #geom_line(data=subset(df.fr.tcge, sample_id==""), size=0.75, alpha=1)
  b1 <- b1 + geom_line(data=tcge_median, aes(x=bin, y=median_combined_centered), size=0.75, alpha=0.5, color="black")
  b1 <- b1 + labs(x="", y="Combined profile\n", color="")
  b1 <- b1 + facet_grid(~arm, switch="x",space="free_x", scales="free_x", labeller=labeller(arm=arm.labels))
  b1 <- b1 + coord_cartesian(xlim = NULL, ylim=c(-.10,.12), expand = TRUE)
  b1 <- b1 + mytheme
  b1
  ggsave(file.path(outdir, paste0("tcge_combined", cancer_type, ".pdf")), b1, width=15, height=3, units="in")
  
  
  return(paste0("Script run for type ", cancer_type))
}

# Filter metadata_df for only 'PE' type
metadata_df_pe <- metadata_df[metadata_df$type == "PE", ]

# Get unique cancer types for PE samples
cancer_types_delfi <- unique(metadata_df_pe$cancer_type)

# Create a table of cancer type counts for PE samples
type_counts <- as.data.frame(table(metadata_df_pe$cancer_type))

## Run the script for each type
for (cancer_type in cancer_types_delfi) {
  if (!is.na(cancer_type)) {
    print(run_script_for_cancer_type(cancer_type))
  }
}



## Do the cancer vs normal script comparison
df.fr.tcge <- df.fr.tcge %>%
  dplyr::mutate(CN_classifier = ifelse(cancer_type_corrected == "normal", "normal", "cancer"))

df.fr.tcge <- df.fr.tcge %>% 
  filter(!is.na(cancer_type_corrected))

df.fr.tcge <- df.fr.tcge %>% 
  mutate(cancer_type = ifelse(cancer_type == "uveal_melanoma", "eye_cancer", cancer_type)) %>% 
  filter(!is.na(cancer_type_corrected))

## Do for cancer vs normal -- added May 2024

## Do it for the cancer types
run_script_for_all_cancer_together <- function(cancer_type) {
  # Filter for cancer_type
  coverage <- df.fr.tcge %>% ungroup() %>% dplyr::group_by(sample_id) %>%
    dplyr::summarize(hqbases_analyzed = 100*sum(nfrags)*2,
                     depth = hqbases_analyzed/(504*5e6)
    )
  
  # Generate tcge median
  tcge_median <- df.fr.tcge %>% dplyr::filter(CN_classifier == "normal") %>%
    dplyr::group_by(bin) %>% 
    dplyr::summarize(sample_id = "median_tcge",
                     seqnames = unique(seqnames),
                     arm = unique(arm),
                     start = unique(start),
                     end = unique(end),
                     gc = unique(gc),
                     median_frag_gc = median(frag_gc, na.rm=TRUE),
                     median_ratio=median(short/long, na.rm=TRUE),
                     sd_ratio=sd(short/long, na.rm=TRUE),
                     median_ratio_corrected=median(short_corrected/long_corrected, na.rm=TRUE),
                     sd_ratio_corrected=sd(short_corrected/long_corrected, na.rm=TRUE),
                     median_ratio_centered=median(ratio_centered, na.rm=TRUE),
                     sd_ratio_centered=sd(ratio_centered, na.rm=TRUE),
                     median_coverage=median(coverage, na.rm=TRUE),
                     sd_coverage=sd(coverage, na.rm=TRUE),
                     median_coverage_corrected=median(coverage_corrected, na.rm=TRUE),
                     sd_coverage_corrected=sd(coverage_corrected, na.rm=TRUE),
                     median_coverage_centered=median(coverage_centered, na.rm=TRUE),
                     sd_coverage_centered=sd(coverage_centered, na.rm=TRUE),
                     median_combined=median(combined, na.rm=TRUE),
                     sd_combined=sd(combined, na.rm=TRUE),
                     median_combined_centered=median(combined_centered, na.rm=TRUE),
                     sd_combined_centered=sd(combined_centered, na.rm=TRUE),
                     median_mode_size=median(mode_size, na.rm=TRUE),
                     median_mean_size=median(mean_size, na.rm=TRUE),
                     sd_mean_size=sd(mean_size, na.rm=TRUE),
                     median_median_size=median(median_size, na.rm=TRUE),
                     sd_median_size=sd(median_size, na.rm=TRUE),
                     median_q25_size=median(q25_size, na.rm=TRUE),
                     median_q75_size=median(q75_size, na.rm=TRUE))
  
  #write.table(tcge_median, file.path(outdir, paste0("tcge.median.hg38", cancer_type, ".txt")), sep = "\t", row.names = FALSE)
  save(tcge_median, file=file.path(outdir, paste0("tcge.median.hg38", cancer_type, ".rda")))
  
  # Generate correlations
  summary_df <- df.fr.tcge %>% ungroup() %>% dplyr::group_by(sample_id) %>%
    dplyr::summarize(
      ratio_cor=cor(ratio, tcge_median$median_ratio, method="pearson", use="complete.obs"),
      ratio_corrected_cor=cor(ratio_corrected, tcge_median$median_ratio_corrected, method="pearson", use="complete.obs"),
      ratio_centered_cor=cor(ratio_centered, tcge_median$median_ratio_centered, method="pearson", use="complete.obs"),
      coverage_cor=cor(coverage, tcge_median$median_coverage, method="pearson", use="complete.obs"),
      coverage_corrected_cor=cor(coverage_corrected, tcge_median$median_coverage_corrected, method="pearson", use="complete.obs"),
      coverage_centered_cor=cor(coverage_centered, tcge_median$median_coverage_centered, method="pearson", use="complete.obs"),
      combined_cor=cor(combined, tcge_median$median_combined, method="pearson", use="complete.obs"),
      combined_centered_cor=cor(combined_centered, tcge_median$median_combined_centered, method="pearson", use="complete.obs"),
      nfrags = sum(nfrags),
      mode_size=unique(mode_size),
      mean_size=unique(mean_size),
      median_size=unique(median_size),
      q25_size=unique(q25_size),
      q75_size=unique(q75_size),
      hqbases_analyzed = 100*sum(nfrags)*2,
      coverage = hqbases_analyzed/(504*5e6)
    )
  
  write.table(summary_df, file.path(outdir, paste0("summary_correlations", cancer_type, ".txt")), sep = "\t", row.names = FALSE)
  
  # Plot profiles
  mytheme <- theme_classic(base_size=12) + theme(
    axis.text.x = element_blank(),
    axis.ticks.x=element_blank(),
    strip.text.x = element_text(size=11),
    strip.text.y = element_text(size=12),
    axis.title.x = element_text(face="bold", size=17),
    axis.title.y = element_text(size=15),
    axis.text.y = element_text(size=15),
    plot.title = element_text(size=15),
    legend.position = "none",
    legend.title = element_text(size=10),
    legend.text = element_text(size=10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background=element_rect(fill="white", color="white"),
    panel.spacing.x=unit(0.1, "lines"))
  
  armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
                 "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
                 "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
                 "19p", "19q","20p","20q","21q","22q")
  df.fr.tcge$arm <- factor(df.fr.tcge$arm, levels=armlevels)
  tcge_median$arm <- factor(tcge_median$arm, levels=armlevels)
  
  arm <- df.fr.tcge %>% dplyr::group_by(arm) %>%
    dplyr::summarize(n=n()) %>%
    dplyr::mutate(arm = as.character(arm))
  small.arms <- setNames(c("", "10q", "", "12q", "", "16q",
                           "", "17q", "", "18q",
                           "", "", "", "",
                           "", ""),
                         c("10p", "10q", "12p", "12q", "16p", "16q",
                           "17p", "17q", "18p", "18q",
                           "19p", "19q", "20p", "20q",
                           "21q", "22q"))
  arm.labels <- setNames(arm$arm, arm$arm)
  arm.labels[names(small.arms)] <- small.arms
  
  # Generate Fragmentation and Coverage plots
  g1 <- ggplot(df.fr.tcge %>% dplyr::filter(CN_classifier == "cancer"), aes(x=bin, y=ratio_centered, group=sample_id, color="red")) + 
    geom_line(linewidth=0.5, alpha=0.5)
  #geom_line(data=subset(df.fr.tcge, sample_id=="", size=0.75, alpha=1))
  g1 <- g1 + geom_line(data=tcge_median, aes(x=bin, y=median_ratio_centered), size=0.75, alpha=0.5, color="black")
  g1 <- g1 + labs(x="", y="Fragmentation profile\n", color="")
  g1 <- g1 + facet_grid(~arm, switch="x",space="free_x", scales="free_x", labeller=labeller(arm=arm.labels))
  g1 <- g1 + coord_cartesian(xlim = NULL, ylim=c(-.10,.12), expand = TRUE)
  g1 <- g1 + mytheme
  g1
  ggsave(file.path(outdir, paste0("tcge_fragment", cancer_type, "vs normal.pdf")), g1, width=15, height=3, units="in")
  
  ## Do colored version 
  col_cancer_type <- c('normal' = '#33a02c',
                       'brain_cancer' = '#1f78b4',
                       'lung_cancer' = '#b2df8a',
                       'prostate_cancer' = '#a6cee3',
                       'blood_cancer' = '#fb9a99',
                       'pancreatic_cancer' = '#e31a1c',
                       'eye_cancer' = '#fdbf6f',
                       'head_and_neck_cancer' = '#ff7f00',
                       'breast_cancer' = '#cab2d6',
                       'colorectal_cancer' = '#6a3d9a',
                       'bladder_cancer' = '#fb6a4a',
                       'renal_cancer' = '#b15928',
                       'lfs_survivor' = '#bdbdbd',
                       'lfs_previvor' = '#969696',
                       'lfs_positive' = '#737373')
  
  g1 <- ggplot(df.fr.tcge %>% dplyr::filter(CN_classifier == "cancer"), 
               aes(x=bin, y=ratio_centered, group=sample_id, color=cancer_type_corrected)) + 
    geom_line(linewidth=0.5, alpha=0.25) +
    geom_line(data=tcge_median, aes(x=bin, y=median_ratio_centered), size=0.75, alpha=0.5, color="black") +
    labs(x="", y="Fragmentation profile\n", color="Cancer Type Corrected") +
    facet_grid(~arm, switch="x", space="free_x", scales="free_x", labeller=labeller(arm=arm.labels)) +
    coord_cartesian(xlim = NULL, ylim=c(-.10,.12), expand = TRUE) +
    mytheme +
    scale_color_manual(values = col_cancer_type)
  g1
  ggsave(file.path(outdir, paste0("tcge_fragment", cancer_type, "vs normal, colored, without normal.pdf")), g1, width=15, height=3, units="in")
  
  g1 <- ggplot(df.fr.tcge, 
               aes(x=bin, y=ratio_centered, group=sample_id, color=cancer_type_corrected)) + 
    geom_line(linewidth=0.5, alpha=0.25) +
    geom_line(data=tcge_median, aes(x=bin, y=median_ratio_centered), size=0.75, alpha=0.5, color="black") +
    labs(x="", y="Fragmentation profile\n", color="Cancer Type Corrected") +
    facet_grid(~arm, switch="x", space="free_x", scales="free_x", labeller=labeller(arm=arm.labels)) +
    coord_cartesian(xlim = NULL, ylim=c(-.10,.12), expand = TRUE) +
    mytheme + 
    scale_color_manual(values = col_cancer_type)
  g1
  ggsave(file.path(outdir, paste0("tcge_fragment", cancer_type, "vs normal, colored, with normal.pdf")), g1, width=15, height=3, units="in")
  
  c1 <- ggplot(df.fr.tcge %>% dplyr::filter(CN_classifier == "cancer"), aes(x=bin, y=coverage_centered, group=sample_id, color="red")) + 
    geom_line(linewidth=0.5, alpha=0.5)
  #geom_line(data=subset(df.fr.tcge, sample_id==""), size=0.75, alpha=1)
  c1 <- c1 + geom_line(data=tcge_median, aes(x=bin, y=median_coverage_centered), size=0.75, alpha=0.5, color="black")
  c1 <- c1 + labs(x="", y="Coverage profile\n", color="")
  c1 <- c1 + facet_grid(~arm, switch="x",space="free_x", scales="free_x", labeller=labeller(arm=arm.labels))
  c1 <- c1 + coord_cartesian(xlim = NULL, ylim=c(-.10,.12), expand = TRUE)
  c1 <- c1 + mytheme
  c1
  ggsave(file.path(outdir, paste0("tcge_coverage", cancer_type, "vs normal.pdf")), c1, width=15, height=3, units="in")
  
  b1 <- ggplot(df.fr.tcge %>% dplyr::filter(CN_classifier == "cancer"), aes(x=bin, y=combined_centered, group=sample_id, color="red")) + 
    geom_line(linewidth=0.5, alpha=0.5)
  #geom_line(data=subset(df.fr.tcge, sample_id==""), size=0.75, alpha=1)
  b1 <- b1 + geom_line(data=tcge_median, aes(x=bin, y=median_combined_centered), size=0.75, alpha=0.5, color="black")
  b1 <- b1 + labs(x="", y="Combined profile\n", color="")
  b1 <- b1 + facet_grid(~arm, switch="x",space="free_x", scales="free_x", labeller=labeller(arm=arm.labels))
  b1 <- b1 + coord_cartesian(xlim = NULL, ylim=c(-.10,.12), expand = TRUE)
  b1 <- b1 + mytheme
  b1
  ggsave(file.path(outdir, paste0("tcge_combined", cancer_type, "vs normal.pdf")), b1, width=15, height=3, units="in")
  
  
  return(paste0("Script run for all ", cancer_type))
}



print(run_script_for_all_cancer_together("cancer"))





### Repeat analysis for validation cohorts
df.fr.tcge <- bind_rows(df.fr.tcge.cancer.validation, df.fr.tcge.normal.validation)

# Filter metadata_df for only 'PE' type
metadata_df_validation <- metadata_df[metadata_df$validation == "1", ]

# Get unique cancer types for PE samples
cancer_types_delfi <- unique(metadata_df_validation$cancer_type)

# Create a table of cancer type counts for PE samples
type_counts <- as.data.frame(table(metadata_df_validation$cancer_type))

## Get subtypes for PE samples
cancer_subtypes <- unique(metadata_df_validation$cancer_subtype)

# First, create a table of counts for each cancer subtype for PE samples
subtype_counts <- table(metadata_df_validation$cancer_subtype)

# Filter for cancer subtypes with at least 3 samples for PE samples
valid_subtypes <- names(subtype_counts)[subtype_counts >= 3]

# Loop over the subtypes and run the script if there are 3 or more occurrences
for (cancer_subtype in valid_subtypes) {
  if (!is.na(cancer_subtype)) {
    print(run_script_for_cancer_subtype(cancer_subtype))
  }
}

## Run the script for each type
for (cancer_type in cancer_types_delfi) {
  if (!is.na(cancer_type)) {
    print(run_script_for_cancer_type(cancer_type))
  }
}


## Do the cancer vs normal script comparison
df.fr.tcge <- df.fr.tcge %>%
  dplyr::mutate(CN_classifier = ifelse(cancer_type_corrected == "normal", "normal", "cancer"))

df.fr.tcge <- df.fr.tcge %>% 
  filter(!is.na(cancer_type_corrected))

df.fr.tcge <- df.fr.tcge %>% 
  mutate(cancer_type = ifelse(cancer_type == "uveal_melanoma", "eye_cancer", cancer_type)) %>% 
  filter(!is.na(cancer_type_corrected))

print(run_script_for_all_cancer_together("cancer"))



### Now need to think of some statistical comparisons between these to see if they have similar distributions across cancer types
## Maybe by comparing the z-scores in the next one
## Need to either add validation to the plot or put seperate

## Clean up 
rm(bins.list)
rm(tib.list)
