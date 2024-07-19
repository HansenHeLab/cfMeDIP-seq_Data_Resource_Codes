# Delfi medians - 5Mb 

# Adapted from:
# file: M4_median.R
# author: Dory Abelman
# adapted from Derek Wong
# Original date: August 16th, 2023
# Updated for new samples Oct-Nov 2023
# Updated May 2024 for final set of adjusted prostate cancer samples

library(dplyr)
library(GenomicRanges)
library(ggplot2)

# First load in the sample ID file 

#metadata_df <- read_csv("TCGE-CFMe-Only-Samples-light.csv")

# Read in files and combine into data frame
filedir <- '~/Documents/Thesis_work/R/M4/Projects/Fragmentation_CfMeDIP_Yong/DELFI_ratios/All_5Mb'
outdir <- '~/Documents/Thesis_work/R/M4/Projects/Fragmentation_CfMeDIP_Yong/Final Data May 2024/DELFI/'
files <- list.files(filedir, pattern = "_5Mb_bins.txt", recursive = TRUE, full.names = TRUE)

dir.create(outdir, showWarnings = FALSE)

bins.list <- lapply(files, read.delim)
tib.list <- lapply(bins.list, as_tibble)

# Get everything together
df.fr.tcge <- tib.list %>%
  bind_rows() %>% dplyr::select(everything())

# Add the sample metadata to the file 
# Adding a new column sample_id with _deduped at the end of sample_name
metadata_df <- read.csv("Updated metadata file/TCGE_cfDM_samples_only_fixed_after_QC_PE.csv")
metadata_df$sample_id <- paste(metadata_df$sequencing_id, "_dedup", sep = "")

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

# Print the missing sample_ids
print(missing_sample_ids)

# Export 
write.table(df.fr.tcge, file.path(outdir, paste0("tcge.median_fr_df.hg38.txt")), sep = "\t", row.names = FALSE)
save(df.fr.tcge, file=file.path(outdir, paste0("tcge.median_fr_df.hg38.rda")))

# Export cancer only 
df.fr.tcge.cancer <- df.fr.tcge %>% filter(cancer_type != "Normal")
df.fr.tcge.normal <- df.fr.tcge %>% filter(cancer_type == "Normal")

# Export 
write.table(df.fr.tcge.cancer, file.path(outdir, paste0("tcge.cancer.median_fr_df.hg38.txt")), sep = "\t", row.names = FALSE)
save(df.fr.tcge.cancer, file=file.path(outdir, paste0("tcge.cancer.median_fr_df.hg38.rda")))
save(df.fr.tcge.normal, file=file.path(outdir, paste0("tcge.normal.median_fr_df.hg38.rda")))


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

write.table(tcge_median, file.path(outdir, paste0("tcge.median.hg38.txt")), sep = "\t", row.names = FALSE)
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
cancer_subtypes <- unique(metadata_df$cancer_subtype)

# First, create a table of counts for each cancer subtype
subtype_counts <- table(metadata_df$cancer_subtype)

# Filter for cancer subtypes with at least 3 samples
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
  
  write.table(tcge_median, file.path(outdir, paste0("tcge.median.hg38", cancer_type, ".txt")), sep = "\t", row.names = FALSE)
  #  save(tcge_median, file=file.path(outdir, paste0("tcge.median.hg38", cancer_type, ".rda")))
  
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

cancer_types_delfi <- unique(metadata_df$cancer_type)
type_counts <- as.data.frame(table(metadata_df$cancer_type))

## Run the script for each type
for (cancer_type in cancer_types_delfi) {
  if (!is.na(cancer_type)) {
    print(run_script_for_cancer_type(cancer_type))
  }
}





## Do the cancer vs normal script comparison
df.fr.tcge <- df.fr.tcge %>%
  dplyr::mutate(CN_classifier = ifelse(cancer_type_corrected == "normal", "normal", "cancer"))
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
  
  write.table(tcge_median, file.path(outdir, paste0("tcge.median.hg38", cancer_type, ".txt")), sep = "\t", row.names = FALSE)
  #  save(tcge_median, file=file.path(outdir, paste0("tcge.median.hg38", cancer_type, ".rda")))
  
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


## Clean up 
rm(bins.list)
rm(tib.list)
