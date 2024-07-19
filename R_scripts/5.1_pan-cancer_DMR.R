rm(list = ls())

## similar script was applied to SE samples !!
setwd("/cluster/projects/tcge/cell_free_epigenomics/processed_data/merge/DMRs/PE")

library(sva)
library(DESeq2)
library(ggrepel)
library(ggplot2)
library(gplots)

#########################################################
## selected samples  Folders with PCA results
##########################################################
## matching s_info and perform DMRs analysis with DESeq2
cnt <- readRDS("/cluster/projects/tcge/cell_free_epigenomics/processed_data/merge/raw_cnt/TCGE_cfDM_samples_only_fixed_PE_auto_bfilt_wCpG_Combat_per_cancer_subtype_ajusted_cnt.RDS")
#norm_cnt <- readRDS("/cluster/projects/tcge/cell_free_epigenomics/processed_data/merge/raw_cnt/TCGE_cfDM_samples_only_fixed_PE_auto_bfilt_wCpG_Combat_per_cancer_subtype_and_DESeq2_norm.RDS")
bin_info <- read.table("/cluster/projects/tcge/cell_free_epigenomics/processed_data/merge/raw_cnt/ALL_auto_bfilt_wCpG_bins.bed")

s_info <- read.csv("../TCGE_cfDM_samples_only_fixed_after_QC_PE.csv")

## project/cancer_type order and colors
s_info$project_id <- factor(s_info$project_id,
                            levels = c("TCGE-CFMe-HBC", "TCGE-CFMe-MCA", "TCGE-CFMe-BCA", "TCGE-CFMe-PRAD",
                                           "TCGE-CFMe-HNSC", "TCGE-CFMe-SCLC", "TCGE-CFMe-UM", "TCGE-CFMe-AML",
                                           "TCGE-CFMe-LFS"))


s_info$cancer_type <- factor(s_info$cancer_type,
                             levels = c("Normal", "Brain Cancer", "Lung Cancer", "Prostate Cancer",
                                            "Blood Cancer", "Pancreatic Cancer", "Eye Cancer",
                                            "Head and Neck Cancer", "Breast Cancer", "Colorectal Cancer",
                                            "Bladder Cancer", "Renal Cancer",
                                            "LFS Survivor", "LFS Previvor", "LFS Positive"))


### all cancer vs normal
## to exclude all LFS samples, which can be used as independent validation cohort

cancer_or_normal <- rep("Cancer", nrow(s_info))
idx_n <- s_info$cancer_type == "Normal"
cancer_or_normal[idx_n] <- "Normal"

idx_LFS <- s_info$project_id == "TCGE-CFMe-LFS"
cancer_or_normal[idx_LFS] <- "LFS"
s_info <- cbind(s_info, cancer_or_normal)

s_info$cancer_or_normal <- factor(s_info$cancer_or_normal,
                                  levels = c("Normal", "Cancer", "LFS"))

##################################
## matching ajusted cnt and s_info
##################################

idx_s <- match(colnames(cnt), s_info$sequencing_id)

## requiring bin with more than 5 CpGs for DMR analysis
idx_bin <- bin_info[, 5] > 5
cnt <- cnt[idx_bin, !is.na(idx_s)]

#norm_cnt <- norm_cnt[, !is.na(idx_s)]
s_info_s <- s_info[idx_s[!is.na(idx_s)], ]

## remove empty levels 
s_info_s$project_id <- droplevels(s_info_s$project_id)   
s_info_s$cancer_type <- droplevels(s_info_s$cancer_type)   
s_info_s$cancer_or_normal <- droplevels(s_info_s$cancer_or_normal)

##############################
### DMRs between cancer types
if(TRUE){
## pair-wise comparison
#dds <- DESeqDataSetFromMatrix(cnt, s_info_s, ~ cancer_type)
#dds <- DESeq(dds)
#saveRDS(dds, "DMRs_Cancer_type_vs_Normal.RDS")

## all cancer samples vs normal sample
dds <- DESeqDataSetFromMatrix(cnt, s_info_s, ~ cancer_or_normal)
dds <- DESeq(dds)
saveRDS(dds, "DMRs_all_Cancer_vs_Normal.RDS")

rm(cnt)
}

## with control of age and sex
##  design formula cannot contain NA !!
#dds <- DESeqDataSetFromMatrix(cnt, s_info_s, ~ age + sex + cancer_type)
#dds <- DESeq(dds)
#saveRDS(dds, "DMRs_Cancer_type_vs_Normal_with_control_age_sex.RDS")

## matching s_info and perform DMRs analysis with DESeq2

#########################################
## pull out all comparison against Normal
#########################################
dds_list <- c("DMRs_all_Cancer_vs_Normal.RDS", "DMRs_Cancer_type_vs_Normal.RDS")
for (dds_file in dds_list)
{

  dds <-readRDS(dds_file)

## cancer types vs normal
if(dds_file == "DMRs_all_Cancer_vs_Normal.RDS"){
  conditions <- levels(s_info_s$cancer_or_normal)
  contrast <- "cancer_or_normal"
} else {
  conditions <- levels(s_info_s$cancer_type)
  contrast <- "cancer_type"
}

conditions             ## Normal at the end
L <- length(conditions)

## all cancer vs normal

for(fc in c(1.5, 2))
{
  fdr_c <- 0.05                     ## FDR cutoff for DMRs
  log2fc_c <- log2(fc)              ## abs log2FC cutoff for DMRs
  DMR_bins <- list()
  k = 1

  for (i in 2:L)
  {
  res <- results(dds, contrast=c(contrast, conditions[i], conditions[1]))
  res_name <- paste0(conditions[i], " vs ",conditions[1])
  print(paste0(res_name, "...."))
  ##############
  ## output DMRs
  {
    idx_up <- res$log2FoldChange > log2fc_c & res$padj < fdr_c      ## Hyper DMRs
    idx_do <- res$log2FoldChange < -log2fc_c & res$padj < fdr_c     ## Hypo DMRs
    idx_up[is.na(idx_up)] <- FALSE
    idx_do[is.na(idx_do)] <- FALSE

    DMR_bins[[k]] <- rownames(res)[idx_up]
    DMR_bins[[k + 1]] <- rownames(res)[idx_do]
    names(DMR_bins)[k : (k + 1)] <- paste0(res_name , c("_Hyper", "_Hypo"))
    k = k + 2
    write.csv(res[idx_up, ], file = paste0("FC_", fc, "_", res_name, "_DESeq_Hyper_DMRs.csv"))
    write.csv(res[idx_do, ], file = paste0("FC_", fc, "_", res_name, "_DESeq_Hypo_DMRs.csv"))

  }

  ###################
  ## DEG volcano plot
  {
    sigGroup <- rep("no-diff", nrow(res))
    sigGroup[idx_up] <- "Hyper"
    sigGroup[idx_do] <- "Hypo"
    #main_title <- paste0(res_name, "_Fold_", fc)
    main_title <- res_name
    sub_title = paste0("Hypo: N= ", sum(idx_do), ";  Hyper: N=", sum(idx_up))

    dat <- data.frame(log2FC = res$log2FoldChange, FDR =  -log10(res$padj), sigGroup)
    g  <- ggplot(dat, aes(log2FC, FDR, col = sigGroup))
    g <- g + geom_point() + scale_color_manual(values = c("#d6604d","#4575b4", "#d9d9d9"))
    g <- g + labs(title = main_title,  subtitle = sub_title, x = "log2(FC)", y = "-log10(FDR)")
    g <- g + geom_hline(yintercept = -log10(0.05), linetype="dashed", color = "red")
    g <- g + geom_vline(xintercept = c(log2fc_c, -log2fc_c), linetype= "dashed", color = "red")
    g <- g  + theme_classic() + theme(legend.position = "none")
    ggsave(file = paste0("FC_", fc, "_", res_name, "_DMRs_Volcano_plot.png"), width = 4, height = 4, units = "in", dpi = 600)
  }

  }

  saveRDS(DMR_bins, paste0("FC_", fc, "_", contrast, "_DESeq_DMRs_bin_IDs.RDS"))
}
}
