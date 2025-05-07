rm(list = ls())

## similar script was applied to SE samples !!
setwd("/cluster/projects/tcge/cell_free_epigenomics/processed_data/merge/DMRs/PE_Specific")

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

## remove LFS 
idx_LFS <- s_info$project_id == "TCGE-CFMe-LFS"
s_info <- s_info[!idx_LFS, ]

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

##############################
### DMRs between cancer types
if(TRUE){
  dds <- DESeqDataSetFromMatrix(cnt, s_info_s, ~ cancer_type - 1)   ## remove intercept
  dds <- DESeq(dds)
  saveRDS(dds, "DMRs_Cancer_type_specific.RDS")
  rm(cnt)
}

## with control of age and sex
##  design formula cannot contain NA !!
#dds <- DESeqDataSetFromMatrix(cnt, s_info_s, ~ age + sex + cancer_type)
#dds <- DESeq(dds)
#saveRDS(dds, "DMRs_Cancer_type_vs_Normal_with_control_age_sex.RDS")

## matching s_info and perform DMRs analysis with DESeq2

################################################################
## pull out all comparison for each cancer type against the rest
## LFS was exclued 
################################################################

# dds <-readRDS("DMRs_Cancer_type_specific.RDS")

## to check what coefficients DESeq estimated
contrast_names <- resultsNames(dds)
L <- length(contrast_names)

for(fc in c(1.5, 2))
{
  fdr_c <- 0.05                     ## FDR cutoff for DMRs
  log2fc_c <- log2(fc)              ## abs log2FC cutoff for DMRs
  DMR_bins <- list()
  k = 1

  for (i in 1:L)
  {
    
  ## each cancer subtype vs all other left types
    
  res <- results(dds, contrast=list(contrast_names[i], contrast_names[-i]),
                      listValues=c(1, -1/(L-1)))
                   
  res_name <- paste0(contrast_names[i], "_specific")
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

  saveRDS(DMR_bins, paste0("FC_", fc, "_DESeq_cancer_specific_DMRs_bin_IDs.RDS"))
}

