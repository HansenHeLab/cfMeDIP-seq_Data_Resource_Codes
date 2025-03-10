rm(list = ls())
setwd("/Users/yong/OneDrive - UHN/Projects/TCGE/cfEpigenomics/Resource/2_batch_norm/DMRs/PE")

library(dplyr)
library(BiocParallel)
library(parallel)
library(doParallel)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(wesanderson)
library(EnvStats)
library(annotatr)
library(gridExtra)
library(grid)

####################
{
## loading permutation and ploting functions
source("/Users/yong/OneDrive - UHN/Projects/TCGE/cfEpigenomics/Resource/2_batch_norm/DMRs/annotations_and_permutaion_functions.R")
  
# load annotation reginos
annot_list <-  readRDS("/Users/yong/OneDrive - UHN/Projects/TCGE/cfEpigenomics/Resource/2_batch_norm/DMRs/List_of_annotated_regions_for_CpGs_Promoter_Enhancer_hg38.RDS")
names(annot_list)

## bin bed infomation 
#bin_bed <- read.table("/Users/yong/OneDrive - UHN/Projects/TCGE/cfEpigenomics/Resource/2_batch_norm/merge/ALL_auto_bfilt_wCpG_bins.bed")
## read in ad granges
bg_bins <- read_regions(con = "/Users/yong/OneDrive - UHN/Projects/TCGE/cfEpigenomics/Resource/2_batch_norm/merge/ALL_auto_bfilt_wCpG_bins.bed",
                        genome = "hg38", rename_name = 'bin_id', rename_score = 'number_of_cpg')
resize(bg_bins, width(bg_bins) + 1, fix="end")   ## adjust ranges 

#bg_bins_reduce <-  reduce(bg_bins)

### bin id for DMRs
DMR_bins <- readRDS("FC_2_cancer_type_DESeq_DMRs_bin_IDs.RDS")
DMR_bins_type <- "DMR_bins_per_cancer_type"

## pan cancer 
#DMR_bins <- readRDS("FC_2_cancer_or_normal_DESeq_DMRs_bin_IDs.RDS")
#DMR_bins_type <- "DMR_bins_panCancer"


}


###########################
## run enrichment analysis
###########################
if(TRUE){
  
  names(DMR_bins)
  
  L <- length(DMR_bins)
  DMR_bins_enrich <- list()
  
  for (i in seq(1, L, 2))
  {
    print(paste0("Processing the prject ", i, " ..."))
    
    ############
    ## hyper DMRs
    hyper_dmr <-  DMR_bins[[i]]
    idx_hyper <- match(hyper_dmr, bg_bins$bin_id)
    hyper_dmr_gr <- bg_bins[idx_hyper, ]
    
    ## permutation test 
    hyper_enrich <- permutation(special = hyper_dmr_gr,
                                back.region = bg_bins,
                                annotations.list = annot_list,
                                annotation.region = names(annot_list)) 
    
    DMR_bins_enrich[[i]] <- hyper_enrich
    
    ############
    ## hypo DMRs
    hypo_dmr <- DMR_bins[[i + 1]]
    idx_hypo <- match(hypo_dmr, bg_bins$bin_id)
    hypo_dmr_gr <- bg_bins[idx_hypo, ]
    
    ## permutation test 
    hypo_enrich <- permutation(special = hypo_dmr_gr,
                               back.region = bg_bins,
                               annotations.list = annot_list,
                               annotation.region = names(annot_list)) 
    
    DMR_bins_enrich[[i + 1]] <- hypo_enrich
  }
  
  names(DMR_bins_enrich) <- names(DMR_bins)
  saveRDS(DMR_bins_enrich,  paste0(DMR_bins_type, "_enrichment_res.RDS"))
}


################################
## plotting enrichment analysis
################################
if(TRUE){
  DMR_bins_enrich <- readRDS(paste0(DMR_bins_type, "_enrichment_res.RDS"))
  
  enrich_pvals  <- data.frame() 
  enrich_z_scores <- data.frame() 
  
  L <- length(DMR_bins_enrich)
  for (i in seq(1, L, 2))
  {
    
    print(paste0("Plotting the prject ", i, " ..."))
    
    ## plotting hyper
    hyper_enrich <- DMR_bins_enrich[[i]]
    plot_name <- paste("DMR_enrichment_", names(DMR_bins_enrich)[i], ".pdf")
    plot_permutation(hyper_enrich, plot_name)
    
    
    ## plotting  hypo
    hypo_enrich <- DMR_bins_enrich[[i + 1]]
    plot_name <- paste("DMR_enrichment_", names(DMR_bins_enrich)[i + 1], ".pdf")
    plot_permutation(hypo_enrich, plot_name)
    
    ###############################
    ## retrieving pvals and z-score 
    N <- length(annot_list)
    hyper_pval <- hypo_pval <- vector()
    hyper_z <- hypo_z <- vector()
    
    for(j in 1:N){
      hyper_pval <- c(hyper_pval, hyper_enrich[[j]]$numOverlaps$pval)
      hypo_pval <- c(hypo_pval, hypo_enrich[[j]]$numOverlaps$pval)
      hyper_z <- c(hyper_z, hyper_enrich[[j]]$numOverlaps$zscore)
      hypo_z <- c(hypo_z, hypo_enrich[[j]]$numOverlaps$zscore)
    }
    enrich_pvals <- rbind(enrich_pvals, hyper_pval, hypo_pval)
    enrich_z_scores <- rbind(enrich_z_scores, hyper_z, hypo_z)
    
  }
  
  rownames(enrich_pvals) <- rownames(enrich_z_scores) <-  names(DMR_bins)
  colnames(enrich_pvals) <- colnames(enrich_z_scores) <- names(annot_list)
  saveRDS(enrich_pvals, paste0(DMR_bins_type, "_enrichment_pval.RDS"))
  saveRDS(enrich_z_scores, paste0(DMR_bins_type, "_enrichment_z_score.RDS"))
  
}
