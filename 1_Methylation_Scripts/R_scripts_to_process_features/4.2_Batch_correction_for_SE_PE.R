
rm(list = ls())
setwd("/cluster/projects/tcge/cell_free_epigenomics/processed_data/merge/raw_cnt")

library(sva)
library(DESeq2)
library(ggrepel)
library(ggplot2)
library(ggfortify)

##########################################################
## function to selected top 10K variable bins based on IQR
## and perform the PCA analysis
## need to omit rows with NA for PCA analysis
##########################################################
iqr_pca <- function(expr)
{
  idx_rc <- colSums(is.na(expr)) == nrow(expr)
  expr <- expr[, !idx_rc]       ## remove samples with NA for all bins
  
  expr <- na.omit(expr)        ##  remove bins with NA record
  iqr <- apply(expr, 1, IQR)  
  names(iqr) <- rownames(expr)
  
  idx_cons <- iqr == 0      ## remove constant bin
  iqr <- iqr[!idx_cons]
  
  iqr_10K <- sort(iqr, decreasing = T)[1: 10000]
  
  idx_ff <- match(rownames(expr), names(iqr_10K))
  
  ## most variable 10K bins for the PCA analysis 
  expr_f <- expr[!is.na(idx_ff), ]
  expr_f_pca <- prcomp(t(expr_f), center = T, scale.=T)
  return(expr_f_pca)
}

#######################
#### ALL samples
#### PE or SE
#######################
## TCGE_cfDM_samples_only_fixed_after_QC_10k         ## random sampling 10k bin for testing 
#name_prefix <- c("TCGE_cfDM_samples_only_fixed_after_QC_10K")  
name_prefix <- c("TCGE_cfDM_samples_only_fixed_SE", "TCGE_cfDM_samples_only_fixed_PE")

L <- length(name_prefix)

for (i in 1:L){

  s_info <- read.csv(paste0("../", name_prefix[i], ".csv"))
  cnt <-  readRDS(paste0(name_prefix[i], "_auto_bfilt_wCpG.RDS"))

  ## matching s_info
  idx_s <- match(colnames(cnt), s_info$sequencing_id)
  cnt <- cnt[, !is.na(idx_s)]
  s_info_s <- s_info[idx_s[!is.na(idx_s)], ]

  ########################################
  ## apply combat + DESEQ2 for all samples
  ## apply two round batch correction
  if(TRUE){
    ## missing age and sex are not considered in coar_mod
    print("Applying round 1 Combat ajustment per cancer subtype ...")
    subtypes <- unique(s_info_s$cancer_subtype)
    N <- length(subtypes)
    cnt_adj <- cnt
    for(k in 1:N)
    {
      print(k)
      idx_ss <- s_info_s$cancer_subtype == subtypes[k]
      s_info_ss <- s_info_s[idx_ss, ]
      
      ## samples in all cnt
      idx_cnt <- match(s_info_ss$sequencing_id, colnames(cnt))
      cnt_ss <- cnt[, idx_cnt]
      
      batch_s = s_info_ss$processing_batch
      print( subtypes[k])
      print(table(batch_s))
      
      ## (sum(is.na(batch_s)) == length(batch_s)) 
      if( length(unique(batch_s)) == 1){
        next;
      } else {
        batch_cnt <- table(batch_s)
        idx_1 <- batch_cnt == 1
        
        if(sum(idx_1) > 1){
          idx_1s <- match(batch_s, names(batch_cnt)[idx_1])
          batch_s[!is.na(idx_1s)] <- paste(names(batch_cnt)[idx_1], collapse = "_")
        } else if (sum(idx_1) == 1) {
          ## 1 sample per batch, mering the nearest batch by names
          idx_1c <- which( batch_cnt == 1)
          ## mergeing to index 
          if(idx_1c == 1) idx_mer = 2 else idx_mer <- idx_1c -1
          idx_1s <- match(batch_s, names(batch_cnt)[c(idx_1c, idx_mer)])
          batch_s[!is.na(idx_1s)] <- paste(names(batch_cnt)[idx_1c], names(batch_cnt)[idx_mer], sep = "_")
        }
        
        ## combat adjusting 
        cnt_ss_adj <- ComBat_seq(cnt_ss, batch = batch_s, group = NULL)
        cnt_adj[, idx_cnt] <- cnt_ss_adj
      }
      
    }
    
    saveRDS(cnt_adj, file =  paste0(name_prefix[i], "_auto_bfilt_wCpG_Combat_per_cancer_subtype_ajusted_cnt.RDS"))
    
    
    ################################################################
    ## round 2 adjustment for SE PE, and reserve the main cancertype
    if (FALSE){
    print("Applying round 2 Combat ajustment for SE vs PE ...")
    batch_s = rep("SE", nrow(s_info_s))
    idx_pe <- s_info_s$project_id != "TCGE-CFMe-MCA"
    batch_s[idx_pe] = "PE"
    
    ## the cancer subtype are confounded ...
    bio_mat <-  model.matrix(~ cancer_type, data = s_info_s)
    cnt_adj <- ComBat_seq(cnt_adj, batch = batch_s, group = NULL, covar_mod = bio_mat)
    saveRDS(cnt_adj, file =  paste0(name_prefix[i], "_auto_bfilt_wCpG_Combat_per_cancer_subtype_and_SE_PE_ajusted_cnt.RDS"))
    }
    
    print("Applying DESeq2 normalization ...")
    ## deseq2 normalization
    dds <- DESeqDataSetFromMatrix(cnt_adj, s_info_s, ~1)
    dds <- estimateSizeFactors(dds)
    cnt_adj_norm <- counts(dds, normalized = T) 
    cnt_adj_norm_10k <- iqr_pca(cnt_adj_norm)
    saveRDS(cnt_adj_norm, file =  paste0(name_prefix[i], "_auto_bfilt_wCpG_Combat_per_cancer_subtype_and_DESeq2_norm.RDS"))
    saveRDS(cnt_adj_norm_10k, file =  paste0(name_prefix[i], "_auto_bfilt_wCpG_Combat_per_cancer_subtype_and_DESeq2_norm_Top10Kbin_PCA.RDS"))
  }
  
  
  print("Done !!")
  
}

