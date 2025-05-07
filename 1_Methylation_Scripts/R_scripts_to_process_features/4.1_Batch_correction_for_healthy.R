
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
#### ALL normal samples
#### PE + SE
#######################
name_prefix <- c("TCGE_cfDM_samples_only_fixed_Normal_SE")
                 #"TCGE_cfDM_samples_only_fixed_Normal", 
                 #"TCGE_cfDM_samples_only_fixed_Normal_PE",
                 #"TCGE_cfDM_samples_only_fixed_AML")

L <- length(name_prefix)

for (i in 1:L){

  s_info <- read.csv(paste0("../", name_prefix[i], ".csv"))
  cnt <-  readRDS(paste0(name_prefix[i], "_auto_bfilt_wCpG.RDS"))

  ## matching s_info
  idx_s <- match(colnames(cnt), s_info$sequencing_id)
  cnt <- cnt[, !is.na(idx_s)]
  s_info_s <- s_info[idx_s[!is.na(idx_s)], ]

  ##################################
  ## apply DESeq2 for normalization 
  ## without design matrix 
  if(TRUE){
    print("Applying DESeq2 normalization ...")
  
    dds <- DESeqDataSetFromMatrix(cnt, s_info_s, ~1)
    dds <- estimateSizeFactors(dds)
    cnt_norm<- counts(dds, normalized = T)
    cnt_norm_10k <- iqr_pca(cnt_norm)
  
    saveRDS(cnt_norm, file =  paste0(name_prefix[i], "_auto_bfilt_wCpG_DESeq2_norm.RDS"))
    saveRDS(cnt_norm_10k, file = paste0(name_prefix[i], "_auto_bfilt_wCpG_DESeq2_norm_Top10Kbin_PCA.RDS")) 
  }  

  ##########################
  ## apply combat + DESEQ2
  if(TRUE){
    print("Applying 1 round Combat + DESeq2 normalization ...")
    
    ## merging batch with single sample
    batch_s = s_info_s$processing_batch
    print(table(batch_s))
   
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
     
    cnt_adj <- ComBat_seq(cnt, batch = batch_s, group=NULL)

    ## deseq2 normalization
   dds <- DESeqDataSetFromMatrix(cnt_adj, s_info_s, ~1)
   dds <- estimateSizeFactors(dds)
   cnt_adj_norm <- counts(dds, normalized = T) 
   cnt_adj_norm_10k <- iqr_pca(cnt_adj_norm)
   saveRDS(cnt_adj_norm, file =  paste0(name_prefix[i], "_auto_bfilt_wCpG_1Round_Combat_DESeq2_norm.RDS"))
   saveRDS(cnt_adj_norm_10k, file =  paste0(name_prefix[i], "_auto_bfilt_wCpG_1Round_Combat_DESeq2_norm_Top10Kbin_PCA.RDS"))
  }
  
  ################################
  ## apply 2 round combat + DESEQ2
  if(i == 1){
    #######################################################
    print("Applying round 1st Combat ajustment per study ...")
    ## merging batch with single sample
    batch_s = s_info_s$processing_batch
    print(table(batch_s))
    
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
    cnt_adj <- ComBat_seq(cnt, batch = batch_s, group=NULL)
    
    ################################################################
    ## round 2 adjustment for SE PE
    print("Applying round 2nd Combat ajustment for SE vs PE ...")
    batch_s = rep("SE", nrow(s_info_s))
    idx_pe <- s_info_s$project_id != "TCGE-CFMe-MCA"
    batch_s[idx_pe] = "PE"

    ## the cancer subtype are confounded ...
    #bio_mat <-  model.matrix(~ cancer_type, data = s_info_s)
    cnt_adj <- ComBat_seq(cnt_adj, batch = batch_s, group = NULL)
    
    #######################
    ## deseq2 normalization
    print("Applying DESeq2 normalization ...")
    dds <- DESeqDataSetFromMatrix(cnt_adj, s_info_s, ~1)
    dds <- estimateSizeFactors(dds)
    cnt_adj_norm <- counts(dds, normalized = T) 
    cnt_adj_norm_10k <- iqr_pca(cnt_adj_norm)
    saveRDS(cnt_adj_norm, file =  paste0(name_prefix[i], "_auto_bfilt_wCpG_2round_Combat_DESeq2_norm.RDS"))
    saveRDS(cnt_adj_norm_10k, file =  paste0(name_prefix[i], "_auto_bfilt_wCpG_2round_Combat_DESeq2_norm_Top10Kbin_PCA.RDS"))
  }
  
  print("All Done !!")
  
}

