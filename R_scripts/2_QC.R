rm(list = ls())
setwd("/Users/yong/OneDrive_UHN/Projects/TCGE/cfEpigenomics/Resource/1_QC/1_QC_cutoffs_fixed")

library(dplyr)
library(ggplot2)
library(ggpubr)
#library(RColorBrewer)
#library(corrplot)
#library(matrixStats)
#library(gplots)   ## for heatmap2

## cell-free samples only sample meta information 
s_info <- read.csv("/Users/yong/OneDrive - UHN/Projects/TCGE/cfEpigenomics/Resource/0_sample_meta_info/meta_info_summary_for_fixed/TCGE_cfDM_samples_only_fixed.csv")

## order by samples types
s_info$project_id <- factor(s_info$project_id, 
                            levels = c("TCGE-CFMe-MCA", "TCGE-CFMe-BCA", "TCGE-CFMe-HNSC", "TCGE-CFMe-PRAD",
                                           "TCGE-CFMe-AML", "TCGE-CFMe-SCLC", "TCGE-CFMe-UM", 
                                           "TCGE-CFMe-HBC", "TCGE-CFMe-LFS"))

project_col <- c('#7fcdbb',  '#984ea3', '#ff7f00',  '#377eb8', '#f781bf', '#807dba','#a65628', '#4daf4a', '#999999')



#########################
## aggregated qc reports
########################
{
qc_path <- "/Users/yong/OneDrive - UHN/Projects/TCGE/cfEpigenomics/Resource/1_QC/0_QC_reports"
qc_reports <- list.files(path = qc_path, pattern = "aggr_qc_report.csv", 
                         recursive = T, full.names = T)

L <- length(qc_reports)
qc <- read.csv(qc_reports[1])

for(i in 2:L)
{
  tmp <- read.csv(qc_reports[i])
  qc  <- rbind(qc, tmp)
}

## mathing with curated samples
idx_s <- match(s_info$sequencing_id, qc$sample)

s_info_s <- s_info[!is.na(idx_s), ]
qc_s <- qc[idx_s[!is.na(idx_s)], ]
qc_f <- cbind(s_info_s, qc_s) 

table(qc_f$project_id)
}

######################
## s_info + QC metircs
## plot all qc metric
######################
{
dat <- qc_f[, -19]                 ## rm redundant columns of group 
col_s <- 19                 ## which column starts

## correct typo for the 
#colnames(dat)[26] <- "coverage_pctCpGw1Read"  ## corrected from the rawfiles

## add cpg coverage 
coverage_pctCpGwReads <- 100 - dat$coverage_pctCpGwoRead
median(coverage_pctCpGwReads)   ## combined SE and PE

idx_pe <- dat$project_id != "TCGE-CFMe-MCA"
median(coverage_pctCpGwReads[idx_pe])

## usable read rates
range()
boxplot(list(dat$usable_reads_depth_pct[!idx_pe], dat$usable_reads_depth_pct[idx_pe]))
t.test(dat$usable_reads_depth_pct[!idx_pe], dat$usable_reads_depth_pct[idx_pe])

## add IP specificity
specificity_pctReadsWCpG <- 100 - dat$coverage_pctReadsWoCpG

dat <- cbind(dat, coverage_pctCpGwReads, specificity_pctReadsWCpG)
write.csv(dat, file = "TCGE_cfDM_samples_only_fixed_with_QC_metrics.csv", row.names = F)

## 
qc_quantiles <- data.frame()
colnames(dat)

## by project 
for(i in col_s:ncol(dat))
{
  name_t <- colnames(dat)[i]

  ####### Quantiles ################
  ## 2.5%, 10%, 25%, 50%, 75%, 97.5%
  ## 
  
  if (i == 27 | i == 28 | i == 29 | i == 30 | i == 32 | i == 33| i == 34 | i == 35 | i == 38 | i == 39){
  ## remove single-end reads for coverage and specificity percentiles
  idx_mca <- dat$project_id == "TCGE-CFMe-MCA"
  quantile_t <- quantile(dat[!idx_mca, i], probs = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.975), na.rm = T)
 
  ## sub tile  with out MCA 
  sub_t <- paste0("Percentiles: 2.5% (Red) = ", round(quantile_t[1], 2), 
                  "; 50% (Green) = ", round(quantile_t[4], 2),
                  "; 97.5% (Red) = ", round(quantile_t[6], 2))
  
  dat_s <- dat[!idx_mca, ]
  dat_s$project_id <- droplevels(dat_s$project_id)
  
  project_col_s <- c('#984ea3', '#ff7f00',  '#377eb8', '#f781bf', '#807dba','#a65628', '#4daf4a', '#999999')
    
  ## per project   
  g <- ggboxplot(dat_s, x = "project_id", y = name_t, color = "project_id", add = "jitter", palette = project_col_s)
  g <- g + geom_hline(yintercept = quantile_t[c(1, 4, 6)], linetype = "dashed", color = c("red", "darkgreen", "red"))
  g <- g + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust= 1))
  g <- g + labs(subtitle = sub_t, x = "") + theme(legend.position = "none")
  ggsave(paste0("QC_per_project_", name_t, ".png"), width = 8.5, height = 4, units = "in", dpi = 300)
  
  ## per cancer type
  g <- ggboxplot(dat_s, x = "cancer_type", y = name_t, fill = "project_id", palette = project_col_s)
  g <- g + geom_hline(yintercept = quantile_t[c(1, 4)], linetype = "dashed", color = c("red", "darkgreen"))
  g <- g + labs(subtitle = sub_t)
  g <- g + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  ggsave(paste0("QC_per_cancer_type_", name_t, ".png"), width = 12, height = 6, units = "in", dpi = 300)
  
  qc_quantiles  <- rbind(qc_quantiles, quantile_t)
  
  
   } else {
  quantile_t <- quantile(dat[, i], probs = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.975), na.rm = T)
  
  ## sub tile 
  sub_t <- paste0("Percentiles: 2.5% (Red) = ", round(quantile_t[1], 2), 
                  "; 50% (Green) = ", round(quantile_t[4], 2),
                  "; 97.5% (Red) = ", round(quantile_t[6], 2))
  
  ## per project   
  g <- ggboxplot(dat, x = "project_id", y = name_t, color = "project_id", add = "jitter", palette = project_col)
  g <- g + geom_hline(yintercept = quantile_t[c(1, 4, 6)], linetype = "dashed", color = c("red", "darkgreen", "red"))
  g <- g + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust= 1))
  g <- g + labs(subtitle = sub_t, x = "") + theme(legend.position = "none")
  ggsave(paste0("QC_per_project_", name_t, ".png"), width = 8.5, height = 4, units = "in", dpi = 300)
  
  ## per cancer type
  g <- ggboxplot(dat, x = "cancer_type", y = name_t, fill = "project_id", palette = project_col)
  g <- g + geom_hline(yintercept = quantile_t[c(1, 4)], linetype = "dashed", color = c("red", "darkgreen"))
  g <- g + labs(subtitle = sub_t)
  g <- g + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  ggsave(paste0("QC_per_cancer_type_", name_t, ".png"), width = 12, height = 6, units = "in", dpi = 300)
  
  qc_quantiles  <- rbind(qc_quantiles, quantile_t)
  
  }
  
  
}

######################
## output QC quantiles 
colnames(qc_quantiles) <- c("Percentile_2.5%", "Percentile_10%", "Percentile_25%", "Percentile_50%", "Percentile_75%", "Percentile_97.5%")
rownames(qc_quantiles) <- colnames(dat)[col_s:ncol(dat)]
write.csv(qc_quantiles, file = "QC_metrics_Percentiles.csv")

}


##########################################
## replot selected metrics for main figure 
#########################################
{

  ## enrichment score GoGe
  g <- ggboxplot(dat, x = "project_id", y = "enrichment_GoGe", color = "project_id", add = "jitter", palette = project_col)
  g <- g + geom_hline(yintercept = c(1.70, 1.95), linetype = "dashed", color = c("red", "darkgreen"))
  g <- g + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust= 1))
  g <- g + labs( x = "") + theme(legend.position = "none")
  ggsave(paste0("QC_per_project_enrichment_GoGe_pub.png"), width = 9, height = 4, units = "in", dpi = 300)
  
  
  ## for legend only 
  if(i == col_s){
    g <- ggbarplot(dat, x = "project_id", y = name_t, fill = "project_id", color = "project_id", palette = project_col)
    g <- g + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust= 1))
    g <- g + labs(subtitle = sub_t, x = "") + theme(legend.position = "right")
    ggsave(paste0("With_legend_QC_per_project_", name_t, ".png"), width = 8.5, height = 6, units = "in", dpi = 300)
  }
  
  
}  


###############
## QC filtering 
###############
{
  idx_k <- dat$enrichment_relH > 3.0 & 
           dat$enrichment_GoGe > 1.7 &
           dat$saturation_maxEstCor > 0.9
  
  sum(idx_k) / length(idx_k)
  
  ## saturation is high for all samples
  #dat$saturation_maxEstCor > qc_quantiles[13, 1] &
  
  s_info_out <- dat[idx_k, ]   
  
  write.csv(s_info_out, "TCGE_cfDM_samples_only_fixed_after_QC.csv", row.names = F)
  
  ## excluded samples
  a <-  dat[!idx_k, ] 
  table(a$project_id)
  
  #TCGE-CFMe-HBC  TCGE-CFMe-MCA  TCGE-CFMe-BCA TCGE-CFMe-PRAD 
  #1             10              9             15 
  #TCGE-CFMe-HNSC TCGE-CFMe-SCLC   TCGE-CFMe-UM  TCGE-CFMe-AML 
  #4              2              6              2 
  #TCGE-CFMe-LFS 
  #51 
  
  ## split SE and PE samples
  idx_se <- s_info_out$project_id == "TCGE-CFMe-MCA"
  write.csv(s_info_out[idx_se, ],  "TCGE_cfDM_samples_only_fixed_after_QC_SE.csv", row.names = F)
  write.csv(s_info_out[!idx_se, ], "TCGE_cfDM_samples_only_fixed_after_QC_PE.csv", row.names = F)
  
  ## Normal samples
  #idx_norm <- s_info_out$cancer_subtype == "Normal"
  idx_norm <- s_info_out$cancer_subtype == "Healthy"
  write.csv(s_info_out[idx_norm, ],  "TCGE_cfDM_samples_only_fixed_after_QC_Normal.csv", row.names = F)
  write.csv(s_info_out[idx_norm & !idx_se, ],  "TCGE_cfDM_samples_only_fixed_after_QC_Norma_PE.csv", row.names = F)
  
  ## AML study with all 3 replicats 
  idx_aml <- s_info_out$project_id  == "TCGE-CFMe-AML"
  write.csv(s_info_out[idx_aml, ],  "TCGE_cfDM_samples_only_fixed_after_QC_AML.csv", row.names = F)

}


########################
## QC metric correaltion
########################
{
library(corrplot)
library(matrixStats)
library(gplots)   ## for heatmap2

#####################################
## for SE  samples after QC filtering 
qc <- read.csv("/Users/yong/OneDrive - UHN/Projects/TCGE/cfEpigenomics/Resource/1_QC/1_QC_cutoffs_fixed/TCGE_cfDM_samples_only_fixed_with_QC_metrics.csv")

idx_se <- qc$project_id == "TCGE-CFMe-MCA"
qc_se <- qc[idx_se, ]

col_s <- c(19:20, 23:26, 31, 36:37)                 ## which column starts
qc_m <- data.matrix(qc_se[, col_s])
cor_m <- cor(qc_m, method = "spearman")

png("TCGE_cfDM_samples_only_fixed_after_QC_SE_Corrplot.png", width = 7, height = 7, units = "in", res = 300)
#corrplot(cor_m) #, order = 'hclust')
#corrplot(cor_m, method = 'square', order = 'FPC', diag = FALSE, tl.col = "Black")
corrplot(cor_m, method = 'square', type = 'upper', order = 'hclust', tl.col = 'black', cl.ratio = 0.2, diag = FALSE)
dev.off()


#####################################
## for PE  samples after QC filtering 

qc_pe <- qc[!idx_se, ]
col_s <- c(19:39)                 ## which column starts
qc_m <- data.matrix(qc_pe[, col_s])
cor_m <- cor(qc_m, method = "spearman")

png("TCGE_cfDM_samples_only_fixed_after_QC_PE_Corrplot.png", width = 7, height = 7, units = "in", res = 300)
#corrplot(cor_m) #, order = 'hclust')
#corrplot(cor_m, method = 'square', order = 'FPC', diag = FALSE)
corrplot(cor_m, method = 'square', type = 'upper', order = 'hclust', tl.col = 'black', cl.ratio = 0.2, diag = FALSE)
dev.off()
}
