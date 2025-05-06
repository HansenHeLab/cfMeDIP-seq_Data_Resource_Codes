rm(list = ls())
library(ggfortify)   ## PCA plot: autoplot
library(umap)
library(dplyr)
library(ggplot2)
library(gplots)

#########################################################
## selected samples  Folders with PCA results
##########################################################
#s_info <- read.csv("/Users/yong/OneDrive - UHN/Projects/TCGE/cfEpigenomics/Resource/1_QC/1_QC_cutoffs_fixed/TCGE_cfDM_samples_only_fixed_after_QC.csv")
s_info <- read.csv("/Users/yong/OneDrive - UHN/Projects/TCGE/cfEpigenomics/Resource/1_QC/1_QC_cutoffs_fixed/TCGE_cfDM_samples_only_fixed_with_QC_metrics.csv")

## project/cancer_type order and colors
s_info$project_id <- factor(s_info$project_id, 
                            levels = c("TCGE-CFMe-MCA", "TCGE-CFMe-BCA", "TCGE-CFMe-HNSC", "TCGE-CFMe-PRAD",
                                       "TCGE-CFMe-AML", "TCGE-CFMe-SCLC", "TCGE-CFMe-UM", 
                                       "TCGE-CFMe-HBC", "TCGE-CFMe-LFS"))

project_col <- c('#7fcdbb',  '#984ea3', '#ff7f00',  '#377eb8', '#f781bf', '#807dba','#a65628', '#4daf4a', '#999999')


s_info$cancer_type <- factor(s_info$cancer_type, 
                             levels = c("Normal", "Brain Cancer", "Lung Cancer", "Prostate Cancer",
                                            "Blood Cancer", "Pancreatic Cancer", "Eye Cancer",
                                            "Head and Neck Cancer", "Breast Cancer", "Colorectal Cancer",
                                            "Bladder Cancer", "Renal Cancer",
                                            "LFS Survivor", "LFS Previvor", "LFS Positive"))

cancer_type_col <- c('#33a02c','#1f78b4','#b2df8a','#a6cee3','#fb9a99',
                         '#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a',
                         '#fb6a4a','#b15928','#bdbdbd','#969696','#737373')

## color for TCGE-CFMe-AML
AML_col <- c('#e7298a','#d95f02','#7570b3')

################
## PCA resutls 
folder <- c("/Users/yong/OneDrive - UHN/Projects/TCGE/cfEpigenomics/Resource/2_batch_norm/merge/raw_cnt",
            "/Users/yong/OneDrive - UHN/Projects/TCGE/cfEpigenomics/Resource/2_batch_norm/merge/raw_cnt_deseq2",
            "/Users/yong/OneDrive - UHN/Projects/TCGE/cfEpigenomics/Resource/2_batch_norm/merge/raw_cnt_combat_deseq2",
            "/Users/yong/OneDrive - UHN/Projects/TCGE/cfEpigenomics/Resource/2_batch_norm/merge/rpkm",
            "/Users/yong/OneDrive - UHN/Projects/TCGE/cfEpigenomics/Resource/2_batch_norm/merge/rms_medestrand",
            "/Users/yong/OneDrive - UHN/Projects/TCGE/cfEpigenomics/Resource/2_batch_norm/merge/beta_qsea")
FL <- length(folder)

title_names <- c("Raw count", "DESeq2 normalized count", "Combat-seq + DESeq2 normalized count",
            "RPKM and/or FPKM", "MEDEStrand estimated rms", "QSEA estimated beta")

for(k in 1:FL)
{

setwd(folder[k])
#########################################################
## read in pca results based on the top 10k variable bins
##########################################################
{
  
files <- list.files(path = ".", pattern = "RDS")
files_names <- tools::file_path_sans_ext(files)
L <- length(files)

for (i in 1:L ){
  ## PCs and match colors
  {
    print(paste0("Processing file ", i, " in the folder ", k, "...."))
    pca <- readRDS(files[i])
    s_pca <- summary(pca)
    name_prefix <- files_names[i]
    
    idx_s <- match(row.names(pca$x), s_info$sequencing_id)
    s_info_s <- s_info[idx_s, ]
    s_info_s$project_id <- droplevels(s_info_s$project_id)   ## remove empty levels
    s_info_s$cancer_type <- droplevels(s_info_s$cancer_type)   ## remove empty levels
    
    ## matching colors for project and cancer type
    idx_p <- match(levels(s_info_s$project_id), levels(s_info$project_id)) 
    project_col_s <- project_col[idx_p]
    
    idx_ct <- match(levels(s_info_s$cancer_type), levels(s_info$cancer_type)) 
    cancer_type_col_s <- cancer_type_col[idx_ct]
  }
  
  if (TRUE){ 
  ## draw PCA plots with different labels
  {
    g <- autoplot(pca, data = s_info_s, label = F, colour = 'project_id')
    g <- g + labs(title = title_names[k])
    g <- g + scale_colour_manual(values = project_col_s) + theme_classic()
    ggsave(paste0(name_prefix, "_labeled_by_project_id.pdf"), width = 5.5, height = 3.5)
    
    g <- autoplot(pca, data = s_info_s, label = F, colour = 'cancer_type')
    g <- g + labs(title = title_names[k])
    g <- g + scale_colour_manual(values = cancer_type_col_s) + theme_classic()
    ggsave(paste0(name_prefix, "_labeled_by_cacer_type.pdf"), width = 6, height = 4)
    
    ## for TCGE-CFMe-AML  only 
    if (i == 1) {
    g <- autoplot(pca, data = s_info_s, label = F, colour = 'group', size = 2)
    g <- g + labs(title = title_names[k])
    g <- g + scale_colour_manual(values = AML_col) + theme_classic()
    ggsave(paste0(name_prefix, "_labeled_by_group.pdf"), width = 5, height = 3.5)
    }
    
    ## Facet by caner subtypes
    g <- autoplot(pca, data = s_info_s, label = F, colour = 'cancer_subtype_abbr', size = 0.5)
    g <- g + facet_wrap(~project_id) + labs(title = title_names[k])
    g <- g + theme_classic() + theme(legend.position = "right")
    ggsave(paste0(name_prefix, "_PCA_cancer_subtype.pdf"), width = 8, height = 4.5)
  }

  ## run umap with 500 PCs 
  {
  ## sum(s_pca$importance[2, 1:500])  checking the PCs to be used for UMAP
  #rownames(pca_10) <- s_info_s$sequencing_id
  npc <- min(500, ncol(pca$x))     ## for the cases that less than 500 PCs
      
  sum(s_pca$importance[2, 1:npc])
  pca_umap <- umap(pca$x[, 1:npc])        

  UMAP_1 <- pca_umap$layout[, 1]
  UMAP_2 <- pca_umap$layout[, 2]
  umap_df <- data.frame(UMAP_1, UMAP_2, s_info_s)
  
  ## UMAP plots
  g <- ggplot(data = umap_df, aes(x = UMAP_1, y = UMAP_2, color = project_id)) + geom_point(size = 0.75)
  g <- g + labs(x = "UMAP_1", y = "UMAP_2", title = title_names[k]) 
  g <- g + scale_colour_manual(values = project_col_s) + theme_classic()
  ggsave(paste0(name_prefix, "_UMAP_project_id.pdf"), width = 5.5, height = 3.5 )
  
  g <- ggplot(data = umap_df, aes(x = UMAP_1, y = UMAP_2, color = cancer_type)) + geom_point(size = 0.75)
  g <- g + labs(x = "UMAP_1", y = "UMAP_2", title = title_names[k]) 
  g <- g + scale_colour_manual(values = cancer_type_col_s) + theme_classic()
  ggsave(paste0(name_prefix, "_UMAP_cancer_type.pdf"), width = 6, height = 4)
  
  if(i == 1){
  g <- ggplot(data = umap_df, aes(x = UMAP_1, y = UMAP_2, color = group)) + geom_point(size = 2)
  g <- g + labs(x = "UMAP_1", y = "UMAP_2", title = title_names[k]) 
  g <- g + scale_colour_manual(values = AML_col) + theme_classic()
  ggsave(paste0(name_prefix, "_UMAP_group.pdf"), width = 5, height = 3.5)
  }
  
  ## Facet by caner subtypes
  g <- ggplot(data = umap_df, aes(x = UMAP_1, y = UMAP_2, color = cancer_subtype_abbr)) + geom_point(size = 0.5)
  g <- g + facet_wrap(~project_id)
  g <- g + labs(x = "UMAP_1", y = "UMAP_2", title = title_names[k]) 
  g <- g + theme_classic() + theme(legend.position = "right")
  ggsave(paste0(name_prefix, "_UMAP_cancer_subtype.pdf"), width = 8, height = 4.5)
  
}
  
  }  
    
  ## correlation_with_QC metrics
  if(i > 1){                ## AML will lead to corr error
    ############################################
    ## draw the proportion of variance top 9 PCs
    #pdf(paste0(name_prefix, "_top9PCs.pdf"), width = 8, height = 4)
    #barplot(s_pca$importance[2, 1:9])
    #dev.off()
    
    PC_var  <- round(s_pca$importance[2, 1:9], 2)
    PC_name <- names(s_pca$importance[2, 1:9])
    dat <- data.frame(PC_var, PC_name)
    g <- ggplot(data = dat, aes(x = PC_name, y = PC_var, fill = PC_name)) 
    g <- g + geom_bar(stat="identity", position=position_dodge())
    #g <- g + geom_text(aes(label = PC_var), position = position_dodge(0.9), size = 1.9, vjust = 0)
    g <- g + scale_y_continuous(position = "left") 
    g <- g + labs(y = "FVE", x  = "")     ## FVE : Fraction of variance explained
    g <- g + scale_fill_manual(values = rev(c('#fcfbfd','#efedf5','#dadaeb','#bcbddc','#9e9ac8','#807dba','#6a51a3','#54278f','#3f007d')))
    g <- g + theme_classic() + theme(legend.position = "none", axis.text.x = element_blank())
    #g <- g + theme_void() + theme(legend.position = "none")
    ggsave(paste0(name_prefix, "_top_9PCs.pdf"), width = 4, height = 1)
    
    #######################################
    ## top 
    #qc_info_s <- s_info_s %>% 
    #  select(raw_reads_depth, usable_reads_depth, saturation_maxEstCor, 
    #         specificity_pctReadsWCpG, coverage_pctCpGwReads,  ## for paired-end only
    #         enrichment_relH, enrichment_GoGe)
    
    qc_info_s <- s_info_s %>% 
      select(raw_reads_depth, usable_reads_depth, saturation_maxEstCor, 
             enrichment_relH, enrichment_GoGe)
    
    scc_pq <- cor(pca$x[, 1:9], qc_info_s, use='complete.obs', method = "spearman")
    pcc_pq <- cor(pca$x[, 1:9], qc_info_s, use='complete.obs', method = "pearson")
    
    ## headmaps
    col_pq <- colorRampPalette(c("blue", "white", "red"))(22)
    
    ## SCC 
    pdf(paste0(name_prefix, '_cor_top_9PCs_QC_SCC.pdf'), width = 8, height = 4)
    heatmap.2(t(scc_pq) , scale = "none", col =  col_pq, trace = "none", Colv = F, Rowv = F, 
              cexRow = 1.2, cexCol =  1.2, density.info = "none", breaks = seq(-1.1, 1.1, 0.1),
              lhei = c(1.5, 3),  lwid = c(1.5, 4), margins=c(4, 12), 
              key.title = "SCC", key.xlab = "", keysize = 0.75) 
    dev.off()
    
    ## PCC
    pdf(paste0(name_prefix, '_cor_top_9PCs_QC_PCC.pdf'), width = 8, height = 4)
    heatmap.2(t(pcc_pq) , scale = "none", col =  col_pq, trace = "none", Colv = F, Rowv = F, 
              cexRow = 1.2, cexCol = 1.2, density.info = "none",  breaks = seq(-1.1, 1.1, 0.1),
              lhei = c(1.5, 3),  lwid = c(1.5, 4), margins=c(4, 12), 
              key.title = "PCC", key.xlab = "", keysize = 0.5) 
    dev.off()
    
  }
  
}

}

}


