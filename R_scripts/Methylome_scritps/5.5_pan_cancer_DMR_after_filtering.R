rm(list = ls())
setwd("/Users/yong/OneDrive - UHN/Projects/TCGE/cfEpigenomics/Resource/2_batch_norm/DMRs_revision/SE_vs_PE_after_filtering")

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(org.Hs.eg.db)
library(annotate)
library(ggplot2)

## bin bed infomation 
bin_bed <- read.table("/Users/yong/OneDrive - UHN/Projects/TCGE/cfEpigenomics/Resource/2_batch_norm/merge/ALL_auto_bfilt_wCpG_bins.bed")

## 
s_info <-  read.csv("/Users/yong/Library/CloudStorage/OneDrive-UHN/Projects/TCGE/cfEpigenomics/Resource/1_QC/1_QC_cutoffs_for_revision/TCGE_cfDM_samples_only_fixed_after_QC_vialA_only.csv")


s_info$cancer_type <- factor(s_info$cancer_type, 
                             levels = rev(c("Healthy", "Brain Cancer", "Lung Cancer", "Prostate Cancer",
                                            "AML", "Pancreatic Cancer", "Uveal Melanoma",
                                            "Head and Neck Cancer", "Breast Cancer", "Colorectal Cancer",
                                            "Bladder Cancer", "Renal Cancer",
                                            "LFS Survivor", "LFS Previvor", "LFS Positive")))

cancer_type_col <- rev(c('#33a02c','#1f78b4','#b2df8a','#a6cee3','#fb9a99',
                         '#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a',
                         '#fb6a4a','#b15928','#bdbdbd','#969696','#737373'))

fc = 2

#################################################################
## identify the pan-cancer DMRs overlapped in SE and PE
#################################################################
{
    ############################
    ## overlapped panCancer DMRs
    {
    pan_hyper_se <- read.table("../SE_after_filtering/Union_Hyper_overlap_with_mergedCancer_vs_healthy_DMRs_binID.txt")$V1
    pan_hypo_se  <- read.table("../SE_after_filtering/Union_Hypo_overlap_with_mergedCancer_vs_healthy_DMRs_binID.txt")$V1
    length(pan_hyper_se)
    length(pan_hypo_se)
    
    
    pan_hyper_pe <- read.table("../PE_after_filtering/Union_Hyper_overlap_with_mergedCancer_vs_healthy_DMRs_binID.txt")$V1
    pan_hypo_pe  <- read.table("../PE_after_filtering/Union_Hypo_overlap_with_mergedCancer_vs_healthy_DMRs_binID.txt")$V1
    length(pan_hyper_pe)
    length(pan_hypo_pe)
    
    
    ## overlaper panCancer DMRs
    pan_hyper <- intersect(pan_hyper_se, pan_hyper_pe)
    pan_hypo <- intersect(pan_hypo_se, pan_hypo_pe)
    length(pan_hyper)
    length(pan_hypo)
    
    ## for venn plot: colors    #dedede,#dedede,#d6604d
    length(pan_hyper_pe) -  length(pan_hyper)
    length(pan_hyper_se) -  length(pan_hyper)
    
    ## for venn plot: colors    #dedede,#dedede,#4575b4
    length(pan_hypo_pe) -  length(pan_hypo)
    length(pan_hypo_se) -  length(pan_hypo)
    
    ## panCancer DMRs IDs
    write.table(pan_hyper, "panCancer_Hyper_DMRs_binID.txt", row.names = F, col.names = F, quote = F)
    write.table(pan_hypo, "panCancer_Hypo_DMRs_binID.txt", row.names = F, col.names = F, quote = F)
    }
    
  ##########################################
  ## overlapped DMRs vs non-overlapped  DMRs
  ## using the combined cancer vs healthy directly, since most of them are overlapped with individual cancer type comparison
  ## checking the hyper only  
  {
    ## before filtering   
    dmr_csv_pe <- read.csv("/Users/yong/Library/CloudStorage/OneDrive-UHN/Projects/TCGE/cfEpigenomics/Resource/2_batch_norm/DMRs_revision/PE/FC_2_Cancer vs Healthy_DESeq_Hyper_DMRs.csv")
    dmr_csv_se <- read.csv("/Users/yong/Library/CloudStorage/OneDrive-UHN/Projects/TCGE/cfEpigenomics/Resource/2_batch_norm/DMRs_revision/SE/FC_2_Cancer vs Healthy_DESeq_Hyper_DMRs.csv")
    
    ## after filtering 
    idx_pe_f <- match(pan_hyper_pe, dmr_csv_pe$X)
    idx_se_f <- match(pan_hyper_se, dmr_csv_se$X)
    
    dmr_csv_pe  <- dmr_csv_pe[idx_pe_f, ]
    dmr_csv_se <- dmr_csv_se[idx_se_f, ]
    dmr_cmb <- rbind(dmr_csv_pe, dmr_csv_se)
    
    dmr_group <- c(rep("PE", nrow(dmr_csv_pe)), rep("SE", nrow(dmr_csv_se)))
    
    dmr_type_pe <- rep("Unique", nrow(dmr_csv_pe))
    idx_pe <- match(pan_hyper, dmr_csv_pe$X)
    dmr_type_pe[idx_pe] <- "Common"
    
    dmr_type_se <- rep("Unique", nrow(dmr_csv_se))
    idx_se <- match(pan_hyper, dmr_csv_se$X)
    dmr_type_se[idx_se] <- "Common"
    
    dmr_type <- c(dmr_type_pe, dmr_type_se)
    
    data <- cbind(dmr_cmb, dmr_type, dmr_group)
    data <- data %>%  mutate(log2baseMean = log2(baseMean), log10p_val = -log10(pvalue))
    
    saveRDS(data, "PE_SE_pan_cancer_hyper_DMRs_combined.RDS")
    
    ## testing 
    idx_pe_ol <- dmr_group == "PE" & dmr_type == "Common"
    idx_pe_sp <- dmr_group == "PE" & dmr_type == "Unique"
    
    idx_se_ol <- dmr_group == "SE" & dmr_type == "Common"
    idx_se_sp <- dmr_group == "SE" & dmr_type == "Unique"
    
    ## boxplot for the baseMean
    wilcox.test(data$log2baseMean[idx_pe_ol], data$log2baseMean[idx_pe_sp], alternative = "greater")
    wilcox.test(data$log2baseMean[idx_se_ol], data$log2baseMean[idx_se_sp], alternative = "greater")

    g <- ggplot(data, aes(x = dmr_group, y = log2baseMean, fill = dmr_type))
    g <- g + geom_boxplot() + labs(y = "log2baseMean", x = "") 
    g <- g + scale_fill_manual(values = c("#d6604d", "#dedede"))
    g <- g + theme_classic() + theme(legend.position = "top")
    ggsave("PE_SE_pan_cancer_hyper_DMRs_log2baseMean.pdf", height = 2.5, width = 3)
    
  
    ## boxplot for the log2FC
    wilcox.test(data$log2FoldChange[idx_pe_ol], data$log2FoldChange[idx_pe_sp], alternative = "greater")
    wilcox.test(data$log2FoldChange[idx_se_ol], data$log2FoldChange[idx_se_sp], alternative = "greater")
    
    g <- ggplot(data, aes(x = dmr_group, y = log2FoldChange, fill = dmr_type))
    g <- g + geom_boxplot()  + labs(y = "log2FoldChange", x = "") 
    g <- g + scale_fill_manual(values = c("#d6604d", "#dedede"))
    g <- g + theme_classic() + theme(legend.position = "top")
    ggsave("PE_SE_pan_cancer_hyper_DMRs_log2FoldChange.pdf",  height = 2.5, width = 3)
    
    ## boxplot for the -log10pval
    wilcox.test(data$log10p_val[idx_pe_ol], data$log10p_val[idx_pe_sp], alternative = "greater")
    wilcox.test(data$log10p_val[idx_se_ol], data$log10p_val[idx_se_sp], alternative = "greater")
    
    g <- ggplot(data, aes(x = dmr_group, y = log10p_val, fill = dmr_type))
    g <- g + geom_boxplot()  + labs(y = "-log10Pvalue", x = "") 
    g <- g + scale_fill_manual(values = c("#d6604d", "#dedede"))
    g <- g + theme_classic() + theme(legend.position = "top")
    ggsave("PE_SE_pan_cancer_hyper_DMRs_log10Pvalue.pdf", height = 2.5, width = 3)
    
    
  }
  
  

    #################################################
    ## panCancer DMRs bed and annotation for panHyper
   {
     ## panCancer DMRs bed
     idx_hyper <- match(pan_hyper, bin_bed$V4)
     idx_hypo <- match(pan_hypo, bin_bed$V4)
     
     hyper_bed <- bin_bed[idx_hyper, ]
     hypo_bed <- bin_bed[idx_hypo, ]
     
     write.table(hyper_bed, file = "panCancer_Hyper_DMRs.bed", 
                 quote = F, col.names = F, row.names = F, sep = "\t")
     
     write.table(hypo_bed , file = "panCancer_Hypo_DMRs.bed", 
                 quote = F, col.names = F, row.names = F, sep = "\t")
     
     ## panCancer hyper peaks annotation
       {
       txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
       #promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
       promoter <- getPromoters(TxDb=txdb, upstream=1500, downstream=1500)         ## to be consistent with missMethyl
         
       {
         files <- c("panCancer_Hyper_DMRs.bed")
         files <- as.list(files)
         names(files) <- "pan-Cancer Hyper DMRs"
         ## distribution around TSS
         tagMatrixList <- lapply(files, getTagMatrix, windows=promoter)
         
         pdf("panCancer_Hyper_DMRs_around_TSS.pdf", width = 5, height = 3.5)
         plotAvgProf(tagMatrixList, xlim=c(-1500, 1500)) + scale_color_manual(values = "#d6604d") + theme_classic() + theme(legend.position = "none") 
         dev.off()
         
         ## peak annotation
         dmrAnnoList <- lapply(files, annotatePeak, TxDb=txdb,
                               tssRegion=c(-1500, 1500), verbose=FALSE)
         
         Anno_col <- c('#d6604d') 
         
         pdf("panCancer_Hyper_DMRs_distribution.pdf", width = 6, height = 1.5 )
         plotAnnoBar(dmrAnnoList, title = "")  + theme_classic() + theme(legend.position = "none")
         dev.off()
         
         ## for legend only 
         pdf("panCancer_Hyper_DMRs_distribution_legend.pdf", width = 5, height = 5 )
         plotAnnoBar(dmrAnnoList, title = "")  + theme_classic() + theme(legend.position = "right")
         dev.off()
         #########################################################
         ## genes functional profiles for for all pan-cancer DMRs
         {
        
         dmrs_genes = lapply(dmrAnnoList, function(i) unique(as.data.frame(i)$geneId))
         
         compKEGG <- compareCluster(geneCluster   = dmrs_genes,
                                    fun           = "enrichKEGG",
                                    pvalueCutoff  = 0.05,
                                    pAdjustMethod = "BH")
         
         write.csv(data.frame(compKEGG), "panCancer_Hyper_DMRs_Genes_KEGG_all.csv", row.names = F)
         
         pdf("panCancer_Hyper_DMRs_Genes_KEGG_top10.pdf", width = 5, height = 6)
         dotplot(compKEGG,  showCategory=10, title = "Pan-cancer hyper DMRs") + labs(x = "") + scale_x_discrete(guide = guide_axis(angle = 45)) 
         dev.off()
         
         #####
         ## GO
         compGO <- compareCluster(geneCluster   = dmrs_genes,
                                  fun           = "enrichGO",
                                  OrgDb='org.Hs.eg.db',
                                  pvalueCutoff  = 0.05,
                                  pAdjustMethod = "BH")
         
         compGO_sim <- clusterProfiler::simplify(compGO )
         write.csv(data.frame(compGO_sim), "panCancer_Hyper_DMRs_Genes_GO_all_simplified.csv", row.names = F)
         
         ## plot by levels
         #  ego_sim_1 <- gofilter(ego_sim, level = 1)
    
         pdf("panCancer_Hyper_DMRs_Genes_GO_top10.pdf", width = 5, height = 6)
         dotplot(compGO_sim, showCategory=10, title = "Pan-cancer hyper DMRs") + labs(x = "") + scale_x_discrete(guide = guide_axis(angle = 45)) 
         dev.off()
         
         }
         
         ###########################################################################
         ## genes functional profiles for for all pan-cancer DMRs in promoter region
         {
           
           ## limited <=1 kb region
           dmr_data = as.data.frame(dmrAnnoList)
           #idx_promoter <- dmr_data$pan.Cancer.Hyper.DMRs.annotation == "Promoter (<=1kb)" | dmr_data$pan.Cancer.Hyper.DMRs.annotation == "Promoter (1-2kb)" | dmr_data$pan.Cancer.Hyper.DMRs.annotation == "Promoter (2-3kb)" 
           #idx_promoter <- dmr_data$pan.Cancer.Hyper.DMRs.annotation == "Promoter (<=1kb)" 
           idx_promoter <- dmr_data$pan.Cancer.Hyper.DMRs.annotation == "Promoter" 
           dmr_promoter <- dmr_data[idx_promoter, ]
           
           nrow(dmr_promoter) 
           nrow(dmr_promoter) / nrow(dmr_data)
           
           #dmrs_genes <- list()
           dmrs_genes[[1]] =  unique(dmr_promoter$pan.Cancer.Hyper.DMRs.geneId)   ## need to be list, 
          
           compKEGG <- compareCluster(geneCluster   = dmrs_genes,
                                      fun           = "enrichKEGG",
                                      pvalueCutoff  = 0.05,
                                      pAdjustMethod = "BH")
           Ranking <- compKEGG_rank <- seq(1, nrow(compKEGG@compareClusterResult), 1)
           write.csv(cbind(Ranking, data.frame(compKEGG)), "panCancer_Hyper_DMRs_Promoter1.5K_Genes_KEGG_all.csv", row.names = F)
           
           pdf("panCancer_Hyper_DMRs_Promoter1.5K_Genes_KEGG_top10.pdf", width = 5, height = 6)
           dotplot(compKEGG, showCategory=10, title = "KEGG pathway enrichmented for pan-cancer DMRs") + labs(x = "") + scale_x_discrete(guide = guide_axis(angle = 45)) 
           dev.off()
           
           ### for 
           compGO <- compareCluster(geneCluster   = dmrs_genes,
                                    fun           = "enrichGO",
                                    OrgDb='org.Hs.eg.db',
                                    pvalueCutoff  = 0.05,
                                    pAdjustMethod = "BH")
           
           compGO_sim <- clusterProfiler::simplify(compGO )
           Ranking <- compGO_sim_rank <- seq(1, nrow(compGO_sim@compareClusterResult), 1)
           
           write.csv(cbind(Ranking, data.frame(compGO_sim)), "panCancer_Hyper_DMRs_Promoter1.5K_Genes_GO_all_simplified.csv", row.names = F)
           
           ## plot by levels
           #  ego_sim_1 <- gofilter(ego_sim, level = 1)
           
           pdf("panCancer_Hyper_DMRs_Promoter1.5K_Genes_GO_top10.pdf", width = 5, height = 6)
           dotplot(compGO_sim, showCategory=10, title = "Pan-cancer hyper DMRs") + labs(x = "") + scale_x_discrete(guide = guide_axis(angle = 45)) 
           dev.off()
         
       }
       
         ###########################################################################
         ## DMRs based GSEA analysis
         {
        
          library(missMethyl)
          
          ## need to be lifeover to hg19
          dmr_bed_hg19 <- read.table("panCancer_Hyper_DMRs_hg19.bed")
           
          hyper_gr_hg19 <- GRanges(seqnames = dmr_bed_hg19$V1,
                              ranges = IRanges(start = dmr_bed_hg19$V2,
                                               end = dmr_bed_hg19$V3,
                                               names = dmr_bed_hg19$V4))
          
          #######
          ## KEGG
          {
          dmrKEGG <- goregion(hyper_gr_hg19, collection="KEGG", array.type="EPIC", 
                              genomic.features = "TSS1500")
          
          idx_sig <- dmrKEGG$P.DE < 0.05
          #idx_sig <- dmrKEGG$FDR < 0.05
          sig_dmrKEGG <- dmrKEGG[idx_sig, ]
          write.csv(sig_dmrKEGG, "panCancer_Hyper_DMRs_TSS1500_Genes_KEGG_all_byMissMehthyl.csv", row.names = F)
          
          
          ## overlap with cluster profiler results
          idx_ol <- match(compKEGG@compareClusterResult$ID, rownames(sig_dmrKEGG))
          
          nrow(compKEGG)
          nrow(sig_dmrKEGG)
          sum(!is.na(idx_ol))
          
          ## draw the overlapped 
          compKEGG_plot <- compKEGG[!is.na(idx_ol)]
          compKEGG_plot  <- merge_result(list(Pan_cancer_hyper_DMRs = compKEGG_plot))  ## dataframe to clusterprofiler object 
          
          pdf("panCancer_Hyper_DMRs_Promoter1.5K_Genes_KEGG_clusterprofiler_MissMehthyl_overlapped_top10.pdf", width = 5, height = 6)
          dotplot(compKEGG_plot, showCategory=10, title = "Top enriched KEGG pathways") + labs(x = "") + scale_x_discrete(guide = guide_axis(angle = 0)) 
          dev.off()
          
          Ranking <- compKEGG_rank[!is.na(idx_ol)]
          
          
          ol_KEGG <- cbind(Ranking, data.frame(compKEGG_plot), sig_dmrKEGG[idx_ol[!is.na(idx_ol)], ])
          write.csv(ol_KEGG, "panCancer_Hyper_DMRs_TSS1500_Genes_KEGG_clusterprofiler_MissMehthyl_overlapped.csv", row.names = F)
          }
          
          ##########
          ## for GO 
          dmrGO <- goregion(hyper_gr_hg19, collection="GO", array.type="EPIC", 
                            genomic.features = "TSS1500")
          
          
          idx_sig <- dmrGO$P.DE < 0.05
          #idx_sig <- dmrGO$FDR < 0.05
          sig_dmrGO <- dmrGO[idx_sig, ]
          write.csv(sig_dmrGO, "panCancer_Hyper_DMRs_TSS1500_Genes_GO_all_byMissMehthyl.csv", row.names = F)
          
          ## overlap with cluster profiler results
          #idx_ol <- match(compGO@compareClusterResult$ID, rownames(sig_dmrGO))  
          #nrow(compGO)
          
          idx_ol <- match(compGO_sim@compareClusterResult$ID, rownames(sig_dmrGO))
          nrow(compGO_sim)
          nrow(sig_dmrGO)
          sum(!is.na(idx_ol))
        
          ## draw the overlapped 
          compGO_plot <- compGO_sim[!is.na(idx_ol)]
          compGO_plot  <- merge_result(list(Pan_cancer_hyper_DMRs = compGO_plot))  ## dataframe to clusterprofiler object 
          
          pdf("panCancer_Hyper_DMRs_Promoter1.5K_Genes_GO_clusterprofiler_MissMehthyl_overlapped_top10.pdf", width = 5, height = 6)
          dotplot(compGO_plot, showCategory=10, title = "Top enriched GO terms") + labs(x = "") + scale_x_discrete(guide = guide_axis(angle = 0)) 
          dev.off()
          
          Ranking <- compGO_sim_rank[!is.na(idx_ol)]
          
          ol_GO <- cbind(Ranking, data.frame(compGO_plot), sig_dmrGO[idx_ol[!is.na(idx_ol)], ])
          write.csv(ol_GO, "panCancer_Hyper_DMRs_TSS1500_Genes_GO_clusterprofiler_MissMehthyl_overlapped.csv", row.names = F)
         }
         
         
     }
     
     
   }
  
}

}



