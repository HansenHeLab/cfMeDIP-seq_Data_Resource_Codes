rm(list = ls())
setwd("/Users/yong/OneDrive - UHN/Projects/TCGE/cfEpigenomics/Resource/2_batch_norm/DMRs/SE_vs_PE")

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(org.Hs.eg.db)
library(annotate)
library(ggplot2)

## bin bed infomation 
bin_bed <- read.table("/Users/yong/OneDrive - UHN/Projects/TCGE/cfEpigenomics/Resource/2_batch_norm/merge/ALL_auto_bfilt_wCpG_bins.bed")


s_info <- read.csv("/Users/yong/OneDrive - UHN/Projects/TCGE/cfEpigenomics/Resource/1_QC/1_QC_cutoffs_fixed/TCGE_cfDM_samples_only_fixed_after_QC.csv")


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
    pan_hyper_se <- read.table("../SE/Union_Hyper_overlap_with_mergedCancer_vs_normal_DMRs_binID.txt")$V1
    pan_hypo_se  <- read.table("../SE/Union_Hypo_overlap_with_mergedCancer_vs_normal_DMRs_binID.txt")$V1
    length(pan_hyper_se)
    length(pan_hypo_se)
    
    
    pan_hyper_pe <- read.table("../PE/Union_Hyper_overlap_with_mergedCancer_vs_normal_DMRs_binID.txt")$V1
    pan_hypo_pe  <- read.table("../PE/Union_Hypo_overlap_with_mergedCancer_vs_normal_DMRs_binID.txt")$V1
    length(pan_hyper_pe)
    length(pan_hypo_pe)
    
    
    length(DMR_bins_pe$`Cancer vs Normal_Hyper`)
    length(DMR_bins_pe$`Cancer vs Normal_Hypo`)
    
    ## overlaper panCancer DMRs
    pan_hyper <- intersect(pan_hyper_se, pan_hyper_pe)
    pan_hypo <- intersect(pan_hypo_se, pan_hypo_pe)
    length(pan_hyper)
    length(pan_hypo)
    
    ## for venn plot 
    length(pan_hyper_pe) -  length(pan_hyper)
    length(pan_hyper_se) -  length(pan_hyper)
    
    length(pan_hypo_pe) -  length(pan_hypo)
    length(pan_hypo_se) -  length(pan_hypo)
    
    ## panCancer DMRs IDs
    write.table(pan_hyper, "panCancer_Hyper_DMRs_binID.txt", row.names = F, col.names = F, quote = F)
    write.table(pan_hypo, "panCancer_Hypo_DMRs_binID.txt", row.names = F, col.names = F, quote = F)
    }
    
    #################################################
    ## panCancer DMRs bed and annotation for panHyper
   {
     ## panCancer DMRs bed
     idx_hyper <- match(pan_hyper, bin_bed$V4)
     idx_hypo <- match(pan_hypo, bin_bed$V4)
     
     write.table(bin_bed[idx_hyper, ], file = "panCancer_Hyper_DMRs.bed", 
                 quote = F, col.names = F, row.names = F, sep = "\t")
     
     write.table(bin_bed[idx_hypo, ], file = "panCancer_Hypo_DMRs.bed", 
                 quote = F, col.names = F, row.names = F, sep = "\t")
     
     ## panCancer hyper peaks annotation
       {
       txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
       promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
      
       {
         files <- c("panCancer_Hyper_DMRs.bed")
         files <- as.list(files)
         names(files) <- "pan-Cancer Hyper DMRs"
         ## distribution around TSS
         tagMatrixList <- lapply(files, getTagMatrix, windows=promoter)
         
         pdf("panCancer_Hyper_DMRs_around_TSS.pdf", width = 6, height = 4)
         plotAvgProf(tagMatrixList, xlim=c(-3000, 3000)) + scale_color_manual(values = "#d6604d") + theme_classic() 
         dev.off()
         
         ## peak annotation
         dmrAnnoList <- lapply(files, annotatePeak, TxDb=txdb,
                               tssRegion=c(-3000, 3000), verbose=FALSE)
         
         Anno_col <- c('#d6604d') 
         
         pdf("panCancer_Hyper_DMRs_distribution.pdf", width = 6, height = 4)
         plotAnnoBar(dmrAnnoList, title = "")  + theme_classic() 
         dev.off()
         
         ###############################################
         ## genes functional profiles for specific DMRs
         dmrs_genes = lapply(dmrAnnoList, function(i) unique(as.data.frame(i)$geneId))
         
         compKEGG <- compareCluster(geneCluster   = dmrs_genes,
                                    fun           = "enrichKEGG",
                                    pvalueCutoff  = 0.05,
                                    pAdjustMethod = "BH")
         
         pdf("panCancer_Hyper_DMRs_Genes_KEGG.pdf", width = 5, height = 6)
         dotplot(compKEGG, showCategory = 2, title = "KEGG pathway enrichmented for pan-cancer DMRs") + labs(x = "") + scale_x_discrete(guide = guide_axis(angle = 45)) 
         dev.off()
         
         compGO <- compareCluster(geneCluster   = dmrs_genes,
                                  fun           = "enrichGO",
                                  OrgDb='org.Hs.eg.db',
                                  pvalueCutoff  = 0.05,
                                  pAdjustMethod = "BH")
         
         pdf("panCancer_Hyper_DMRs_Genes_GO.pdf", width = 5, height = 6)
         dotplot(compGO, showCategory = 1, title = "GO enrichment analysis for cancer specific DMRs associated genes") + labs(x = "") + scale_x_discrete(guide = guide_axis(angle = 45)) 
         dev.off()
         
       }
       
     }
     
     
   }
  
}




####################################################################
## cancer type specific DMRs: cancer type vs the rest 
## Check overlapping of DMRs for lung and blood cancer  in PE can SE
####################################################################
{
  
  ###################################################
  ## hyper in all 
  if(TRUE){
  #DMR_bins <- readRDS("FC_1.5_DESeq_DMRs_bin_IDs.RDS"
  
    
  ## using FC_2 specific results 
  DMR_bins_pe <- readRDS("../PE/PE_Specific/FC_2_DESeq_cancer_specific_DMRs_bin_IDs.RDS")
  DMR_bins_se <- readRDS("../SE/SE_Specific/FC_2_DESeq_cancer_specific_DMRs_bin_IDs.RDS")  
    
  names(DMR_bins_pe)
  names(DMR_bins_se)
  
  #################################################
  ## check lung and blood overlapping in SE and PE
  ## using intersect specific DMRs
  length(DMR_bins_pe$cancer_typeLung.Cancer_specific_Hyper)
  length(DMR_bins_se$cancer_typeLung.Cancer_specific_Hyper)
  #lung_hyper <- unique(c(DMR_bins_pe$cancer_typeLung.Cancer_specific_Hyper, DMR_bins_se$cancer_typeLung.Cancer_specific_Hyper)) 
  lung_hyper <- intersect(DMR_bins_pe$cancer_typeLung.Cancer_specific_Hyper, DMR_bins_se$cancer_typeLung.Cancer_specific_Hyper)
  length(lung_hyper)
  
 
  length(DMR_bins_pe$cancer_typeLung.Cancer_specific_Hypo)
  length(DMR_bins_se$cancer_typeLung.Cancer_specific_Hypo)
  #lung_hypo <- unique(c(DMR_bins_pe$cancer_typeLung.Cancer_specific_Hypo, DMR_bins_se$cancer_typeLung.Cancer_specific_Hypo)) 
  lung_hypo <- intersect(DMR_bins_pe$cancer_typeLung.Cancer_specific_Hypo, DMR_bins_se$cancer_typeLung.Cancer_specific_Hypo)
  length(lung_hypo)
  
  ## blood cancer 
  length(DMR_bins_pe$cancer_typeBlood.Cancer_specific_Hyper)
  length(DMR_bins_se$cancer_typeBlood.Cancer_specific_Hyper)
  # blood_hyper <- unique(c(DMR_bins_pe$cancer_typeBlood.Cancer_specific_Hyper, DMR_bins_se$cancer_typeBlood.Cancer_specific_Hyper)) 
  blood_hyper <- intersect(DMR_bins_pe$cancer_typeBlood.Cancer_specific_Hyper, DMR_bins_se$cancer_typeBlood.Cancer_specific_Hyper)
  length(blood_hyper)
  
  length(DMR_bins_pe$cancer_typeLung.Cancer_specific_Hypo)
  length(DMR_bins_se$cancer_typeLung.Cancer_specific_Hypo)
  #blood_hypo <- unique(c(DMR_bins_pe$cancer_typeLung.Cancer_specific_Hypo, DMR_bins_se$cancer_typeLung.Cancer_specific_Hypo)) 
  blood_hypo <- intersect(DMR_bins_pe$cancer_typeLung.Cancer_specific_Hypo, DMR_bins_se$cancer_typeLung.Cancer_specific_Hypo)
  length(blood_hypo)
  
  
  
  ###################################
  ## merging pe and se specific DMRs
  names(DMR_bins_pe)
  names(DMR_bins_se)
  
  DMR_bins <- c(DMR_bins_pe[-c(1:2)], DMR_bins_se[-c(1:6)])
  names(DMR_bins)
  
  DMR_bins[[3]] <- lung_hyper
  DMR_bins[[4]] <- lung_hypo
  DMR_bins[[7]] <- blood_hyper
  DMR_bins[[8]] <- blood_hypo
  
  L <- length(DMR_bins)
  hyper_dmr <- hypo_dmr <- vector()
  
  for (i in seq(1, L, 2))
  {
    hyper_dmr <- c(hyper_dmr, DMR_bins[[i]])
    hypo_dmr <- c(hypo_dmr, DMR_bins[[i + 1]])
  }
  
  hyper_dmr_cnt <- table(hyper_dmr)
  hypo_dmr_cnt <- table(hypo_dmr)
  
  ########
  ## hyper
  pdf(paste0("FC_", fc, "_Cancer_Specific_DMRs_Hyper_frequency.pdf"), width = 5, height = 3)
  barplot(table(hyper_dmr_cnt), col = "#d6604d", ylim = c(0, (max(table(hyper_dmr_cnt)) + 0.1 * max(table(hyper_dmr_cnt))) ),
          ylab = "Number of Specific Hyper DMRs", xlab = "Frequency of DMRs identified from comparison of cancer type vs the rest")
  dev.off()
  
  #######
  ## hypo
  pdf(paste0("FC_", fc, "_Cancer_Specific_DMRs_Hypo_frequency.pdf"), width = 5, height = 3)
  barplot(table(hypo_dmr_cnt), col = "#4575b4", ylim = c(0, (max(table(hypo_dmr_cnt)) + 0.1 * max(table(hypo_dmr_cnt))) ),
          ylab = "Number of Hypo DMRs", xlab = "Frequency of DMRs identified from comparison of each cancer type vs the rest")
  dev.off()
 
  }

  
  ######################################################
  ##### cancer-specific hyper and hypo  peaks annotation
  #### remove specific DMRs with frequcncy > 1
  if(TRUE){
    hyper_u <- names(hyper_dmr_cnt)[hyper_dmr_cnt  == 1]
    hypo_u <- names(hypo_dmr_cnt)[hypo_dmr_cnt  == 1]
    
    L <- length(DMR_bins)
    hyper_spec_cnt <- hypo_spec_cnt <- vector()
    k = 1
    
    for (i in seq(1, L, 2))
    {
      
      ## cancer type
      ct <-  strsplit(names(DMR_bins)[i], "_")[[1]][2]
      ct <- gsub("type", "", ct)
      ct <- gsub("\\.", " ", ct)
      
      ## matched hyper dmr bins
      idx_bin <- match(DMR_bins[[i]], hyper_u)    
      hyper_bin <- DMR_bins[[i]][!is.na(idx_bin)]
      hyper_spec_cnt[k] <- length(hyper_bin) 
        
      idx_bed <- match(hyper_bin, bin_bed$V4)
      write.table(bin_bed[idx_bed, ], file = paste0(ct, " specific Hyper DMRs.bed"), 
                  quote = F, col.names = F, row.names = F, sep = "\t")
      

      ## ## matched hypo dmr bins
      idx_bin <- match(DMR_bins[[i+1]], hypo_u)    
      hypo_bin <- DMR_bins[[i+1]][!is.na(idx_bin)]
      hypo_spec_cnt[k] <- length(hypo_bin) 
      
      idx_bed <- match(hypo_bin, bin_bed$V4)
      write.table(bin_bed[idx_bed, ], file = paste0(ct, " specific Hypo DMRs.bed"), 
                  quote = F, col.names = F, row.names = F, sep = "\t")
      
      names(hyper_spec_cnt)[k] <-  names(hypo_spec_cnt)[k] <- ct
      k <- k + 1
    }
    
    ## updating blood and eye cancer 
    names(hyper_spec_cnt)[4:5] <- c("AML", "Uveal Melanoma")
    names(hypo_spec_cnt)[4:5] <- c("AML", "Uveal Melanoma")
    
    ###########################
    ## number of specific DMRs
    {
      cancer_type <- names(hyper_spec_cnt)
      
      hyper_spec_cnt_plot <- hyper_spec_cnt
      hyper_spec_cnt_plot[ hyper_spec_cnt_plot > 10000] <- 10000
   
      dat_hyper <- data.frame(cancer_type, hyper_spec_cnt, hyper_spec_cnt_plot)
      write.csv(dat_hyper, "Cancer_type_specific_Hyper_DMRs_after_filtering.csv")
      
      dat_hyper$cancer_type <- factor(dat_hyper$cancer_type,
                                levels = names(sort(hyper_spec_cnt)))
      
      ## matching colors for project and cancer type
      idx_ct <- match(levels(dat_hyper$cancer_type), levels(s_info$cancer_type)) 
      cancer_type_col_s <- cancer_type_col[idx_ct]
      
      ## hyper specific
      g <- ggplot(data = dat_hyper, aes(x = cancer_type, y = hyper_spec_cnt_plot, fill = cancer_type)) 
      g <- g + geom_bar(stat="identity", position=position_dodge())
      g <- g + geom_text(aes(label = hyper_spec_cnt), position = position_dodge(0.9), size=3, hjust = 0.75)
      g <- g + labs(y = "Number of specific Hyper DMRs", x  = "")
      g <- g + scale_fill_manual(values = cancer_type_col_s)
      g <- g + coord_flip() + theme_classic()
      g <- g + theme_classic() + theme(legend.position = "none")
      ggsave("Hyper_cancer_type_specific_DMR_count.pdf", width = 5, height = 4)
      
  
      ## hypo specific
      cancer_type <- names(hypo_spec_cnt)
      
      hypo_spec_cnt_plot <- hypo_spec_cnt
      hypo_spec_cnt_plot[hypo_spec_cnt_plot > 8000] <- 8000
      
      dat_hypo <- data.frame(cancer_type, hypo_spec_cnt, hypo_spec_cnt_plot)
      
      write.csv(dat_hypo, "Cancer_type_specific_Hypo_DMRs_after_filtering.csv")
      
      ## remove 0 cancer types
      idx_rm <- dat_hypo$hypo_spec_cnt == 0
      dat_hypo <- dat_hypo[!idx_rm, ]
      
      dat_hypo$cancer_type <- factor(dat_hypo$cancer_type,
                                levels = names(sort(hypo_spec_cnt)))
      
      g <- ggplot(data = dat_hypo, aes(x = cancer_type, y = hypo_spec_cnt_plot, fill = cancer_type)) 
      g <- g + geom_bar(stat="identity", position=position_dodge())
      g <- g + geom_text(aes(label = hypo_spec_cnt), position = position_dodge(0.9), size=3, hjust = 0.75)
      g <- g + labs(y = "Number of specific Hypo DMRs", x  = "")
      g <- g + scale_fill_manual(values = cancer_type_col_s)
      g <- g + coord_flip() + theme_classic()
      g <- g + theme_classic() + theme(legend.position = "none")
      ggsave("Hypo_cancer_type_specific_DMRs_count.pdf", width = 5, height = 4)
      
    }
    
  }

  
  
  ######################################
  ##   peaks annotations
  ######################################
  if(TRUE){
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
    
    ## for hyper
    peak_type <- "Hyper"
    {
    files <- paste(names(sort(hyper_spec_cnt, decreasing = T)), "specific", peak_type, "DMRs.bed")
    files <- as.list(files)
    names(files) <- names(sort(hyper_spec_cnt, decreasing = T))
    
    ## distribution around TSS
    tagMatrixList <- lapply(files, getTagMatrix, windows=promoter)
    
    idx_col <- match(names(tagMatrixList), levels(s_info$cancer_type))
    cancer_type_col_s <- cancer_type_col[idx_col]
    
    pdf(paste0(peak_type, "_cancer_type_specific_DMRs_around_TSS.pdf"), width = 6, height = 2)
    plotAvgProf(tagMatrixList, xlim=c(-3000, 3000)) + scale_color_manual(values = cancer_type_col_s) + theme_classic() 
    dev.off()
    
    ## peak annotation
    dmrAnnoList <- lapply(files, annotatePeak, TxDb=txdb,
                           tssRegion=c(-3000, 3000), verbose=FALSE)
    
    Anno_col <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#b15928')
    
    pdf(paste0(peak_type, "_cancer_type_specific_DMRs_annotation.pdf"), width = 6, height = 4)
    plotAnnoBar(dmrAnnoList, title = "") +  scale_fill_manual(values = Anno_col)  + theme_classic() 
    dev.off()
    
    ###############################################
    ## genes functional profiles for specific DMRs
    dmrs_genes = lapply(dmrAnnoList, function(i) unique(as.data.frame(i)$geneId))
    
    compKEGG <- compareCluster(geneCluster   = dmrs_genes,
                               fun           = "enrichKEGG",
                               pvalueCutoff  = 0.05,
                               pAdjustMethod = "BH")
  
    pdf(paste0(peak_type, "_cancer_type_specific_DMRs_annotation_genes_KEGG.pdf"), width = 10, height = 6)
    dotplot(compKEGG, showCategory = 2, title = "KEGG pathway enrichment analysis for cancer specific DMRs associated genes") + labs(x = "") + scale_x_discrete(guide = guide_axis(angle = 45)) 
    dev.off()
    
    compGO <- compareCluster(geneCluster   = dmrs_genes,
                             fun           = "enrichGO",
                             OrgDb='org.Hs.eg.db',
                             pvalueCutoff  = 0.05,
                             pAdjustMethod = "BH")
    
    pdf(paste0(peak_type, "_cancer_type_specific_DMRs_annotation_genes_GO.pdf"), width = 10, height = 6)
    dotplot(compGO, showCategory = 1, title = "GO enrichment analysis for cancer specific DMRs associated genes") + labs(x = "") + scale_x_discrete(guide = guide_axis(angle = 45)) 
    dev.off()
    
    }
    
    ## for hyper
    peak_type <- "Hypo"
    {
      
      ## Lung Cancer and Blood Cancer  hype DMRS == 0
      files <- paste(names(sort(hypo_spec_cnt, decreasing = T)[1:9]), "specific", peak_type, "DMRs.bed")
      files <- as.list(files)
      names(files) <- names(sort(hypo_spec_cnt, decreasing = T)[1:9])
      
      ## distribution around TSS
      tagMatrixList <- lapply(files, getTagMatrix, windows=promoter)
      
      idx_col <- match(names(tagMatrixList), levels(s_info$cancer_type))
      cancer_type_col_s <- cancer_type_col[idx_col]
      
      pdf(paste0(peak_type, "_cancer_type_specific_DMRs_count_around_TSS.pdf"), width = 6, height = 2)
      plotAvgProf(tagMatrixList, xlim=c(-3000, 3000)) + scale_color_manual(values = cancer_type_col_s) + theme_classic() 
      dev.off()
      
      ## peak annotation
      dmrAnnoList <- lapply(files, annotatePeak, TxDb=txdb,
                            tssRegion=c(-3000, 3000), verbose=FALSE)
      
      Anno_col <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#b15928')
      
      pdf(paste0(peak_type, "_cancer_type_specific_DMRs_annotation.pdf"), width = 6, height = 4)
      plotAnnoBar(dmrAnnoList, title = "") +  scale_fill_manual(values = Anno_col)  + theme_classic() 
      dev.off()
      
      ###############################################
      ## genes functional profiles for specific DMRs
      dmrs_genes = lapply(dmrAnnoList, function(i) unique(as.data.frame(i)$geneId))
      
      compKEGG <- compareCluster(geneCluster   = dmrs_genes,
                                 fun           = "enrichKEGG",
                                 pvalueCutoff  = 0.05,
                                 pAdjustMethod = "BH")
      
      pdf(paste0(peak_type, "_cancer_type_specific_DMRs_annotation_genes_KEGG.pdf"), width = 10, height = 6)
      dotplot(compKEGG, showCategory = 2, title = "KEGG pathway enrichment analysis for cancer specific DMRs associated genes") + labs(x = "") + scale_x_discrete(guide = guide_axis(angle = 45)) 
      dev.off()
      
      compGO <- compareCluster(geneCluster   = dmrs_genes,
                               fun           = "enrichGO",
                               OrgDb='org.Hs.eg.db',
                               pvalueCutoff  = 0.05,
                               pAdjustMethod = "BH")
      
      pdf(paste0(peak_type, "_cancer_type_specific_DMRs_annotation_genes_GO.pdf"), width = 10, height = 6)
      dotplot(compGO, showCategory = 2, title = "GO enrichment analysis for cancer specific DMRs associated genes") + labs(x = "") + scale_x_discrete(guide = guide_axis(angle = 45)) 
      dev.off()
      
    }
    
    
    
  }
  

  
  
  
}


