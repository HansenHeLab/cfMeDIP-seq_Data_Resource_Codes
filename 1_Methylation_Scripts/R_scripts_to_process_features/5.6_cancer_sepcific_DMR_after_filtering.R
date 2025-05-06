rm(list = ls())
setwd("/Users/yong/Library/CloudStorage/OneDrive-UHN/Projects/TCGE/cfEpigenomics/Resource/2_batch_norm/DMRs_revision/SE_vs_PE_after_filtering/cancer_specific")

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(org.Hs.eg.db)
library(annotate)
library(ggplot2)

{
fc = 2  ## fold change
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
    #DMR_bins_pe <- readRDS("/Users/yong/Library/CloudStorage/OneDrive-UHN/Projects/TCGE/cfEpigenomics/Resource/2_batch_norm/DMRs_revision/DMR_filtering/PE_cancer_specific_DMRs_after_PBL_depletion_bin_IDs.RDS")
    #DMR_bins_se <- readRDS("/Users/yong/Library/CloudStorage/OneDrive-UHN/Projects/TCGE/cfEpigenomics/Resource/2_batch_norm/DMRs_revision/DMR_filtering/SE_cancer_specific_DMRs_after_PBL_depletion_bin_IDs.RDS")  
    
    DMR_bins_pe <- readRDS("/Users/yong/Library/CloudStorage/OneDrive-UHN/Projects/TCGE/cfEpigenomics/Resource/2_batch_norm/DMRs_revision/DMR_filtering/PE_cancer_specific_DMRs_after_blood_cell_age_sex_filtered_bin_IDs.RDS")
    DMR_bins_se <- readRDS("/Users/yong/Library/CloudStorage/OneDrive-UHN/Projects/TCGE/cfEpigenomics/Resource/2_batch_norm/DMRs_revision/DMR_filtering/SE_cancer_specific_DMRs_after_blood_cell_age_sex_filtered_bin_IDs.RDS")  
    
    
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
    
    ## AML
    length(DMR_bins_pe$cancer_typeAML_specific_Hyper)
    length(DMR_bins_se$cancer_typeAML_specific_Hyper)
    # blood_hyper <- unique(c(DMR_bins_pe$cancer_typeBlood.Cancer_specific_Hyper, DMR_bins_se$cancer_typeBlood.Cancer_specific_Hyper)) 
    AML_hyper <- intersect(DMR_bins_pe$cancer_typeAML_specific_Hyper, DMR_bins_se$cancer_typeAML_specific_Hyper)
    length(AML_hyper)
    
    length(DMR_bins_pe$cancer_typeAML_specific_Hypo)
    length(DMR_bins_se$cancer_typeAML_specific_Hypo)
    #blood_hypo <- unique(c(DMR_bins_pe$cancer_typeLung.Cancer_specific_Hypo, DMR_bins_se$cancer_typeLung.Cancer_specific_Hypo)) 
    AML_hypo <- intersect(DMR_bins_pe$cancer_typeAML_specific_Hypo, DMR_bins_se$cancer_typeAML_specific_Hypo)
    length(AML_hypo)
    
    
    ###################################
    ## merging pe and se specific DMRs
    names(DMR_bins_pe)
    names(DMR_bins_se)
    
    DMR_bins <- c(DMR_bins_pe[-c(1:2)], DMR_bins_se[-c(1:6)])
    names(DMR_bins)
    
    DMR_bins[[3]] <- lung_hyper
    DMR_bins[[4]] <- lung_hypo
    DMR_bins[[7]] <- AML_hyper
    DMR_bins[[8]] <- AML_hypo
    
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
    length(hyper_dmr_cnt)
    length(hypo_dmr_cnt)
    
    hyper_u <- names(hyper_dmr_cnt)[hyper_dmr_cnt  == 1]
    hypo_u <- names(hypo_dmr_cnt)[hypo_dmr_cnt  == 1]
    
    length(hyper_u) / length(hyper_dmr_cnt)
    length(hypo_u) / length(hypo_dmr_cnt)
    
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
    
    
    ###########################
    ## number of specific DMRs
    {
      cancer_type <- names(hyper_spec_cnt)
      
      hyper_spec_cnt_plot <- hyper_spec_cnt
      hyper_spec_cnt_plot[ hyper_spec_cnt_plot > 2500] <- 2500
      
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
      ggsave("Hyper_cancer_type_specific_DMR_count.pdf", width = 4, height = 4)
      
      
      ## hypo specific
      cancer_type <- names(hypo_spec_cnt)
      
      hypo_spec_cnt_plot <- hypo_spec_cnt
      hypo_spec_cnt_plot[hypo_spec_cnt_plot > 1500] <- 1500
      
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
      ggsave("Hypo_cancer_type_specific_DMRs_count.pdf", width = 4, height = 4)
      
    }
    
  }
  
  ######################################
  ##   peaks annotations
  ######################################
  if(TRUE){
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    promoter <- getPromoters(TxDb=txdb, upstream=1500, downstream=1500)
    
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
      
      pdf(paste0(peak_type, "_cancer_type_specific_DMRs_around_TSS.pdf"), width = 5, height = 2)
      plotAvgProf(tagMatrixList, xlim=c(-1500, 1500)) + scale_color_manual(values = cancer_type_col_s) + theme_classic() 
      dev.off()
      
      ## at leat 100 DMRs
      {
      tagMatrixList_100  <- lapply(files[1:8], getTagMatrix, windows=promoter)
      idx_col_100 <- match(names(tagMatrixList_100), levels(s_info$cancer_type))
      cancer_type_col_s_100 <- cancer_type_col[idx_col_100]
      
      pdf(paste0(peak_type, "_cancer_type_specific_DMRs_around_TSS_100.pdf"), width = 5, height = 2)
      plotAvgProf(tagMatrixList_100, xlim=c(-1500, 1500)) + scale_color_manual(values = cancer_type_col_s_100) + theme_classic() 
      dev.off()
      }
      
      ## peak annotation
      dmrAnnoList <- lapply(files, annotatePeak, TxDb=txdb,
                            tssRegion=c(-1500, 1500), verbose=FALSE)
      
      Anno_col <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#b15928')
      
      pdf(paste0(peak_type, "_cancer_type_specific_DMRs_annotation.pdf"), width = 5, height = 4)
      plotAnnoBar(dmrAnnoList, title = "") +  scale_fill_manual(values = Anno_col)  + theme_classic() 
      dev.off()
      
      ###############################################
      ## genes functional profiles for all specific DMRs
      {
      dmrs_genes = lapply(dmrAnnoList, function(i) unique(as.data.frame(i)$geneId))
      
      compKEGG <- compareCluster(geneCluster   = dmrs_genes,
                                 fun           = "enrichKEGG",
                                 pvalueCutoff  = 0.05,
                                 pAdjustMethod = "BH")
      
      pdf(paste0(peak_type, "_cancer_type_specific_DMRs_annotation_genes_KEGG.pdf"), width = 10, height = 6)
      dotplot(compKEGG, title = "KEGG pathway enrichment analysis for cancer specific DMRs associated genes") + labs(x = "") + scale_x_discrete(guide = guide_axis(angle = 45)) 
      dev.off()
      
      compGO <- compareCluster(geneCluster   = dmrs_genes,
                               fun           = "enrichGO",
                               OrgDb='org.Hs.eg.db',
                               pvalueCutoff  = 0.05,
                               pAdjustMethod = "BH")
      
      pdf(paste0(peak_type, "_cancer_type_specific_DMRs_annotation_genes_GO.pdf"), width = 10, height = 6)
      dotplot(compGO,  title = "GO enrichment analysis for cancer specific DMRs associated genes") + labs(x = "") + scale_x_discrete(guide = guide_axis(angle = 45)) 
      dev.off()
      }
      
      ###############################################
      ## genes functional profiles for specific DMRs with promoter region 
      # dmrs_genes = lapply(dmrAnnoList, function(i) unique(as.data.frame(i)$geneId))
      {
        dmrs_genes <- list()
        N <- length(dmrAnnoList)
        frac_p <- vector()
        for(k in 1:N)
        {
          dmr_data = as.data.frame(dmrAnnoList[[k]])
          #idx_promoter <- dmr_data[, 9] == "Promoter (<=1kb)" | dmr_data[, 9] == "Promoter (1-2kb)" | dmr_data[, 9] == "Promoter (2-3kb)" 
          idx_promoter <- dmr_data[, 9] == "Promoter"
          dmr_promoter <- dmr_data[idx_promoter, ]
          
          frac_p[k] <- (nrow(dmr_promoter) / nrow(dmr_data))
          
          dmrs_genes[[k]] =  unique(dmr_promoter[, 15])   ## need to be list, 
        }
        
        names(dmrs_genes) <- names(frac_p)<- names(dmrAnnoList)
        median(frac_p)
        mean(frac_p)
        
        compKEGG <- compareCluster(geneCluster   = dmrs_genes,
                                   fun           = "enrichKEGG",
                                   pvalueCutoff  = 0.05,
                                   pAdjustMethod = "BH")
        
        write.csv(data.frame(compKEGG), "Cancer_Specific_Hyper_DMRs_Promoter1.5K_Genes_KEGG_all.csv", row.names = F)
        
        
        pdf(paste0(peak_type, "_cancer_type_specific_DMRs_Promoter1.5K_genes_KEGG.pdf"), width = 10, height = 6)
        dotplot(compKEGG, showCategory=5, title = "Top enriched KEGG pathways for cancer-specific hyper-DMRs ") + labs(x = "") + scale_x_discrete(guide = guide_axis(angle = 45)) 
        dev.off()
        
        compGO <- compareCluster(geneCluster   = dmrs_genes,
                                 fun           = "enrichGO",
                                 OrgDb='org.Hs.eg.db',
                                 pvalueCutoff  = 0.05,
                                 pAdjustMethod = "BH")
        
        compGO_sim <- clusterProfiler::simplify(compGO )
        
        
        write.csv(data.frame(compGO_sim), "Cancer_Specific_Hyper_DMRs_Promoter1.5K_Genes_GO_all.csv", row.names = F)
        
        pdf(paste0(peak_type, "_cancer_type_specific_DMRs_Promoter1.5_genes_GO.pdf"), width = 10, height = 6)
        dotplot(compGO_sim,  showCategory=5,  title = "Top enriched GO terms for cancer-specific hyper-DMRs") + labs(x = "") + scale_x_discrete(guide = guide_axis(angle = 45)) 
        dev.off()
      }
      
      ###########################################################################
      ## DMRs based GSEA analysis
      {
        
        library(missMethyl)
        miss_KEGG <- miss_GO <- data.frame()
        ol_KEGG <- ol_GO <- list()
        ## list bed files 
        
        files <- paste(names(sort(hyper_spec_cnt, decreasing = T)), "specific", peak_type, "DMRs.bed")
        #files <- as.list(files)
        names(files) <- names(sort(hyper_spec_cnt, decreasing = T))
        
        ## need to be lifeover to hg19
        library(rtracklayer)
        path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
        lift_ch <- import.chain(path)     ## liftover chain
        
        L <- length(files)
        
        for(i in 1:L)  
        {
          print(i)
          dmr_bed_hg38 <- read.table(files[i])
          dmr_bed_gr_hg38 <- GRanges(seqnames = dmr_bed_hg38$V1,
                                     ranges = IRanges(start = dmr_bed_hg38$V2,
                                                     end = dmr_bed_hg38$V3,
                                                     names = dmr_bed_hg38$V4))
          seqlevelsStyle(dmr_bed_gr_hg38) = "UCSC"  # necessary
          dmr_bed_gr_hg19 = liftOver(dmr_bed_gr_hg38, lift_ch)
          class(dmr_bed_gr_hg19)
          
          dmr_bed_gr_hg19 = unlist(dmr_bed_gr_hg19)
          
          #######
          ## KEGG
          {
            dmrKEGG <- goregion(dmr_bed_gr_hg19, collection="KEGG", array.type="EPIC", 
                                genomic.features = "TSS1500")
            
            #idx_sig <- dmrKEGG$FDR < 0.05
            idx_sig <- dmrKEGG$P.DE < 0.05
            sig_dmrKEGG <- dmrKEGG[idx_sig, ]
            
            Cluster <- rep(names(files)[i], nrow(sig_dmrKEGG))
            
            if(length(Cluster) > 0)
            {
              term_id <- rownames(sig_dmrKEGG)
              sig_dmrKEGG <- cbind(Cluster, term_id, sig_dmrKEGG)
              miss_KEGG <- rbind(miss_KEGG, sig_dmrKEGG)
              
              ## ## overlap with cluster profiler results
              idx_cmp <- compKEGG@compareClusterResult$Cluster == names(files)[i]
              idx_ol <- !is.na(match(compKEGG@compareClusterResult$ID, rownames(sig_dmrKEGG)))
              
              idx_ss <- idx_cmp & idx_ol
              ol_KEGG[[i]] <- compKEGG[idx_ss]
              names(ol_KEGG)[i] <- names(files)[i]
    
            }
          
          }
          
          #######
          ## GO
          {
            dmrGO <- goregion(dmr_bed_gr_hg19, collection="GO", array.type="EPIC", 
                                genomic.features = "TSS1500")
            
            #idx_sig <- dmrGO$FDR < 0.05
            idx_sig <- dmrGO$P.DE < 0.05
            sig_dmrGO <- dmrGO[idx_sig, ]
            
            Cluster <- rep(names(files)[i], nrow(sig_dmrGO))
            
            if(length(Cluster) > 0)
            {
              term_id <- rownames(sig_dmrGO)
              sig_dmrGO <- cbind(Cluster, term_id, sig_dmrGO)
              miss_GO <- rbind(miss_GO, sig_dmrGO)
              
              ## ## overlap with cluster profiler results
              idx_cmp <- compGO_sim@compareClusterResult$Cluster == names(files)[i]
              idx_ol <- !is.na(match(compGO_sim@compareClusterResult$ID, rownames(sig_dmrGO)))
              
              idx_ss <- idx_cmp & idx_ol
              ol_GO[[i]] <- compGO_sim[idx_ss]
              names(ol_GO)[i] <- names(files)[i]
              
            }
            
          }
          
          
          
        
        }
        
        ###############
        ## output  KEGG
        write.csv(miss_KEGG, "Cancer_specific_Hyper_DMRs_TSS1500_Genes_KEGG_all_byMissMehthyl.csv", row.names = F)
        
        ## ploting overlaps
        compKEGG_plot_ol  <- merge_result(ol_KEGG)    ## dataframe to clusterprofiler object 
        
        pdf("Cancer_specific_Hyper_DMRs_Promoter1.5K_Genes_KEGG_clusterprofiler_MissMehthyl_overlapped_top10.pdf", width = 7.5, height = 6)
        dotplot(compKEGG_plot_ol , showCategory=5, title = "Top enriched KEGG pathways") + labs(x = "") + scale_x_discrete(guide = guide_axis(angle = 0)) 
        dev.off()
        
        # ol_KEGG <- data.frame(compKEGG_plot_ol)
        write.csv(data.frame(compKEGG_plot_ol), "Cancer_specific_Hyper_DMRs_TSS1500_Genes_KEGG_clusterprofiler_MissMehthyl_overlapped.csv", row.names = F)
        
        
        ##########
        ## out GO
        write.csv(miss_GO, "Cancer_specific_Hyper_DMRs_TSS1500_Genes_GO_all_byMissMehthyl.csv", row.names = F)
        
        ## ploting overlaps
        compGO_plot_ol  <- merge_result(ol_GO)    ## dataframe to clusterprofiler object 
        
        pdf("Cancer_specific_Hyper_DMRs_Promoter1.5K_Genes_GO_clusterprofiler_MissMehthyl_overlapped_top10.pdf", width = 8.5, height = 6)
        dotplot(compGO_plot_ol , showCategory=5, title = "Top enriched GO terms") + labs(x = "") + scale_x_discrete(guide = guide_axis(angle = 0)) 
        dev.off()
        
        # ol_GO <- data.frame(compGO_plot_ol)
        write.csv(data.frame(compGO_plot_ol), "Cancer_specific_Hyper_DMRs_TSS1500_Genes_GO_clusterprofiler_MissMehthyl_overlapped.csv", row.names = F)
        
        
           
       }
    }       
       
    
    ## for hypo
    peak_type <- "Hypo"
    {
      files <- paste(names(sort(hypo_spec_cnt, decreasing = T)), "specific", peak_type, "DMRs.bed")
      files <- as.list(files)
      names(files) <- names(sort(hypo_spec_cnt, decreasing = T))
      
      ## distribution around TSS
      tagMatrixList <- lapply(files, getTagMatrix, windows=promoter)
      
      idx_col <- match(names(tagMatrixList), levels(s_info$cancer_type))
      cancer_type_col_s <- cancer_type_col[idx_col]
      
      pdf(paste0(peak_type, "_cancer_type_specific_DMRs_around_TSS.pdf"), width = 5, height = 2)
      plotAvgProf(tagMatrixList, xlim=c(-1500, 1500)) + scale_color_manual(values = cancer_type_col_s) + theme_classic() 
      dev.off()
      
      ## at leat 100 DMRs
      {
        tagMatrixList_100  <- lapply(files[1:6], getTagMatrix, windows=promoter)
        idx_col_100 <- match(names(tagMatrixList_100), levels(s_info$cancer_type))
        cancer_type_col_s_100 <- cancer_type_col[idx_col_100]
       
         pdf(paste0(peak_type, "_cancer_type_specific_DMRs_around_TSS_100.pdf"), width = 5, height = 2)
        plotAvgProf(tagMatrixList_100, xlim=c(-1500, 1500)) + scale_color_manual(values = cancer_type_col_s_100) + theme_classic() 
        dev.off()
      }
      
      ## peak annotation
      dmrAnnoList <- lapply(files, annotatePeak, TxDb=txdb,
                            tssRegion=c(-1500, 1500), verbose=FALSE)
      
      Anno_col <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#b15928')
      
      pdf(paste0(peak_type, "_cancer_type_specific_DMRs_annotation.pdf"), width = 5, height = 4)
      plotAnnoBar(dmrAnnoList, title = "") +  scale_fill_manual(values = Anno_col)  + theme_classic() 
      dev.off()
      
      
      
      
      ###############################################
      ## genes functional profiles for all specific DMRs
      {
        dmrs_genes = lapply(dmrAnnoList, function(i) unique(as.data.frame(i)$geneId))
        
        compKEGG <- compareCluster(geneCluster   = dmrs_genes,
                                   fun           = "enrichKEGG",
                                   pvalueCutoff  = 0.05,
                                   pAdjustMethod = "BH")
        
        pdf(paste0(peak_type, "_cancer_type_specific_DMRs_annotation_genes_KEGG.pdf"), width = 10, height = 6)
        dotplot(compKEGG, title = "KEGG pathway enrichment analysis for cancer specific DMRs associated genes") + labs(x = "") + scale_x_discrete(guide = guide_axis(angle = 45)) 
        dev.off()
        
        compGO <- compareCluster(geneCluster   = dmrs_genes,
                                 fun           = "enrichGO",
                                 OrgDb='org.Hs.eg.db',
                                 pvalueCutoff  = 0.05,
                                 pAdjustMethod = "BH")
        
        pdf(paste0(peak_type, "_cancer_type_specific_DMRs_annotation_genes_GO.pdf"), width = 10, height = 6)
        dotplot(compGO,  title = "GO enrichment analysis for cancer specific DMRs associated genes") + labs(x = "") + scale_x_discrete(guide = guide_axis(angle = 45)) 
        dev.off()
      }
      
      ###############################################
      ## genes functional profiles for specific DMRs with promoter region 
      # dmrs_genes = lapply(dmrAnnoList, function(i) unique(as.data.frame(i)$geneId))
      {
        dmrs_genes <- list()
        N <- length(dmrAnnoList)
        frac_p <- vector()
        for(k in 1:N)
        {
          dmr_data = as.data.frame(dmrAnnoList[[k]])
          #idx_promoter <- dmr_data[, 9] == "Promoter (<=1kb)" | dmr_data[, 9] == "Promoter (1-2kb)" | dmr_data[, 9] == "Promoter (2-3kb)" 
          idx_promoter <- dmr_data[, 9] == "Promoter"
          dmr_promoter <- dmr_data[idx_promoter, ]
          
          frac_p[k] <- (nrow(dmr_promoter) / nrow(dmr_data))
          
          dmrs_genes[[k]] =  unique(dmr_promoter[, 15])   ## need to be list, 
        }
        
        names(dmrs_genes) <- names(frac_p)<- names(dmrAnnoList)
        median(frac_p)
        mean(frac_p)
        
        compKEGG <- compareCluster(geneCluster   = dmrs_genes,
                                   fun           = "enrichKEGG",
                                   pvalueCutoff  = 0.05,
                                   pAdjustMethod = "BH")
        
        write.csv(data.frame(compKEGG), "Cancer_Specific_Hypo_DMRs_Promoter1.5K_Genes_KEGG_all.csv", row.names = F)
        
        
        pdf(paste0(peak_type, "_cancer_type_specific_DMRs_Promoter1.5K_genes_KEGG.pdf"), width = 10, height = 6)
        dotplot(compKEGG, showCategory=5, title = "Top enriched KEGG pathways for cancer-specific hypo-DMRs ") + labs(x = "") + scale_x_discrete(guide = guide_axis(angle = 45)) 
        dev.off()
        
        compGO <- compareCluster(geneCluster   = dmrs_genes,
                                 fun           = "enrichGO",
                                 OrgDb='org.Hs.eg.db',
                                 pvalueCutoff  = 0.05,
                                 pAdjustMethod = "BH")
        
        compGO_sim <- clusterProfiler::simplify(compGO )
        
        
        write.csv(data.frame(compGO_sim), "Cancer_Specific_Hypo_DMRs_Promoter1.5K_Genes_GO_all.csv", row.names = F)
        
        pdf(paste0(peak_type, "_cancer_type_specific_DMRs_Promoter1.5_genes_GO.pdf"), width = 10, height = 6)
        dotplot(compGO_sim,  showCategory=5,  title = "Top enriched GO terms for cancer-specific hypo-DMRs") + labs(x = "") + scale_x_discrete(guide = guide_axis(angle = 45)) 
        dev.off()
      }
      
      ###########################################################################
      ## DMRs based GSEA analysis
      {
        
        library(missMethyl)
        miss_KEGG <- miss_GO <- data.frame()
        ol_KEGG <- ol_GO <- list()
        ## list bed files 
        
        files <- paste(names(sort(hypo_spec_cnt, decreasing = T)), "specific", peak_type, "DMRs.bed")
        #files <- as.list(files)
        names(files) <- names(sort(hypo_spec_cnt, decreasing = T))
        
        ## need to be lifeover to hg19
        library(rtracklayer)
        path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
        lift_ch <- import.chain(path)     ## liftover chain
        
        L <- length(files)
        
        for(i in 1:L)
        {
          print(i)
          dmr_bed_hg38 <- read.table(files[i])
          dmr_bed_gr_hg38 <- GRanges(seqnames = dmr_bed_hg38$V1,
                                     ranges = IRanges(start = dmr_bed_hg38$V2,
                                                      end = dmr_bed_hg38$V3,
                                                      names = dmr_bed_hg38$V4))
          seqlevelsStyle(dmr_bed_gr_hg38) = "UCSC"  # necessary
          dmr_bed_gr_hg19 = liftOver(dmr_bed_gr_hg38, lift_ch)
          class(dmr_bed_gr_hg19)
          
          dmr_bed_gr_hg19 = unlist(dmr_bed_gr_hg19)
          
          #######
          ## KEGG
          {
            dmrKEGG <- goregion(dmr_bed_gr_hg19, collection="KEGG", array.type="EPIC", 
                                genomic.features = "TSS1500")
            
            #idx_sig <- dmrKEGG$FDR < 0.05
            idx_sig <- dmrKEGG$P.DE < 0.05
            sig_dmrKEGG <- dmrKEGG[idx_sig, ]
            
            Cluster <- rep(names(files)[i], nrow(sig_dmrKEGG))
            
            if(length(Cluster) > 0)
            {
              term_id <- rownames(sig_dmrKEGG)
              sig_dmrKEGG <- cbind(Cluster, term_id, sig_dmrKEGG)
              miss_KEGG <- rbind(miss_KEGG, sig_dmrKEGG)
              
              ## ## overlap with cluster profiler results
              idx_cmp <- compKEGG@compareClusterResult$Cluster == names(files)[i]
              idx_ol <- !is.na(match(compKEGG@compareClusterResult$ID, rownames(sig_dmrKEGG)))
              
              idx_ss <- idx_cmp & idx_ol
              ol_KEGG[[i]] <- compKEGG[idx_ss]
              names(ol_KEGG)[i] <- names(files)[i]
              
            }
            
           
          }
          
          #######
          ## GO
          {
            dmrGO <- goregion(dmr_bed_gr_hg19, collection="GO", array.type="EPIC", 
                              genomic.features = "TSS1500")
            
            #idx_sig <- dmrGO$FDR < 0.05
            idx_sig <- dmrGO$P.DE < 0.05
            sig_dmrGO <- dmrGO[idx_sig, ]
            
            Cluster <- rep(names(files)[i], nrow(sig_dmrGO))
            
            if(length(Cluster) > 0)
            {
              term_id <- rownames(sig_dmrGO)
              sig_dmrGO <- cbind(Cluster, term_id, sig_dmrGO)
              miss_GO <- rbind(miss_GO, sig_dmrGO)
              
              ## ## overlap with cluster profiler results
              idx_cmp <- compGO_sim@compareClusterResult$Cluster == names(files)[i]
              idx_ol <- !is.na(match(compGO_sim@compareClusterResult$ID, rownames(sig_dmrGO)))
              
              idx_ss <- idx_cmp & idx_ol
              ol_GO[[i]] <- compGO_sim[idx_ss]
              names(ol_GO)[i] <- names(files)[i]
              
            }
            
          }
          
        }
        
        ###############
        ## output  KEGG
        write.csv(miss_KEGG, "Cancer_specific_Hypo_DMRs_TSS1500_Genes_KEGG_all_byMissMehthyl.csv", row.names = F)
        
        ## ploting overlaps
        compKEGG_plot_ol  <- merge_result(ol_KEGG)    ## dataframe to clusterprofiler object 
        
        pdf("Cancer_specific_Hypo_DMRs_Promoter1.5K_Genes_KEGG_clusterprofiler_MissMehthyl_overlapped_top10.pdf", width = 7, height = 6)
        dotplot(compKEGG_plot_ol , showCategory=5, title = "Top enriched KEGG pathways") + labs(x = "") + scale_x_discrete(guide = guide_axis(angle = 0)) 
        dev.off()
        
        # ol_KEGG <- data.frame(compKEGG_plot_ol)
        write.csv(data.frame(compKEGG_plot_ol), "Cancer_specific_Hypo_DMRs_TSS1500_Genes_KEGG_clusterprofiler_MissMehthyl_overlapped.csv", row.names = F)
        
        
        ##########
        ## out GO
        write.csv(miss_GO, "Cancer_specific_Hypo_DMRs_TSS1500_Genes_GO_all_byMissMehthyl.csv", row.names = F)
        
        ## ploting overlaps
        compGO_plot_ol  <- merge_result(ol_GO)    ## dataframe to clusterprofiler object 
        
        pdf("Cancer_specific_Hypo_DMRs_Promoter1.5K_Genes_GO_clusterprofiler_MissMehthyl_overlapped_top10.pdf", width = 10.5, height = 6)
        dotplot(compGO_plot_ol , showCategory=5, title = "Top enriched GO pathways") + labs(x = "") + scale_x_discrete(guide = guide_axis(angle = 0)) 
        dev.off()
        
        # ol_GO <- data.frame(compGO_plot_ol)
        write.csv(data.frame(compGO_plot_ol), "Cancer_specific_Hypo_DMRs_TSS1500_Genes_GO_clusterprofiler_MissMehthyl_overlapped.csv", row.names = F)
        
        
      }
    }  
  }
  
}


##############################################
## Check overlapping of DMRs using all samples
##############################################
if(FALSE){
  library("ggVennDiagram")
  
  ### differential peaks overlapping per matched cancer types
  ### with SE and PE reuslt
  DMR_bins_spe <- readRDS("../ALL/FC_2_cancer_type_DESeq_DMRs_bin_IDs.RDS")
  
  idx_m <- match(names(DMR_bins), names(DMR_bins_spe))
  DMR_bins_spe_m <- DMR_bins_spe[idx_m]
  
  L <- length(DMR_bins_spe_m)
  
  for (i in seq(1, L, 2)){
    dmrs_cmp <- list(Hyper = DMR_bins_spe_m[[i]],  Hyper_SE = DMR_bins[[i]],
                     Hypo_SE = DMR_bins[[i + 1]], Hypo = DMR_bins_spe_m[[i + 1]])
    
    g_title <- strsplit(names(DMR_bins)[i], "_")[[1]][1]
    ggVennDiagram(dmrs_cmp) + ggtitle(g_title) + theme(plot.title = element_text(hjust = 0.5, vjust = -100), legend.position = "none")
    ggsave(paste0("DMRs of ", g_title, " compaired to SE study.pdf"), width = 5, height = 5) 
  }
  
}