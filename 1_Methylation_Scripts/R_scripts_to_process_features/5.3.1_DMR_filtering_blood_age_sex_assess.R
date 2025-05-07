rm(list = ls())
setwd("/Users/yong/Library/CloudStorage/OneDrive-UHN/Projects/TCGE/cfEpigenomics/Resource/2_batch_norm/DMRs_revision/DMR_filtering")

library(ggplot2)
library(GenomicFeatures)

## FC = 2 was used
dmr_files <- c("/Users/yong/Library/CloudStorage/OneDrive-UHN/Projects/TCGE/cfEpigenomics/Resource/2_batch_norm/DMRs_revision/PE/FC_2_cancer_type_DESeq_DMRs_bin_IDs.RDS",
               "/Users/yong/Library/CloudStorage/OneDrive-UHN/Projects/TCGE/cfEpigenomics/Resource/2_batch_norm/DMRs_revision/PE/FC_2_cancer_or_healthy_DESeq_DMRs_bin_IDs.RDS",
               "/Users/yong/Library/CloudStorage/OneDrive-UHN/Projects/TCGE/cfEpigenomics/Resource/2_batch_norm/DMRs_revision/SE/FC_2_cancer_type_DESeq_DMRs_bin_IDs.RDS",
               "/Users/yong/Library/CloudStorage/OneDrive-UHN/Projects/TCGE/cfEpigenomics/Resource/2_batch_norm/DMRs_revision/SE/FC_2_cancer_or_healthy_DESeq_DMRs_bin_IDs.RDS",
               "/Users/yong/Library/CloudStorage/OneDrive-UHN/Projects/TCGE/cfEpigenomics/Resource/2_batch_norm/DMRs_revision/PE/PE_Specific/FC_2_DESeq_cancer_specific_DMRs_bin_IDs.RDS",
               "/Users/yong/Library/CloudStorage/OneDrive-UHN/Projects/TCGE/cfEpigenomics/Resource/2_batch_norm/DMRs_revision/SE/SE_Specific/FC_2_DESeq_cancer_specific_DMRs_bin_IDs.RDS")


#####################################################
## for SE and PE DMRs after considering mixed effects
## filtering base on healthy PBL signals
#####################################################
{
## bins lowly methylated in healthy PBL retained 
## using medetrand rms 0.1 
pbl <- read.csv("/Users/yong/Library/CloudStorage/OneDrive-UHN/Projects/TCGE/cfEpigenomics/Resource/2_batch_norm/merge/PBL/Healthy_PBL_bins_median_rms_smaller_than_0.1.csv")
bin_keep <- as.character(pbl$x)

DMR_PBL_list <- c("PE_DMRs_per_cancer_type_after_PBL_depletion",
             "PE_DMRs_cancer_vs_healthy_after_PBL_depletion",
             "SE_DMRs_per_cancer_type_after_PBL_depletion",
             "SE_DMRs_cancer_vs_healthy_after_PBL_depletion",
             "PE_cancer_specific_DMRs_after_PBL_depletion",
             "SE_cancer_specific_DMRs_after_PBL_depletion")

for(j in 1:6)
{
dmr <- readRDS(dmr_files[j])
tt <- DMR_PBL_list[j]

dmr_keep <- list()
count <- vector()  ## overlap and non-overlapped between dmr and bin_keep

k = 1
L <- length(dmr)
for(i in 1:L)
{
  
  dmr_keep[[i]] <- intersect(dmr[[i]], bin_keep)
  
  count[k] <- length(dmr_keep[[i]])
  count[k + 1] <- length(dmr[[i]]) - length(dmr_keep[[i]])
  
  k <- k + 2
}

names(dmr_keep) <- names(dmr)
saveRDS(dmr_keep, paste0(tt, "_bin_IDs.RDS") )

dmr_group <- rep(names(dmr), each = 2)
ol_group <- rep(c("Retained", "PBL_filtered"), length(dmr))

data <- data.frame(count, ol_group, dmr_group)

## for count 
g <- ggplot(data, aes(y=count, x = dmr_group, fill = ol_group)) 
g <- g + geom_bar(position="stack", stat="identity") + labs(x = "", y = "Number of DMRs", title = tt)
g <- g + scale_fill_manual(values = c("gray", "darkcyan"))
g <- g + theme_classic() + coord_flip() + theme(legend.position = "bottom") 
ggsave(paste0(tt, "_count.pdf"), height = 10, width = 7)

## for percentage 
g <- ggplot(data, aes(y=count, x = dmr_group, fill = ol_group)) 
g <- g + geom_bar(position="fill", stat="identity") + labs(x = "", y = "Fraction of DMRs", title = tt)
g <- g + scale_fill_manual(values = c("gray", "darkcyan"))
g <- g + theme_classic() + coord_flip() + theme(legend.position = "bottom") 
ggsave(paste0(tt, "_percentages.pdf"), height = 10, width = 7)


## output summary data
dmr_count <- count[seq(1, length(count), 2)] + count[seq(2, length(count), 2)]
PBL_retained <- count[seq(1, length(count), 2)]
PBL_filtered <- count[seq(2, length(count), 2)]
PBL_retained_frac <- count[seq(1, length(count), 2)] / dmr_count 
PBL_filtered_frac <- count[seq(2, length(count), 2)] / dmr_count 

out <- data.frame(dmr_count, PBL_retained, PBL_filtered, PBL_retained_frac, PBL_filtered_frac)
rownames(out) <- dmr_group[seq(1, length(count), 2)]
write.csv(out, paste0(tt, "_summary.csv"))

}

}


#####################################################
## for SE and PE DMRs after considering mixed effects
## checking the overlapping with cfDNA tissure diconvolutoin reference
#####################################################
{
ref <- read.csv("/Users/yong/Library/CloudStorage/OneDrive-UHN/Projects/TCGE/cfEpigenomics/Resource/4_tissue_deconvolution/Moss_et_al_2018_hg38_clean.csv")
dim(ref)
## cell types
cells <- unique(ref$name)
unique(ref$name)
## CpG associated with belowing cells
#[1] "Monocytes EPIC"            
#[2] "B-cells EPIC"              
#[3] "CD4T-cells EPIC"           
#[4] "NK-cells EPIC"             
#[5] "CD8T-cells EPIC"           
#[6] "Neutrophils EPIC"          
#[7] "Erythrocyte progenitors"   
#[15] "Vascular endothelial cells"

idx_c <- match(ref$name, cells[c(1:7, 15)])
ref_s <- ref[!is.na(idx_c), ]

## extend to 300 bin as well 
## start and end have been extend to 100bp 
ref_s_gr <-  GRanges(seqnames = ref_s$seqnames,
                      ranges = IRanges(start = ref_s$start - 100,
                                       end = ref_s$end + 99,
                                      names = ref_s$acc))

## 
tt_list <- c("PE_DMRs_per_cancer_type_overlap_blood_cells_specific_CpGs_ext300",
             "PE_DMRs_cancer_vs_healthy_overlap_blood_cells_specific_CpGs_ext300",
             "SE_DMRs_per_cancer_type_overlap_blood_cells_specific_CpGs_ext300",
             "SE_DMRs_cancer_vs_healthy_overlap_blood_cells_specific_CpGs_ext300",
             "PE_cancer_specific_DMRs_overlap_blood_cells_specific_CpGs_ext300",
             "SE_cancer_specific_DMRs_overlap_blood_cells_specific_CpGs_ext300")

bed_ref <- read.table("/Users/yong/Library/CloudStorage/OneDrive-UHN/Projects/TCGE/cfEpigenomics/Resource/2_batch_norm/merge/ALL_auto_bfilt_wCpG_bins.bed")

for(j in 1:6)
{
  
  tt <- tt_list[j]
  
  ## before PBL filtering 
  {
  dmr <- readRDS(dmr_files[j])

  out <- data.frame()
  
  k = 1
  L <- length(dmr)
  
  for(i in 1:L)
  {
    
    idx_id <- match(dmr[[i]], bed_ref$V4)
    
    dmr_name <- names(dmr)[i]
    
    dmr_bed <- bed_ref[idx_id, ]
    dmr_bed_gr <- GRanges(seqnames = dmr_bed$V1,
                          ranges = IRanges(start = dmr_bed$V2,
                                           end = dmr_bed$V3,
                                           names = dmr_bed$V4))
    
    ## checking overlappping 
    ol <- findOverlaps(dmr_bed_gr, ref_s_gr)
    
    ## 
    dmr_n <-  length(dmr_bed_gr)
    dmr_ol <- length(unique(ol@from))
    dmr_ol_f <- dmr_ol / dmr_n
    ref_s_n <- length(ref_s_gr)
    ref_s_ol <- length(unique(ol@to))
    
    #n_bin <- 7445098    # ALL_auto_bfilt_wCpG_bins.bed
    #hyper_p <- phyper(dmr_ol, dmr_n,  n_bin , ref_s_n)
    
    cmb <- cbind(dmr_n, dmr_ol, dmr_ol_f, ref_s_n , ref_s_ol)
    rownames(cmb) <- dmr_name
    out <- rbind(out, cmb)
  }
  }
  
  
  ## after PBL filtering 
  {
    ## after PBL filtering 
    dmr <- readRDS(paste0(DMR_PBL_list[j], "_bin_IDs.RDS"))
    
    out_pbl <- data.frame()
    
    k = 1
    L <- length(dmr)
    
    for(i in 1:L)
    {
      
      idx_id <- match(dmr[[i]], bed_ref$V4)
      
      dmr_name <- names(dmr)[i]
      
      dmr_bed <- bed_ref[idx_id, ]
      dmr_bed_gr <- GRanges(seqnames = dmr_bed$V1,
                            ranges = IRanges(start = dmr_bed$V2,
                                             end = dmr_bed$V3,
                                             names = dmr_bed$V4))
      
      ## checking overlappping 
      ol <- findOverlaps(dmr_bed_gr, ref_s_gr)
      
      ## 
      dmr_n <-  length(dmr_bed_gr)
      dmr_ol <- length(unique(ol@from))
      dmr_ol_f <- dmr_ol / dmr_n
      ref_s_n <- length(ref_s_gr)
      ref_s_ol <- length(unique(ol@to))
      
      #n_bin <- 7445098    # ALL_auto_bfilt_wCpG_bins.bed
      #hyper_p <- phyper(dmr_ol, dmr_n,  n_bin , ref_s_n)
      
      cmb <- cbind(dmr_n, dmr_ol, dmr_ol_f, ref_s_n , ref_s_ol)
      rownames(cmb) <- dmr_name
      out_pbl <- rbind(out_pbl, cmb)
    }
    
    colnames(out_pbl) <- paste0("PBL_filtered_", colnames(out_pbl))
  }
  

  PBL_filtered_dmr_ol_over_dmr_n <- out_pbl$PBL_filtered_dmr_ol / out$dmr_n
  out_cmb <- cbind(out, out_pbl, PBL_filtered_dmr_ol_over_dmr_n)
  
  print(tt)
  print("fraction_of_dmrs_ol_blood_cells_specific_CpGs:")
  print(summary(out_cmb$dmr_ol_f))
  print("fraction_of_dmrs_after_pbL_filitering_ol_blood_cells_specific_CpGs:")
  print(summary(out_cmb$PBL_filtered_dmr_ol_over_dmr_n))
  
  write.csv(out_cmb, paste0(tt, "_before_and_after_PBL_filtering.csv"))
  
  
}

}



#######################################
## DMRs after filtering for age and sex
#######################################
DMR_PBL_list <- c("PE_DMRs_per_cancer_type_after_PBL_depletion",
                  "PE_DMRs_cancer_vs_healthy_after_PBL_depletion",
                  "SE_DMRs_per_cancer_type_after_PBL_depletion",
                  "SE_DMRs_cancer_vs_healthy_after_PBL_depletion",
                  "PE_cancer_specific_DMRs_after_PBL_depletion",
                  "SE_cancer_specific_DMRs_after_PBL_depletion")

bed_ref <- read.table("/Users/yong/Library/CloudStorage/OneDrive-UHN/Projects/TCGE/cfEpigenomics/Resource/2_batch_norm/merge/ALL_auto_bfilt_wCpG_bins.bed")


#####################################################
## for SE and PE DMRs after considering mixed effects
## checking the overlapping with age associated CpGs 
#####################################################
{
  
  # ref <- read.table("/Users/yong/Library/CloudStorage/OneDrive-UHN/Projects/TCGE/cfEpigenomics/Resource/2_batch_norm/methylclock/DNA_methylclocks_associated_CpGs_hg38.bed")
  ref <- read.table("/Users/yong/Library/CloudStorage/OneDrive-UHN/Projects/TCGE/cfEpigenomics/Resource/2_batch_norm/age_CpGs/Cor_age/norm/Cor0.6_combined_age_associated_regions_in_hg38.bed")  ##have extended to 300bp
  dim(ref)
  
  ## extend to 300 bin as well 
  ref_s_gr <-  GRanges(seqnames = ref$V1,
                       ranges = IRanges(start = ref$V2,    
                                        end = ref$V3,
                                        names = ref$V4))
  
  ## 
  tt_list <- c("PE_DMRs_per_cancer_type_overlap_age_associted_CpGs_ext300",
               "PE_DMRs_cancer_vs_healthy_overlap_age_associted_CpGs_ext300",
               "SE_DMRs_per_cancer_type_overlap_age_associted_CpGs_ext300",
               "SE_DMRs_cancer_vs_healthy_overlap_age_associted_CpGs_ext300",
               "PE_cancer_specific_DMRs_overlap_age_associted_CpGs_ext300",
               "SE_cancer_specific_DMRs_overlap__age_associted__CpGs_ext300")
  
  
  for(j in 1:6)
  {
    
    tt <- tt_list[j]
    
    ## before PBL filtering 
    {
      dmr <- readRDS(dmr_files[j])
      
      out <- data.frame()
      
      k = 1
      L <- length(dmr)
      
      for(i in 1:L)
      {
        
        idx_id <- match(dmr[[i]], bed_ref$V4)
        
        dmr_name <- names(dmr)[i]
        
        dmr_bed <- bed_ref[idx_id, ]
        dmr_bed_gr <- GRanges(seqnames = dmr_bed$V1,
                              ranges = IRanges(start = dmr_bed$V2,
                                               end = dmr_bed$V3,
                                               names = dmr_bed$V4))
        
        ## checking overlappping 
        ol <- findOverlaps(dmr_bed_gr, ref_s_gr)
        
        ## 
        dmr_n <-  length(dmr_bed_gr)
        dmr_ol <- length(unique(ol@from))
        dmr_ol_f <- dmr_ol / dmr_n
        ref_s_n <- length(ref_s_gr)
        ref_s_ol <- length(unique(ol@to))
        
        #n_bin <- 7445098    # ALL_auto_bfilt_wCpG_bins.bed
        #hyper_p <- phyper(dmr_ol, dmr_n,  n_bin , ref_s_n)
        
        cmb <- cbind(dmr_n, dmr_ol, dmr_ol_f, ref_s_n , ref_s_ol)
        rownames(cmb) <- dmr_name
        out <- rbind(out, cmb)
      }
    }
    
    
    ## after PBL filtering 
    {
      ## after PBL filtering 
      dmr <- readRDS(paste0(DMR_PBL_list[j], "_bin_IDs.RDS"))
      
      out_pbl <- data.frame()
      
      k = 1
      L <- length(dmr)
      
      for(i in 1:L)
      {
        
        idx_id <- match(dmr[[i]], bed_ref$V4)
        
        dmr_name <- names(dmr)[i]
        
        dmr_bed <- bed_ref[idx_id, ]
        dmr_bed_gr <- GRanges(seqnames = dmr_bed$V1,
                              ranges = IRanges(start = dmr_bed$V2,
                                               end = dmr_bed$V3,
                                               names = dmr_bed$V4))
        
        ## checking overlappping 
        ol <- findOverlaps(dmr_bed_gr, ref_s_gr)
        
        ## 
        dmr_n <-  length(dmr_bed_gr)
        dmr_ol <- length(unique(ol@from))
        dmr_ol_f <- dmr_ol / dmr_n
        ref_s_n <- length(ref_s_gr)
        ref_s_ol <- length(unique(ol@to))
        
        #n_bin <- 7445098    # ALL_auto_bfilt_wCpG_bins.bed
        #hyper_p <- phyper(dmr_ol, dmr_n,  n_bin , ref_s_n)
        
        cmb <- cbind(dmr_n, dmr_ol, dmr_ol_f, ref_s_n , ref_s_ol)
        rownames(cmb) <- dmr_name
        out_pbl <- rbind(out_pbl, cmb)
      }
      
      colnames(out_pbl) <- paste0("PBL_filtered_", colnames(out_pbl))
    }
    
    
    PBL_filtered_dmr_ol_over_dmr_n <- out_pbl$PBL_filtered_dmr_ol / out$dmr_n
    out_cmb <- cbind(out, out_pbl, PBL_filtered_dmr_ol_over_dmr_n)
    
    print(tt)
    print("fraction_of_dmrs_ol_age_associated_CpGs:")
    print(summary(out_cmb$dmr_ol_f))
    print("fraction_of_dmrs_after_pbL_filitering_ol__age_associated_CpGs:")
    print(summary(out_cmb$PBL_filtered_dmr_ol_over_dmr_n))
    
    write.csv(out_cmb, paste0(tt, "_before_and_after_PBL_filtering.csv"))
    
    
  }
  
}


#####################################################
## for SE and PE DMRs after considering mixed effects
## checking the overlapping with sex associated CpGs 
#####################################################
{
  
  #ref <- read.table("/Users/yong/Library/CloudStorage/OneDrive-UHN/Projects/TCGE/cfEpigenomics/Resource/2_batch_norm/sex_CpGs/Sex_associated_CpGs_in_autosomes_hg38.bed")
  #dim(ref)
  
  #ref_s_gr <-  GRanges(seqnames = ref$V1,
  #                     ranges = IRanges(start = ref$V2 - 150,
  #                                      end = ref$V3 + 149,
  #                                      names = ref$V4))
  
  ## have extended to 300 bin as well 
  ref <- read.table("/Users/yong/Library/CloudStorage/OneDrive-UHN/Projects/TCGE/cfEpigenomics/Resource/2_batch_norm/sex_CpGs/sex_DMRs/Combined_sex_associated_regions_in_autosomes_hg38.bed")
  ref_s_gr <-  GRanges(seqnames = ref$V1,
                       ranges = IRanges(start = ref$V2,
                                        end = ref$V3,
                                        names = ref$V4))
  
  ## 
  tt_list <- c("PE_DMRs_per_cancer_type_overlap_sex_associted_CpGs_ext300",
               "PE_DMRs_cancer_vs_healthy_overlap_sex_associted_CpGs_ext300",
               "SE_DMRs_per_cancer_type_overlap_sex_associted_CpGs_ext300",
               "SE_DMRs_cancer_vs_healthy_overlap_sex_associted_CpGs_ext300",
               "PE_cancer_specific_DMRs_overlap_sex_associted_CpGs_ext300",
               "SE_cancer_specific_DMRs_overlap__sex_associted__CpGs_ext300")
  
  
  for(j in 1:6)
  {
    
    tt <- tt_list[j]
    
    ## before PBL filtering 
    {
      dmr <- readRDS(dmr_files[j])
      
      out <- data.frame()
      
      k = 1
      L <- length(dmr)
      
      for(i in 1:L)
      {
        
        idx_id <- match(dmr[[i]], bed_ref$V4)
        
        dmr_name <- names(dmr)[i]
        
        dmr_bed <- bed_ref[idx_id, ]
        dmr_bed_gr <- GRanges(seqnames = dmr_bed$V1,
                              ranges = IRanges(start = dmr_bed$V2,
                                               end = dmr_bed$V3,
                                               names = dmr_bed$V4))
        
        ## checking overlappping 
        ol <- findOverlaps(dmr_bed_gr, ref_s_gr)
        
        ## 
        dmr_n <-  length(dmr_bed_gr)
        dmr_ol <- length(unique(ol@from))
        dmr_ol_f <- dmr_ol / dmr_n
        ref_s_n <- length(ref_s_gr)
        ref_s_ol <- length(unique(ol@to))
        
        #n_bin <- 7445098    # ALL_auto_bfilt_wCpG_bins.bed
        #hyper_p <- phyper(dmr_ol, dmr_n,  n_bin , ref_s_n)
        
        cmb <- cbind(dmr_n, dmr_ol, dmr_ol_f, ref_s_n , ref_s_ol)
        rownames(cmb) <- dmr_name
        out <- rbind(out, cmb)
      }
    }
    
    
    ## after PBL filtering 
    {
      ## after PBL filtering 
      dmr <- readRDS(paste0(DMR_PBL_list[j], "_bin_IDs.RDS"))
      
      out_pbl <- data.frame()
      
      k = 1
      L <- length(dmr)
      
      for(i in 1:L)
      {
        
        idx_id <- match(dmr[[i]], bed_ref$V4)
        
        dmr_name <- names(dmr)[i]
        
        dmr_bed <- bed_ref[idx_id, ]
        dmr_bed_gr <- GRanges(seqnames = dmr_bed$V1,
                              ranges = IRanges(start = dmr_bed$V2,
                                               end = dmr_bed$V3,
                                               names = dmr_bed$V4))
        
        ## checking overlappping 
        ol <- findOverlaps(dmr_bed_gr, ref_s_gr)
        
        ## 
        dmr_n <-  length(dmr_bed_gr)
        dmr_ol <- length(unique(ol@from))
        dmr_ol_f <- dmr_ol / dmr_n
        ref_s_n <- length(ref_s_gr)
        ref_s_ol <- length(unique(ol@to))
        
        #n_bin <- 7445098    # ALL_auto_bfilt_wCpG_bins.bed
        #hyper_p <- phyper(dmr_ol, dmr_n,  n_bin , ref_s_n)
        
        cmb <- cbind(dmr_n, dmr_ol, dmr_ol_f, ref_s_n , ref_s_ol)
        rownames(cmb) <- dmr_name
        out_pbl <- rbind(out_pbl, cmb)
      }
      
      colnames(out_pbl) <- paste0("PBL_filtered_", colnames(out_pbl))
    }
    
    
    PBL_filtered_dmr_ol_over_dmr_n <- out_pbl$PBL_filtered_dmr_ol / out$dmr_n
    out_cmb <- cbind(out, out_pbl, PBL_filtered_dmr_ol_over_dmr_n)
    
    print(tt)
    print("fraction_of_dmrs_ol_sex_associated_CpGs:")
    print(summary(out_cmb$dmr_ol_f))
    print("fraction_of_dmrs_after_pbL_filitering_ol__sex_associated_CpGs:")
    print(summary(out_cmb$PBL_filtered_dmr_ol_over_dmr_n))
    
    write.csv(out_cmb, paste0(tt, "_before_and_after_PBL_filtering.csv"))
    
    
  }
  
}


