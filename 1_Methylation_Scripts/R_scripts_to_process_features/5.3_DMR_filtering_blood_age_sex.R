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


dmr_filterd_list <- c("PE_DMRs_per_cancer_type_after_blood_cell_age_sex_filtered",
                      "PE_DMRs_cancer_vs_healthy_after_blood_cell_age_sex_filtered",
                      "SE_DMRs_per_cancer_type_after_blood_cell_age_sex_filtered",
                      "SE_DMRs_cancer_vs_healthy_after_blood_cell_age_sex_filtered",
                      "PE_cancer_specific_DMRs_after_blood_cell_age_sex_filtered",
                      "SE_cancer_specific_DMRs_after_blood_cell_age_sex_filtered")

## dmr bin info
bed_ref <- read.table("/Users/yong/Library/CloudStorage/OneDrive-UHN/Projects/TCGE/cfEpigenomics/Resource/2_batch_norm/merge/ALL_auto_bfilt_wCpG_bins.bed")

#####################################################
## for SE and PE DMRs after considering mixed effects
## filtering base on healthy PBL signals
#####################################################
{
  ################################################
  ## bins lowly methylated in healthy PBL retained 
  ## using medetrand rms 0.1 
  pbl <- read.csv("/Users/yong/Library/CloudStorage/OneDrive-UHN/Projects/TCGE/cfEpigenomics/Resource/2_batch_norm/merge/PBL/Healthy_PBL_bins_median_rms_smaller_than_0.1.csv")
  bin_keep <- as.character(pbl$x)

  ##############################
  ## blood cell specific regions
  blood_ref <- read.csv("/Users/yong/Library/CloudStorage/OneDrive-UHN/Projects/TCGE/cfEpigenomics/Resource/4_tissue_deconvolution/Moss_et_al_2018_hg38_clean.csv")
  
  ## cell types
  cells <- unique(blood_ref$name)
  unique(blood_ref$name)
  ## CpG associated with belowing cells
  #[1] "Monocytes EPIC"            
  #[2] "B-cells EPIC"              
  #[3] "CD4T-cells EPIC"           
  #[4] "NK-cells EPIC"             
  #[5] "CD8T-cells EPIC"           
  #[6] "Neutrophils EPIC"          
  #[7] "Erythrocyte progenitors"   
  #[15] "Vascular endothelial cells"
  
  idx_c <- match(blood_ref$name, cells[c(1:7, 15)])
  blood_ref_s <- blood_ref[!is.na(idx_c), ]
  dim(blood_ref_s)
  
  ## extend to 300 bin as well 
  blood_ref_s_gr <-  GRanges(seqnames = blood_ref_s$seqnames,
                       ranges = IRanges(start = blood_ref_s$start - 100,
                                        end = blood_ref_s$end + 99,
                                        names = blood_ref_s$acc))
  
  ## for AML specifically 
  ## only exclued B-cells, T cell, Nk cells and vascular endothelia cells
  idx_aml <- match(blood_ref$name, cells[c(2:5, 15)])
  aml_rm_ref <-  blood_ref[!is.na(idx_aml), ]
  aml_rm_ref_gr <-  GRanges(seqnames = aml_rm_ref$seqnames,
                             ranges = IRanges(start = aml_rm_ref$start - 100,
                                              end = aml_rm_ref$end + 99,
                                              names = aml_rm_ref$acc))
  
  #########################
  ## age associated regions 
  age_ref <- read.table("/Users/yong/Library/CloudStorage/OneDrive-UHN/Projects/TCGE/cfEpigenomics/Resource/2_batch_norm/age_CpGs/Cor_age/norm/Cor0.6_combined_age_associated_regions_in_hg38.bed")
  dim(age_ref)
  
  ## extend to 300 bin as well 
  age_ref_s_gr <-  GRanges(seqnames = age_ref$V1,
                       ranges = IRanges(start = age_ref$V2,
                                        end = age_ref$V3,
                                        names = age_ref$V4))
  
  #########################
  ## sex associated regions 
  sex_ref <- read.table("/Users/yong/Library/CloudStorage/OneDrive-UHN/Projects/TCGE/cfEpigenomics/Resource/2_batch_norm/sex_CpGs/sex_DMRs/Combined_sex_associated_regions_in_autosomes_hg38.bed")
  dim(sex_ref)
  
  ## extend to 300 bin as well 
  sex_ref_s_gr <-  GRanges(seqnames = sex_ref$V1,
                       ranges = IRanges(start = sex_ref$V2,
                                        end = sex_ref$V3,
                                        names = sex_ref$V4))
  
  
  
  ## combined regions to filter out 
  rm_gr <- c(blood_ref_s_gr, age_ref_s_gr, sex_ref_s_gr)
  
  ## for aml specifically 
  aml_rm_gr <- c(aml_rm_ref_gr, age_ref_s_gr, sex_ref_s_gr)
 
  bed_ref_gr <- GRanges(seqnames = bed_ref$V1,
                        ranges = IRanges(start = bed_ref$V2,
                                         end = bed_ref$V3,
                                         names = bed_ref$V4))
  
  ####################################################################
  ## bins to be removed based on blood, age and sex associated regions 
  ol <- findOverlaps(bed_ref_gr, rm_gr)
  bin_rm <- unique(names(bed_ref_gr)[ol@from])
  saveRDS(bin_rm, "Associated_with_blood_age_sex_bin_IDs.RDS")
  
  ## for AML 
  ol_aml <- findOverlaps(bed_ref_gr, aml_rm_gr)
  aml_bin_rm <- unique(names(bed_ref_gr)[ol_aml@from])
  saveRDS(aml_bin_rm, "Associated_with_blood_age_sex_bin_IDs_for_AML.RDS")
  
  
   
  for(j in 1:6)
  {
    print(j)
    dmr <- readRDS(dmr_files[j])
    tt <- dmr_filterd_list[j]

    dmr_keep <- list()
    L <- length(dmr)
    
    for(i in 1:L)
    {
      ## PBL depeltion 
      #print(names(dmr)[i])
      keep_tmp <- intersect(dmr[[i]], bin_keep)
      
      if((names(dmr)[i] == "AML vs Healthy_Hyper") |  (names(dmr)[i] == "AML vs Healthy_Hypo" | 
         (names(dmr)[i] == "cancer_typeAML_specific_Hyper") | (names(dmr)[i] == "cancer_typeAML_specific_Hypo"))){
        idx_rm <- match(keep_tmp, aml_bin_rm)
        dmr_keep[[i]] <- keep_tmp[is.na(idx_rm)]
      } else {
        idx_rm <- match(keep_tmp, bin_rm)
        dmr_keep[[i]] <- keep_tmp[is.na(idx_rm)]
      }
    }
    
    names(dmr_keep) <- names(dmr)
    saveRDS(dmr_keep, paste0(tt, "_bin_IDs.RDS") )

}

}
