rm(list = ls())
setwd("/Users/yong/OneDrive - UHN/Projects/TCGE/cfEpigenomics/Resource/0_sample_meta_info/meta_info_summary_for_revision")

library(ggplot2)
library(dplyr)
library(microViz)
library(RColorBrewer)

#s_info <- read.csv("../TCGE_CFMe_Samples_MetaInfo_Master_Table_light_fixed_202410.csv")
## fill with NA for the missing elements
s_info_all <- read.csv("../TCGE_CFMe_Samples_MetaInfo_Master_Table_202411_all.csv",
                       na.strings=c("NA","NaN", ""))

col_rm <- c("cfDNA_concentration_ngPml", "antibody",
            "sequencing_platform", "Input_DNA_ng")

## healthy PBL sample 
s_info_pbl <- s_info_all %>% 
  filter(sample_type == 18) %>% 
  dplyr::select(!col_rm)
write.csv(s_info_pbl, "TCGE-CFMe-HNSC_PBL_Healthy.csv", row.names = F)  

## cell free only 
s_info_cf <- s_info_all %>% 
  filter(Analyte == "cfDM") %>% 
  dplyr::select(!col_rm)

length(unique(s_info_cf$project_id))   # 11
length(unique(s_info_cf$cancer_type))  # 19

### validation only 
idx_val <- match(s_info_cf$project_id, c("TCGE-CFMe-INSPIRE", "TCGE-CFMe-HCC"))
s_info_val <- s_info_cf[!is.na(idx_val), ]

write.csv(s_info_val, "TCGE_cfDM_samples_only_fixed_validation.csv", row.names = F)  

##############################
## originally included samples
##############################
{
  s_info <- s_info_cf[is.na(idx_val), ] 
  write.csv(s_info, file = "TCGE_cfDM_samples_only_fixed.csv", row.names = F)
  
  ## order by samples types: rev for coord_flip()
  s_info$project_id <- factor(s_info$project_id, 
                              levels = rev(c("TCGE-CFMe-MCA", "TCGE-CFMe-BCA", "TCGE-CFMe-HNSC", "TCGE-CFMe-PRAD",
                                             "TCGE-CFMe-AML", "TCGE-CFMe-SCLC", "TCGE-CFMe-UM", 
                                             "TCGE-CFMe-HBC", "TCGE-CFMe-LFS")))
  
  project_col <- rev(c('#7fcdbb',  '#984ea3', '#ff7f00',  '#377eb8', '#f781bf', '#807dba','#a65628', '#4daf4a', '#999999'))
  
  
  ## since 20240624 Normal to Healthy; blood cancer --> AML;  Eye cancer -> Uveal Melanoma in the TCGE_CFMe_Samples_MetaInfo_Master_Table_light_fixed_20231018.csv
  
  s_info$cancer_type <- factor(s_info$cancer_type, 
                               levels = rev(c("Healthy", "Brain Cancer", "Lung Cancer", "Prostate Cancer",
                                              "AML", "Pancreatic Cancer", "Uveal Melanoma",
                                              "Head and Neck Cancer", "Breast Cancer", "Colorectal Cancer",
                                              "Bladder Cancer", "Renal Cancer",
                                              "LFS Survivor", "LFS Previvor", "LFS Positive")))
  
  
  cancer_type_col <- rev(c('#33a02c','#1f78b4','#b2df8a','#a6cee3','#fb9a99',
                           '#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a',
                           '#fb6a4a','#b15928','#bdbdbd','#969696','#737373'))
  
  
  ## all paired-end samples
  idx_mca <- s_info$project_id == "TCGE-CFMe-MCA"
  write.csv(s_info[idx_mca, ], file = "TCGE_cfDM_samples_only_fixed_SE.csv",  row.names = F)
  write.csv(s_info[!idx_mca, ], file = "TCGE_cfDM_samples_only_fixed_PE.csv",  row.names = F)
  
  ## cf sample health control only
  idx_normal <- s_info$cancer_type == "Healthy"
  write.csv(s_info[idx_normal, ], file = "TCGE_cfDM_samples_only_fixed_Normal.csv",  row.names = F)
  write.csv(s_info[idx_normal & !idx_mca, ], file = "TCGE_cfDM_samples_only_fixed_Normal_PE.csv",  row.names = F)
  write.csv(s_info[idx_normal & idx_mca, ], file = "TCGE_cfDM_samples_only_fixed_Normal_SE.csv",  row.names = F)
  
  
  idx_aml <- s_info$project_id == "TCGE-CFMe-AML"
  write.csv(s_info[idx_aml, ], file = "TCGE_cfDM_samples_only_fixed_AML.csv",  row.names = F)
  
  
  ## cell-free samples per cancer type 
  dat <- s_info %>% count(cancer_type) 
  g <- ggplot(data = dat, aes(x = cancer_type, y = n, fill = cancer_type)) + geom_bar(stat="identity")
  g <- g + geom_text(aes(label = n), size = 3.5) +  coord_flip() + theme_classic()
  g <- g + scale_fill_manual(values = cancer_type_col)
  ggsave("cf_samples_per_cancer_type.pdf", width = 5, height = 7)
  
  ## cell-free samples per cancer subtype 
  dat <- s_info %>% count(cancer_subtype) 
  g <- ggplot(data = dat, aes(x = cancer_subtype, y = n)) + geom_bar(stat="identity", fill = "steelblue")
  g <- g + geom_text(aes(label = n), size = 3.5) +  coord_flip() + theme_classic()
  ggsave("cf_samples_per_cancer_subtype.pdf", width = 7, height = 7)
  
  
  ## cell-free samples per project per cancer type
  dat <- s_info %>% 
    group_by(project_id) %>% 
    count(cancer_type)
  write.csv(dat, "cf_samples_per_project_per_cancer_type.csv", row.names = F)
  
  ## group by project
  g <- ggplot(data = dat, aes(x = project_id, y = n, fill = cancer_type)) 
  g <- g + geom_bar(stat="identity", position=position_dodge())
  g <- g + geom_text(aes(label = n), position = position_dodge(0.9), size=3)
  g <- g + coord_flip() + theme_classic()
  g <- g + scale_fill_manual(values = cancer_type_col)
  ggsave("cf_samples_per_project_per_cancer_type.pdf", width = 7, height = 7)
  
  ## group by cancer types
  g <- ggplot(data = dat, aes(x = cancer_type, y = n, fill = project_id)) 
  g <- g + geom_bar(stat="identity", position=position_dodge())
  g <- g + geom_text(aes(label = n), position = position_dodge(0.9), size=3)
  g <- g + labs(y = "Number of Samples", x = "") + coord_flip() + theme_classic()
  g <- g + scale_fill_manual(values = project_col)
  ggsave("cf_samples_per_cancer_per_project_type.pdf", width = 6, height = 6)
  
  ## cancer_subtype group by cancer types
  dat <- s_info %>% 
    group_by(cancer_type) %>% 
    count(cancer_subtype) 
  
  g <- ggplot(data = dat, aes(x = cancer_type, y = n, fill = cancer_subtype)) 
  g <- g + geom_bar(stat="identity", position=position_dodge())
  g <- g + geom_text(aes(label = n), position = position_dodge(0.9), size=3)
  g <- g + coord_flip() + theme_classic()
  ggsave("cf_samples_cancer_subtype_per_cancer.pdf", width = 12, height = 14)
  
  
  ## summary table per project and cancer type 
  s_info_tmp <- s_info
  s_info_tmp$project_id <- factor(s_info_tmp$project_id, 
                                  levels = c("TCGE-CFMe-MCA", "TCGE-CFMe-BCA", "TCGE-CFMe-HNSC", "TCGE-CFMe-PRAD",
                                             "TCGE-CFMe-AML", "TCGE-CFMe-SCLC", "TCGE-CFMe-UM", 
                                             "TCGE-CFMe-HBC", "TCGE-CFMe-LFS"))
  
  cnt_sex <- s_info_tmp %>% 
    group_by(project_id, cancer_type) %>% 
    count(sex) 
  
  cnt_stage <- s_info_tmp %>% mutate(na = !is.na(cancer_grade_or_stage)) %>% 
    group_by(project_id, cancer_type) %>% 
    summarise(sum(na))
  
  cnt_meta <- s_info_tmp  %>% mutate(na = !is.na(metastasis)) %>%  
    group_by(project_id, cancer_type) %>% 
    summarise(sum(na))
  
  cnt_relapse <- s_info_tmp  %>% mutate(na = !is.na(relapse)) %>%  
    group_by(project_id, cancer_type) %>% 
    summarise(sum(na))
  
  sample_cnt <- s_info_tmp %>% group_by(project_id, cancer_type) %>% dplyr::count()
  
  cnt_age <- s_info_tmp %>% mutate(age_na = !is.na(age)) %>% group_by(project_id, cancer_type) %>% summarise(sum(age_na, na.rm = T))
  age_mean <- s_info_tmp %>% group_by(project_id, cancer_type) %>% summarise(mean(age, na.rm = T))
  age_sd <- s_info_tmp %>% group_by(project_id, cancer_type) %>% summarise(sd(age, na.rm = T))
  
  
  
  
  
  
  ##################
  ## unique patients 
  ##################
  s_info  <- s_info %>% 
    mutate(TSS_Participant = paste(TSS, Participant, sep = "_")) 
  
  
  tt <- paste0("# unique participants : ", length(unique(s_info$TSS_Participant)))
  
  pdf("cf_unique_participants_frequncey.pdf", width = 4, height = 3)
  hist(table(s_info$TSS_Participant), main = tt, col = "steelblue", xlab = "Samples per participant")
  dev.off()
  
  ## by ggplot
  cnt <- data.frame(table(s_info$TSS_Participant))
  dat_cnt <- cnt %>% 
    group_by(Freq) %>% 
    count() 
  dat_cnt$Freq <- as.factor(dat_cnt$Freq)
  
  ## unique participants
  length(cnt[, 2])
  sum(cnt[, 2] > 1)
  
  
  g <- ggplot(data = dat_cnt, aes(x = Freq, y = n, fill = Freq)) 
  g <- g + geom_bar(stat="identity", position=position_dodge())
  g <- g + geom_text(aes(label = n), position = position_dodge(0.9), size=3, vjust =  -0.25)
  g <- g + labs(y = "Number of Participants", x  = "Samples per participant")
  g <- g + scale_fill_manual(values = rev(c('#f7fcfd','#e0ecf4','#bfd3e6','#9ebcda','#8c96c6','#8c6bb1','#88419d','#810f7c')))
  g <- g + theme_classic() + theme(legend.position = "none")
  
  ggsave("cf_unique_participants_frequncey_pub.pdf", width = 4, height = 3)
  
  
  ## summary of sample tumor grad or stage info
  sum(s_info$cancer_type == "Healthy")
  
  sum(s_info$cancer_type != "Healthy")
  sum(!is.na(s_info$cancer_grade_or_stage)) 
  sum(!is.na(s_info$relapse))
  sum(!is.na(s_info$metastasis))
  
  
  ################################
  ## remove duplicated participant
  s_info_u <- s_info %>% filter(vial == "A") %>% 
    dplyr::select(!c("TSS_Participant"))
  write.csv(s_info_u, file = "TCGE_cfDM_samples_only_fixed_vialA_only.csv", row.names = F)
  
  ## cell-free samples sex
  dat_ug <- s_info_u %>% count(sex) 
  dat_ug$sex[3] <- "No record"
  dat_ug$sex <- factor(dat_ug$sex, levels = c("Male", "Female", "No record"))
  
  g <- ggplot(data = dat_ug, aes(x = sex, y = n, fill = sex)) + geom_bar(stat="identity")
  g <- g + geom_text(aes(label = n), size = 3.25, vjust =  -0.25) 
  g <- g + labs(y = "Number of Participants") 
  g <- g + scale_fill_manual(values = c('#0868ac', '#df65b0', '#969696'))
  g <- g + theme_classic() + theme(legend.position = "none")
  ggsave("cf_unique_participants_per_sex.pdf", width = 2.5, height = 3)
  
  ###########
  ## pie plot
  library(scales)
  
  ## balnk_theme
  blank_theme <- theme_minimal()+
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.border = element_blank(),
      panel.grid=element_blank(),
      axis.ticks = element_blank(),
      plot.title=element_text(size=14, face="bold")
    )
  
  g <- ggplot(dat_ug,  aes(x = "", y = n, fill = sex)) + geom_bar(width = 1, stat = "identity")
  g <- g + coord_polar("y") + blank_theme  + theme(axis.text.x=element_blank()) 
  g <- g + scale_fill_manual(values = c('#0868ac', '#df65b0', '#969696'))
  #g <- g + geom_text(aes(y = n/3 + c(0, cumsum(n)[-length(n)]), label = n), size=5)   ##will resort 
  g <- g + theme(legend.position = "top")
  ggsave("cf_unique_participants_per_sex_Pie_chart.pdf", width = 3, height = 3)
  
  
  #############################
  ## samples ages distribution 
  ## participate with multiple samples
  
  summary(s_info_u$age)
  age_cdf <- ecdf(s_info_u$age)
  tt <- paste0("Samples N = ", length(s_info_u$age), "; Age NA: ", sum(is.na(s_info_u$age)))
  
  pdf("cf_samples_age_CDF.pdf", width = 4, height = 3)
  par(mar = c(4, 3, 1, 1), mgp = c(2, 1, 0))
  plot(age_cdf, sub = tt, main = "", xlab = "Age")
  abline(h = 0.2, lty = 2, col = "red")
  dev.off()
  
  age <- s_info_u$age
  age_group <- vector()
  age_group[age <= 20] <- "(0, 20]"
  age_group[age > 20 & age <= 40] <- "(20, 40]"
  age_group[age > 40 & age <= 60] <- "(40, 60]"
  age_group[age > 60 & age <= 80] <- "(60, 80]"
  age_group[age > 80 & age <= 100] <- "(80, 100]"
  age_group[is.na(age_group)] <- "No record"
  
  dat_age <- data.frame(age, age_group)
  dat_age_cnt <- dat_age %>% 
    group_by(age_group) %>% 
    count() 
  
  g <- ggplot(data = dat_age_cnt, aes(x = age_group, y = n, fill = age_group)) 
  g <- g + geom_bar(stat="identity", position=position_dodge())
  g <- g + geom_text(aes(label = n), position = position_dodge(0.9), size=3, vjust =  -0.25)
  g <- g + labs(y = "Number of Participants", x  = "Age Ranges")
  g <- g + scale_fill_manual(values = c('#ccebc5','#a8ddb5','#7bccc4','#4eb3d3','#2b8cbe','#969696'))
  g <- g + theme_classic() + theme(legend.position = "none")
  ggsave("cf_samples_age_frequncey_pub.pdf", width = 4, height = 3)
  
  ## summary of sample tumor grad or stage info unisamples
  sum(s_info_u$cancer_type == "Healthy")
  
  sum(s_info_u$cancer_type != "Healthy")
  sum(!is.na(s_info_u$cancer_grade_or_stage)) 
  sum(!is.na(s_info_u$relapse))
  sum(!is.na(s_info_u$metastasis))
  
  
  
}

######################
## validation samples
######################
{
  s_info <- s_info_val
  
  ## order by samples types: rev for coord_flip()
  s_info$project_id <- factor(s_info$project_id, 
                              levels = c("TCGE-CFMe-HCC", "TCGE-CFMe-INSPIRE"))
  
  project_col <- c('#7fcdbb',  '#984ea3')
  
  
  ## since 20240624 Normal to Healthy; blood cancer --> AML;  Eye cancer -> Uveal Melanoma in the TCGE_CFMe_Samples_MetaInfo_Master_Table_light_fixed_20231018.csv
  
  s_info$cancer_type <- factor(s_info$cancer_type, 
                               levels = rev(c("Healthy", "Melanoma","Head and Neck Cancer", "Breast Cancer", 
                                              "Ovarian Cancer", "Liver Cancer", "Mixed Cancer")))
  
  
  
  cancer_type_col <- rev(c('#33a02c','#1f78b4','#b2df8a','#a6cee3','#fb9a99',
                           '#e31a1c','#fdbf6f','#ff7f00','#cab2d6'))
  
  
  ## cell-free samples per cancer type 
  dat <- s_info %>% count(cancer_type) 
  g <- ggplot(data = dat, aes(x = cancer_type, y = n, fill = cancer_type)) + geom_bar(stat="identity")
  g <- g + geom_text(aes(label = n), size = 3.5) +  coord_flip() + theme_classic()
  g <- g + scale_fill_manual(values = cancer_type_col)
  ggsave("Validation_cf_samples_per_cancer_type.pdf", width = 5, height = 7)
  
  ## cell-free samples per cancer subtype 
  dat <- s_info %>% count(cancer_subtype) 
  g <- ggplot(data = dat, aes(x = cancer_subtype, y = n)) + geom_bar(stat="identity", fill = "steelblue")
  g <- g + geom_text(aes(label = n), size = 3.5) +  coord_flip() + theme_classic()
  ggsave("Validation_cf_samples_per_cancer_subtype.pdf", width = 7, height = 7)
  
  
  ## cell-free samples per project per cancer type
  dat <- s_info %>% 
    group_by(project_id) %>% 
    count(cancer_type)
  write.csv(dat, "Validation_cf_samples_per_project_per_cancer_type.csv", row.names = F)
  
  ## group by project
  g <- ggplot(data = dat, aes(x = project_id, y = n, fill = cancer_type)) 
  g <- g + geom_bar(stat="identity", position=position_dodge())
  g <- g + geom_text(aes(label = n), position = position_dodge(0.9), size=3)
  g <- g + coord_flip() + theme_classic()
  g <- g + scale_fill_manual(values = cancer_type_col)
  ggsave("Validation_cf_samples_per_project_per_cancer_type.pdf", width = 7, height = 7)
  
  ## group by cancer types
  g <- ggplot(data = dat, aes(x = cancer_type, y = n, fill = project_id)) 
  g <- g + geom_bar(stat="identity", position=position_dodge())
  g <- g + geom_text(aes(label = n), position = position_dodge(0.9), size=3)
  g <- g + labs(y = "Number of Samples", x = "") + coord_flip() + theme_classic()
  g <- g + scale_fill_manual(values = project_col)
  ggsave("Validation_cf_samples_per_cancer_per_project_type.pdf", width = 6, height = 6)
  
  ## cancer_subtype group by cancer types
  dat <- s_info %>% 
    group_by(cancer_type) %>% 
    count(cancer_subtype) 
  
  g <- ggplot(data = dat, aes(x = cancer_type, y = n, fill = cancer_subtype)) 
  g <- g + geom_bar(stat="identity", position=position_dodge())
  g <- g + geom_text(aes(label = n), position = position_dodge(0.9), size=3)
  g <- g + coord_flip() + theme_classic()
  ggsave("Validation_cf_samples_cancer_subtype_per_cancer.pdf", width = 12, height = 14)
  
  
  ##################
  ## unique patients 
  ##################
  
  
  s_info  <- s_info %>% 
    mutate(TSS_Participant = paste(TSS, Participant, sep = "_")) 
  
  #a <- unique(s_info$TSS_Participant)
  #b <- unique(s_info$TSS_Participant[s_info$vial == "A"])
  #idx <- match(a, b)  
  #a[is.na(idx)]
  ## "10_4"  "10_10" "10_65" "10_68"  samples are not at baseline time point 
  
  tt <- paste0("# unique participants : ", length(unique(s_info$TSS_Participant)))
  
  pdf("Validation_cf_unique_participants_frequncey.pdf", width = 4, height = 3)
  hist(table(s_info$TSS_Participant), main = tt, col = "steelblue", xlab = "Samples per participant")
  dev.off()
  
  ## by ggplot
  cnt <- data.frame(table(s_info$TSS_Participant))
  dat_cnt <- cnt %>% 
    group_by(Freq) %>% 
    count() 
  dat_cnt$Freq <- as.factor(dat_cnt$Freq)
  
  ## unique participants
  length(cnt[, 2])
  sum(cnt[, 2] > 1)
  
  
  g <- ggplot(data = dat_cnt, aes(x = Freq, y = n, fill = Freq)) 
  g <- g + geom_bar(stat="identity", position=position_dodge())
  g <- g + geom_text(aes(label = n), position = position_dodge(0.9), size=3, vjust =  -0.25)
  g <- g + labs(y = "Number of Participants", x  = "Samples per participant")
  g <- g + scale_fill_manual(values = rev(c('#f7fcfd', '#f7fcfd','#e0ecf4','#bfd3e6','#9ebcda','#8c96c6','#8c6bb1','#88419d','#810f7c' )))
  g <- g + theme_classic() + theme(legend.position = "none")
  
  ggsave("Validation_cf_unique_participants_frequncey_pub.pdf", width = 4, height = 3)
  
  
  
  ################################
  ## remove duplicated participant
  s_info_u <- s_info %>% filter(vial == "A") %>% 
    select(!c("TSS_Participant"))
  write.csv(s_info_u, file = "TCGE_cfDM_samples_only_fixed_validation_vialA_only.csv", row.names = F)
  
  ## cell-free samples sex
  dat_ug <- s_info_u %>% count(sex) 
  dat_ug$sex <- factor(dat_ug$sex, levels = c("Male", "Female"))
  
  g <- ggplot(data = dat_ug, aes(x = sex, y = n, fill = sex)) + geom_bar(stat="identity")
  g <- g + geom_text(aes(label = n), size = 3.25, vjust =  -0.25) 
  g <- g + labs(y = "Number of Participants") 
  g <- g + scale_fill_manual(values = c('#0868ac', '#df65b0', '#969696'))
  g <- g + theme_classic() + theme(legend.position = "none")
  ggsave("Validation_cf_unique_participants_per_sex.pdf", width = 2.5, height = 3)
  
  ###########
  ## pie plot
  library(scales)
  
  ## balnk_theme
  blank_theme <- theme_minimal()+
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.border = element_blank(),
      panel.grid=element_blank(),
      axis.ticks = element_blank(),
      plot.title=element_text(size=14, face="bold")
    )
  
  g <- ggplot(dat_ug,  aes(x = "", y = n, fill = sex)) + geom_bar(width = 1, stat = "identity")
  g <- g + coord_polar("y") + blank_theme  + theme(axis.text.x=element_blank()) 
  g <- g + scale_fill_manual(values = c('#0868ac', '#df65b0', '#969696'))
  #g <- g + geom_text(aes(y = n/3 + c(0, cumsum(n)[-length(n)]), label = n), size=5)   ##will resort 
  g <- g + theme(legend.position = "top")
  ggsave("Validation_cf_unique_participants_per_sex_Pie_chart.pdf", width = 3, height = 3)
  
  
  #############################
  ## samples ages distribution 
  ## participate with multiple samples
  
  summary(s_info_u$age)
  age_cdf <- ecdf(s_info_u$age)
  tt <- paste0("Samples N = ", length(s_info$age), "; Age NA: ", sum(is.na(s_info_u$age)))
  
  pdf("Validation_cf_samples_age_CDF.pdf", width = 4, height = 3)
  par(mar = c(4, 3, 1, 1), mgp = c(2, 1, 0))
  plot(age_cdf, sub = tt, main = "", xlab = "Age")
  abline(h = 0.2, lty = 2, col = "red")
  dev.off()
  
  age <- s_info_u$age
  age_group <- vector()
  age_group[age <= 20] <- "(0, 20]"
  age_group[age > 20 & age <= 40] <- "(20, 40]"
  age_group[age > 40 & age <= 60] <- "(40, 60]"
  age_group[age > 60 & age <= 80] <- "(60, 80]"
  age_group[age > 80 & age <= 100] <- "(80, 100]"
  age_group[is.na(age_group)] <- "No record"
  
  dat_age <- data.frame(age, age_group)
  dat_age_cnt <- dat_age %>% 
    group_by(age_group) %>% 
    count() 
  
  g <- ggplot(data = dat_age_cnt, aes(x = age_group, y = n, fill = age_group)) 
  g <- g + geom_bar(stat="identity", position=position_dodge())
  g <- g + geom_text(aes(label = n), position = position_dodge(0.9), size=3, vjust =  -0.25)
  g <- g + labs(y = "Number of Participants", x  = "Age Ranges")
  g <- g + scale_fill_manual(values = c('#ccebc5','#a8ddb5','#7bccc4','#4eb3d3','#2b8cbe','#969696'))
  g <- g + theme_classic() + theme(legend.position = "none")
  ggsave("Validation_cf_samples_age_frequncey_pub.pdf", width = 4, height = 3)
  
}


