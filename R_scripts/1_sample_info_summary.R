rm(list = ls())
setwd("/Users/yong/OneDrive - UHN/Projects/TCGE/cfEpigenomics/Resource/0_sample_meta_info/meta_info_summary_for_fixed")

library(ggplot2)
library(dplyr)
library(microViz)
library(RColorBrewer)

s_info <- read.csv("../TCGE_CFMe_Samples_MetaInfo_Master_Table_light_fixed_20231018.csv")

length(unique(s_info$project_id))   # 9
length(unique(s_info$cancer_type))  # 15

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


##########################
## cell-free samples only
##########################
{
## cf sample only
dat_cf <- s_info %>% filter(Analyte == "cfDM")
write.csv(dat_cf, file = "TCGE_cfDM_samples_only_fixed.csv", row.names = F)

## all paired-end samples
idx_mca <- dat_cf$project_id == "TCGE-CFMe-MCA"
write.csv(dat_cf[idx_mca, ], file = "TCGE_cfDM_samples_only_fixed_SE.csv",  row.names = F)
write.csv(dat_cf[!idx_mca, ], file = "TCGE_cfDM_samples_only_fixed_PE.csv",  row.names = F)

## cf sample health control only
idx_normal <- dat_cf$cancer_type == "Normal"
write.csv(dat_cf[idx_normal, ], file = "TCGE_cfDM_samples_only_fixed_Normal.csv",  row.names = F)
write.csv(dat_cf[idx_normal & !idx_mca, ], file = "TCGE_cfDM_samples_only_fixed_Normal_PE.csv",  row.names = F)
write.csv(dat_cf[idx_normal & idx_mca, ], file = "TCGE_cfDM_samples_only_fixed_Normal_SE.csv",  row.names = F)


idx_aml <- dat_cf$project_id == "TCGE-CFMe-AML"
write.csv(dat_cf[idx_aml, ], file = "TCGE_cfDM_samples_only_fixed_AML.csv",  row.names = F)

## cell-free samples per cancer type 
dat <- dat_cf %>% count(cancer_type) 
g <- ggplot(data = dat, aes(x = cancer_type, y = n, fill = cancer_type)) + geom_bar(stat="identity")
g <- g + geom_text(aes(label = n), size = 3.5) +  coord_flip() + theme_classic()
g <- g + scale_fill_manual(values = cancer_type_col)
ggsave("cf_samples_per_cancer_type.pdf", width = 5, height = 7)

## cell-free samples per cancer subtype 
dat <- dat_cf %>% count(cancer_subtype) 
g <- ggplot(data = dat, aes(x = cancer_subtype, y = n)) + geom_bar(stat="identity", fill = "steelblue")
g <- g + geom_text(aes(label = n), size = 3.5) +  coord_flip() + theme_classic()
ggsave("cf_samples_per_cancer_subtype.pdf", width = 7, height = 7)


## cell-free samples per project per cancer type
dat <- dat_cf %>% 
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
dat <- dat_cf %>% 
  group_by(cancer_type) %>% 
  count(cancer_subtype) 

g <- ggplot(data = dat, aes(x = cancer_type, y = n, fill = cancer_subtype)) 
g <- g + geom_bar(stat="identity", position=position_dodge())
g <- g + geom_text(aes(label = n), position = position_dodge(0.9), size=3)
g <- g + coord_flip() + theme_classic()
ggsave("cf_samples_cancer_subtype_per_cancer.pdf", width = 12, height = 14)



##################
## unique patients 
##################
dat_cf <- dat_cf %>% 
  mutate(TSS_Participant = paste(TSS, Participant, sep = "_")) 

tt <- paste0("# unique participants : ", length(unique(dat$TSS_Participant)))
pdf("cf_unique_participants_frequncey.pdf", width = 4, height = 3)
hist(table(dat_cf$TSS_Participant), main = tt, col = "steelblue", xlab = "Samples per participant")
dev.off()

## by ggplot
cnt <- data.frame(table(dat_cf$TSS_Participant))
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

#############################
## samples ages distribution 
## participate with multiple samples

summary(dat_cf$age)
age_cdf <- ecdf(dat_cf$age)
tt <- paste0("Samples N = ", length(dat_cf$age), "; Age NA: ", sum(is.na(dat_cf$age)))

pdf("cf_samples_age_CDF.pdf", width = 4, height = 3)
par(mar = c(4, 3, 1, 1), mgp = c(2, 1, 0))
plot(age_cdf, sub = tt, main = "", xlab = "Age")
abline(h = 0.2, lty = 2, col = "red")
dev.off()

age <- dat_cf$age
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
g <- g + labs(y = "Number of samples", x  = "Age Ranges")
g <- g + scale_fill_manual(values = c('#ccebc5','#a8ddb5','#7bccc4','#4eb3d3','#2b8cbe','#969696'))
g <- g + theme_classic() + theme(legend.position = "none")
ggsave("cf_samples_age_frequncey_pub.pdf", width = 4, height = 3)


################################
## remove duplicated participant
dat_cfu <- dat_cf %>% distinct(TSS_Participant, .keep_all = TRUE)

## cell-free samples sex
dat_ug <- dat_cfu %>% count(sex) 
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


}

