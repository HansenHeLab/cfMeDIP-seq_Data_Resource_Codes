set.seed(12345)
library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(NMF)
library(RColorBrewer)
library(caret)

directory <- "Final Data Dec 2024/Insert_size/"
ifelse(!dir.exists(directory),dir.create(directory),FALSE)
setwd(directory)

TCGE_cfDM_samples_only_fixed_Dory <- read.csv("metadata_updated_Dec2024.csv")
insert_matrix_df_for_Althaf <- readRDS("Final Data Dec 2024/Insert_size/Without validation - other plots/insert_matrix_df_for_Althaf_updated.rds")

## Annotate the metadata df
TCGE_cfDM_samples_only_fixed_Dory$sample_id <- paste(TCGE_cfDM_samples_only_fixed_Dory$sequencing_id, "_dedup", sep = "")
TCGE_cfDM_samples_only_fixed_Dory$cancer_type_corrected_updated <- tolower(gsub(" ", "_", TCGE_cfDM_samples_only_fixed_Dory$cancer_type))
TCGE_cfDM_samples_only_fixed_Dory$cancer_type_corrected_updated <- as.character(TCGE_cfDM_samples_only_fixed_Dory$cancer_type_corrected_updated)
TCGE_cfDM_samples_only_fixed_Dory$cancer_type_corrected_updated <- gsub('blood_cancer', 'Aml', TCGE_cfDM_samples_only_fixed_Dory$cancer_type_corrected_updated)
TCGE_cfDM_samples_only_fixed_Dory$cancer_type_corrected_updated <- gsub('normal', 'healthy', TCGE_cfDM_samples_only_fixed_Dory$cancer_type_corrected_updated)
TCGE_cfDM_samples_only_fixed_Dory$cancer_type_title_case <- tools::toTitleCase(TCGE_cfDM_samples_only_fixed_Dory$cancer_type_corrected_updated)
TCGE_cfDM_samples_only_fixed_Dory$cancer_type_title_case  <- gsub("Lfs", "LFS", TCGE_cfDM_samples_only_fixed_Dory$cancer_type_title_case)
TCGE_cfDM_samples_only_fixed_Dory$cancer_type_title_case  <- gsub("Aml", "AML", TCGE_cfDM_samples_only_fixed_Dory$cancer_type_title_case)

mono <- insert_matrix_df_for_Althaf[,66:241]
rownames(mono) <- rownames(insert_matrix_df_for_Althaf)
mono_normalized <- mono/rowSums(mono)

TCGE_cfDM_samples_only_fixed_Dory$is_sample_id <- paste0(TCGE_cfDM_samples_only_fixed_Dory$sequencing_id,"_dedup")

TCGE_cfDM_samples_only_fixed_Dory <- TCGE_cfDM_samples_only_fixed_Dory %>% filter(is_sample_id %in% rownames(mono_normalized))

healthy_TCGE <- TCGE_cfDM_samples_only_fixed_Dory %>% filter(cancer_type == "Normal")
LFS_TCGE <- TCGE_cfDM_samples_only_fixed_Dory %>% filter(cancer_type == "LFS Positive" | cancer_type == "LFS Previvor" | cancer_type == "LFS Survivor")
cancer_TCGE <- TCGE_cfDM_samples_only_fixed_Dory %>% filter(cancer_type != "Normal" & cancer_type != "LFS Positive" & cancer_type != "LFS Previvor" & cancer_type != "LFS Survivor")

LFS_mono_normalized_is <- mono_normalized[LFS_TCGE$is_sample_id,]
healthy_mono_normalized_is <- mono_normalized[healthy_TCGE$is_sample_id,]
cancer_mono_normalized_is <- mono_normalized[cancer_TCGE$is_sample_id,]

healthy_df <- data.frame(ID = rownames(healthy_mono_normalized_is), Classes = healthy_TCGE$project_id)
healthy_samples_training <- createDataPartition(healthy_df$Classes, p = 0.7)

cancer_df <- data.frame(ID = rownames(cancer_mono_normalized_is), Classes = cancer_TCGE$cancer_type)
cancer_samples_training <- createDataPartition(cancer_df$Classes, p = 0.7)

healthy_mono_normalized_is_training <- healthy_mono_normalized_is[healthy_samples_training$Resample1,]
healthy_mono_normalized_is_testing <- healthy_mono_normalized_is[!(rownames(healthy_mono_normalized_is) %in% rownames(healthy_mono_normalized_is_training)),]
healthy_TCGE_training <- healthy_TCGE[healthy_samples_training$Resample1,]
healthy_TCGE_testing <- healthy_TCGE[!(rownames(healthy_mono_normalized_is) %in% rownames(healthy_mono_normalized_is_training)),]

cancer_mono_normalized_is_training <- cancer_mono_normalized_is[cancer_samples_training$Resample1,]
cancer_mono_normalized_is_testing <- cancer_mono_normalized_is[!(rownames(cancer_mono_normalized_is) %in% rownames(cancer_mono_normalized_is_training)),]
cancer_TCGE_training <- cancer_TCGE[cancer_samples_training$Resample1,]
cancer_TCGE_testing <- cancer_TCGE[!(rownames(cancer_mono_normalized_is) %in% rownames(cancer_mono_normalized_is_training)),]

saveRDS(healthy_mono_normalized_is_training, "healthy_mono_normalized_is_training.rds")
saveRDS(healthy_mono_normalized_is_testing, "healthy_mono_normalized_is_testing.rds")
saveRDS(healthy_TCGE_training, "healthy_TCGE_training.rds")
saveRDS(healthy_TCGE_testing, "healthy_TCGE_testing.rds")

saveRDS(cancer_mono_normalized_is_training, "cancer_mono_normalized_is_training.rds")
saveRDS(cancer_mono_normalized_is_testing, "cancer_mono_normalized_is_testing.rds")
saveRDS(cancer_TCGE_training, "cancer_TCGE_training.rds")
saveRDS(cancer_TCGE_testing, "cancer_TCGE_testing.rds")

saveRDS(LFS_mono_normalized_is, "LFS_mono_normalized_is.rds")
saveRDS(LFS_TCGE, "LFS_TCGE.rds")

all_training_is <- rbind(healthy_mono_normalized_is_training,cancer_mono_normalized_is_training)
all_training_TCGE <- rbind(healthy_TCGE_training,cancer_TCGE_training)

all_testing_is <- rbind(healthy_mono_normalized_is_testing,cancer_mono_normalized_is_testing)
all_testing_TCGE <- rbind(healthy_TCGE_testing,cancer_TCGE_testing)

NMF <- nmf(x = all_training_is,rank = 2, nrun = 20)

saveRDS(NMF, "NMF_result.rds")

W <-NMF@fit@W
H <- NMF@fit@H
H_rank1 <- H[1,]/sum(H[1,])
H_rank2 <- H[2,]/sum(H[2,])
H_scale <- as.matrix(data.frame(Rank1 = H_rank2, Rank2 = H_rank1))
colnames(H_scale) <- c("Signature 1 (normal)", "Signature 2 (cancer)")
H_melt <- melt(H_scale)

cbp1 <- c('#33a02c','#1f78b4','#b2df8a','#a6cee3','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#fb6a4a','#b15928','#bdbdbd','#969696','#737373')
mytheme <- theme(plot.title = element_text(hjust = 0.5, size = 13), 
                     axis.line = element_line(colour = "black"),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.border = element_rect(fill = NA),
                     panel.background = element_blank(),
                     legend.position = "bottom",
                     legend.key = element_blank(),
                     legend.title = element_text(size = 12),
                     legend.text = element_text(size = 12),
                     strip.background = element_blank(),
                     strip.text = element_text(size = 13),
                     axis.text = element_text(size = 13),
                     axis.title = element_text(size = 13))

NMF_plot <- ggplot(H_melt, aes(Var1, value)) + geom_line(aes(colour = Var2)) +
  xlab("Fragment Length") + ylab("Relative Frequency") + labs(colour = "NMF Signatures") + coord_cartesian(xlim = c(75,250)) +  
  scale_color_manual(values = cbp1[c(1,2)]) + mytheme 

png(file = "Relative_Frequency_NMF_All_fragments.png",
    width = 500, height = 450, units = "px")
print(NMF_plot)
dev.off()

H_diff <- log2(H_scale[,2]/H_scale[,1])

saveRDS(H_diff, "Log2SignatureRatio.rds")

Vessies_Ref <- read_csv("Final Data Dec 2024/Insert_size/vessies_reference_set.csv", col_names = FALSE)
Vessies_Ref_50_275 <- Vessies_Ref[75:250,]

Compare_NMF_Vessies <- data.frame(NMF_log2_ratio = H_diff, Vessies_Ref_Val = Vessies_Ref_50_275$X1)

colnames(Compare_NMF_Vessies) <- c("Log2(Signature Ratio)", "Vessies Per Fragment\nFragmentation Score")

Compare_NMF_Vessies_melt <- melt(as.matrix(Compare_NMF_Vessies))
Compare_NMF_plot <- ggplot(Compare_NMF_Vessies_melt, aes(Var1, value)) + geom_line(aes(colour = Var2)) +
  xlab("Fragment Length") + ylab("Value") + labs(colour = "Reference") + coord_cartesian(xlim = c(75,250)) +  
  scale_color_manual(values = cbp1[c(2,1)]) + mytheme 

png(file = "Compare_NMF_to_Vessies_Fragment_Score_All_fragments.png",
    width = 500, height = 450, units = "px")
print(Compare_NMF_plot)
dev.off()

cancer_type_col <- c('#33a02c','#1f78b4','#b2df8a','#a6cee3','#fb9a99',
                     
                     '#fdbf6f','#ff7f00')

all_training_is_FS <- t(t(all_training_is)*H_diff)
all_training_FS <- data.frame(Sample = rownames(all_training_is), FS = rowSums(t(t(all_training_is)*H_diff)), W_cancer_sig = W[,1], w_normal_sig = W[,2])
all_training_FS <- cbind(all_training_FS,all_training_TCGE)
norm_median <- median(all_training_FS$FS[all_training_FS$cancer_type == "Normal"])
all_training_FS$cancer_type <- factor(all_training_FS$cancer_type,
                                        
                                        levels = c("Normal", "Brain Cancer", "Lung Cancer", "Prostate Cancer",
                                                   
                                                   "AML", "Eye Cancer",
                                                   
                                                   "Head and Neck Cancer"))

Training_FS_plot <- ggplot(all_training_FS, aes(y = FS, x = cancer_type, fill = cancer_type)) + geom_boxplot(outlier.shape = NA) + 
  geom_jitter(pch = 16, alpha = 0.5) + 
  geom_hline(yintercept = norm_median, linetype = "dashed", size = 0.5) +
  xlab("Cancer Type") + ylab("Patient Level Fragmentation Score") + labs(fill = "Cancer Type") + 
  scale_fill_manual(values = cancer_type_col) + mytheme + theme(axis.text.x = element_text(angle =90))
png(file = "Training_Data_Patient_Level_Fragmentation_Score.png",
    width = 940, height = 630, units = "px")
print(Training_FS_plot)
dev.off()


norm_median <- median(all_training_FS$W_cancer_sig[all_training_FS$cancer_type == "Normal"])
Training_W_Cancer_plot <- ggplot(all_training_FS, aes(y = W_cancer_sig, x = cancer_type, fill = cancer_type)) + geom_boxplot(outlier.shape = NA) + 
  geom_jitter(pch = 16, alpha = 0.5) + 
  geom_hline(yintercept = norm_median, linetype = "dashed", size = 0.5) + xlab("Cancer Type") + ylab("Signature 1 (cancer) Weight") + labs(fill = "Cancer Type") + 
  scale_fill_manual(values = cancer_type_col) + mytheme + theme(axis.text.x = element_text(angle =90))
png(file = "Training_W_Cancer_Signature_Weights.png",
    width = 940, height = 630, units = "px")
print(Training_W_Cancer_plot)
dev.off()


all_testing_FS <- data.frame(Sample = rownames(all_testing_is), FS = rowSums(t(t(all_testing_is)*H_diff)))
all_testing_FS <- cbind(all_testing_FS,all_testing_TCGE)
all_testing_FS$cancer_type <- factor(all_testing_FS$cancer_type,
                                      
                                      levels = c("Normal", "Brain Cancer", "Lung Cancer", "Prostate Cancer",
                                                 
                                                 "AML", "Eye Cancer",
                                                 
                                                 "Head and Neck Cancer"))
norm_median <- median(all_testing_FS$FS[all_testing_FS$cancer_type == "Normal"])

Testing_FS_plot <- ggplot(all_testing_FS, aes(y = FS, x = cancer_type, fill = cancer_type)) + geom_boxplot(outlier.shape = NA) + 
  geom_jitter(pch = 16, alpha = 0.5) + 
  geom_hline(yintercept = norm_median, linetype = "dashed", size = 0.5) + xlab("Cancer Type") + ylab("Patient Level Fragmentation Score") + labs(fill = "Cancer Type") + 
  scale_fill_manual(values = cancer_type_col) + mytheme + theme(axis.text.x = element_text(angle =90))
png(file = "Testing_Data_Patient_Level_Fragmentation_Score.png",
    width = 940, height = 630, units = "px")
print(Testing_FS_plot)
dev.off()

LFS_FS <- data.frame(Sample = rownames(LFS_mono_normalized_is), FS = rowSums(t(t(LFS_mono_normalized_is)*H_diff)))
LFS_FS <- cbind(LFS_FS,LFS_TCGE)
LFS_FS$cancer_type <- factor(LFS_FS$cancer_type,
                             
                             levels = c("LFS Survivor", "LFS Previvor", "LFS Positive"))
LFS_type_col <- c('#bdbdbd','#969696','#737373')
LFS_FS_plot <- ggplot(LFS_FS, aes(y = FS, x = cancer_type, fill = cancer_type)) + geom_boxplot(outlier.shape = NA) + 
  geom_jitter(pch = 16, alpha = 0.5) + 
  geom_hline(yintercept = norm_median, linetype = "dashed", size = 0.5) + xlab("LFS Class") + ylab("Patient Level Fragmentation Score") + labs(fill = "LFS Class") + 
  scale_fill_manual(values = LFS_type_col) + mytheme + theme(axis.text.x = element_text(angle =90))
png(file = "LFS_Data_Patient_Level_Fragmentation_Score.png",
    width = 550, height = 630, units = "px")
print(LFS_FS_plot)
dev.off()

test_FS <- rbind(all_testing_FS,LFS_FS)
test_FS$cancer_type <- factor(test_FS$cancer_type,
                             
                             levels = c("Normal", "Brain Cancer", "Lung Cancer", "Prostate Cancer",
                                        
                                        "AML", "Eye Cancer",
                                        
                                        "Head and Neck Cancer", "LFS Survivor", "LFS Previvor", "LFS Positive"))
all_type_col <- c('#33a02c','#1f78b4','#b2df8a','#a6cee3','#fb9a99',
                  
                  '#fdbf6f','#ff7f00','#bdbdbd','#969696','#737373')
test_FS_plot <- ggplot(test_FS, aes(y = FS, x = cancer_type, fill = cancer_type)) + geom_boxplot(outlier.shape = NA) + 
  geom_jitter(pch = 16, alpha = 0.5) + 
  geom_hline(yintercept = norm_median, linetype = "dashed", size = 0.5) + xlab("Cancer Type") + ylab("Patient Level Fragmentation Score") + labs(fill = "Cancer Type") + 
  scale_fill_manual(values = all_type_col) + mytheme + theme(axis.text.x = element_text(angle =90))
png(file = "All_Testing_Data_Patient_Level_Fragmentation_Score.png",
    width = 1080, height = 630, units = "px")
print(test_FS_plot)
dev.off()


### Added by Dory 
# Find common columns between the two datasets
common_columns <- intersect(names(all_training_FS[, c(1, 2, 5:44)]), names(test_FS))

# Subset both datasets to only include these common columns
aligned_training_FS <- all_training_FS[, common_columns, drop = FALSE]
aligned_test_FS <- test_FS[, common_columns, drop = FALSE]

# Combine the datasets by rows
all_FS <- rbind(aligned_training_FS, aligned_test_FS)


#all_FS <- rbind(all_training_FS[,c(1,2,5:44)],test_FS)

saveRDS(all_FS, "All_samples_FS.rds")

all_FS$cancer_type <- factor(all_FS$cancer_type,
                             
                             levels = c("Normal", "Brain Cancer", "Lung Cancer", "Prostate Cancer",
                                        
                                        "AML", "Eye Cancer",
                                        
                                        "Head and Neck Cancer", "LFS Survivor", "LFS Previvor", "LFS Positive"))
all_type_col <- c('#33a02c','#1f78b4','#b2df8a','#a6cee3','#fb9a99',
                  
                  '#fdbf6f','#ff7f00','#bdbdbd','#969696','#737373')
norm_median <- median(all_FS$FS[all_FS$cancer_type == "Normal"])
all_FS_plot <- ggplot(all_FS, aes(y = FS, x = cancer_type, fill = cancer_type)) + geom_boxplot(outlier.shape = NA) + 
  geom_jitter(pch = 16, alpha = 0.5) + 
  geom_hline(yintercept = norm_median, linetype = "dashed", size = 0.5) + xlab("Cancer Type") + ylab("Patient Level Fragmentation Score") + labs(fill = "Cancer Type") + 
  scale_fill_manual(values = all_type_col) + mytheme + theme(axis.text.x = element_text(angle =90))
png(file = "All_Data_Patient_Level_Fragmentation_Score.png",
    width = 1080, height = 630, units = "px")
print(all_FS_plot)
dev.off()

g <- ggplotGrob(all_FS_plot)$grobs
Cancer_legend_plot <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
dev.off()
lheight <- sum(Cancer_legend_plot$height)


all_plots <- grid.arrange(NMF_plot, 
                          Training_W_Cancer_plot+theme(legend.position = "none",axis.title.x=element_blank(),
                                                       axis.text.x=element_blank(),
                                                       axis.ticks.x=element_blank()), 
                          Compare_NMF_plot,
                          Training_FS_plot+theme(legend.position = "none",axis.title.x=element_blank(),
                                                 axis.text.x=element_blank(),
                                                 axis.ticks.x=element_blank()), 
                          test_FS_plot+theme(legend.position = "none",axis.title.x=element_blank(),
                                             axis.text.x=element_blank(),
                                             axis.ticks.x=element_blank()), 
                          all_FS_plot+theme(legend.position = "none",axis.title.x=element_blank(),
                                            axis.text.x=element_blank(),
                                            axis.ticks.x=element_blank()),
                          Cancer_legend_plot,
                          layout_matrix = rbind(c(1,2),c(3,4),c(5,5),c(6,6),c(7,7)), heights = unit.c((unit(1, "npc") - lheight)/4,
                                                                                                      (unit(1, "npc") - lheight)/4,
                                                                                                      (unit(1, "npc") - lheight)/4,
                                                                                                      (unit(1, "npc") - lheight)/4,
                                                                                                      lheight))
png(file = "All_Plots.png",
    width = 1280, height = 1440, units = "px")
grid.arrange(all_plots)
dev.off()

ggsave(file.path(directory, paste("NMF fragmentation scores, Feb 2024 all frags.pdf", sep = "")), all_plots, width = 11.5, height = 13, dpi = 500, units = "in")
c("Normal" = '#33a02c', "Brain Cancer" = '#1f78b4', "Lung Cancer" = '#b2df8a', "Prostate Cancer" = '#a6cee3',
  
  "AML" = '#fb9a99', "Eye Cancer" = '#fdbf6f',
  
  "Head and Neck\nCancer" = '#ff7f00')

g <- ggplotGrob(Training_FS_plot+theme(legend.position = "bottom"))$grobs
training_legend_plot <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
dev.off()
lt_height <- sum(training_legend_plot$height)
lt_width <- sum(training_legend_plot$width)

new_figure_plots <- grid.arrange(NMF_plot+theme(legend.byrow = F), 
                                 Training_W_Cancer_plot+theme(legend.position = "none",axis.title.x=element_blank(),
                                                              axis.text.x=element_blank(),
                                                              axis.ticks.x=element_blank()), 
                                 Compare_NMF_plot+theme(legend.byrow = F),
                                 Training_FS_plot+theme(legend.position = "none",axis.title.x=element_blank(),
                                                        axis.text.x=element_blank(),
                                                        axis.ticks.x=element_blank()), 
                                 training_legend_plot,
                                 all_FS_plot+theme(legend.position = "none",axis.title.x=element_blank(),
                                                   axis.text.x=element_blank(),
                                                   axis.ticks.x=element_blank()),
                                 Cancer_legend_plot,
                                 layout_matrix = rbind(c(NA,NA,NA,NA,NA,NA,NA,NA,NA),
                                                       c(1,1,1,1,2,2,2,2,2)
                                                       ,c(3,3,3,3,4,4,4,4,4),
                                                       c(NA,NA,NA,NA,5,5,5,5,5),
                                                       c(NA,NA,NA,NA,NA,NA,NA,NA,NA),
                                                       c(6,6,6,6,6,6,6,6,6),
                                                       c(7,7,7,7,7,7,7,7,7)),
                                 heights = unit.c(lheight/2,
                                                  (unit(1, "npc") - (2.5*lheight+lt_height))/4,
                                                  (unit(1, "npc") - (2.5*lheight+lt_height))/4, 
                                                  lt_height,
                                                  lheight,
                                                  (unit(1, "npc") - (2.5*lheight+lt_height))/2,
                                                  lheight))


png(file = "new_figure_plots2.png",
    width = 1280, height = 1440, units = "px")
grid.arrange(new_figure_plots)
dev.off()

pdf(file = "NMF fragmentation scores, May 2024 all frags.pdf",
    width = 13.5, height = 14.4)
grid.arrange(new_figure_plots)
dev.off()

ggsave(file.path(directory, paste("NMF fragmentation scores, Dec 2024 all frags.pdf", sep = "")), new_figure_plots, width = 11.5, height = 13, dpi = 1000, units = "in")
