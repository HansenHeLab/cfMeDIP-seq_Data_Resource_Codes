args = commandArgs(trailingOnly=TRUE)

## assigning commandArgs
sample_file  =  args[1]

### loading libraries ####
#library(genefilter)   ## requied for sva
library(sva)
library(DESeq2)
library(umap)
library(dplyr)
library(ggfortify)   ## PCA plot: autoplot
library(ggrepel)
library(ggplot2)
library(gplots)
library(tools)

name_prefix  =  tools::file_path_sans_ext(basename(sample_file))

s_info  =  read.csv(sample_file)

##########################################################
## function to selected top 10K variable bins based on IQR
## and perform the PCA analysis
## need to omit rows with NA for PCA analysis
##########################################################
iqr_pca <- function(expr)
{
  idx_rc <- colSums(is.na(expr)) == nrow(expr)
  expr <- expr[, !idx_rc]       ## remove samples with NA for all bins

  expr <- na.omit(expr)        ##  remove bins with NA record
  iqr <- apply(expr, 1, IQR)
  names(iqr) <- rownames(expr)

  idx_cons <- iqr == 0      ## remove constant bin
  iqr <- iqr[!idx_cons]

  iqr_10K <- sort(iqr, decreasing = T)[1: 10000]

  idx_ff <- match(rownames(expr), names(iqr_10K))

  ## most variable 10K bins for the PCA analysis
  expr_f <- expr[!is.na(idx_ff), ]
  expr_f_pca <- prcomp(t(expr_f), center = T, scale.=T)
  return(expr_f_pca)
}

######################
## merge all data sets
#####################

if(FALSE){
  files <- list.files(path = ".", pattern = "txt.gz")

  #############################
  ## merge samples from studies
  ## chek.names = F to ensure id start with nubmers or with hypen to be reamined
  print(files[1])
  aggr <- read.table(gzfile(files[1]), header = T, comment.char = "*", check.names = F)

  for (i in 2: length(files))
  {
    print(files[i])
    tmp <- read.table(gzfile(files[i]), header = T, comment.char = "*", check.names =F)
    tmp <- tmp[, -c(1:6)]
    aggr <- cbind(aggr, tmp)
  }

  ### extract bin
  idx_b <- aggr[, 5] >= 1         ## require the bin contain at least 1 CpG sites
  saveRDS(aggr[idx_b, ], file = "ALL_auto_bfilt_wCpG.RDS")

  write.table(aggr[idx_b, 1:6], "ALL_auto_bfilt_wCpG_bins.bed", row.names = F,
              col.names = F, quote = F, sep = "\t")

}

#########################################
## pull out samples from provided s_info
## and perform PCA analysis
#########################################
{
print("Loading pre-aggregated data ....")

aggr <- readRDS("ALL_auto_bfilt_wCpG.RDS")

## split bin and samples
bin_info <- aggr[, 1:6]
# write.table(bin_info, "ALL_auto_bfilt_wCpG_bins.bed", row.names = F, col.names = F, quote = F, sep = "\t")

quant <- data.matrix(aggr[, -(1:6)])
rownames(quant) <- bin_info$bin_id

## quantified samples info
## cnt might contain samples were excluded in meta table
idx_s <- match(colnames(quant), s_info$sequencing_id)
quant_s <- quant[, !is.na(idx_s)]
s_info_s <- s_info[idx_s[!is.na(idx_s)], ]
rm(quant)
saveRDS(quant_s, file = paste0(name_prefix, "_auto_bfilt_wCpG.RDS"))

print("Extracted quants saved !!!")



## run PCA analysis
print("Performing the PCA analysis ...")
quant_s10k <- iqr_pca(quant_s)
saveRDS(quant_s10k, file = paste0(name_prefix, "_auto_bfilt_wCpG_Top10Kbin_PCA.RDS"))

print("PCA Done !!")
}
