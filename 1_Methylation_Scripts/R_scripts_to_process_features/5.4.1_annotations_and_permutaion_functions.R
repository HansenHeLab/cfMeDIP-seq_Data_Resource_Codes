######################################
## annotation regions using "annotatr"
#####################################
if(FALSE){
  library(annotatr)
  
  ## hg38
  {
  ## regions of interest
  annot_regions <- c("hg38_cpg_islands", "hg38_cpg_shores", "hg38_cpg_shelves", "hg38_cpg_inter", 
                     "hg38_genes_promoters", "hg38_enhancers_fantom")
  
  annot_regions_names <- c("CpG island", "CpG shore", "CpG shelf", "CpG sea", 
                     "Promoter", "Enhancer")
  
  annot_regions_list <- list()
  
  
  for(i in 1:length(annot_regions))
  {
    ## retrieving annotations
    annot_regions_list[[i]] <- build_annotations(genome = 'hg38', annotations = annot_regions[i])
  }
  
  names(annot_regions_list) <- annot_regions_names
  saveRDS(annot_regions_list, "List_of_annotated_regions_for_CpGs_Promoter_Enhancer_hg38.RDS")
  }
  
  ## hg19
  {
    ## regions of interest
    annot_regions <- c("hg19_cpg_islands", "hg19_cpg_shores", "hg19_cpg_shelves", "hg19_cpg_inter", 
                       "hg19_genes_promoters", "hg19_enhancers_fantom")
    
    annot_regions_names <- c("CpG island", "CpG shore", "CpG shelf", "CpG sea", 
                             "Promoter", "Enhancer")
    
    annot_regions_list <- list()
    
    
    for(i in 1:length(annot_regions))
    {
      ## retrieving annotations
      annot_regions_list[[i]] <- build_annotations(genome = 'hg19', annotations = annot_regions[i])
    }
    
    names(annot_regions_list) <- annot_regions_names
    saveRDS(annot_regions_list, "List_of_annotated_regions_for_CpGs_Promoter_Enhancer_hg19.RDS")
  }
  
}



################
## ggplot2 theme
################
{
p_theme<- theme_bw()+theme(
  text=element_text(family="Helvetica"),
  axis.title.x =element_text(color="black",size=14,family="Helvetica") ,
  axis.title.y =element_text(color="black",size=14,family="Helvetica") ,
  axis.text.x =element_text(color="black",size=12,family="Helvetica") ,
  axis.text.y =element_text(color="black",size=12,family="Helvetica") ,
  legend.text =element_text(color="black",size=10,family="Helvetica"),
  legend.title=element_text(color="black",size=12,family="Helvetica"),
  legend.background = element_blank(),
  panel.border = element_blank(),panel.grid.major = element_blank(),
  panel.grid.minor =element_blank(),axis.line=element_line(colour = "black",size=0.4))
}

####################
## permutation test 
###################
{
library(regioneR)

permutation <- function(special=NULL,back.region=NULL,annotations.list=NULL,annotation.region=NULL){
  Pvalue <- list(NULL)
  for(annote.type in annotation.region ){
    pt.100 <- permTest(A= special, ntimes=1000, randomize.function=resampleRegions,
                       evaluate.function=numOverlaps, B=annotations.list[[annote.type]], verbose=TRUE,
                       count.once=TRUE, per.chromosome=FALSE, non.overlapping=FALSE,universe=back.region)
    
    Pvalue[[annote.type]] <- pt.100
    print("#=====================")
    print(annote.type);
    print(pt.100$numOverlaps$pval)
  }
  Pvalue[[1]] <- NULL
  return(Pvalue)
}

}

###########################
## plot permutation results
###########################
{
library(gridExtra)
library(grid)
  
plot_permutation <- function(Pvalue = NULL, plot_name){
  #Pvalue <- pvalue.dmrs.hyper
  pvalthres <- 0.05
  p.value.table <- data.frame(NULL)
  for (i in 1:length(Pvalue)){
    pt.100  <- Pvalue[[i]]
    xcoords <- pt.100$numOverlaps$permuted
    xcoords <- xcoords[order(xcoords)]
    aux.big <- qnorm((1-pvalthres),mean=mean(xcoords,na.rm=TRUE),sd=sd(xcoords,na.rm=TRUE))
    aux.less<- qnorm(pvalthres,mean=mean(xcoords,na.rm=TRUE),sd=sd(xcoords,na.rm=TRUE))
    overlap.time <- c(xcoords,pt.100$numOverlaps$observed)
    z_score <-  as.numeric(scale(overlap.time,center = TRUE,scale = TRUE))
    type <- c(rep("permuted",1000),"observed")
    temp <- data.frame(numOP= overlap.time,Zscore=z_score,type=type)
    temp$anno.region <- names(Pvalue)[i]
    temp$pval        <- pt.100$numOverlaps$pval
    temp$alternative <- pt.100$numOverlaps$alternative
    temp$z_score     <- pt.100$numOverlaps$zscore
    p.value.table    <- rbind(p.value.table,temp)
  }
  table(p.value.table$alternative)
  summary(p.value.table$pval)
  
  ## setup order
  #p.value.table$tanno.region <- factor(p.value.table$anno.region, levels = c("CpG island", 
  #                                    "CpG shore", "CpG shelf", "CpG sea", "Promoter", "Enhancer"))
                                
  observed <-subset(p.value.table, type=="observed")
  observed$color <- "black"
  observed$color[which(observed$pval<0.05&observed$alternative=="greater")] <-"red"
  observed$color[which(observed$pval<0.05&observed$alternative=="less")] <-"blue"
  observed$anno.region <- factor(observed$anno.region, levels = c("CpG island", 
                                                                       "CpG shore", "CpG shelf", "CpG sea", "Promoter", "Enhancer"))
  ## sorting by z-scores
  #p.value.table$anno.region <- factor(p.value.table$anno.region, levels = observed$tanno.region[order(observed$z_score,decreasing = TRUE)])
  permuted <- subset(p.value.table, type=="permuted")
  permuted$anno.region <- factor(permuted$anno.region, levels = c("CpG island", 
                                                                  "CpG shore", "CpG shelf", "CpG sea", "Promoter", "Enhancer"))
  
  p1 <- ggplot(permuted, aes(x=anno.region, y=Zscore, color=anno.region)) 
  p1 <- p1 + geom_boxplot(outlier.shape = NA)
  p1 <- p1 + labs(title = NULL, x = NULL, y = "Z-score")
  p1 <- p1 + p_theme + theme(axis.text.x = element_text(angle = 90,vjust=0.6),
                          legend.position  =c(0.1,0.9))+guides(color=FALSE)
  p1 <- p1 + scale_color_grey(start=0.1,end=0.4)
  p1 <- p1 + geom_point(data=observed ,col=  observed$color,shape=16,size=3)
  
  #observed$anno.region <- factor(observed$anno.region,levels = observed$anno.region[order(observed$z_score,decreasing = TRUE)])
  p2 <- ggplot(observed, aes(x=anno.region, y=numOP, fill=anno.region))
  p2 <- p2 + geom_bar(stat="identity", position=position_dodge())
  p2 <- p2 + p_theme + theme(axis.text.x = element_blank())
  p2 <- p2 + labs(x=NULL, y="No. of overlapping", title=NULL) 
  p2 <- p2 + guides(fill=FALSE) + scale_fill_grey(start=0.4, end=0.8)
  
  ## combine p1 and p2
  g <- grid.arrange(p2, p1)
  ggsave(plot_name, width = 6, height = 5, g)
}

}
