library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(pheatmap)
library(Seurat)
library(gridExtra)
library(ggpubr)
library(grid)
library(reshape2)

######
# Code used to create correlation plots
#
# SCTransform normalized gene matrix is taken for genes that are either 
# HVG or LVG in both genotypes. Then the Pearson correlations for each 
# gene pair is calculated independently in each genotype. Then, the 
# difference between the absolute values of the T21 and Euploid matrices
# is calculated and visualized. 
#####


fc_thresh <- 1.1
columns <- c(0,1,2) #tissue
f <- "DS_scrnaseq_5"

gene_list <- function(fc_thresh, columns, dx, file){
  genes <- list()
  for (j in 1:length(columns)){
    
    if(dx == "T21"){
      out_file <- paste(file,"/down_cluster_",columns[j],".csv",sep="")
    }
    else{
      out_file <- paste(file,"/WT_cluster_",columns[j],".csv",sep="")
    }
    
    data <- read.csv(out_file)
    fc <- list()
    for (i in 1:length(data$Unnamed..0)){
      if (data[i,6] <1 ){
        fc[i] <- (-1/data[i,6])
      }
      else{
        fc[i] <- data[i,6]
      }
    }
    data$fc <- fc
    df2 <- as.data.frame(lapply(data, unlist))
    df2 <- df2[order(df2$fc, decreasing = TRUE), ]
    
    if(dx == "T21"){
      print("T21")
      data_up <- df2[df2$fc >fc_thresh,]
      lb <- -1*fc_thresh
      data_down <- df2[df2$fc < lb,]
      out1 <- paste("T21_HVG_", columns[j], sep="")
      out2 <- paste("T21_LVG_", columns[j], sep="")
    }
    else{
      data_up <- df2[df2$fc >fc_thresh,]
      lb <- -1*fc_thresh
      data_down <- df2[df2$fc < lb,]
      out1 <- paste("Eu_HVG_", columns[j], sep="")
      out2 <- paste("Eu_LVG_", columns[j], sep="")
    }
    
    genes[[out1]] <- data_up$Gene
    genes[[out2]] <- data_down$Gene
  }
  return(genes)
}

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((cormat))
  hc <- hclust(dd, method="average")
  cormat <-cormat[hc$order, hc$order]
}

genes_t21 <- gene_list(fc_thresh, columns, 'T21', f)
genes_Eu <- gene_list(fc_thresh, columns, "Eu", f)


obj <- readRDS("seurat_objs/ds_scrnaseq_SCT_v1.rds")


co <- list()
for (i in 1:length(columns)){
  c <- columns[i]
  obj1 <- obj[,obj$cca_clusters==c]
  
  objDS <- obj1[,obj1$dx == 'T21']
  objeu <- obj1[,obj1$dx=='Euploid']
  #objDS <- obj1[,obj1$dx == 'DS']
  #objeu <- obj1[,obj1$dx=='WT']
  
  var_list_t21 <- paste("T21_HVG_",c,sep="")
  #var_list_t21_L <- paste("T21_LVG_",c,sep="")
  var_list_eu <- paste("Eu_HVG_",c,sep="")
  #var_list_eu_L <- paste("Eu_LVG_",c,sep="")
  
  genes_overlap <- Reduce(intersect,list(genes_t21[[var_list_t21]], genes_Eu[[var_list_eu]]))
  

  sct_counts_DS <- as.matrix(objDS@assays[["SCT"]]@data)
  sct_counts_DS <- sct_counts_DS[rownames(sct_counts_DS) %in% genes_overlap,]
  sct_counts_DS <- t(sct_counts_DS)
  res <- cor(sct_counts_DS)
  

  sct_counts_eu <- as.matrix(objeu@assays[["SCT"]]@data)
  sct_counts_eu <- sct_counts_eu[rownames(sct_counts_eu) %in% genes_overlap,]
  sct_counts_eu <- t(sct_counts_eu)
  res2 <- cor(sct_counts_eu)
  
  res3 <- abs(res) - abs(res2)
  

  co[[paste("clusMean_",c,sep="")]] <- mean(res3)

  
  col<- colorRampPalette(c("blue", "white", "red"))(20)
  out_png <- paste("corr_take7_HVG_gene_gene_",f,"_" ,c, ".png",sep="")
  
  print(c)
  
  p<- ggplot(data = xx, aes(x=Var1, y=Var2, fill=value)) +
   geom_tile() + scale_fill_gradient2(low = "darkblue", high = "darkred",
                         mid = "white", midpoint = 0,
                         name="Change in\nCorrelation") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())

  png(file=out_png, width=5.5,height=5,units="in",res=1500)
  grid.draw(p)
  dev.off()

  rm(res,res2)
  res3<- data.frame(row=rownames(res3)[row(res3)[upper.tri(res3)]],
                    col=colnames(res3)[col(res3)[upper.tri(res3)]],
                    corr=res3[upper.tri(res3)])

  res3 <- res3[order(res3$corr), ]
  res3$rank <-  rank(res3$corr)

  co[[paste(c,"_dec", sep="")]] <- (nrow(res3[res3$corr<0,]) / nrow(res3)) * 100  #%dec
  co[[paste(c,"_inc", sep="")]] <- (nrow(res3[res3$corr>0,]) / nrow(res3)) * 100  #%inc
  co[[paste(c,"_size", sep="")]] <- length(genes_overlap)


  out_pdf <- paste("gene_gene_pairHVG_",c,"_" ,f,".pdf",sep="")
  p <- ggplot(data=res3, aes(x=rank, y=corr)) +
    geom_line() +
    labs(x = 'Gene-Gene Pairs', y = "Change in Correlation Strength", title = "Progenitor 0") +
    theme(
      # LABELS APPEARANCE
      axis.title.x = element_text(size=14, face="bold", colour = "black"),
      axis.title.y = element_text(size=14, face="bold", colour = "black"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size=12, face="bold", colour = "black"), # bold
      strip.text.x = element_text(size = 10, face="bold", colour = "black" ),
      strip.text.y = element_text(size = 10, face="bold", colour = "black"),
      axis.line.x = element_line(color="black", size = 0.3),
      axis.line.y = element_line(color="black", size = 0.3),
      panel.border = element_rect(colour = "black", fill=NA, size=0.3),
      plot.title = element_text(hjust=0.5, face='bold')
    )


  pdf(file=out_pdf, width= 6, height = 3)
  grid.draw(p)
  dev.off()

}

saveRDS(co, "corr_HVG.rds")

## LVG
co <- list()
for (i in 1:length(columns)){
  c <- columns[i]
  obj1 <- obj[,obj$cca_clusters==c]
  
  objDS <- obj1[,obj1$dx == 'T21']
  objeu <- obj1[,obj1$dx=='Euploid']
  #objDS <- obj1[,obj1$dx == 'DS']
  #objeu <- obj1[,obj1$dx=='WT']
  
  #var_list_t21 <- paste("T21_HVG_",c,sep="")
  var_list_t21_L <- paste("T21_LVG_",c,sep="")
  #var_list_eu <- paste("Eu_HVG_",c,sep="")
  var_list_eu_L <- paste("Eu_LVG_",c,sep="")
  
  genes_overlap <- Reduce(intersect,list(genes_t21[[var_list_t21_L]], genes_Eu[[var_list_eu_L]]))
  

  sct_counts_DS <- as.matrix(objDS@assays[["SCT"]]@data)
  sct_counts_DS <- sct_counts_DS[rownames(sct_counts_DS) %in% genes_overlap,]
  sct_counts_DS <- t(sct_counts_DS)
  res <- cor(sct_counts_DS)

  sct_counts_eu <- as.matrix(objeu@assays[["SCT"]]@data)
  sct_counts_eu <- sct_counts_eu[rownames(sct_counts_eu) %in% genes_overlap,]
  sct_counts_eu <- t(sct_counts_eu)
  res2 <- cor(sct_counts_eu)
  
  res3 <- abs(res) - abs(res2)
  
  co[[paste("clusMean_",c,sep="")]] <- mean(res3)

  
  col<- colorRampPalette(c("blue", "white", "red"))(20)
  out_png <- paste("corr_take4LVG_gene_gene_",f,"_" ,c, ".png",sep="")

  setwd("~/BROAD/Figures")
  print(c)

  #heatmap(x = res3, col = col, symm = TRUE, labRow = FALSE, labCol = FALSE)

  p<- ggplot(data = melt(reorder_cormat(res3)), aes(x=Var1, y=Var2, fill=value)) +
    geom_tile() + scale_fill_gradient2(low = "blue", high = "red",
                                       mid = "white", midpoint = 0,
                                       name="Change in\nCorrelation") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())

  png(file=out_png, width=5.5,height=5,units="in",res=1500)
  grid.draw(p)
  dev.off()


  rm(res,res2)
  res3<- data.frame(row=rownames(res3)[row(res3)[upper.tri(res3)]],
                    col=colnames(res3)[col(res3)[upper.tri(res3)]],
                    corr=res3[upper.tri(res3)])

  res3 <- res3[order(res3$corr), ]
  res3$rank <-  rank(res3$corr)

  co[[paste(c,"_dec", sep="")]] <- (nrow(res3[res3$corr<0,]) / nrow(res3)) * 100  #%dec
  co[[paste(c,"_inc", sep="")]] <- (nrow(res3[res3$corr>0,]) / nrow(res3)) * 100  #%inc
  co[[paste(c,"_size", sep="")]] <- length(genes_overlap)
  #co[[paste("clus_",c,sep="")]] <- res3

  out_pdf <- paste("gene_gene_pairLVG_",c,"_" ,f,".pdf",sep="")

  p<- ggplot(data=res3, aes(x=rank, y=corr)) +
    geom_line() +
    labs(x = 'Gene-Gene Pairs', y = "Change in Correlation Strength", title = "Progenitor 0") +
    theme(
      # LABELS APPEARANCE
      axis.title.x = element_text(size=14, face="bold", colour = "black"),
      axis.title.y = element_text(size=14, face="bold", colour = "black"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size=12, face="bold", colour = "black"), # bold
      strip.text.x = element_text(size = 10, face="bold", colour = "black" ),
      strip.text.y = element_text(size = 10, face="bold", colour = "black"),
      axis.line.x = element_line(color="black", size = 0.3),
      axis.line.y = element_line(color="black", size = 0.3),
      panel.border = element_rect(colour = "black", fill=NA, size=0.3),
      plot.title = element_text(hjust=0.5, face='bold')
    )

  pdf(file=out_pdf, width= 6, height = 3)
  grid.draw(p)
  dev.off()
  
}
saveRDS(co, "corr_LVG.rds")
