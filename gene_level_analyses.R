library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(pheatmap)
library(clusterProfiler) # for PEA analysis
library('org.Hs.eg.db')
library(DOSE)
library(enrichplot) # for visualisations
library(ggupset) # for visualisations
library("org.Hs.eg.db", character.only = TRUE)
library('biomaRt')
library(plotly)
library(AnnotationDbi)
library(clusterProfiler)
library(Seurat)

##################### 
#GO analysis

fc_thresh <- 1.1
columns <- c(0,1,2) #Progenitors qiu and jenny
columns <- c(0,5,9)  #ExN 120
columns <- c(0,8,14)  #CPNs CHD8_3m_r1
columns <- c(0,1,11) #tissue

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
    df2 <- as.data.frame(lapply(data, unlist))
    
    if(dx == "T21"){
      print("T21")
      df2 <- df2[order(df2$DS_resid, decreasing = TRUE), ]
      data_up <- df2[df2$DS_resid >fc_thresh,]
      lb <- 1 - (fc_thresh - 1)
      data_down <- df2[df2$DS_resid < lb,]
      out1 <- paste("T21_HVG_", columns[j], sep="")
      out2 <- paste("T21_LVG_", columns[j], sep="")
    }
    else{
      df2 <- df2[order(df2$Nor_resid, decreasing = TRUE), ]
      data_up <- df2[df2$Nor_resid >fc_thresh,]
      lb <- 1 - (fc_thresh - 1)
      data_down <- df2[df2$Nor_resid < lb,]
      out1 <- paste("Eu_HVG_", columns[j], sep="")
      out2 <- paste("Eu_LVG_", columns[j], sep="")
    }
    
    genes[[out1]] <- data_up$Gene
    genes[[out2]] <- data_down$Gene
  }
  return(genes)
}

genes_t21 <- gene_list(fc_thresh, columns, 'T21', "DS_scrnaseq_5")
genes_Eu <- gene_list(fc_thresh, columns, "Eu", "DS_scrnaseq_5")

genes_eu_overlap_HVG <- Reduce(intersect,list(genes_Eu[[1]],genes_Eu[[3]],genes_Eu[[5]]))
genes_eu_overlap_LVG <- Reduce(intersect,list(genes_Eu[[2]],genes_Eu[[4]],genes_Eu[[6]]))
genes_t21_overlap_HVG <- Reduce(intersect,list(genes_t21[[1]],genes_t21[[3]],genes_t21[[5]]))
genes_t21_overlap_LVG <- Reduce(intersect,list(genes_t21[[2]],genes_t21[[4]],genes_t21[[6]]))

genes_eu_all_HVG <- unique(c(genes_Eu[[1]],genes_Eu[[3]],genes_Eu[[5]]))
genes_eu_all_LVG <- unique(c(genes_Eu[[2]],genes_Eu[[4]],genes_Eu[[6]]))
genes_t21_all_HVG <- unique(c(genes_t21[[1]],genes_t21[[3]],genes_t21[[5]]))
genes_t21_all_LVG <- unique(c(genes_t21[[2]],genes_t21[[4]],genes_t21[[6]]))


egn = bitr(genes_eu_overlap_HVG, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
eg2 = bitr(genes_eu_overlap_LVG, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
eg3 = bitr(genes_t21_overlap_HVG, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
eg4 = bitr(genes_t21_overlap_LVG, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

egnA = bitr(genes_eu_all_HVG, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
eg2A = bitr(genes_eu_all_LVG, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
eg3A = bitr(genes_t21_all_HVG, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
eg4A = bitr(genes_t21_all_LVG, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

genes_all <- c(egnA$ENTREZID, eg2A$ENTREZID, eg3A$ENTREZID, eg4A$ENTREZID)


variable <- rep("HVG", length(egnA$ENTREZID))
variable2 <- rep("LVG", length(eg2A$ENTREZID))
variable3 <- rep("HVG", length(eg3A$ENTREZID))
variable4 <- rep("LVG", length(eg4A$ENTREZID))
variable_all <- c(variable, variable2, variable3, variable4)

#T21
dx <- rep("Euploid", length(egnA$ENTREZID)+length(eg2A$ENTREZID))
dx2 <- rep("T21", length(eg3A$ENTREZID)+length(eg4A$ENTREZID))
dx_all <- c(dx, dx2)

#CHD8
dx <- rep("WT", length(egnA$ENTREZID)+length(eg2A$ENTREZID))
dx2 <- rep("MUT", length(eg3A$ENTREZID)+length(eg4A$ENTREZID))
dx_all <- c(dx, dx2)


mydf_all <- data.frame(Entrez=genes_all, Variability=variable_all, DX=dx_all)

formula_res_all <- compareCluster(Entrez~Variability+DX, data=mydf_all,
                                  fun = "enrichGO", pvalueCutoff = 0.05,
                                  OrgDb = org.Hs.eg.db, ont='ALL')
#dotplot(formula_res, x="Variability") + facet_grid(~DX)
bp2_all <- clusterProfiler::simplify(formula_res_all, cutoff=0.7, by="p.adjust", select_fun=min)
saveRDS(bp2_all, "ds_scrnaseq_plot_data_all.rds")

#bp2_all <- readRDS("Qiu_plot_data_all.rds")


pdf(file="ds_scrnaseq_GO.pdf", width=7, height=10.5)
dotplot(bp2_all, x="Variability", font.size = 15) + facet_grid(~DX)
dev.off()

#######################




####################### 
# Distributions

obj <- readRDS("seurat_objs/ds_scrnaseq_SCT_v1.rds")
DefaultAssay(obj)<-"RNA"
obj<-NormalizeData(obj)
obj<-ScaleData(obj)
obj1 <- obj[,obj$cca_clusters==0]
obj1 <- obj1[,obj1$dx=='T21']
counts <- as.matrix(obj1@assays[["RNA"]]@layers[["data"]])
DefaultAssay(obj1) <- "RNA"
rownames(counts) <- rownames(obj1)
hvg_counts <- counts[row.names(counts) %in% genes_t21[["T21_HVG_0"]], ]
lvg_counts <- counts[row.names(counts) %in% genes_t21[["T21_LVG_0"]], ]
hvg_counts <- as.data.frame(hvg_counts)
hvg_counts$means <- rowMeans(hvg_counts)
hvg_counts_bin1 <- hvg_counts[hvg_counts$means<1, ]
hvg_counts_bin1 <- hvg_counts_bin1[hvg_counts_bin1$means>0.25, ]
hvg_counts_bin2 <- hvg_counts[hvg_counts$means>=1, ]
hvg_counts_bin0 <- hvg_counts[hvg_counts$means<=.25, ]

hvg_counts_bin0$means <- NULL
hvg_counts_bin1$means <- NULL
hvg_counts_bin2$means <- NULL

lvg_counts <- as.data.frame(lvg_counts)
lvg_counts$means <- rowMeans(lvg_counts)
lvg_counts_bin1 <- lvg_counts[lvg_counts$means<1, ]
lvg_counts_bin1 <- lvg_counts_bin1[lvg_counts_bin1$means>0.25, ]
lvg_counts_bin2 <- lvg_counts[lvg_counts$means>=1, ]
lvg_counts_bin0 <- lvg_counts[lvg_counts$means<=.25, ]

lvg_counts_bin0$means <- NULL
lvg_counts_bin1$means <- NULL
lvg_counts_bin2$means <- NULL

library(ggplot2)
library(reshape2)


pdf(file="DS_HVG_lowbin.pdf", width=5, height=3)
ggplot(melt(t(hvg_counts_bin0)),aes(x=value, fill=Var2)) + geom_density(alpha=0.25, linewidth=0.05)+theme(legend.position="none") +  
  xlab("Normalized Gene Expression") + xlim(0, 1) + ggtitle("HVG <0.25 Expression") + 
  theme(
    # LABELS APPEARANCE
    axis.title.x = element_text(size=14, face="bold", colour = "black"),    
    axis.title.y = element_text(size=14, face="bold", colour = "black"),    
    axis.text.x = element_text(size=12, face="bold", colour = "black"), 
    axis.text.y = element_text(size=12, face="bold", colour = "black"), # bold
    strip.text.x = element_text(size = 10, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 10, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3)
  ) + ylab("Density")
dev.off()

pdf(file="DS_HVG_midbin.pdf", width=5, height=3)
ggplot(melt(t(hvg_counts_bin1)),aes(x=value, fill=Var2)) + geom_density(alpha=0.25)+theme(legend.position="none") +  
  xlab("Normalized Gene Expression") + xlim(0, 5) + ggtitle("HVG 1>x>0.25 Expression") +
  theme(
    # LABELS APPEARANCE
    axis.title.x = element_text(size=14, face="bold", colour = "black"),    
    axis.title.y = element_text(size=14, face="bold", colour = "black"),    
    axis.text.x = element_text(size=12, face="bold", colour = "black"), 
    axis.text.y = element_text(size=12, face="bold", colour = "black"), # bold
    strip.text.x = element_text(size = 10, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 10, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3)
  ) + ylab("Density")
dev.off()

pdf(file="DS_HVG_highbin.pdf", width=5, height=3)
ggplot(melt(t(hvg_counts_bin2)),aes(x=value, fill=Var2)) + geom_density(alpha=0.25)+theme(legend.position="none") +  
  xlab("Normalized Gene Expression") + xlim(0,5) + ggtitle("HVG >1 Expression") + ylim(0, 3) +
  theme(
    # LABELS APPEARANCE
    axis.title.x = element_text(size=14, face="bold", colour = "black"),    
    axis.title.y = element_text(size=14, face="bold", colour = "black"),    
    axis.text.x = element_text(size=12, face="bold", colour = "black"), 
    axis.text.y = element_text(size=12, face="bold", colour = "black"), # bold
    strip.text.x = element_text(size = 10, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 10, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3)
  ) + ylab("Density")
dev.off()

pdf(file="DS_LVG_lowbin.pdf", width=5, height=3)
ggplot(melt(t(lvg_counts_bin0)),aes(x=value, fill=Var2)) + geom_density(alpha=0.25)+theme(legend.position="none") +  
  xlab("Normalized Gene Expression") + xlim(0, 1) + ggtitle("LVG <0.25 Expression") +
  theme(
    # LABELS APPEARANCE
    axis.title.x = element_text(size=14, face="bold", colour = "black"),    
    axis.title.y = element_text(size=14, face="bold", colour = "black"),    
    axis.text.x = element_text(size=12, face="bold", colour = "black"), 
    axis.text.y = element_text(size=12, face="bold", colour = "black"), # bold
    strip.text.x = element_text(size = 10, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 10, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3)
  ) + ylab("Density")
dev.off()
pdf(file="DS_LVG_midbin.pdf", width=5, height=3)
ggplot(melt(t(lvg_counts_bin1)),aes(x=value, fill=Var2)) + geom_density(alpha=0.25)+theme(legend.position="none") +  
  xlab("Normalized Gene Expression") + xlim(0, 5) + ggtitle("LVG 1>x>0.25 Expression")  + ylim(0, 6) +
  theme(
    # LABELS APPEARANCE
    axis.title.x = element_text(size=14, face="bold", colour = "black"),    
    axis.title.y = element_text(size=14, face="bold", colour = "black"),    
    axis.text.x = element_text(size=12, face="bold", colour = "black"), 
    axis.text.y = element_text(size=12, face="bold", colour = "black"), # bold
    strip.text.x = element_text(size = 10, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 10, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3)
  ) + ylab("Density")
dev.off()
pdf(file="DS_LVG_highbin.pdf", width=5, height=3)
ggplot(melt(t(lvg_counts_bin2)),aes(x=value, fill=Var2)) + geom_density(alpha=0.25)+theme(legend.position="none") +  
  xlab("Normalized Gene Expression") + xlim(0,5) + ggtitle("LVG >1 Expression") + ylim(0, 3) +
  theme(
    # LABELS APPEARANCE
    axis.title.x = element_text(size=14, face="bold", colour = "black"),    
    axis.title.y = element_text(size=14, face="bold", colour = "black"),    
    axis.text.x = element_text(size=12, face="bold", colour = "black"), 
    axis.text.y = element_text(size=12, face="bold", colour = "black"), # bold
    strip.text.x = element_text(size = 10, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 10, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3)
  ) + ylab("Density")
dev.off()
####################


####################
# Constraint analysis

gnomad = read.delim("gnomad.v4.1.constraint_metrics.tsv", header = TRUE, sep = "\t")
gnomad <- gnomad[gnomad$canonical=='true',]
gnomad <- gnomad[!duplicated(gnomad$gene), ]
gnomad <- dplyr::select(gnomad, gene, lof.oe_ci.upper)
fc_thresh <- 1.5

gene_list <- function(column, file){
  out_file <- paste(file,"/down_cluster_",column,".csv",sep="")
  data <- read.csv(out_file)
  data <- data[data$Mean_DS >0.7, ]
  data <- data[order(data$DS_resid),]
  LVG <- data[data$DS_resid<1,]
  HVG <- data[data$DS_resid>1,]
  gene_list <- list()
  if (length(LVG$Gene) >150){
    lvg_genes <- head(LVG, 150)
  }
  else{
    lvg_genes <- LVG
  }
  
  if (length(HVG$Gene) >150){
    hvg_genes <- top_n(HVG, 150, DS_resid)
  }
  else{
    hvg_genes <- HVG
  }
  gene_list[["LVG"]] <- lvg_genes$Gene
  gene_list[["HVG"]] <- hvg_genes$Gene
  return(gene_list)
  
}

genes_t21 <- gene_list(2, 'CHD8/120_5')

t21_LVG <- gnomad[gnomad$gene %in% genes_t21[['LVG']],] 
t21_HVG <- gnomad[gnomad$gene %in% genes_t21[['HVG']],] 

dx <- rep("LVG", length(t21_LVG$lof.oe_ci.upper))
dx2 <- rep("HVG", length(t21_HVG$lof.oe_ci.upper))
dx_ov <- c(dx, dx2)
euploid <- data.frame(Variability=dx_ov, LOEUF=rbind(t21_LVG,t21_HVG)$lof.oe_ci.upper)

library(ggplot2)
library(ggpubr)
library(rstatix)
bxp <- ggboxplot(
  euploid, x = "Variability", y = "LOEUF", 
  ylab = "LOEUF", xlab = "Variability"
)
bxp
stat.test <- euploid %>% 
  wilcox_test(LOEUF ~ Variability) %>%
  add_significance()
stat.test
euploid %>% wilcox_effsize(LOEUF ~ Variability)
stat.test <- stat.test %>% add_xy_position(x = "Variability")


pdf(file="ds_scrnaseq_constraint.pdf", width=4, height=3)
bxp + 
  stat_pvalue_manual(stat.test, tip.length = 0) +
  labs(subtitle = get_test_label(stat.test, detailed = TRUE))
dev.off()
#######################


######################
# EnrichR Analysis

library(enrichR)
gene_list <- function(column, file){
  out_file <- paste(file,"/down_cluster_",column,".csv",sep="")
  data <- read.csv(out_file)
  data <- data[data$Mean_DS >0.7, ]
  data <- data[order(data$DS_resid),]
  LVG <- data[data$DS_resid<1,]
  HVG <- data[data$DS_resid>1,]
  gene_list <- list()
  if (length(LVG$Gene) >150){
    lvg_genes <- head(LVG, 150)
  }
  else{
    lvg_genes <- LVG
  }
  
  if (length(HVG$Gene) >150){
    hvg_genes <- top_n(HVG, 150, DS_resid)
  }
  else{
    hvg_genes <- HVG
  }
  gene_list[["LVG"]] <- lvg_genes$Gene
  gene_list[["HVG"]] <- hvg_genes$Gene
  return(gene_list)
  
}

websiteLive <- getOption("enrichR.live")
dbs <- c("Epigenomics_Roadmap_HM_ChIP-seq", "ENCODE_TF_ChIP-seq_2015", "ENCODE_Histone_Modifications_2015")

datasets <- c("DS_scrnaseq_5", "Qiu_5", "tissue_5", 'CHD8/CHD8_GM_3m_r1_5', 'CHD8/120_5')
cluster <- c(2,2,0,8,0)
hvg_tfs <- list()
lvg_tfs <- list()
hvg_hist <- list()
lvg_hist <- list()

genH <- list()
genL <- list()

for(i in 1:length(datasets)){
  print(datasets[i])
  genes_t21 <- gene_list(cluster[i], datasets[i])
  genH[[i]] <- genes_t21[["HVG"]]
  genL[[i]] <- genes_t21[["LVG"]]
  
  print(length(genes_t21[["HVG"]]))
  if (websiteLive) {
    enriched_HVG <- enrichr(genes_t21[["HVG"]], dbs)
    #print(genes_t21[["HVG"]])
  }
  if (websiteLive) {
    enriched_LVG <- enrichr(genes_t21[["LVG"]], dbs)
  }
  
  enriched_HVG[[1]] <- enriched_HVG[[1]][order(enriched_HVG[[1]]$Term),]
  enriched_HVG[[2]] <- enriched_HVG[[2]][order(enriched_HVG[[2]]$Term),]
  enriched_LVG[[1]] <- enriched_LVG[[1]][order(enriched_LVG[[1]]$Term),]
  enriched_LVG[[2]] <- enriched_LVG[[2]][order(enriched_LVG[[2]]$Term),]
  print(enriched_HVG[[2]])
  
  enriched_HVG[[1]] <- dplyr::select(enriched_HVG[[1]], Term, Adjusted.P.value)
  enriched_HVG[[2]] <- dplyr::select(enriched_HVG[[2]], Term, Adjusted.P.value)
  enriched_LVG[[1]] <- dplyr::select(enriched_LVG[[1]], Term, Adjusted.P.value)
  enriched_LVG[[2]] <- dplyr::select(enriched_LVG[[2]], Term, Adjusted.P.value)
  
  hvg_hist[[i]]<-enriched_HVG[[1]]
  lvg_hist[[i]]<-enriched_LVG[[1]]
  
  hvg_tfs[[i]]<-enriched_HVG[[2]]
  lvg_tfs[[i]]<-enriched_LVG[[2]]
  
}
col.names <- c("DS LVG","Qiu LVG", "Tissue LVG", "CHD8 LVG","Villa LVG",
               "DS HVG","Qiu HVG", "Tissue HVG", "CHD8 HVG","Villa HVG")

data <- c(lvg_hist, hvg_hist)
for(i in 1:length(col.names)){
  if(i==1 | i==2){
    merged <-  merge(data[[1]], data[[2]], by = "Term")
  }
  else{
    merged <- merge(merged, data[[i]], by = "Term")
  }
}



merged_subset <- merged[grep("Fetal Brain", merged$Term), ]
merged_subset2 <- merged[grep("H1 Derived Neuronal Progenitor Cultured Cells", merged$Term), ]
merged_subset3 <- merged[grep("H9", merged$Term), ]
merged_subset4 <- merged[grep("iPS DF 19.11", merged$Term), ]
merged_subset5 <- merged[grep("iPS DF 6.9", merged$Term), ]
merged_subset6 <- merged[grep("Hippocampus", merged$Term), ]

merged_all <- rbind(merged_subset, merged_subset2)
merged_all <- rbind(merged_all, merged_subset3)
merged_all <- rbind(merged_all, merged_subset4)
merged_all <- rbind(merged_all, merged_subset5)
merged_all <- rbind(merged_all, merged_subset6)
merged_all <- merged_all[order(merged_all$Term),]

merged_all2 <- merged_all[grep("27ac", merged_all$Term), ]
merged_all3 <- merged_all[grep("27me3", merged_all$Term), ]
merged_all4 <- merged_all[grep("K4me3", merged_all$Term), ]
merged_all6 <- merged_all[grep("K4me1", merged_all$Term), ]
merged_all7 <- merged_all[grep("36me3", merged_all$Term), ]

merged_all5 <- rbind(merged_all2, merged_all4)
merged_all5 <- rbind(merged_all5, merged_all3)
merged_all5 <- rbind(merged_all5, merged_all6)
merged_all5 <- rbind(merged_all5, merged_all7)
#merged_all5 <- merged_all5[order(merged_all5$Term),]

rownames(merged_all5) <- merged_all5$Term
merged_all5$Term <- NULL
colnames(merged_all5) <- col.names

merged_all5 <- -1*log10(merged_all5[,1:ncol(merged_all5)])

library(gplots)
library(pheatmap)  
library(RColorBrewer)


pdf(file="histone_modifications_matrix.pdf",width=10, height=8)
temp<- pheatmap(as.matrix(merged_all5), cluster_rows=FALSE, cluster_cols=FALSE,
                fontsize_row= 15,
                labels_col = col.names,
                fontsize_col = 12, angle_col= 90, main = "Histone Modifications T21",
                legend_breaks = c(1,3,5,7,9,11,13,15,17,max(merged_all5)),
                legend_labels = c('1','3','5','7','9','11','13','15','17',"-log(p)\n"),
)


pdf(file="histone_modifications_matrix.pdf",width=9, height=8)
grid.draw(temp)
dev.off()



data <- c(lvg_tfs, hvg_tfs)
for(i in 1:length(col.names)){
  if(i==1 | i==2){
    merged <-  merge(data[[1]], data[[2]], by = "Term")
  }
  else{
    merged <- merge(merged, data[[i]], by = "Term")
  }
}

rownames(merged) <- merged$Term
merged$Term <- NULL
colnames(merged) <- col.names

rowMin(as.data.frame(merged))
merged <- as.data.frame(merged)
merged$min <- do.call(pmin, merged)
merged_filter <- merged[merged$min<1e-10,]

merged_filter$tfs <- str_split_fixed(row.names(merged_filter), " ",3)[,1]
#write.csv(merged, "TF_allDatasets_500.csv")
all_tfs <- unique(merged_filter$tfs)

for(i in 1:length(all_tfs)){
  tf <- all_tfs[i]
  db <- merged_filter[merged_filter$tfs == tf,]
  if (nrow(db)>1){
    db <- db[which.min(db$min),]
  }
  if(i==1){
    tf_db <- db
  }
  else{
    tf_db <- rbind(tf_db, db)
  }
  
}

#view(merged_filter)
#

rownames(tf_db) <- tf_db$tfs
tf_db$tfs <- NULL
tf_db$min <- NULL
tf_db <- -1*log10(tf_db[,1:ncol(tf_db)])

temp <- pheatmap(as.matrix(tf_db), cluster_rows=FALSE, cluster_cols=FALSE,
                 fontsize_row= 10,
                 labels_col = col.names,
                 fontsize_col = 12, angle_col= 90, main = "Enriched Transcripton Factors",
                 legend_breaks = c(1,10,20,30,40,50,max(tf_db)),
                 legend_labels = c('1','10','20','30','40','50',"-log(p)\n"),
)

pdf(file="tf_modifications_matrix.pdf",width=9, height=8)
grid.draw(temp)
dev.off()

tfs_want <- c("KAT2A", "TAF1", "YY1","FOXM1")
write.csv(tf_db[rownames(tf_db) %in% tfs_want, ], "Example_tfs.csv")

############################


###########################
# Half_life Analysis

setwd("~/BROAD/Variable_gene_lists")
hl <- read.csv("half_life.csv", header=TRUE, sep =",", 
               quote="", comment.char="", row.names=1)

#data <- read.csv("DS_scrnaseq_5/WT_cluster_2.csv")
# data_hl <- merge(data, hl, by.x = 'Gene', by.y = 'Gene.name')

images <- list()
cluster <- c(2,2,0,8,0)
datasets <- c("DS_scrnaseq_5", "Qiu_5", "tissue_5", 'CHD8/CHD8_GM_3m_r1_5', 'CHD8/120_5')
for (i in 1:length(datasets)){
  file <- datasets[i]
  out_file <- paste(file,"/down_cluster_",cluster[i],".csv",sep="")
  data <- read.csv(out_file)
  data_hl <- merge(data, hl, by.x = 'Gene', by.y = 'Gene.name')
  
  p <- ggplot(data_hl, aes(x=half.life..PC1., y=DS_resid)) +
    geom_point(size=1, color='#ADBCE6')+
    theme_bw() + scale_y_log10() +
    labs(x = "mRNA Half Life", y = expression("Residual " [T21]), title = file) +
    theme(
      # LABELS APPEARANCE
      axis.title.x = element_text(size=14, face="bold", colour = "black"),    
      axis.title.y = element_text(size=14, face="bold", colour = "black"),    
      axis.text.x = element_text(size=12, face="bold", colour = "black"), 
      axis.text.y = element_text(size=12, face="bold", colour = "black"), # bold
      strip.text.x = element_text(size = 10, face="bold", colour = "black" ),
      strip.text.y = element_text(size = 10, face="bold", colour = "black"),
      axis.line.x = element_line(color="black", size = 0.3),
      axis.line.y = element_line(color="black", size = 0.3),
      panel.border = element_rect(colour = "black", fill=NA, size=0.3),
      plot.title = element_text(hjust=0.5, face='bold')
    ) +geom_smooth(method = "lm", formula = y ~ x) +stat_cor() 
  images[[i]] <- p
}
images_w <- list()
for (i in 1:length(datasets)){
  file <- datasets[i]
  out_file <- paste(file,"/WT_cluster_",cluster[i],".csv",sep="")
  data <- read.csv(out_file)
  data_hl <- merge(data, hl, by.x = 'Gene', by.y = 'Gene.name')
  
  p <- ggplot(data_hl, aes(x=half.life..PC1., y=Nor_resid)) +
    geom_point(size=1, color='#ADBCE6')+
    theme_bw() + scale_y_log10() +
    labs(x = "mRNA Half Life", y = expression("Residual " [WT]), title = file) +
    theme(
      # LABELS APPEARANCE
      axis.title.x = element_text(size=14, face="bold", colour = "black"),    
      axis.title.y = element_text(size=14, face="bold", colour = "black"),    
      axis.text.x = element_text(size=12, face="bold", colour = "black"), 
      axis.text.y = element_text(size=12, face="bold", colour = "black"), # bold
      strip.text.x = element_text(size = 10, face="bold", colour = "black" ),
      strip.text.y = element_text(size = 10, face="bold", colour = "black"),
      axis.line.x = element_line(color="black", size = 0.3),
      axis.line.y = element_line(color="black", size = 0.3),
      panel.border = element_rect(colour = "black", fill=NA, size=0.3),
      plot.title = element_text(hjust=0.5, face='bold')
    ) +geom_smooth(method = "lm", formula = y ~ x) +stat_cor()
  images_w[[i]] <- p
}


library(ggplot2)
library(gridExtra)
library(ggpubr)

ggarrange(images[[1]],images[[2]],images[[3]],images[[4]],images[[5]],
          images_w[[1]], images_w[[2]], images_w[[3]], images_w[[4]], images_w[[5]],
          ncol = 5, nrow = 2)
#####################

#####################
# Differential Gene Expression vs Differential Variance

pdf(file="ds_scrnaseq_cluster2_varVsExp.pdf", width=5, height=5)
ggplot(res, aes(x=avg_log2FC, y=resid)) +
  geom_point(size=1, color=keyvals)+ geom_hline(yintercept=1,linetype="dashed") +geom_hline(yintercept=-1,linetype="dashed") +
  geom_vline(xintercept = 1.5, linetype="dashed")+geom_vline(xintercept = -1.5, linetype="dashed") +
  xlim(-8,8) + ylim(-25,25) + theme_bw() +
  labs(x = expression(Differential~Gene~Expression~(log[2]~FC)), y = expression(atop("Differential Variability", paste("Residual " [T21], " - Residual " [Euploid])))) +
  theme(
    # LABELS APPEARANCE
    axis.title.x = element_text(size=14, face="bold", colour = "black"),    
    axis.title.y = element_text(size=14, face="bold", colour = "black"),    
    axis.text.x = element_text(size=12, face="bold", colour = "black"), 
    axis.text.y = element_text(size=12, face="bold", colour = "black"), # bold
    strip.text.x = element_text(size = 10, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 10, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3)
  ) +geom_smooth(method = "lm", formula = y ~ x) +stat_cor()
dev.off()


keyvals <- ifelse(res$avg_log2FC < -1.5 & res$resid < -1 , 'red',
                  ifelse(res$avg_log2FC < -1.5 & res$resid > 1, 'red',
                         ifelse(res$avg_log2FC > 1.5 & res$resid < -1, 'red',
                                ifelse(res$avg_log2FC > 1.5 & res$resid > 1, 'red',
                                       'gray'))))
#keyvals[is.na(keyvals)] <- 'gray'
names(keyvals)[keyvals == 'red'] <- 'Downregulated, More Variable'
names(keyvals)[keyvals == 'gray'] <- 'NS'
names(keyvals)[keyvals == 'red'] <- 'Upregulated, Less Variable'
names(keyvals)[keyvals == 'red'] <- 'Downregulated, Less Variable'
names(keyvals)[keyvals == 'red'] <- 'Upregulated, more variable'



## venn diagram overlap figure
genes_t21_qiu <- genes_t21
genes_Eu_qiu <- genes_Eu
genes_t21_DS <- genes_t21
genes_Eu_DS <- genes_Eu

genes_qiu_HVG <- Reduce(intersect,list(genes_t21_qiu[[1]],genes_t21_qiu[[3]],genes_t21_qiu[[5]]))
genes_qiu_LVG <- Reduce(intersect,list(genes_t21_qiu[[2]],genes_t21_qiu[[4]],genes_t21_qiu[[6]]))
genes_DS_HVG <- Reduce(intersect,list(genes_t21_DS[[1]],genes_t21_DS[[3]],genes_t21_DS[[5]]))
genes_DS_LVG <- Reduce(intersect,list(genes_t21_DS[[2]],genes_t21_DS[[4]],genes_t21_DS[[6]]))

HVG_list <- list(
  "Qiu et al" = genes_qiu_HVG,
  "Upadhya et al" = genes_DS_HVG
)

LVG_list <- list(
  "Qiu et al" = genes_qiu_LVG,
  "Upadhya et al" = genes_DS_LVG
)


library(VennDiagram) 
temp <- venn.diagram(HVG_list, fill = c("#0073C2FF", "#CD534CFF"), 
                     alpha = c(0.5, 0.5), lwd =2, filename=NULL,
                     cat.fontface = "bold",fontface = "bold",cat.fontfamily = "sans",
                     fontfamily = "sans",cat.pos = c(-27, 27),
                     cat.dist = c(0.055, 0.055),)

temp_LVG <- venn.diagram(LVG_list, fill = c("#0073C2FF", "#CD534CFF"), 
                         alpha = c(0.5, 0.5), lwd =2, filename=NULL,
                         cat.fontface = "bold",fontface = "bold",cat.fontfamily = "sans",
                         fontfamily = "sans",cat.pos = c(-27, 27),
                         cat.dist = c(0.055, 0.055),)


pdf(file="venn_HVG_overlap.pdf")
grid.draw(temp)
dev.off()

pdf(file="venn_LVG_overlap.pdf")
grid.draw(temp_LVG)
dev.off()



genes_t21_120 <- genes_t21
genes_Eu_120 <- genes_Eu
genes_t21_gm <- genes_t21
genes_Eu_gm <- genes_Eu

genes_120_HVG <- Reduce(intersect,list(genes_t21_120[[1]],genes_t21_120[[3]],genes_t21_120[[5]]))
genes_120_LVG <- Reduce(intersect,list(genes_t21_120[[2]],genes_t21_120[[4]],genes_t21_120[[6]]))
genes_GM_HVG <- Reduce(intersect,list(genes_t21_gm[[1]],genes_t21_gm[[3]],genes_t21_gm[[5]]))
genes_GM_LVG <- Reduce(intersect,list(genes_t21_gm[[2]],genes_t21_gm[[4]],genes_t21_gm[[6]]))

HVG_list <- list(
  "Villa et al" = genes_120_HVG,
  "Paulsen et al" = genes_GM_HVG
)

LVG_list <- list(
  "Villa et al" = genes_120_LVG,
  "Paulsen et al" = genes_GM_LVG
)


temp <- venn.diagram(HVG_list, fill = c("#0073C2FF", "#CD534CFF"), 
                     alpha = c(0.5, 0.5), lwd =2, filename=NULL,
                     cat.fontface = "bold",fontface = "bold",cat.fontfamily = "sans",
                     fontfamily = "sans",cat.pos = c(-27, 27),
                     cat.dist = c(0.055, 0.055),)

temp_LVG <- venn.diagram(LVG_list, fill = c("#0073C2FF", "#CD534CFF"), 
                         alpha = c(0.5, 0.5), lwd =2, filename=NULL,
                         cat.fontface = "bold",fontface = "bold",cat.fontfamily = "sans",
                         fontfamily = "sans",cat.pos = c(-27, 27),
                         cat.dist = c(0.055, 0.055),)


pdf(file="venn_HVG_overlap_CHD8.pdf")
grid.draw(temp)
dev.off()

pdf(file="venn_LVG_overlap_CHD8.pdf")
grid.draw(temp_LVG)
dev.off()

#####################
# Supplemental
# Plotting CV2 vs Mean

ct <- c("Progenitor 1", "Cycling Progenitor", "Progenitor 2",
        "Neuronal Progenitor", "Unknown", "Gliogenic Progenitor")
exp_plot <- list()
for(i in 1:6){
  out_file <- paste("DS_scrnaseq_5/all/down_cluster_",i-1,".csv",sep="")
  data <- read.csv(out_file)
  keyvals <- ifelse(data$DS_resid < 0.9, 'blue',
                    ifelse(data$DS_resid > 1.1, 'red', 'gray'))
  p<-ggplot(data, aes(x=Mean_DS, y=CV2_DS)) +
    geom_point(size=1, color=keyvals) +
    xlim(-8,8) + ylim(-25,25) + theme_bw() +
    scale_y_log10() + scale_x_log10(labels = scales::comma) +
    labs(x = 'log(Mean Expression)', y = expression('log('~CV^{2}~')'), title = ct[i]) +
    theme(
      # LABELS APPEARANCE
      axis.title.x = element_text(size=14, face="bold", colour = "black"),    
      axis.title.y = element_text(size=14, face="bold", colour = "black"),    
      axis.text.x = element_text(size=12, face="bold", colour = "black"), 
      axis.text.y = element_text(size=12, face="bold", colour = "black"), # bold
      strip.text.x = element_text(size = 10, face="bold", colour = "black" ),
      strip.text.y = element_text(size = 10, face="bold", colour = "black"),
      axis.line.x = element_line(color="black", size = 0.3),
      axis.line.y = element_line(color="black", size = 0.3),
      panel.border = element_rect(colour = "black", fill=NA, size=0.3),
      plot.title = element_text(hjust=0.5, face='bold')
    ) + geom_line(aes(Mean_DS, CV2_DS_Exp))
  exp_plot[[i]] <- p
}
exp_plot_wt <- list()
for(i in 1:6){
  out_file <- paste("DS_scrnaseq_5/all/WT_cluster_",i-1,".csv",sep="")
  data <- read.csv(out_file)
  keyvals_wt <- ifelse(data$Nor_resid < 0.9, 'blue',
                       ifelse(data$Nor_resid > 1.1, 'red', 'gray'))
  
  p<-ggplot(data, aes(x=Mean_Nor, y=CV2_Nor)) +
    geom_point(size=1, color=keyvals_wt) +
    xlim(-8,8) + ylim(-25,25) + theme_bw() +
    scale_y_log10() + scale_x_log10(labels = scales::comma) +
    labs(x = 'log(Mean Expression)', y = expression('log('~CV^{2}~')'), title = ct[i]) +
    theme(
      # LABELS APPEARANCE
      axis.title.x = element_text(size=14, face="bold", colour = "black"),    
      axis.title.y = element_text(size=14, face="bold", colour = "black"),    
      axis.text.x = element_text(size=12, face="bold", colour = "black"), 
      axis.text.y = element_text(size=12, face="bold", colour = "black"), # bold
      strip.text.x = element_text(size = 10, face="bold", colour = "black" ),
      strip.text.y = element_text(size = 10, face="bold", colour = "black"),
      axis.line.x = element_line(color="black", size = 0.3),
      axis.line.y = element_line(color="black", size = 0.3),
      panel.border = element_rect(colour = "black", fill=NA, size=0.3),
      plot.title = element_text(hjust=0.5, face='bold')
    ) + geom_line(aes(Mean_Nor, CV2_Nor_Exp))
  exp_plot_wt[[i]] <- p
}


ggarrange(exp_plot[[1]],exp_plot[[2]],exp_plot[[3]],exp_plot[[4]],exp_plot[[5]],exp_plot[[6]],
          exp_plot_wt[[1]], exp_plot_wt[[2]], exp_plot_wt[[3]], exp_plot_wt[[4]], exp_plot_wt[[5]],exp_plot_wt[[6]],
          ncol = 6, nrow = 2)


