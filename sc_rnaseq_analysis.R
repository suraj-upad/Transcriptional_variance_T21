library(Seurat)
library(ggplot2)
library(SeuratData)
library(SeuratWrappers)
library(Azimuth)
library(patchwork)
library(dplyr)
library(DoubletFinder)
options(future.globals.maxSize = 1e10)
options(Seurat.object.assay.version = "v5")

#pull in cell-cycle genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

#load in data
ds.data <- Read10X(data.dir = "DS_scrnaseq/DS1")
euploid.data <- Read10X(data.dir = "DS_scrnaseq/DS2U")

#create initial v3 seurat objects and annotate Mitochondrial genes and genotype
ds_obj <- CreateSeuratObject(counts = ds.data, min.cells = 3, min.features = 200)
ds_obj$dx <- "T21"
ds_obj[["percent.mt"]] <- PercentageFeatureSet(ds_obj, pattern = "^MT-")

euploid_obj <- CreateSeuratObject(counts = euploid.data, min.cells = 3, min.features = 200)
euploid_obj$dx <- "Euploid"
euploid_obj[["percent.mt"]] <- PercentageFeatureSet(euploid_obj, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(ds_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# filter out low quality cells
ds_obj <- subset(ds_obj, subset =  nCount_RNA < 300000  & nFeature_RNA > 500 & percent.mt < 20)
euploid_obj <- subset(euploid_obj, subset =  nCount_RNA < 200000  & nFeature_RNA > 500 & percent.mt < 20)

# create list of seurat objects
seu_obj <- list()
seu_obj[["ds"]] <- ds_obj
seu_obj[["euploid"]] <- euploid_obj

estimates <- c(0.08, 0.065) #doublet formation rates
doublet_finder_use <- TRUE


# loop through seurat object list, normalize using SCTransform
# remove doublets using doublet finder
for (i in 1:length(seu_obj)) {
  print(i)
  seu_obj[[i]] <- SCTransform(seu_obj[[i]], verbose = FALSE)

  if (doublet_finder_use == TRUE){
    seu_obj[[i]] <- RunPCA(seu_obj[[i]])
    seu_obj[[i]] <- RunUMAP(seu_obj[[i]], dims = 1:10)

    #pk
    sweep.list <- paramSweep(seu_obj[[i]], PCs = 1:10, sct = TRUE)
    sweep.stats <- summarizeSweep(sweep.list)
    bcmvn <- find.pK(sweep.stats)
    bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
    optimal.pk <- bcmvn.max$pK
    optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
    
    annotations <- seu_obj[[i]]@meta.data$seurat_clusters
    homotypic.prop <- modelHomotypic(annotations)           
    nExp_poi <- round(estimates[i]*nrow(seu_obj[[i]]@meta.data))  
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    
    seu_obj[[i]] <- doubletFinder(seu = seu_obj[[i]], 
                                    PCs = 1:10, 
                                    pK = optimal.pk,
                                    nExp = nExp_poi.adj, sct = TRUE)
    metadata <- seu_obj[[i]]@meta.data
    colnames(metadata)[9] <- "doublet_finder"
    seu_obj[[i]]@meta.data <- metadata 
    
    # subset and save
    singlets <- subset(seu_obj[[i]], doublet_finder == "Singlet")
    seu_obj[[i]] <- singlets
    remove(singlets)
    
  }
  
}

# merge objects
obj <- merge(seu_obj[["ds"]], seu_obj[["normal"]])
rm(seu_obj)

# Assign cell-cycle metrics and regress out cell-cycle by running SCTransform.
DefaultAssay(obj) <- "SCT"
obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
obj$CC.Difference <- obj$S.Score - obj$G2M.Score
obj <- SCTransform(obj, assay = 'RNA', vars.to.regress = c("CC.Difference"), verbose = FALSE,  vst.flavor = "v2")
obj <- RunPCA(obj)

#integrate using CCA, RPCA, and Harmony
obj <- IntegrateLayers(
  object = obj, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca", normalization.method="SCT",
  verbose = FALSE
)

obj <- IntegrateLayers(
  object = obj, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",normalization.method="SCT",
  verbose = FALSE
)
obj <- IntegrateLayers(
  object = obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",normalization.method="SCT",
  verbose = FALSE
)

# Find clusters and UMAP coordinates for each integration
DefaultAssay(obj) <- "SCT"
obj <- FindNeighbors(obj, reduction = "integrated.cca", dims = 1:30)
obj <- FindClusters(obj, resolution = .2, cluster.name = "cca_clusters")
obj <- RunUMAP(obj, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")

obj <- FindNeighbors(obj, reduction = "integrated.rpca", dims = 1:30)
obj <- FindClusters(obj, resolution = .2, cluster.name = "rpca_clusters")
obj <- RunUMAP(obj, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
obj <- FindClusters(obj, resolution = .25, cluster.name = "harmony_clusters")
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")



# Plot UMAP of each integration, CCA selected for our data
DimPlot(obj, reduction = "umap.harmony", split.by = "dx", label=TRUE)
DimPlot(obj, reduction = "umap.cca", split.by = "dx", label=TRUE)
DimPlot(obj, reduction = "umap.rpca",  split.by = "dx", label=TRUE) 


# visualize gene markers for cell-type annotation
features <- c("PAX6", "SOX2", 'VIM', 'NEB', 'DNM3', 'AEBP1', 'SOX2-OT', 'LINC01579', 'CFAP47', 'CFAP54', 'ARMC3', 'WDR49','PIF1','CDCA3', 'DCX', 'NEUROD1', 'STMN2', 'SLC1A4')
DotPlot(obj, features=features)
FeaturePlot(obj, features=features)


obj <- PrepSCTFindMarkers(obj)
Idents(obj) <- "cca_clusters"
obj.markers <- FindAllMarkers(obj, only.pos = TRUE)


obj.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(obj, features = top10$gene) + NoLegend()

pdf(file="ds_marker_pred.pdf", width=12, height=8)
DoHeatmap(obj, features = top10$gene) + NoLegend()
dev.off()



# Assigning cell-type labels
new.cluster.ids <- c("Progenitor 1", "Cycling Progenitor", "Progenitor 2",
                     "Neuronal Progenitor", "Unknown", "Gliogenic Progenitor")
names(new.cluster.ids) <- levels(obj)
obj <- RenameIdents(obj, new.cluster.ids)

# plotting CCA integrated UMAP with cell-type labels
pdf(file="ds_npcs_UMAP.pdf", width=16, height=6)
DimPlot(obj, reduction = "umap.cca", split.by = "dx", label=TRUE)
dev.off()

# Differential Gene Expression analysis
# Compare genotypes within each cluster
obj$celltype.dx <- paste(obj$cca_clusters, obj$dx, sep = "_")
Idents(obj) <- "celltype.dx"

DefaultAssay(obj) <- "RNA"
clus0 <- FindMarkers(obj, ident.1 = "0_T21", ident.2 = "0_Euploid", verbose = FALSE, min.pct = 0.25)
clus1 <- FindMarkers(obj, ident.1 = "1_T21", ident.2 = "1_Euploid", verbose = FALSE, min.pct = 0.25)
clus2 <- FindMarkers(obj, ident.1 = "2_T21", ident.2 = "2_Euploid", verbose = FALSE, min.pct = 0.25)
clus3 <- FindMarkers(obj, ident.1 = "3_T21", ident.2 = "3_Euploid", verbose = FALSE, min.pct = 0.25)
clus4 <- FindMarkers(obj, ident.1 = "4_T21", ident.2 = "4_Euploid", verbose = FALSE, min.pct = 0.25)
clus5 <- FindMarkers(obj, ident.1 = "5_T21", ident.2 = "5_Euploid", verbose = FALSE, min.pct = 0.25)


####################################
# Cell-type variability analysis
# Distance to the Medoid analysis from Liu et al (2023)
# 
# Following code from Liu et al.:
# github: https://github.com/jiayiliujiayi/scRNA_Seq-Differential_Variability_Analysis/
# 
# Paper reference: 
# Liu, J., Kreimer, A. & Li, W. V. Differential variability analysis of
# single-cell gene expression data. Brief. Bioinform. 24, (2023).

DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
obj <- NormalizeData(obj, normalization.method = "RC", scale.factor = 1000)
obj <- ScaleData(obj, cell.type = NULL)
#pca
DefaultAssay(obj) <- "RNA"
all.genes <- rownames(obj)
obj <- RunPCA(obj, features = all.genes, npcs = min(200, ncol(obj) - 1), verbose = F)
#dense
coord_pca <- Embeddings(obj, reduction = "pca")
coord_pca <- coord_pca[, 1:min(200, ncol(coord_pca))]

source("densvis/densne/densne.R")
out_densne <- run_densne(coord_pca, no_dims=2, perplexity=50,
                         theta=0.5,
                         randseed=-1,
                         verbose=FALSE, use_pca=FALSE,
                         max_iter=1000, dens_frac=0.3, # as in the paper
                         final_dens=TRUE)
densne_embeddings = out_densne[[1]]
colnames(densne_embeddings) <- paste0("denSNE_", 1:2)
rownames(densne_embeddings) <- colnames(obj)
obj[["denSNE"]] <- CreateDimReducObject(embeddings = densne_embeddings,
                                        key = "denSNE_", 
                                        assay = DefaultAssay(obj))

# visualize DeSNE clustering
DimPlot(obj, reduction = "denSNE", group.by = c("dx", "cca_clusters"), combine = TRUE, label=TRUE)
DimPlot(obj, reduction = "denSNE", split.by = 'dx')


write.csv(densne_embeddings, 'dense_ds_scrnaseq_NPCs.csv')




####################################
# DESeq2 normalization for gene-level variance analysis

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = obj@assays[["RNA"]]@layers[["counts"]], colData = obj@meta.data, design = ~ dx)
dds <- estimateSizeFactors(dds, type='poscounts')   #poscounts to handle genes with some zeros
normalized_counts <- counts(dds, normalized=TRUE)

write.csv(normalized_counts, "Deseq2_norm_DS_scrnaseq.csv")

write.csv(obj@meta.data, 'meta_data_DS_scrnaseq_NPCS_2.csv')
write.csv(rownames(obj), 'genes_rows_DS_scrnaseq_NPCs.csv')
