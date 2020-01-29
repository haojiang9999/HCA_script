### Seurat3 colon adult cell filtered
dir.create("Seurat3_analysis")
setwd("Seurat3_analysis/")
getwd()
## Input
#Expression matrix
filePath <- "/stor/jianghao/GEO/Normal_tissues/GSE103239_Tang_digestive_tract/GSE103154_adult_tissues/R_GSE103154_Tang_adult_colon_normal"
file <- paste0(filePath,"/GSE103254.adult.tpm.cellFiltered.rds")
GSE103254.adult.tpm <- readRDS(file)
exprMat <- GSE103254.adult.tpm$GSE103254.adult.cellFlitered.tpm
cellInfo <- GSE103254.adult.tpm$cellInfo.cellFiltered
#### Using Seurat do cluster jobs
# load data set
### Setup the Seurat Object
library(dplyr)
library(Seurat)
Seurat.adult.colon.1464.1464 <- CreateSeuratObject(counts = exprMat, 
                                         project = "GSE103254.adult.tpm.cellFiltered",
                                         min.cells = 3, min.features = 200)
Seurat.adult.colon.1464.1464
### Standard pre-processing workflow
# Visualize QC metrics as a violin plot
VlnPlot(Seurat.adult.colon.1464, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
#plot1 <- FeatureScatter(Seurat.adult.colon.1464, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(Seurat.adult.colon.1464, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#CombinePlots(plots = list(plot1, plot2))
### Normalizing the data
Seurat.adult.colon.1464 <- NormalizeData(Seurat.adult.colon.1464, normalization.method = "LogNormalize", scale.factor = 10000)
### Identification of highly variable features (feature selection)
Seurat.adult.colon.1464 <- FindVariableFeatures(Seurat.adult.colon.1464, 
                                           selection.method = "mean.var.plot", 
                                           nfeatures = 2000)
#Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Seurat.adult.colon.1464), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Seurat.adult.colon.1464)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
### Scaling the data
all.genes <- rownames(Seurat.adult.colon.1464)
Seurat.adult.colon.1464 <- ScaleData(Seurat.adult.colon.1464, features = all.genes)
### Perform linear dimensional reduction
Seurat.adult.colon.1464 <- RunPCA(Seurat.adult.colon.1464, features = VariableFeatures(object = Seurat.adult.colon.1464))

# Examine and visualize PCA results a few different ways
print(Seurat.adult.colon.1464[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Seurat.adult.colon.1464, dims = 1:2, reduction = "pca")
DimPlot(Seurat.adult.colon.1464, reduction = "pca")
DimHeatmap(Seurat.adult.colon.1464, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(Seurat.adult.colon.1464, dims = 1:15, cells = 500, balanced = TRUE)
### Determine the ‘dimensionality’ of the dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Seurat.adult.colon.1464 <- JackStraw(Seurat.adult.colon.1464, num.replicate = 100)
Seurat.adult.colon.1464 <- ScoreJackStraw(Seurat.adult.colon.1464, dims = 1:20)

### Cluster the cells
Seurat.adult.colon.1464 <- FindNeighbors(Seurat.adult.colon.1464, dims = 1:10)
Seurat.adult.colon.1464 <- FindClusters(Seurat.adult.colon.1464, resolution = 0.5)
head(Idents(Seurat.adult.colon.1464), 5)

### Run non-linear dimensional reduction (UMAP/tSNE)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
Seurat.adult.colon.1464 <- RunUMAP(Seurat.adult.colon.1464, dims = 1:20)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Seurat.adult.colon.1464, reduction = "umap")
### Finding differentially expressed features (cluster biomarkers)
# find all markers of cluster 1
cluster1.markers <- FindMarkers(Seurat.adult.colon.1464, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(Seurat.adult.colon.1464, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)
# find markers for every cluster compared to all remaining cells, report only the positive ones
Seurat.adult.colon.1464.markers <- FindAllMarkers(Seurat.adult.colon.1464, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Seurat.adult.colon.1464.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
VlnPlot(Seurat.adult.colon.1464, features = c("MS4A1", "CD79A"))
# you can plot raw counts as well
VlnPlot(Seurat.adult.colon.1464, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
FeaturePlot(Seurat.adult.colon.1464, features = c("CHGA", "CHGB", "PYY", 
                                             "SCG2", "CA1", "FABP1", 
                                             "MS4A12", "SI", 
                                             "CDH1"))
FeaturePlot(Seurat.adult.colon.1464, features = c("CHGA", "PYY", "TPH1", 
                                                  "MUC2", "TFF3", "FCGBP", 
                                                  "CA1", "CA2", "FABP1",
                                                  "PLA2G2A", "OLFM4", "MKI67",
                                                  "CD79A", "VIM", "EPCAM"))
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
#### Assigning cell type identity to clusters
new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

saveRDS(pbmc, file = "../output/pbmc3k_final.rds")


########### some thing
Seurat.adult.colon.1464@meta.data$orig.ident <- "adult.colon"
Seurat.adult.colon.1464@meta.data$seurat_clusters
table(Seurat.adult.colon.1464@meta.data$seurat_clusters)