#### Using Seurat do cluster jobs
# load data set
GSE103254.adult.tpm.dataset <- readRDS("GSE103254.adult.tpm.dataset.rds")
df.GSE103254.adult.tpm <- GSE103254.adult.tpm.dataset$df.adult.colon.tpm
pd <- GSE103254.adult.tpm.dataset$phenotype
### Setup the Seurat Object
library(dplyr)
library(Seurat)
Seurat.adult.colon <- CreateSeuratObject(counts = df.GSE103254.adult.tpm, 
                                         project = "GSE103254.adult.tpm",
                                         min.cells = 3, min.features = 200)
Seurat.adult.colon
### Standard pre-processing workflow
# Visualize QC metrics as a violin plot
VlnPlot(Seurat.adult.colon, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
#plot1 <- FeatureScatter(Seurat.adult.colon, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(Seurat.adult.colon, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#CombinePlots(plots = list(plot1, plot2))
### Normalizing the data
Seurat.adult.colon <- NormalizeData(Seurat.adult.colon, normalization.method = "LogNormalize", scale.factor = 10000)
### Identification of highly variable features (feature selection)
Seurat.adult.colon <- FindVariableFeatures(Seurat.adult.colon, 
                                           selection.method = "mean.var.plot", 
                                           nfeatures = 2000)
#Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Seurat.adult.colon), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Seurat.adult.colon)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
### Scaling the data
all.genes <- rownames(Seurat.adult.colon)
Seurat.adult.colon <- ScaleData(Seurat.adult.colon, features = all.genes)
### Perform linear dimensional reduction
Seurat.adult.colon <- RunPCA(Seurat.adult.colon, features = VariableFeatures(object = Seurat.adult.colon))

# Examine and visualize PCA results a few different ways
print(Seurat.adult.colon[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Seurat.adult.colon, dims = 1:2, reduction = "pca")
DimPlot(Seurat.adult.colon, reduction = "pca")
DimHeatmap(Seurat.adult.colon, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(Seurat.adult.colon, dims = 1:15, cells = 500, balanced = TRUE)
### Determine the ‘dimensionality’ of the dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Seurat.adult.colon <- JackStraw(Seurat.adult.colon, num.replicate = 100)
Seurat.adult.colon <- ScoreJackStraw(Seurat.adult.colon, dims = 1:20)

### Cluster the cells
Seurat.adult.colon <- FindNeighbors(Seurat.adult.colon, dims = 1:10)
Seurat.adult.colon <- FindClusters(Seurat.adult.colon, resolution = 0.5)
head(Idents(Seurat.adult.colon), 5)

### Run non-linear dimensional reduction (UMAP/tSNE)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
Seurat.adult.colon <- RunUMAP(Seurat.adult.colon, dims = 1:20)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Seurat.adult.colon, reduction = "umap")
### Finding differentially expressed features (cluster biomarkers)
# find all markers of cluster 1
cluster1.markers <- FindMarkers(Seurat.adult.colon, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(Seurat.adult.colon, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)
# find markers for every cluster compared to all remaining cells, report only the positive ones
Seurat.adult.colon.markers <- FindAllMarkers(Seurat.adult.colon, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Seurat.adult.colon.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
VlnPlot(Seurat.adult.colon, features = c("MS4A1", "CD79A"))
# you can plot raw counts as well
VlnPlot(Seurat.adult.colon, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
FeaturePlot(Seurat.adult.colon, features = c("CHGA", "CHGB", "PYY", 
                                             "SCG2", "CA1", "FABP1", 
                                             "MS4A12", "SI", 
                                              "CDH1"))
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
Seurat.adult.colon@meta.data$orig.ident <- "adult.colon"
Seurat.adult.colon@meta.data$seurat_clusters
table(Seurat.adult.colon@meta.data$seurat_clusters)
