#### 3.Basic_Molecular_analysis_heatmap_version.R
cutree.res <- Cluster.20200201.V7.Tumor$cutree.res
### Mutation geens
## 1.Read data
COAD_TCGA_PanCan33_mc3_nonsilentGene_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/TCGA_Pan_Cancer/somatic_mutation_SNP_INDEL/MC3_Gene_level_non_silent_mutation/COAD_TCGA_PanCan33_mc3_nonsilentGene_dataset.rds")
# COAD_TCGA_PanCan33_mc3_nonsilentGene_dataset$mc3.nonsilentGene.metadata
COAD.mc3.nonsilentGene.xena <- COAD_TCGA_PanCan33_mc3_nonsilentGene_dataset$COAD.mc3.nonsilentGene.xena

geneList <- c("APC", "TP53", "KRAS", "BRAF", "MLH1")
Gene.Mut <-t(COAD.mc3.nonsilentGene.xena[geneList,])
Gene.Mut <- as.data.frame(Gene.Mut)
Gene.Mut$rownames <- rownames(Gene.Mut)
MergeTable.Gene.Mut<- dplyr::left_join(Cluster.df, Gene.Mut, by = "rownames")
library(pheatmap)
test <- MergeTable.Gene.Mut
test[is.na(test)] <-2 
rownames(test) <- rownames(MergeTable.Gene.Mut)
anno_col <- test[,1:2]
rownames(anno_col) <- rownames(test)
pheatmap(t(test[,4:8]), cluster_cols = hclust.Res, annotation_col = anno_col)




