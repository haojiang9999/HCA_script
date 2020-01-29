#### Heatmap plot of candidata lncRNA
library(pheatmap)
# construct matrix for plotting
# GeneNames 
geneNmae <- as.character(epi_stem_DE_TvN.lnc.NameFiltered_fc_0.2$X)
# expression data
NM.fpkm.exp.stemTA <- RCA_epithelial.fpkm.exp.stemTA.pQ[geneNmae, samNM.stemTA]
# sample names had value in expression data
samTumor.stemTA <- samTumor.stemTA[samTumor.stemTA %in% colnames(RCA_epithelial.fpkm.exp.stemTA.pQ)]
tumor.fpkm.exp.stemTA <- RCA_epithelial.fpkm.exp.stemTA.pQ[geneNmae, samTumor.stemTA]
NM.num <- length(samNM.stemTA)
tumor.num <- length(samTumor.stemTA)
data_subset <- as.matrix(cbind(NM.fpkm.exp.stemTA,tumor.fpkm.exp.stemTA))

# create heatmap using pheatmap
my_sample_col <- data.frame(sample = rep(c("Normal", "Tumor"), c(NM.num,tumor.num)))
row.names(my_sample_col) <- colnames(data_subset)
data_subset_scaled <- scale(t(log2(data_subset)))
data_subset_scaled <- scale(t(data_subset_scaled))
pheatmap(data_subset_scaled[,1:NM.num],scale = c("none"),
         cluster_cols = T)

pheatmap(data_subset_scaled[,NM.num+1:tumor.num],scale = c("none"),
         cluster_cols = T)

