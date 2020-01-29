# test for pheatmap plot
col.test <- samples.anno.df[,2,drop = F]
col.test <- sapply(test, as.character)
col.test <- as.data.frame(test)
rownames(col.test)<- samples
library(pheatmap)
pheatmap(t(pmat), annotation_col = col.test,
         show_colnames = F,cutree_cols = 10 )
