#### 5.1.Surpervised_Methyl_Plot.R
COAD.methy.set
COAD.methy.blue <- COAD.methy.set[ , COAD.methy.set$dynamicColors == "blue"]
pheatmap::pheatmap(COAD.methy.blue,annotation_col = pData(COAD.methy.set)[,1:2],
                   main = "Methylation_Blue VS orthers_rescaled",annotation_colors = ann_colors)




