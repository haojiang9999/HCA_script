data(sysdata, envir=environment())
point_cex=1
cluster_color_labels=NULL
c = dynamicColors;
pr = prcomp(t(scale(fpkm_for_clust))); 
pc_projection = as.data.frame(pr$x);
cell_projection = pc_projection[,1:2];
pch_to_use = 21;
cex_to_use = point_cex;

png("GSE97693_globel_RCAplot_PCA_RCA_clusters.png")
main_name = "PCA of cell clusters in RCA space";
par(mar=c(1,1,1,1)*7)  
plot(cell_projection[,1],cell_projection[,2],
     type="p",pch = pch_to_use,
     xlab = colnames(cell_projection)[1],ylab = colnames(cell_projection)[2],
     col = c,lwd = 1,bg = c,
     main = main_name, cex = cex_to_use, cex.main = 2,
     font.axis = 1, font.lab = 1, font.main = 1
);
dev.off();


if (!is.null(cluster_color_labels)){
  
#  png("RCAplot_PCA_external_labels.png")
  main_name = "PCA of cell clusters in RCA space";
  #par(mar=c(1,1,1,1)*7)  
  plot(cell_projection[,1],cell_projection[,2],
       type="p",pch = pch_to_use,
       xlab = colnames(cell_projection)[1],ylab = colnames(cell_projection)[2],
       col = cluster_color_labels,lwd = 1,bg = cluster_color_labels,
       main = main_name, cex = cex_to_use, cex.main = 2,
       font.axis = 1, font.lab = 1, font.main = 1
  );
 # dev.off();
library(gplots)
png("GSE97693_globel_RCAplot_heatmap_RCA_clusters.png",width = 3000, height = 3000, units = "px", pointsize = 15);
  color_scheme =     colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                        "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(100);
  heatmap.2(as.matrix(fpkm_for_clust),
            col=color_scheme,
            Colv=as.dendrogram(cellTree),
            ColSideColors=c,
            scale= "column", 
            margins=c(5,20),
            trace="none",
            key = TRUE,
            keysize = 0.5,
            cexCol = 1,cexRow =1,
            labCol = "")  
  dev.off();
  