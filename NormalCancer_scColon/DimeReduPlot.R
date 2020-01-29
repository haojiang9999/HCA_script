#### Dimension reductjion Plot
# Input was a matrix or table
# scaled mx maybe better
DimeReduPlot <- function(mx , color, tiltle = "Something", print = T){
  ### Check the packages
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  if (!require("ggplot2")) 
    BiocManager::install("ggplot2")
  if (!require("Rtsne")) 
    BiocManager::install("Rtsne")
  if (!require("umap")) 
    BiocManager::install("umap")
  if (!require("ggpubr")) 
    BiocManager::install("ggpubr")
  require(ggplot2)
  # mx = scale(Cor.merged)
  #### PCA plot
  pca <- prcomp(t(mx)); 
  pca_out <- as.data.frame(pca$x)
  p1<-ggplot(pca_out,aes(x=PC1,y=PC2,color = color ))
  p1<-p+geom_point() + labs(title = tiltle, subtitle = "PCA plot")
  
  ### t-SNE plot
  require(Rtsne)
  #t(test_for_clust_power)
  TSNE <- Rtsne(t(mx), dims = 2, initial_dims = 50, perplexity = 50)
  tSNE_out <- as.matrix(TSNE$Y)
  tsne_plot <- data.frame(x = tSNE_out[,1], y =tSNE_out[,2], col = color)
  p2 <- ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=col)) + labs(title = tiltle, subtitle = "t-SNE plot")

  #### UMAP plot
  require(umap)
  UMAP <- umap(t(mx))
  UMAP_out <- UMAP$layout
  umap_plot <- data.frame(x =UMAP_out[,1] , y =UMAP_out[,2], col = color)
  p3 <- ggplot(umap_plot) + geom_point(aes(x=x, y=y, color=col))+ labs(title = tiltle, subtitle = "UMAP plot")
  #### Plot or Save the results
  require(ggpubr)
  if(print){
  print(p1)
  print(p2)
  print(p3)}else{
    splots <- list()
    splots[[1]] <- p1
    splots[[2]] <- p2
    splots[[3]] <- p3
    res.arr <- ggarrange(plotlist = splots, nrow =3)
    ggexport(res.arr, filename = paste0(tiltle, "_PCA_tSNE_UMAP_Plots.pdf"),width = 15, height = 50)
  }
}
