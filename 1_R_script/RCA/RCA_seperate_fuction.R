### prepare input data
##input data were dataframe
#input <- df.TPM
#reference <- refPanel.filtered_of_fetal


### Function calculate corelation Matrix
corMatrix <- function(input, reference, method){
  geneList <- intersect(rownames(input),rownames(reference))
  df_combine = cbind(input[geneList,],reference[geneList,])
  cor_matrix = cor(df_combine,method = method)
  cor_matrix = as.data.frame(cor_matrix)
  cell_type_matrix = cor_matrix[(dim(input)[2]+1):dim(cor_matrix)[2],1:dim(input)[2]]
}
# f
#test <- corMatrix(log10(df.FPKM+1),refPanel.filtered_of_fetal, method = "pearson")
#test <- corMatrix(df.TPM,refPanel.filtered_of_fetal, method = "pearson")
# convert the cor value
#test_for_clust_power = abs(test)^power * sign(test)

#### Generate cell clusters
RCA.cluster <- function(matrix, deepSplit_wgcna=1, min_group_Size_wgcna=5){
  #### 1: reading input ####
  if (!require(flashClust)) install.packages("flashClust",repos = "http://cran.us.r-project.org") 
  require(flashClust)
  if (!require(WGCNA)){
    source("http://bioconductor.org/biocLite.R");
    biocLite(c("impute", "GO.db", "preprocessCore")); 
    install.packages("WGCNA");
  }
  require(WGCNA)  
  
  #### 2: choosing method ####
  #  if (method == "hclust"){
  matrix[is.na(matrix)]<-0
  d = as.dist(1-cor(matrix,method="pearson"));
  cellTree = flashClust(d,method = "average");     
  dynamicGroups = cutreeDynamic(dendro = cellTree,distM = as.matrix(d),deepSplit = deepSplit_wgcna,
                                pamStage = FALSE,minClusterSize= min_group_Size_wgcna);
  
  dynamicColors = labels2colors(dynamicGroups);
  group_labels = matrix(dynamicGroups,nrow = length(names(matrix)),ncol = 1,list(names(matrix),c("groupLabel")),byrow=FALSE)
  group_labels = as.data.frame(group_labels)
  group_labels_color = cbind(group_labels,dynamicColors);
  clusterResault <- list()
  clusterResault[[1]] <- cellTree
  clusterResault[[2]] <- group_labels_color
  return(clusterResault)
}

#clusterResault <- RCA.cluster(test_for_clust_power)
#clusterResault[[1]]
### Function to plot PCA figure
#matrix <- test_for_clust_power
RCA.PCA <- function(matrix, color){
  ## PCA analysis
  pr = prcomp(t(scale(matrix))); 
  pc_projection = as.data.frame(pr$x);
  cell_projection = pc_projection[,1:2];
  ############## figure plot ################
  ## plot PCA
  plot(cell_projection[,1],cell_projection[,2],
       type="p",pch = 21,
       xlab = colnames(cell_projection)[1],ylab = colnames(cell_projection)[2],
       col = as.character(color),lwd = 1,bg = as.character(color),
       main = "PCA of cell clusters in RCA space", cex = 1, cex.main = 2,
       font.axis = 1, font.lab = 1, font.main = 1);
}
#RCA.PCA(test_for_clust_power, color = clusterResault[[2]]$dynamicColors)

### Function to plot heatmap figure
RCA.heatmap <- function(matrix, cellTree, color){
  #### 1: reading input ####
  if (!require(gplots)) install.packages("gplots",repos = "http://cran.us.r-project.org") 
  require(gplots)
  color_scheme = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                    "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(100);
  heatmap.2(as.matrix(matrix),
            col=color_scheme,
            Colv=as.dendrogram(cellTree),
            ColSideColors=as.vector(color),
            scale= "column", 
            margins=c(5,20),
            trace="none",
            key = TRUE,
            keysize = 0.5,
            cexCol = 1,cexRow =1,
            labCol = "") 
}
#RCA.heatmap(test_for_clust_power, cellTree = clusterResault[[1]], color = as.character(group_labels_color$dynamicColors))
