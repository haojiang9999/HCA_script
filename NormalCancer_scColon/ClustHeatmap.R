#### Function for WGCNA cluster and heat map plot
ClustHeatmap <- function(mx,Pheno.merged, title,...){
  # mx = Cor.Res.CV4000$Cor.merged
  ## Step1 Cluster
  source("/data8t_4/JH/MyJobs/1_R_script/RCA/RCA_seperate_fuction.R")
  mx.cluster.res <- RCA.cluster(matrix = mx, deepSplit_wgcna=1, min_group_Size_wgcna=5)
  cellTree <- mx.cluster.res[[1]]
  # add cluster results to annotation
  group_labels_color <- mx.cluster.res[[2]]
  Pheno.merged.clus <- cbind(Pheno.merged,group_labels_color$groupLabel)
  Pheno.merged.clus <- Pheno.merged.clus[,-1]
  ## Step2 Heatmap plot
  source("/data8t_4/JH/MyJobs/1_R_script/R_Plot/heatmap.JH.R")
  #table(is.na(mx))
  heatmap.JH(scale(mx), main = paste0(title,"_heatmap"),
             filename = paste0(title,"_heatmap.pdf"),width = 15, height = 10,
             cluster_cols= cellTree, annotation_col = Pheno.merged.clus,show_rownames = T, 
             show_colnames = F,...)
  
}
