## 3.Projectjion results analysis
#### 1.Multiple dimension deduction plot of the data
source("/data8t_4/JH/MyJobs/1_R_script/R_Plot/DimeReduPlot.R")
#### Its sacled

DimeReduPlot(mx = scale(Cor.Res.CV8000$Cor.merged), color = Pheno.merged$Cell_info, 
             tiltle = "Cor.Res.CV8000_Cell_info", print = F)
DimeReduPlot(mx = scale(Cor.Res.CV4000$Cor.merged), color = Pheno.merged$Cell_info, 
             tiltle = "Cor.Res.CV4000_Cell_info", print = F)
DimeReduPlot(mx = scale(Cor.Res.CV2000$Cor.merged), color = Pheno.merged$Cell_info, 
             tiltle = "Cor.Res.CV2000_Cell_info", print = F)
DimeReduPlot(mx = scale(Cor.Res.CV1500$Cor.merged), color = Pheno.merged$Cell_info, 
             tiltle = "Cor.Res.CV1500_Cell_info", print = F)
DimeReduPlot(mx = scale(Cor.Res.CV1000$Cor.merged), color = Pheno.merged$Cell_info, 
             tiltle = "Cor.Res.CV1000_Cell_info", print = F)

### Cell Type

DimeReduPlot(mx = scale(Cor.Res.CV8000$Cor.merged), color = Pheno.merged$CellType, 
             tiltle = "Cor.Res.CV8000_Cell_Type", print = F)
DimeReduPlot(mx = scale(Cor.Res.CV4000$Cor.merged), color = Pheno.merged$CellType, 
             tiltle = "Cor.Res.CV4000_Cell_Type", print = F)
DimeReduPlot(mx = scale(Cor.Res.CV2000$Cor.merged), color = Pheno.merged$CellType, 
             tiltle = "Cor.Res.CV2000_Cell_Type", print = F)
DimeReduPlot(mx = scale(Cor.Res.CV1500$Cor.merged), color = Pheno.merged$CellType, 
             tiltle = "Cor.Res.CV1500_Cell_Type", print = F)
DimeReduPlot(mx = scale(Cor.Res.CV1000$Cor.merged), color = Pheno.merged$CellType, 
             tiltle = "Cor.Res.CV1000_Cell_Type", print = F)

#### 2. Cluster and heatmap plot
### cluster using WGCNA
## Step1 Cluster
source("/data8t_4/JH/MyJobs/1_R_script/RCA/RCA_seperate_fuction.R")
mx.cluster.res <- RCA.cluster(matrix = mx, deepSplit_wgcna=1, min_group_Size_wgcna=5)
cellTree <- mx.cluster.res[[1]]
# add cluster results to annotation
group_labels_color <- mx.cluster.res[[2]]
Pheno.merged.clus <- cbind(Pheno.merged,group_labels_color$groupLabel)
Pheno.merged.clus <- Pheno.merged.clus[,-1]
## Step2 Heatmap plot
source("/data8t_4/JH/MyJobs/1_R_script/NormalCancer/ClustHeatmap.R")
ClustHeatmap( Cor.Res.CV8000$Cor.merged, Pheno.merged= Pheno.merged, title = "Cor.Res.CV8000",scale = c("column"))
ClustHeatmap( Cor.Res.CV4000$Cor.merged, Pheno.merged= Pheno.merged, title = "Cor.Res.CV4000",scale = c("column"))
ClustHeatmap( Cor.Res.CV2000$Cor.merged, Pheno.merged= Pheno.merged, title = "Cor.Res.CV2000",scale = c("column"))
ClustHeatmap( Cor.Res.CV1500$Cor.merged, Pheno.merged= Pheno.merged, title = "Cor.Res.CV1500",scale = c("column"))
ClustHeatmap( Cor.Res.CV1000$Cor.merged, Pheno.merged= Pheno.merged, title = "Cor.Res.CV1000",scale = c("column"))

ClustHeatmap( Cor.Res.CV2000$Cor.merged, Pheno.merged= Pheno.merged, title = "Cor.Res.CV2000_Unscale")
ClustHeatmap( Cor.Res.CV1500$Cor.merged, Pheno.merged= Pheno.merged, title = "Cor.Res.CV1500_Unscale")
ClustHeatmap( Cor.Res.CV1000$Cor.merged, Pheno.merged= Pheno.merged, title = "Cor.Res.CV1000_Unscale")


