#### Distance Cluster
#### 1.Multiple dimension deduction plot of the data
source("/data8t_4/JH/MyJobs/1_R_script/R_Plot/DimeReduPlot.R")

DimeReduPlot(mx = Min_max_norm.8000, color = Pheno.merged$Cell_info, 
             tiltle = "Cor.Res.CV8000_Cell_info", print = F)
DimeReduPlot(mx = Min_max_norm.4000, color = Pheno.merged$Cell_info, 
             tiltle = "Cor.Res.CV4000_Cell_info", print = F)
DimeReduPlot(mx = Min_max_norm.2000, color = Pheno.merged$Cell_info, 
             tiltle = "Cor.Res.CV2000_Cell_info", print = F)
DimeReduPlot(mx = Min_max_norm.1500, color = Pheno.merged$Cell_info, 
             tiltle = "Cor.Res.CV1500_Cell_info", print = F,
             check_duplicates = FALSE)
DimeReduPlot(mx = Min_max_norm.1000, color = Pheno.merged$Cell_info, 
             tiltle = "Cor.Res.CV1000_Cell_info", print = F, 
             check_duplicates = FALSE)
### Cell Type

DimeReduPlot(mx = Min_max_norm.8000, color = Pheno.merged$CellType, 
             tiltle = "Cor.Res.CV8000_Cell_Type", print = F)
DimeReduPlot(mx = Min_max_norm.4000, color = Pheno.merged$CellType, 
             tiltle = "Cor.Res.CV4000_Cell_Type", print = F)
DimeReduPlot(mx = Min_max_norm.2000, color = Pheno.merged$CellType, 
             tiltle = "Cor.Res.CV2000_Cell_Type", print = F)
DimeReduPlot(mx = Min_max_norm.1500, color = Pheno.merged$CellType, 
             tiltle = "Cor.Res.CV1500_Cell_Type", print = F,
             check_duplicates = FALSE)
DimeReduPlot(mx = Min_max_norm.1000, color = Pheno.merged$CellType, 
             tiltle = "Cor.Res.CV1000_Cell_Type", print = F,
             check_duplicates = FALSE)


#### 2. Cluster and heatmap plot
source("/data8t_4/JH/MyJobs/1_R_script/R_Plot/heatmap.JH.R")
heatmap.JH(Min_max_norm.8000,show_colnames = F,
           annotation_col = Pheno.merged)

heatmap.JH(Min_max_norm.4000,show_colnames = F,
           annotation_col = Pheno.merged)
heatmap.JH(Min_max_norm.2000,show_colnames = F,
           annotation_col = Pheno.merged)
### replace the NA into 0
Min_max_norm.1500[is.na(Min_max_norm.1500)] <- 0
heatmap.JH(Min_max_norm.1500,show_colnames = F,
           annotation_col = Pheno.merged)
Min_max_norm.1000[is.na(Min_max_norm.1000)] <- 0 # replace the NA into 0
heatmap.JH(Min_max_norm.1000,show_colnames = F,
           annotation_col = Pheno.merged)


