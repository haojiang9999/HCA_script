##### Seperate cancer cells into different sub groups
### Step1.Keep Cancer calls only
Pheno.Tumor <- Pheno.merged[Pheno.merged$CellType == "Tumor cells",]
TumorCells <- rownames(Pheno.Tumor)
Tumor.norm.1500 <- Min_max_norm.1500[,TumorCells]
Tumor.norm.2000 <- Min_max_norm.2000[,TumorCells]
### Step2.Check the heatmap of tumor cell cluster
#### Distance Cluster
## 1.Multiple dimension deduction plot of the data
source("/data8t_4/JH/MyJobs/1_R_script/R_Plot/DimeReduPlot.R")
DimeReduPlot(mx = Tumor.norm.2000, color = Pheno.Tumor$Cell_info, 
             tiltle = "Cor.Res.CV2000_Cell_info", print = T)
DimeReduPlot(mx = Tumor.norm.1500, color = Pheno.Tumor$Cell_info, 
             tiltle = "Cor.Res.CV1500_Cell_info", print = T,
             check_duplicates = FALSE)
## 2.Cluster and heatmap plot
source("/data8t_4/JH/MyJobs/1_R_script/R_Plot/heatmap.JH.R")
heatmap.JH(Tumor.norm.2000,show_colnames = F,
           annotation_col = Pheno.Tumor, filename = "Tumor_heatmap_2000.pdf", width = 15)
### replace the NA into 0
Tumor.norm.1500[is.na(Tumor.norm.1500)] <- 0  # Remove NA values
heatmap.JH(Tumor.norm.1500,show_colnames = F,
           annotation_col = Pheno.Tumor,filename = "Tumor_heatmap_1500.pdf", width = 15)
## 3.Re-cluster using SC3
#### Cluster using SC3
source("/data8t_4/JH/MyJobs/1_R_script/NormalCancer/JH_SC3_cluster.R")
SC3_cluster_2000 <- JH_SC3_cluster(Tumor.norm.2000,Pheno.Tumor,ks=2:4)
hc_2000 <- SC3_cluster_2000$`3`$hc
#SC3_cluster_1500 <- JH_SC3_cluster(Tumor.norm.1500,Pheno.Tumor,ks=2:4)
#hc_1500 <- SC3_cluster_1500$`3`$hc
## Using SC3_cluster_2000 first
source("/data8t_4/JH/MyJobs/1_R_script/R_Plot/heatmap.JH.R")
heatmap.JH(Tumor.norm.2000,show_colnames = F,
           annotation_col = Pheno.merged[,c(1,2)], cluster_cols = hc_2000,filename = "Tumor_heatmap_2000_SC3.pdf", width = 15)
### replace the NA into 0
#heatmap.JH(Tumor.norm.1500,show_colnames = F,
#           annotation_col = Pheno.merged[,c(1,2)], cluster_cols = hc_1500,filename = "Tumor_heatmap_1500_SC3.pdf", width = 15)

############### Using SC3_cluster_2000 group cancer cells #################
# 1.Check the hierarchical structure
plot(as.dendrogram(hc_2000))           
d1 <- cut(as.dendrogram(hc_2000), h=16) # seperate into one samll(15) one big(1425)
d1.1 <- d1$lower[[1]]
d1.2 <- d1$lower[[2]]
plot(d1.1)
plot(d1.2)
d2 <- cut(d1.2, h=13)
d2.1 <- d2$lower[[1]]
d2.2 <- d2$lower[[2]]
plot(d2.1)
plot(d2.2)
# Cut d2.1 into two groups
d3 <- cut(d2.1, h=8)
d3.1<- d3$lower[[1]]
d3.2<- d3$lower[[2]]
plot(d3.1)
plot(d3.2)
library(ggdendro)
#ggdendro::dendro_data	
d3.1_data <- dendro_data(d3.1)
d3.2_data <- dendro_data(d3.2)
Group1 <- d3.1_data$labels
Group2 <- d3.2_data$labels
# Cut d2.2 into four groups
d4 <- cut(d2.2, h=10)
d4.1<- d4$lower[[1]]
# Group3
plot(d4.1)
d4.1_data <- dendro_data(d4.1)
Group3 <- d4.1_data$labels
# Group 4,5,6
d4.2<- d4$lower[[2]]
plot(d4.2)
d5 <- cut(d4.2, h=8)
# Group4
d5.1 <- d5$lower[[1]]
plot(d5.1)
d5.1_data <- dendro_data(d5.1)
Group4 <- d5.1_data$labels
# Group5,6
d5.2 <- d5$lower[[2]]
plot(d5.2)
d6 <- cut(d5.2, h=6)
# Group5
d6.1 <- d6$lower[[1]]
plot(d6.1)
d6.1_data <- dendro_data(d6.1)
Group5 <- d6.1_data$labels
# Group6
d6.2 <- d6$lower[[2]]
plot(d6.2)
d6.2_data <- dendro_data(d6.2)
Group6 <- d6.2_data$labels
# 2.Find cell names in each groups
## These cells were ordered by SC3 cluster
G1_cells <- as.character(Group1$label)
G2_cells <- as.character(Group2$label)
G3_cells <- as.character(Group3$label)
G4_cells <- as.character(Group4$label)
G5_cells <- as.character(Group5$label)
G6_cells <- as.character(Group6$label)
