# RCA classify single cells
# load 
# data transformation of data.frame log10
pseudo_count = 1;
log_fpkm_temp = df.FPKM;
log_fpkm_temp[log_fpkm_temp<=pseudo_count] = pseudo_count;
log_fpkm_temp = log10(log_fpkm_temp);
fpkm_transform = log_fpkm_temp;
################ featureConstruct ##############      
#### 1: reading input #### 
# load reference data
load("/data8t_4/JH/MyJobs/Colon_SC_Project/GSE97693_Tang_colorectal_cancer/R_RNA_data/sysdata.rda")
fpkm_temp = fpkm_transform;
#### 2: choosing method ####
## "ColonEpitheliumPanel"
ColonEpiPanel <- sysdata$ColonEpiPanel
# modify the gene names
temp = rownames(ColonEpiPanel)
temp1 = gsub("^.*?_","",temp);
temp2 <- strsplit(temp1,"_ENS");
temp3 <- paste("",lapply(temp2,"[[",1),sep=""); # converted gene names in sysdata$ColonEpiPanel
temp4 = intersect(temp3,rownames(data1))        # gene names in input data and  sysdata$ColonEpiPanel
temp5 = temp[match(temp4,temp3)]; # find the orignal gene names of commone genes 
ColonEpiPanel = ColonEpiPanel[temp5,];
cbind(rownames(ColonEpiPanel),temp3[match(temp4,temp3)])
# convert gene names of ColonEpiPanel
rownames(ColonEpiPanel) <- temp3[match(temp4,temp3)]
# select features(genes) in ColonEpiPanel
fc = apply(ColonEpiPanel,1,function(x) x-median(x)); # normalized by median
fs = fc>1.5;                  # selecte genes over-expreasioned in these cell types??
fs1 = rownames(ColonEpiPanel[apply(fs,1,function(x) sum(x))>0,]);
# gene list of transformed gene expression data 
gl_intersect = intersect(rownames(fpkm_temp),fs1)
data_combine0 = cbind(fpkm_temp[gl_intersect,],ColonEpiPanel[gl_intersect,]);
data_combine = data_combine0;
cor_matrix = cor(data_combine,method = "pearson");
cor_matrix = as.data.frame(cor_matrix);
individual_cell_type_prediction_matrix = cor_matrix[(dim(fpkm_temp)[2]+1):dim(cor_matrix)[2],1:dim(fpkm_temp)[2]];
fpkm_for_clust_0 = individual_cell_type_prediction_matrix;
fpkm_for_clust_power = abs(fpkm_for_clust_0)^power * sign(fpkm_for_clust_0);
fpkm_for_clust = fpkm_for_clust_power;
#### I dont known if these was right or not??
rownames(fpkm_for_clust) <- c("black","Enterocyte type 1B","Stem/TA","Goblet B","Enterocyte type 2",
                              "Goblet C","Enterocyte type 1A","Goblet A","Non-stem 2")

#### Generate cell clusters#########
#### 1: reading input ####  
fpkm_temp = fpkm_for_clust
deepSplit_wgcna=1 
min_group_Size_wgcna=5

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
d = as.dist(1-cor(fpkm_temp,method="pearson"));
cellTree = flashClust(d,method = "average");     
dynamicGroups = cutreeDynamic(dendro = cellTree,distM = as.matrix(d),deepSplit = deepSplit_wgcna,
                              pamStage = FALSE,minClusterSize= min_group_Size_wgcna);

dynamicColors = labels2colors(dynamicGroups);
group_labels = matrix(dynamicGroups,nrow = length(names(fpkm_temp)),ncol = 1,list(names(fpkm_temp),c("groupLabel")),byrow=FALSE)
group_labels = as.data.frame(group_labels)
group_labels_color = cbind(group_labels,dynamicColors);
#  }
############## figure plot ################
## plot PCA
plot(cell_projection[,1],cell_projection[,2],
     type="p",pch = pch_to_use,
     xlab = colnames(cell_projection)[1],ylab = colnames(cell_projection)[2],
     col = c,lwd = 1,bg = c,
     main = main_name, cex = cex_to_use, cex.main = 2,
     font.axis = 1, font.lab = 1, font.main = 1
);

## plot heatmap
library(gplots)
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
## heatmap of gene expression 
markerGenes <- read.table("markerGenes.txt")
markerGenes <- markerGenes[,1]
markerGenes <- as.character(markerGenes)
markerGenes <- rev(markerGenes)
expMatrix <- df.FPKM[markerGenes,]

heatmap.2(as.matrix(expMatrix),
          col=color_scheme,
          Rowv = F,
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
