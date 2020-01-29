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
## "GlobalPanel"
data=fpkm_temp;
data11=list();
power = 4
for (i in 1:length(sysdata[[1]])){
  data1 = sysdata[[1]][[i]];
  d1 = as.dist(1-cor(data1));
  t1 = hclust(d1,method="average");
  temp = rownames(data);
#  temp1 = gsub("^.*?_","",temp);
#  temp2 <- strsplit(temp1,"_ENS");
#  temp3 <- paste("",lapply(temp2,"[[",1),sep="");
  temp4 = intersect(temp,rownames(data1));  # for the gene name duplicated issue they find intercect names first
  temp5 = temp[match(temp4,temp3)];
  data2 = data1[temp4,];
  data4 = data[temp5,,drop=FALSE];
  data3 = data2;
  data3[data3<=(sysdata$at)[i]] = (sysdata$at)[i];
  data5 = cbind(data4,data3);
  data6 = cor(data5,method = "pearson");
  data6 = as.data.frame(data6);
  data7 = data6[(dim(data)[2]+1):dim(data6)[2],1:dim(data)[2]];
  data8 = data7;
  #data9 = abs(data8)^(sysdata$sp) * sign(data8);
  data9 = abs(data8)^(power) * sign(data8);
  data10 = scale(data9,center=TRUE,scale=TRUE);
  data11[[i]]=data10;
} 
## combines two matrix  
data12 = as.data.frame(t(cbind(as.data.frame(t(data11[[1]])),as.data.frame(t(data11[[2]])))));    
fpkm_for_clust = data12;  
  

#### Generate cell clusters#########
#### 1: reading input ####  
fpkm_temp = data10;
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
    
