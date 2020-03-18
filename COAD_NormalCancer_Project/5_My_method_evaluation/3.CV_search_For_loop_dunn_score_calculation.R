#### 3.For_loop_dunn_score_calculation.R
TopNs <- seq(1000,6000,by = 500)

for(i in TopNs){
  ### 2.Find genes had high coefficience-varience
  TopN = i ###########################################
  source("/data8t_4/JH/MyJobs/1_R_script/FUN_TopCV.R")
  scReference.list.TopN <- lapply(scReference.list.V8.GI, function(x){
    TopCV(x, TopN = TopN, MARGIN = 1)
  })
  
  ### 3.log10(x+1) transformed scReference.list.TopN
  log10.scReference.list.TopN <- lapply(scReference.list.TopN, function(x){
    log10.x <- log10(x+1)
    return(log10.x)
  })
  
  
  ### 4.Correlation calculation 
  source("/data8t_4/JH/MyJobs/1_R_script/NormalCancer/refCorMerge.R")
  Cor.Res.TopN <- refCorMerge(COAD_tpm.Log10.expList, scReference.list.TopN)
  Cor.merged.TopN <- Cor.Res.TopN$Cor.merged
  COAD.pheno <- COAD.pheno
  #Cor.merged.TopN
  
  ### 5.Correlation normalization steps
  ## 1)The first Normalization
  #normal the correlation across samples
  Trans.Rang1.TopN<- base::apply(Cor.merged.TopN, 1, function(x){
    (x-min(x))/(max(x)-min(x))
  })
  Trans.Rang1.TopN <- t(Trans.Rang1.TopN)
  
  ## 2)The second Normalization
  # Normalize the correlation across each cell type cluster; It's the weight of each cell type cluster
  cellTypesName <- rownames(Trans.Rang1.TopN)
  #grep("Adult", cellTypesName)
  Trans.Rang1.TopN.Normalized <- Trans.Rang1.TopN
  Trans.Rang1.TopN.Normalized[grep("Adult", cellTypesName),]<- base::apply(Trans.Rang1.TopN.Normalized[grep("Adult", cellTypesName),], 2, function(x){
    (x-min(x))/(max(x)-min(x))
  })
  Trans.Rang1.TopN.Normalized[grep("Fetal", cellTypesName),]<- base::apply(Trans.Rang1.TopN.Normalized[grep("Fetal", cellTypesName),], 2, function(x){
    (x-min(x))/(max(x)-min(x))
  })
  
  ### 6.Distance calculation
  Trans.Rang1.TopN.Normalized
  TopN.Normalized.dist <- dist(t(Trans.Rang1.TopN.Normalized), method="euclidean")
  
  ### 7.Calculate the Dunn
  library(clValid)
  #cbind(rownames(COAD.pheno),colnames(Trans.Rang1.TopN.Normalized))
  NTcluster <- COAD.pheno$sampleTypes
  NTcluster <- as.numeric(NTcluster)
  #table(NTcluster)
  names(NTcluster) <- rownames(COAD.pheno)
  dunn.score <- dunn(TopN.Normalized.dist,NTcluster)
  print(dunn.score)
}






