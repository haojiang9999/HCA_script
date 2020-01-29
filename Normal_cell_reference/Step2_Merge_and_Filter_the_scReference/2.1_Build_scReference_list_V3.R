### 2.1_Build_scReference_list_V3.R
#### Step1.Find common genes
GeneNames.Tang.Adult.colon <- rownames(Tang.Adult.colon.ref)
GeneNames.Tang.Fetal.GI <- rownames(Tang.Fetal.GI.ref)
GeneNames.Tang.Normal.embryo <- rownames(Tang.Normal.embryo.ref)
GeneNames.Lanner.Preim.embryo <- rownames(Lanner.Preim.embryo.ref)
#### find the common genes
GeneNames.common <- Reduce(intersect, list(GeneNames.Tang.Adult.colon,GeneNames.Tang.Fetal.GI,
                       GeneNames.Tang.Normal.embryo,GeneNames.Lanner.Preim.embryo))
#### Step2.Build_scReference_list_V3
scReference.list.ComGenes <- list(Tang.Adult.colon = Tang.Adult.colon.ref[GeneNames.common,],
                         Tang.Fetal.GI = Tang.Fetal.GI.ref[GeneNames.common,],
                         Tang.Normal.embryo = Tang.Normal.embryo.ref[GeneNames.common,],
                         Lanner.Preim.embryo = Lanner.Preim.embryo.ref[GeneNames.common,])
dim(scReference.list.ComGenes$Tang.Fetal.GI)
#### Step3.log10(x+1) transformed
scReference.list.ComGenes.log10<- lapply(scReference.list.ComGenes, function(x){
  log10.x <- log10(x+1)
  return(log10.x)
})

### Step4.Filter gene by high varience
source("/data8t_4/JH/MyJobs/1_R_script/FUN_TopCV.R")
scReference.list.ComGenes.log10.CV8000 <- lapply(scReference.list.ComGenes.log10, function(x){
  TopCV(x, TopN = 8000, MARGIN = 1)
})
scReference.list.ComGenes.log10.CV4000 <- lapply(scReference.list.ComGenes.log10, function(x){
  TopCV(x, TopN = 4000, MARGIN = 1)
})
scReference.list.ComGenes.log10.CV2000 <- lapply(scReference.list.ComGenes.log10, function(x){
  TopCV(x, TopN = 2000, MARGIN = 1)
})
scReference.list.ComGenes.log10.CV1500 <- lapply(scReference.list.ComGenes.log10, function(x){
  TopCV(x, TopN = 1500, MARGIN = 1)
})
scReference.list.ComGenes.log10.CV1000 <- lapply(scReference.list.ComGenes.log10, function(x){
  TopCV(x, TopN = 1000, MARGIN = 1)
})


### Step5.Build list and save
scReference.ComGenes.log10.CV <- list(scReference.list.ComGenes.log10 = scReference.list.ComGenes.log10, 
                             scReference.list.ComGenes.log10.CV8000 = scReference.list.ComGenes.log10.CV8000,
                             scReference.list.ComGenes.log10.CV4000 = scReference.list.ComGenes.log10.CV4000, 
                             scReference.list.ComGenes.log10.CV2000 = scReference.list.ComGenes.log10.CV2000,
                             scReference.list.ComGenes.log10.CV1500 = scReference.list.ComGenes.log10.CV1500, 
                             scReference.list.ComGenes.log10.CV1000 = scReference.list.ComGenes.log10.CV1000)
saveRDS(scReference.ComGenes.log10.CV, file = "2020_1_19_scReferenceV3.ComGenes.log10.CV.rds")
