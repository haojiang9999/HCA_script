### 2.1_Build_scReference_list_V5_Subset_Features.R
### This vesrion was to select reference data that I think better represent the nature 
# of the cancer cells
## Remove 
# 1.Tang and Lanner embryo data 
# 2.Adult_Immune immuno cell from Tang.Adult.colon.ref
#### Build_scReference_list
scReference.list.V5.subset <- list(Tang.Adult.colon = Tang.Adult.colon.ref[,-which(names(Tang.Adult.colon.ref) %in% "Adult_Immune")],
                         Tang.Fetal.GI = Tang.Fetal.GI.ref,
                         Zemin.CRC.Tcell = Zemin.CRC.Tcell.ref,
                         Zemin.CRC.Tcell.Cluster = Zemin.CRC.Tcell.Cluster.ref)
### Step1.Find genes had high coefficience-varience
source("/data8t_4/JH/MyJobs/1_R_script/FUN_TopCV.R")
scReference.list.CV8000 <- lapply(scReference.list.V5.subset, function(x){
  TopCV(x, TopN = 8000, MARGIN = 1)
})
scReference.list.CV4500 <- lapply(scReference.list.V5.subset, function(x){
  TopCV(x, TopN = 4500, MARGIN = 1)
})
scReference.list.CV4000 <- lapply(scReference.list.V5.subset, function(x){
  TopCV(x, TopN = 4000, MARGIN = 1)
})
scReference.list.CV3500 <- lapply(scReference.list.V5.subset, function(x){
  TopCV(x, TopN = 3500, MARGIN = 1)
})
scReference.list.CV3000 <- lapply(scReference.list.V5.subset, function(x){
  TopCV(x, TopN = 3000, MARGIN = 1)
})
scReference.list.CV2500 <- lapply(scReference.list.V5.subset, function(x){
  TopCV(x, TopN = 2500, MARGIN = 1)
})
scReference.list.CV2000 <- lapply(scReference.list.V5.subset, function(x){
  TopCV(x, TopN = 2000, MARGIN = 1)
})
#### Step2.log10(x+1) transformed
log10.scReference.list.CV.8000 <- lapply(scReference.list.CV8000, function(x){
  log10.x <- log10(x+1)
  return(log10.x)
})

log10.scReference.list.CV.4500 <- lapply(scReference.list.CV4500, function(x){
  log10.x <- log10(x+1)
  return(log10.x)
})
log10.scReference.list.CV.4000 <- lapply(scReference.list.CV4000, function(x){
  log10.x <- log10(x+1)
  return(log10.x)
})
log10.scReference.list.CV.3500 <- lapply(scReference.list.CV3500, function(x){
  log10.x <- log10(x+1)
  return(log10.x)
})
log10.scReference.list.CV.3000 <- lapply(scReference.list.CV3000, function(x){
  log10.x <- log10(x+1)
  return(log10.x)
})
log10.scReference.list.CV.2500 <- lapply(scReference.list.CV2500, function(x){
  log10.x <- log10(x+1)
  return(log10.x)
})
log10.scReference.list.CV.2000 <- lapply(scReference.list.CV2000, function(x){
  log10.x <- log10(x+1)
  return(log10.x)
})
#### Step3.Save the reference list
scReference.log10.CV.V5.Subset.Features <- list(scReference.list.log10.CV8000 = log10.scReference.list.CV.8000,
                                                scReference.list.log10.CV4500 = log10.scReference.list.CV.4500,
                                                scReference.list.log10.CV4000 = log10.scReference.list.CV.4000,
                                                scReference.list.log10.CV3500 = log10.scReference.list.CV.3500,
                                                scReference.list.log10.CV3000 = log10.scReference.list.CV.3000,
                                                scReference.list.log10.CV2500 = log10.scReference.list.CV.2500,
                                                scReference.list.log10.CV2000 = log10.scReference.list.CV.2000
                                                )
saveRDS(scReference.log10.CV.V5.Subset.Features, file = "2020_1_29_scReferenceV5.Subset.Features.log10.CV.rds")
