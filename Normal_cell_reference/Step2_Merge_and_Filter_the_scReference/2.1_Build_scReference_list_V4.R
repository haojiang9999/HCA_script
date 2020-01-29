### 2.1_Build_scReference_list_V4.R
#### Build_scReference_list_Version2
scReference.list <- list(Tang.Adult.colon = Tang.Adult.colon.ref,
                         Tang.Fetal.GI = Tang.Fetal.GI.ref,
                         Tang.Normal.embryo = Tang.Normal.embryo.ref,
                         Lanner.Preim.embryo = Lanner.Preim.embryo.ref,
                         Zemin.CRC.Tcell = Zemin.CRC.Tcell.ref,
                         Zemin.CRC.Tcell.Cluster = Zemin.CRC.Tcell.Cluster.ref)

### Step1.Find genes had high coefficience-varience
source("/data8t_4/JH/MyJobs/1_R_script/FUN_TopCV.R")
scReference.list.CV8000 <- lapply(scReference.list, function(x){
  TopCV(x, TopN = 8000, MARGIN = 1)
})
scReference.list.CV4000 <- lapply(scReference.list, function(x){
  TopCV(x, TopN = 4000, MARGIN = 1)
})
scReference.list.CV3000 <- lapply(scReference.list, function(x){
  TopCV(x, TopN = 3000, MARGIN = 1)
})
scReference.list.CV2500 <- lapply(scReference.list, function(x){
  TopCV(x, TopN = 2500, MARGIN = 1)
})
scReference.list.CV2000 <- lapply(scReference.list, function(x){
  TopCV(x, TopN = 2000, MARGIN = 1)
})
scReference.list.CV1500 <- lapply(scReference.list, function(x){
  TopCV(x, TopN = 1500, MARGIN = 1)
})
scReference.list.CV1000 <- lapply(scReference.list, function(x){
  TopCV(x, TopN = 1000, MARGIN = 1)
})


#### Step2.log10(x+1) transformed
log10.scReference.list <- lapply(scReference.list, function(x){
  log10.x <- log10(x+1)
  return(log10.x)
})
log10.scReference.list.CV.8000 <- lapply(scReference.list.CV8000, function(x){
  log10.x <- log10(x+1)
  return(log10.x)
})

log10.scReference.list.CV.4000 <- lapply(scReference.list.CV4000, function(x){
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

log10.scReference.list.CV.1500 <- lapply(scReference.list.CV1500, function(x){
  log10.x <- log10(x+1)
  return(log10.x)
})

log10.scReference.list.CV.1000 <- lapply(scReference.list.CV1000, function(x){
  log10.x <- log10(x+1)
  return(log10.x)
})
## Check the transformation
#head(scReference.list$Lanner.Preim.embryo)
#head(scReference.list.log10$Lanner.Preim.embryo)
#log10(2.5759387)

scReference.log10.CV.V4 <- list(scReference.list.log10 = log10.scReference.list, 
                             scReference.list.log10.CV8000 = log10.scReference.list.CV.8000,
                             scReference.list.log10.CV4000 = log10.scReference.list.CV.4000,
                             scReference.list.log10.CV3000 = log10.scReference.list.CV.3000,
                             scReference.list.log10.CV2500 = log10.scReference.list.CV.2500,
                             scReference.list.log10.CV2000 = log10.scReference.list.CV.2000,
                             scReference.list.log10.CV1500 = log10.scReference.list.CV.1500, 
                             scReference.list.log10.CV1000 = log10.scReference.list.CV.1000)
saveRDS(scReference.log10.CV.V4, file = "2020_1_19_scReferenceV5.log10.CV.rds")

