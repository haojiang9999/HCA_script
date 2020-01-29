
#### Build_scReference_list_Version2
scReference.list <- list(Tang.Adult.colon = Tang.Adult.colon.ref,
                         Tang.Fetal.GI = Tang.Fetal.GI.ref,
                         Tang.Normal.embryo = Tang.Normal.embryo.ref,
                         Lanner.Preim.embryo = Lanner.Preim.embryo.ref)
#### Step1.log10(x+1) transformed
scReference.list.log10<- lapply(scReference.list, function(x){
  log10.x <- log10(x+1)
  return(log10.x)
})
## Check the transformation
#head(scReference.list$Lanner.Preim.embryo)
#head(scReference.list.log10$Lanner.Preim.embryo)
#log10(2.5759387)
### Step2.Filter gene by high varience
source("/data8t_4/JH/MyJobs/1_R_script/FUN_TopCV.R")
scReference.list.log10.CV8000 <- lapply(scReference.list.log10, function(x){
  TopCV(x, TopN = 8000, MARGIN = 1)
})
scReference.list.log10.CV4000 <- lapply(scReference.list.log10, function(x){
  TopCV(x, TopN = 4000, MARGIN = 1)
})
scReference.list.log10.CV2000 <- lapply(scReference.list.log10, function(x){
  TopCV(x, TopN = 2000, MARGIN = 1)
})
scReference.list.log10.CV1500 <- lapply(scReference.list.log10, function(x){
  TopCV(x, TopN = 1500, MARGIN = 1)
})
scReference.list.log10.CV1000 <- lapply(scReference.list.log10, function(x){
  TopCV(x, TopN = 1000, MARGIN = 1)
})


scReference.log10.CV <- list(scReference.list.log10 = scReference.list.log10, 
                             scReference.list.log10.CV8000 = scReference.list.log10.CV8000,
                             scReference.list.log10.CV4000 = scReference.list.log10.CV4000, 
                             scReference.list.log10.CV2000 = scReference.list.log10.CV2000,
                             scReference.list.log10.CV1500 = scReference.list.log10.CV1500, 
                             scReference.list.log10.CV1000 = scReference.list.log10.CV1000)
saveRDS(scReference.log10.CV, file = "2020_1_19_scReference.log10.CV.rds")








