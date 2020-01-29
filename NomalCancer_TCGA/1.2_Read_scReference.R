## Read scReference 
#### Reading data
### 1.Read the reference scRNA data
###
scReference.V1 <- readRDS("/data8t_4/JH/MyJobs/Normal_cell_reference/Step2_Merge_and_Filter_the_scReference/2019_11_16_scReference.V1.rds")
scReference.list <- scReference.V1$scReference.list

### 2.Reference data transformation
#### Reference transformation 
#### log2(x+1) transformed the same as COAD data
log2.scReference.list.CV.8000 <- lapply(scReference.V1[["scReference.list.CV.8000"]], function(x){
  log2.x <- log2(x + 1)
  return(log2.x)
})

log2.scReference.list.CV.4000 <- lapply(scReference.V1[["scReference.list.CV.4000"]], function(x){
  log2.x <- log2(x + 1)
  return(log2.x)
})
log2.scReference.list.CV.2000 <- lapply(scReference.V1[["scReference.list.CV.2000"]], function(x){
  log2.x <- log2(x + 1)
  return(log2.x)
})

log2.scReference.list.CV.1500 <- lapply(scReference.V1[["scReference.list.CV.1500"]], function(x){
  log2.x <- log2(x + 1)
  return(log2.x)
})

log2.scReference.list.CV.1000 <- lapply(scReference.V1[["scReference.list.CV.1000"]], function(x){
  log2.x <- log2(x + 1)
  return(log2.x)
})

