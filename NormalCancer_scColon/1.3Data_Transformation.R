#### Data transformation ####
#### log2 (M + 1)
log2.Tang.colon.cancer.FPKM.500 <- log2(Tang.colon.cancer.FPKM.500 + 1)
log2.Tang.colon.cancer.TPM.688 <- log2(Tang.colon.cancer.TPM.688 + 1)
log2.NM.epithelial.COUNT.normal.uniGene <- log2(NM.epithelial.COUNT.normal.uniGene + 1)
log2.Tumor.epithelial.COUNT.normal.uniGene <- log2(Tumor.epithelial.COUNT.normal.uniGene + 1)
## build a list
Log2.expList <- list(log2.Tang.colon.cancer.FPKM.500 = log2.Tang.colon.cancer.FPKM.500,
     log2.Tang.colon.cancer.TPM.688 = log2.Tang.colon.cancer.TPM.688,
     log2.NM.epithelial.COUNT.normal.uniGene = log2.NM.epithelial.COUNT.normal.uniGene,
     log2.Tumor.epithelial.COUNT.normal.uniGene = log2.Tumor.epithelial.COUNT.normal.uniGene)


#### Reference transformation
#### log2 (M + 1)
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




