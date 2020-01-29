##### Construct cancer cell group exprtession matrix
### Using CV2000 get cancer cell groups 
# 
### Step1.Construct expression matrix for each cell group
## group1
G1_cells
Cancer.Log2.expList <- list(log2.Tang.colon.cancer.FPKM.500 = log2.Tang.colon.cancer.FPKM.500,
                     log2.Tang.colon.cancer.TPM.688 = log2.Tang.colon.cancer.TPM.688,
                     log2.Tumor.epithelial.COUNT.normal.uniGene = log2.Tumor.epithelial.COUNT.normal.uniGene)
G1.exp<- lapply(Cancer.Log2.expList, function(x){
  x[,colnames(x) %in% G1_cells]
})

G2.exp<- lapply(Cancer.Log2.expList, function(x){
  x[,colnames(x) %in% G2_cells]
})
G3.exp<- lapply(Cancer.Log2.expList, function(x){
  x[,colnames(x) %in% G3_cells]
})
G4.exp<- lapply(Cancer.Log2.expList, function(x){
  x[,colnames(x) %in% G4_cells]
})
G5.exp<- lapply(Cancer.Log2.expList, function(x){
  x[,colnames(x) %in% G5_cells]
})
G6.exp<- lapply(Cancer.Log2.expList, function(x){
  x[,colnames(x) %in% G6_cells]
})

Cancer.group.list <- list(G1.exp = G1.exp, G2.exp = G2.exp, G3.exp = G3.exp,
     G4.exp = G4.exp, G5.exp = G5.exp, G6.exp = G6.exp)
saveRDS(Cancer.group.list, file = "Cancer.group.list.exp.rds")



