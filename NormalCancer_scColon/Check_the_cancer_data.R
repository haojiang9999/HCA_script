#Check the cancer data
Tang.colon.cancer.FPKM.500
Tang.colon.cancer.FPKM.500[1:10,1:10]
colnames(Tang.colon.cancer.FPKM.500)
Tang.colon.cancer.TPM.688[1:10,1:10]
RCA.Tumor.epithelial.COUNT.normal[1:10,1:10]
### How many Non-zero genes in each sample
summary(apply(Tang.colon.cancer.FPKM.500, MARGIN = 2 ,FUN = function(x){
  sum(x != 0)
}))
summary(apply(Tang.colon.cancer.TPM.688, MARGIN = 2 ,FUN = function(x){
  sum(x != 0)
}))
summary(apply(RCA.NM.epithelial.COUNT.normal, MARGIN = 2 ,FUN = function(x){
  sum(x != 0)
}))
summary(apply(RCA.Tumor.epithelial.COUNT.normal, MARGIN = 2 ,FUN = function(x){
  sum(x != 0)
}))















