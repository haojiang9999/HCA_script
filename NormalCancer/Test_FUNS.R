# test for PreForEnAnno
test.PreForEnAnno <- PreForEnAnno(test2, GeneBackGround_ENTREZID)
# subset example
test.PreForEnAnno.sub <- test.PreForEnAnno[1:10]

time.anno <- system.time(test <- mclapply(test.PreForEnAnno.sub, function(x){
  EnAnno(x$geneMerge, x$geneBackGround)
}, mc.cores = mc.cores))
