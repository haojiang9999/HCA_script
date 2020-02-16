#### 2.GSA_COAD_RSEM_normal_count_Log2(x+1).R
library(GSA)
help(GSA)
### 1.Group information
Cluster.20200201.V7.Tumor <- readRDS("/data8t_4/JH/MyJobs/NormalCancer_TCGA_V2/Cluster.20200201.V7.Tumor.rds")
cutree.res <- Cluster.20200201.V7.Tumor$cutree.res
dynamicColors <- Cluster.20200201.V7.Tumor$dynamicColors
table(cutree.res)
table(dynamicColors)
#### 2.Construct GSA data
## V1 blue vs other
x = as.matrix(COAD.RSEM.gene.norm.count.Log2[,names(cutree.res)])
y = dynamicColors
y[y == "brown"] <- 2
y[y != 2] <- 1
genenames = rownames(COAD.RSEM.gene.norm.count.Log2)
genesets = genesets.cms.list
geneset.names =  genesets.cms.names
GSA.COAD.log2.counts.obj<-GSA(x,y, genenames=genenames, genesets=genesets,  
                              resp.type="Two class unpaired", nperms=1000)
GSA.listsets(GSA.COAD.log2.counts.obj, geneset.names=geneset.names,FDRcut=.5)
GSA.COAD.log2.counts.obj$pvalues.lo
GSA.listsets(GSA.COAD.log2.counts.obj, geneset.names=geneset.names,FDRcut=.5)
GSA.plot(GSA.COAD.log2.counts.obj, fac=1, FDRcut = 1)

##### Filter genes all >0 and expressed in at least 3 samples
x2 <- x[rowSums(x) > 3,]
x2 <- x2[rowSums(x2>0) >5,]
genenames2 <- rownames(x2)
dim(x2)
GSA.COAD.log2.counts.obj2<-GSA(x2,y, genenames=genenames2, genesets=genesets,  
                              resp.type="Two class unpaired", nperms=1000)
GSA.listsets(GSA.COAD.log2.counts.obj2, geneset.names=geneset.names,FDRcut=.5)
