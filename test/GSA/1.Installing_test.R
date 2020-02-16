#### 1.Installing.R
install.packages("GSA")
library(GSA)
help(GSA)
## example test
set.seed(100)
x<-matrix(rnorm(1000*20),ncol=20)
dd<-sample(1:1000,size=100)
### add the difference
u<-matrix(2*rnorm(100),ncol=10,nrow=100)
x[dd,11:20]<-x[dd,11:20]+u
## Group information
y<-c(rep(1,10),rep(2,10))

genenames=paste("g",1:1000,sep="")
#create some random gene sets
genesets=vector("list",50)
for(i in 1:50){
  genesets[[i]]=paste("g",sample(1:1000,size=30),sep="")
}
geneset.names=paste("set",as.character(1:50),sep="")
GSA.obj<-GSA(x,y, genenames=genenames, genesets=genesets,  resp.type="Two class unpaired", nperms=100)


GSA.listsets(GSA.obj, geneset.names=geneset.names,FDRcut=.5)








