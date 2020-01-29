#### Ordered by cluster resaults plotChromosomeHeatmap plot
#pre-set
colo=NULL
colo2=NULL
expThresh=0.4
thresh=1
gene_pos = gene_pos
windowsize=121 # genes
chr=T
mat <- suva_expr
#cbind(colnames(mat),Tang_TPM_688.anno[,"sampleName"])
normal= which(Tang_TPM_688.anno[,"cellSites"] == "NC")
tumor=which(Tang_TPM_688.anno[,"cellSites"] != "NC")
plotcells = c(normal,tumor)
#Create average of reference (for all non-tumor cells, even if not plotted)
ref=rowMeans(mat[,colnames(mat)[normal]])
#Define matrix to be plotted (remove not plotted normal and tumor cells)
normal=intersect(normal,plotcells)
tumor=setdiff(plotcells,normal)
gexp=mat[,c(normal,tumor)] ## In here separate normal and cancer cells
#Filter matrices
ref=ref[which(ref>expThresh)]
gexp=gexp[names(ref),]
gexp=gexp[which(rowMeans(gexp)>expThresh),]
ref=ref[rownames(gexp)]
#Get sorted gene positions
#gp=getGenePositions(rownames(gexp))
gp=gene_pos
gp=gp[which(gp[,"chromosome_name"] %in% 1:22),]
gp=gp[order(as.numeric(gp[,"chromosome_name"]),as.numeric(gp[,"start_position"])),]
#Reduce matrix to genes with known positions
gexp=gexp[intersect(gp[,2],rownames(gexp)),]
ref=ref[rownames(gexp)]
#Calculate ratios
rat=gexp-ref
#Smoothed expression
print("Generating distance Matrix")
d=apply (rat,2,function (x) zoo::rollapply(x,mean,width=windowsize,align="center"))
#Center in each cell and set to boundaries
d=apply(d,2,function(x) x-mean(x))
#Set min/max 
d[which(d>thresh)]=thresh
d[which(d<(-thresh))]=(-thresh)
# Order by chromosome
# And order the tumor samples by AnnoRes cluster
normal=1:length(normal)
tumor=(length(normal)+1):ncol(d)
#dmat=as.dist(pdist(t(d[,tumor]))) # cluster tumors by chromosome region
#hc = hclust(dmat)
#### In here insert my group resaults ####
samples.anno.df.ordered <- samples.anno.df[order(samples.anno.df$samples.group),]
samples.ordered <- samples.anno.df.ordered$samples

#cellOrder = c(normal,tumor[hc$order])
d=d[,samples.ordered]

####### plot #########
dsize=800/length(plotcells)
i=0
#dev.off()
plot.new()
apply(d,2,function(x){
  res=x
  #Counter increase
  i <<- i + 1
  map = squash::makecmap(c(-thresh,thresh), colFn = squash::darkbluered,n=256)
  #map = squash::makecmap(c(max(-thresh,mind),min(thresh,maxd)), colFn = squash::darkbluered,n=256)
  if (i==1){
    plot(c(1:length(res)),rep(1,length(res)),col=squash::cmap(res, map = map),cex=dsize,pch=".",xaxt='n',yaxt='n',ann=FALSE,ylim=c(0,ncol(d)))
    if (!is.null(colo)){
      points((length(res)+1):(length(res)+20),rep(1,20),col="white",pch=".",cex=dsize)
      points((length(res)+21):(length(res)+110),rep(1,90),col=colo[i],pch=".",cex=dsize)
    }
    if (!is.null(colo2)){
      points((length(res)+111):(length(res)+160),rep(1,50),col="white",pch=".",cex=dsize)
      points((length(res)+161):(length(res)+250),rep(1,90),col=colo2[i],pch=".",cex=dsize)
    }
    squash::hkey(map, title="",stretch = 0.1)
  }
  else{
    points(c(1:length(res)),rep(i,length(res)),col=squash::cmap(res, map = map),pch=".",cex=dsize)
    if (!is.null(colo)){
      points((length(res)+1):(length(res)+20),rep(i,20),col="white",pch=".",cex=dsize)
      points((length(res)+21):(length(res)+110),rep(i,90),col=colo[i],pch=".",cex=dsize)
    }
    if (!is.null(colo2)){
      points((length(res)+111):(length(res)+160),rep(i,50),col="white",pch=".",cex=dsize)
      points((length(res)+161):(length(res)+250),rep(i,90),col=colo2[i],pch=".",cex=dsize)
    }
  }
})

bps=c(0)
for (i in 1:21){
  gp=gp[which(gp[,2] %in% rownames(gexp)),]
  bp1=which(gp[,"chromosome_name"]==i+1)[1]-windowsize/2
  bps=c(bps,bp1)
  abline(v=bp1,lty=16)
}
axis(1, at=bps[2:23]-((bps[2:23]-bps[1:22])/2), labels=1:22, las=2,col="grey")
bps=round(bps,0)
abline(h=length(normal),lty=16)
if (retMat){
  chrpos=c()
  for (i in 2:22){chrpos=c(chrpos,rep(i-1,bps[i]-bps[i-1]))};chrpos=c(chrpos,rep(22,nrow(d)-bps[22]))
  return(rbind(chrpos,t(d)))
}
