### modify of plotHistogram ###
#Pre-set
pmat <- l
expmat <- suva_expr
zscoreThreshold=4
patients=NULL
celltypes=NULL
clusters=2
#function (pmat,expmat,clusters,zscoreThreshold=4,patients=NULL,celltypes=NULL){
  t=zscoreThreshold
  pmat=scale(pmat)
  if (max(pmat)>t){
    pmat[which(pmat>t)]=t
    pmat[which(pmat<(-t))]=(-t)
  }
  else {
    mx=min(max(pmat),abs(min(pmat)))
    sc=t/mx
    pmat=pmat*sc
    pmat[which(pmat>t)]=t
    pmat[which(pmat<(-t))]=(-t)
  }
  if(is.null(patients) & is.null(celltypes)){
    p=pheatmap::pheatmap(t(pmat),cluster_rows=F, cutree_cols = clusters, col=squash::bluered(100),gaps_col=50,show_colnames = F,clustering_distance_cols="euclidean")
  }
  else if (is.null(celltypes)){
    patientcolors =data.frame(patients)
    rownames(patientcolors)=colnames(expmat)
    pmat=pmat[colnames(expmat),]
    p=pheatmap::pheatmap(t(pmat),cluster_rows=F,silent=F, cutree_cols = clusters, col=squash::bluered(100),gaps_col=50,annotation=patientcolors,show_colnames = F,clustering_distance_cols="euclidean")
  }
  else if(is.null(patients)){
    patientcolors =data.frame(celltypes)
    rownames(patientcolors)=colnames(expmat)
    pmat=pmat[colnames(expmat),]
    p=pheatmap::pheatmap(t(pmat),cluster_rows=F, cutree_cols = clusters, col=squash::bluered(100),gaps_col=50,annotation=patientcolors,show_colnames = F,clustering_distance_cols="euclidean")
  }
  else {
    patientcolors =data.frame(celltypes)
    patientcolors=cbind(patientcolors,patients)
    rownames(patientcolors)=colnames(expmat)
    pmat=pmat[colnames(expmat),]
    p=pheatmap::pheatmap(t(pmat),cluster_rows=F, cutree_cols = clusters, col=squash::bluered(100),gaps_col=50,annotation=patientcolors,show_colnames = F,clustering_distance_cols="euclidean")
  }
  ord=unique(cutree(p$tree_col, k = clusters)[p$tree_col[["order"]]])
  numb=table(cutree(p$tree_col, k = clusters))[ord]
  n=length(numb)
  grid::grid.text(expression(bold("Cluster ID \n(left to right)")),x=rep(0.92),y=c(n*0.03+0.03),gp=grid::gpar(fontsize=8, col="grey"))
  grid::grid.text(ord,x=rep(0.92,length(numb)),y=seq(n*0.03, 0.03, -0.03),gp=grid::gpar(fontsize=8, col="grey"))
  return(cutree(p$tree_col, k = clusters))
}