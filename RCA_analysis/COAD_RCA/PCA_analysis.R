color = as.factor(samplType)
matrix = test_for_clust_power
pc1 = 3
pc2 = 4


## PCA analysis
pr = prcomp(t(scale(matrix))); 
pc_projection = as.data.frame(pr$x);
cell_projection = pc_projection;
############## figure plot ################
## plot PCA

plot(cell_projection[,pc1],cell_projection[,pc2],
     type="p",pch = 21,
     xlab = colnames(cell_projection)[1],ylab = colnames(cell_projection)[2],
     col = as.character(color),lwd = 1,bg = as.character(color),
     main = "PCA of cell clusters in RCA space", cex = 1, cex.main = 2,
     font.axis = 1, font.lab = 1, font.main = 1);




