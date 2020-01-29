### PCA plot
install.packages("ggfortify")
library("ggfortify")
autoplot(prcomp(t(test_for_clust_power)), data = COAD.pheno.exp,colour = "treatment_outcome_first_course" )
autoplot(prcomp(t(test_for_clust_power)),colour = samplType )
### T-SNE plot
library(ggplot2)
tsne_out <- Rtsne(t(test_for_clust_power), dims = 2, initial_dims = 50, perplexity = 50)
col = COAD.pheno.exp$ajcc_pathologic_tumor_stage
tsne_plot <- data.frame(x = tsne_out$Y[,1], y = tsne_out$Y[,2], 
                        col = col # change phenotypes
                        )
ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=col))

