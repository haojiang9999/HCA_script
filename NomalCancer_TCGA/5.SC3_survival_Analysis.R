##### TCGA Survaival analysis base on SC3 cluster
### Step1.Cut the tree
### Using hcluster
test<- heatmap.JH(Cor.tumor.cv20,show_colnames = F,
           annotation_col = Pheno.merged.tumor)
plot(as.dendrogram(test$tree_col))  
COAD.cluster <- test$tree_col
COAD.cluster <- cutree(COAD.cluster, h = 0.6)

table(COAD.cluster)
heatmap.JH(Cor.tumor.cv20,show_colnames = F,
           annotation_col = cbind(Pheno.merged.tumor,as.character(COAD.cluster)))
### Using SC3 resualt
hcTumor.cv20
heatmap.JH(Cor.tumor.cv10,show_colnames = F,
           annotation_col = Pheno.merged.tumor, cluster_cols = hcTumor.cv10)
plot(as.dendrogram(hcTumor.cv10)) 
saveRDS(hcTumor.cv10, file = "hcTumor.cv10.rds")
### cut the tree
COAD.cluster <- cutree(hcTumor.cv10, h = 8)
table(COAD.cluster)
COAD.sur.df <- cbind(Pheno.merged.tumor,as.character(COAD.cluster))
heatmap.JH(Cor.tumor.cv10,show_colnames = F,cluster_cols = hcTumor.cv10,
           annotation_col =COAD.sur.df,
           scale = "column")

heatmap.JH(expr_test1[,x],show_colnames = F,cluster_cols = hcTumor.cv10,
           annotation_col =COAD.sur.df,
           scale = "column")
colnames(expr_test1)
x <- colnames(Cor.tumor.cv10)
cbind(colnames(expr_test1),colnames(Cor.tumor.cv10))
saveRDS(COAD.sur.df, file = "COAD.sur.df.clusterRes.rds")
### Step2.Build the survival table
sampleID <- names(COAD.cluster)
surTime <- COAD.pheno[sampleID,c("OS","OS.time")]
COAD.sur.df <- cbind(COAD.cluster,surTime)
head(surTime)

### Step3.Survival analysis
library(survminer)
require("survival")
fit <- survfit(Surv(OS.time, OS) ~ COAD.cluster, data = COAD.sur.df)
# Drawing curves
ggsurvplot(fit,pval = T,pval.method = T)

table(COAD.cluster)
COAD.pheno




