#### kmeans analysis
kmeangroup <- kmeans(t(test_for_clust_power),5)
plot(tSNEdata.used[,1],tSNEdata.used[,2], col = as.factor(kmeangroup$cluster)) ## OK
### survival analysis
library(survminer)
require("survival")
COAD.sur.df.more <- cbind(COAD.sur.df, kmeangroup = kmeangroup$cluster)
fit <- survfit(Surv(OS.time, OS) ~ kmeangroup, data = COAD.sur.df.more)
# Drawing curves
ggsurvplot(fit)

#### hcluster


#### RCA.cluster
clusterResault <- RCA.cluster(test_for_clust_power, deepSplit_wgcna=0.00001, min_group_Size_wgcna=1)
table(clusterResault[[2]]$dynamicColors)
plot(tSNEdata.used[,1],tSNEdata.used[,2], col = as.factor(clusterResault[[2]]$dynamicColors)) ## OK
###
#### cluster by sample type
plot(tSNEdata.used[,1],tSNEdata.used[,2], col = as.factor(samplType)) ## OK
plot(tSNEdata.used[,1],tSNEdata.used[,2], col = as.factor(tumor_stage)) ## OK
plot(tSNEdata.used[,1],tSNEdata.used[,2], col = as.factor(COAD.pheno.exp$treatment_outcome_first_course)) ## OK

table(COAD.pheno.exp$treatment_outcome_first_course)

# RCA
RCA.PCA(test_for_clust_power, color = as.factor(samplType))
RCA.PCA(test_for_clust_power, color = as.factor(tumor_stage))
RCA.PCA(test_for_clust_power, color = as.factor(COAD.pheno.exp$treatment_outcome_first_course))

