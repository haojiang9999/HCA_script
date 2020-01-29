#### COAD RNA-seq result
COAD_tpm<- COAD_UCSC_Toil_tpm_dataset$COAD_RSEM_gene_tpm
head(COAD_tpm)
#### t-SNE
library(Rtsne)
t(COAD_tpm)
TSNE <- Rtsne(t(COAD_tpm), dims = 2, initial_dims = 50, perplexity = 50)
tSNEdata.used <- as.matrix(TSNE$Y)
plot(tSNEdata.used[,1],tSNEdata.used[,2]) ## OK
plot(tSNEdata.used[,1],tSNEdata.used[,2], col = as.factor(samplType)) ## OK

#### PCA
RCA.PCA(COAD_tpm, color = as.factor(samplType))
RCA.PCA(COAD_tpm, color = as.factor(tumor_stage))
