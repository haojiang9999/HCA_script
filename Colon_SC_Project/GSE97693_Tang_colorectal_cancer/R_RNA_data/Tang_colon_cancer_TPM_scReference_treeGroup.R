#### construct a single cell reference set for scBio deconvolution
# acoording to the cluster resaults
table(clusterResault[[2]]$dynamicColors)
table(clusterResault[[2]]$dynamicColors %in% c("turquoise","blue","brown"))
choosedSamplesIndex <- clusterResault[[2]]$dynamicColors %in% c("turquoise","blue","brown")
df.TPM.ref <- df.TPM[,choosedSamplesIndex]
clusterResault.ref <- clusterResault[[2]]$dynamicColors[choosedSamplesIndex]
pca_out.ref <- pca_out[choosedSamplesIndex,]
# construct reference detaset
Tang_GSE97693_scRef_set <- list(df.TPM.ref = df.TPM.ref,
                                clusterResault.ref = clusterResault.ref,
                                pca_out.ref = pca_out.ref)
saveRDS(Tang_GSE97693_scRef_set, file = "Tang_GSE97693_scRef_set_three_group.rds")
