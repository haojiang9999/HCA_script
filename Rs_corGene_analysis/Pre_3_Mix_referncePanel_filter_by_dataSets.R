####calculate CV for different datasets ####
### sseperate dtasets
GI_ref <- referPanel_of_mix[,1:40]
embryo_ref <- referPanel_of_mix[,41:44]
### fatal and adults data sets
# Coefficient Variance function
CV <- function(x){
  (sd(x)/mean(x))*100
}
## Step1 filtered low expression genes
# genes 
#GI_ref[summary(rowSums(GI_ref)) == 4,]
#GI_ref[rowSums(GI_ref) < 5,]
#summary(t(GI_ref[rowSums(GI_ref) < 5,]))
## Step2 calculate CV for each genes
GI_genesCV <- apply(GI_ref, 1, CV)
summary(GI_genesCV)
## Step3 find the genes with largest cv value Top 5000
GI_5000_genes <- head(sort(GI_genesCV,decreasing = T), 5000)
GI_ref[c("OR10G2", "LOC100287010", "CRNN"),] # check top tree gene expression
tail(GI_5000_genes)

########### embryo data sets ########
## Step1 filtered low expression genes
table(rowSums(embryo_ref) > 1) # how many low expression genes
embryo_ref_filter1 <- embryo_ref[rowSums(embryo_ref) > 1,]
## Step2 calculate CV for each genes
embryo_genesCV <- apply(embryo_ref_filter1, 1, CV)
summary(embryo_genesCV)
## Step3 find the genes with largest cv value Top 5000
embryo_5000_genes <- head(sort(embryo_genesCV,decreasing = T), 5000)
# Check top genes expression
embryo_ref_filter1[c("LINC01208", "RFPL2", "KLK11", "CXCL1","CXCL6","CXorf66"),]


###### re-combine the reference panel #########
filtered_geneList <- unique(c(names(GI_5000_genes),names(embryo_5000_genes)))
Filtered_referPanel_of_mix <- referPanel_of_mix[filtered_geneList,]
# save data sets
saveRDS(Filtered_referPanel_of_mix, file = "Filtered_referPanel_of_mix_V1.rds")




