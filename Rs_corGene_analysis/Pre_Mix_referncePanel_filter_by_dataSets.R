#### reference panel filter
referPanel_of_mix
dim(referPanel_of_mix)
# calculate CV of each genes
CV <- function(x){
  (sd(x)/mean(x))*100
}
genesCV <- apply(referPanel_of_mix, 1, CV)
summary(genesCV)
tail(sort(genesCV, decreasing = T))
sort(genesCV, decreasing = T)[1:100]
sort(genesCV, decreasing = T)[1:100]
colSums(referPanel_of_mix)
grepl("Embryo",colnames(referPanel_of_mix))


### find most variable genes
CV <- function(x){
  (sd(x)/mean(x))*100
}

GI_genesCV <- apply(GI_ref, 1, CV)
summary(GI_genesCV)
GI_8000_genes <- head(sort(GI_genesCV,decreasing = T), 8000)

test <- log10(GI_ref)
RCAvariableGenes <- function(x){
  (x/median(x))
}
(test[1,])/(median(test[1,]))

test_2 <- apply(test[1,], 1,RCAvariableGenes)
head(test_2)

