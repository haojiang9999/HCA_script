#### Build separate scReference table
#### scReference filter
### 1.Using genes that common in different dataset
GeneNames.common
scReference.list <- list(Tang.Adult.colon = Tang.Adult.colon.ref[GeneNames.common,],
                         Tang.Fetal.GI = Tang.Fetal.GI.ref[GeneNames.common,],
                         Tang.Normal.embryo = Tang.Normal.embryo.ref[GeneNames.common,],
                         Lanner.Preim.embryo = Lanner.Preim.embryo.ref[GeneNames.common,])
### 2.Filter gene by high varience
TopN = 8000
scReference.list.CV.8000 <- lapply(scReference.list, FUN = function(ref){
  # Genes Coefficient-Variance(CV)
  # ref = scReference.list[[1]]
  genesCV <- apply(ref, 1, function(x) {sd(x) / mean(x)})
  Genes.TopN <- names(head(sort(genesCV, decreasing = T), TopN))
  ref[Genes.TopN,]
})
scReference.list.CV.8000[[1]]
colnames(scReference.list.CV.8000[[1]])
rownames(scReference.list.CV.8000[[1]])
scReference.list[[1]]

TopN = 4000
scReference.list.CV.4000 <- lapply(scReference.list, FUN = function(ref){
  # Genes Coefficient-Variance(CV)
  # ref = scReference.list[[1]]
  genesCV <- apply(ref, 1, function(x) {sd(x) / mean(x)})
  Genes.TopN <- names(head(sort(genesCV, decreasing = T), TopN))
  ref[Genes.TopN,]
})

TopN = 2000
scReference.list.CV.2000 <- lapply(scReference.list, FUN = function(ref){
  # Genes Coefficient-Variance(CV)
  # ref = scReference.list[[1]]
  genesCV <- apply(ref, 1, function(x) {sd(x) / mean(x)})
  Genes.TopN <- names(head(sort(genesCV, decreasing = T), TopN))
  ref[Genes.TopN,]
})

TopN = 1500
scReference.list.CV.1500 <- lapply(scReference.list, FUN = function(ref){
  # Genes Coefficient-Variance(CV)
  # ref = scReference.list[[1]]
  genesCV <- apply(ref, 1, function(x) {sd(x) / mean(x)})
  Genes.TopN <- names(head(sort(genesCV, decreasing = T), TopN))
  ref[Genes.TopN,]
})


TopN = 1000
scReference.list.CV.1000 <- lapply(scReference.list, FUN = function(ref){
  # Genes Coefficient-Variance(CV)
  # ref = scReference.list[[1]]
  genesCV <- apply(ref, 1, function(x) {sd(x) / mean(x)})
  Genes.TopN <- names(head(sort(genesCV, decreasing = T), TopN))
  ref[Genes.TopN,]
})

scReference.V1 <- list(scReference.list = scReference.list, scReference.list.CV.8000 = scReference.list.CV.8000,
     scReference.list.CV.4000 = scReference.list.CV.4000, scReference.list.CV.2000 = scReference.list.CV.2000,
     scReference.list.CV.1500 = scReference.list.CV.1500, scReference.list.CV.1000 = scReference.list.CV.1000)
saveRDS(scReference.V1, file = "2019_11_16_scReference.V1.rds")




