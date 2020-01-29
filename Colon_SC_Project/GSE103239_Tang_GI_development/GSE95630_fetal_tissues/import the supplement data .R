# import the supplement data 
Esophagus_KNN <- read.csv("Esophagus_KNN_Result.csv")
LIntes_KNN <- read.csv("L-Intes_KNN_Result.csv")
SIntes_KNN <- read.csv("S-Intes_KNN_Result.csv")
Stomach_KNN <- read.csv("Stomach_KNN_Result.csv")
Feathers_organs <- read.csv("Table_1_Features_of_fetal_digestive_organs.csv")
head(Feathers_organs)
Feathers_organs$Cell.Type
table(Stomach_KNN$Cluster)
cellTypeCluster <- list(Esophagus_KNN = Esophagus_KNN,
                        Stomach_KNN = Stomach_KNN,
                        SIntes_KNN = SIntes_KNN,
                        LIntes_KNN = LIntes_KNN,
                        Feathers_organs =Feathers_organs)
saveRDS(cellTypeCluster, file = "GSE95630_Tang_fetal_tissues_cellTypeCluster.rds")








