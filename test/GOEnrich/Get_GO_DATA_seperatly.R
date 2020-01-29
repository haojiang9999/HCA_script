library(clusterProfiler)
library(GOSemSim)
library(org.Hs.eg.db)
library(GO.db)
library(tidyr)

######## 1. GO_CC term ######
keyType = "ENTREZID"
OrgDb = "org.Hs.eg.db"
ont="CC"
GO_DATA_CC <- get_GO_data(OrgDb, ont, keyType)

saveRDS(GO_DATA_CC, file = "GO_DATA_CC.rds")
######## 2. GO_MF term ######
keyType = "ENTREZID"
OrgDb = "org.Hs.eg.db"
ont="MF"
GO_DATA_MF <- get_GO_data(OrgDb, ont = "MF", keyType)
saveRDS(GO_DATA_MF, file = "GO_DATA_MF.rds")



######## 3. GO_BP term ######
keyType = "ENTREZID"
OrgDb = "org.Hs.eg.db"
ont="BP"
GO_DATA_BP <- get_GO_data(OrgDb, ont, keyType)
saveRDS(GO_DATA_BP, file = "GO_DATA_BP.rds")

