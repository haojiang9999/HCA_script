#### GOenrich_JH ####
### load data
GO_DATA_BP <- readRDS("/stor/jianghao/Database/GO_DATA/GO_DATA_BP.rds")
GO_DATA_CC <- readRDS("/stor/jianghao/Database/GO_DATA/GO_DATA_CC.rds")
GO_DATA_MF <- readRDS("/stor/jianghao/Database/GO_DATA/GO_DATA_MF.rds")
### load clusterProfiler_JH
library(devtools)
packPath <- "/data8t_4/JH/MyJobs/test/R_package_process/clusterProfiler_JH"
load_all(packPath)
enrichGO_JH(gene = inputID,
            # universe      = background,
            OrgDb         = "org.Hs.eg.db",
            ont           = "CC",
            pAdjustMethod = "BH",
            pvalueCutoff  = 0.01,
            qvalueCutoff  = 0.05,
            readable      = TRUE)
