#### 1.Laoding_data.R 
# 1)Read COAD expression data
TCGA_COAD_RNAseqV2_normalized_log2_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/TCGA_Hub/TCGA_Colon_Cancer_COAD/TCGA_COAD_RNAseqV2_normalized_log2_dataset.rds")
COAD.HiSeqV2.log2 <- TCGA_COAD_RNAseqV2_normalized_log2_dataset$COAD.HiSeqV2.log2
COAD.pheno <- TCGA_COAD_RNAseqV2_normalized_log2_dataset$COAD.pheno
TCGA_COAD_RNAseqV2_normalized_log2_dataset$COAD.HiSeqV2.metadata
# 2)Read cluster resaults
Cluster.20200201.V7.Tumor <- readRDS("/data8t_4/JH/MyJobs/NormalCancer_TCGA_V2/Cluster.20200201.V7.Tumor.rds")
cutree.res <- Cluster.20200201.V7.Tumor$cutree.res
dynamicColors <- Cluster.20200201.V7.Tumor$dynamicColors
Cluster.df <- cbind(cutree.res,dynamicColors) 
Cluster.df <- as.data.frame(Cluster.df)
Cluster.df$rownames <- rownames(Cluster.df)
hclust.Res <- Cluster.20200201.V7.Tumor$hclust.Res
# 3)Read GeneSetCollection object gmt download from MSigDB
MSigDB_gmt_V7_symbols_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/MSigDB/MSigDB_gmt_V7_symbols_dataset.rds")
MSigDB_gmt_V7_symbols <- MSigDB_gmt_V7_symbols_dataset$symbols.gmt.list
MSigDB_V7_metadata <- MSigDB_gmt_V7_symbols_dataset$MSigDB_V7_metadata
c2.all.v7.0.symbols.gmt <- MSigDB_gmt_V7_symbols$c2.all.v7.0.symbols.gmt
c2.all.v7.0.symbols.gmt
description(c2.all.v7.0.symbols.gmt[[1]])





