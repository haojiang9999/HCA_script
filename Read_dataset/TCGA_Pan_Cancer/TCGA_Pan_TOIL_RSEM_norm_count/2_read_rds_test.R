ACC<- readRDS("ACC_UCSC_Toil_norm_count_dataset.rds")
ACC.RSEM.gene.norm_count <- ACC$ACC.RSEM.gene.norm_count_round
ACC.RSEM.gene.norm_count[1:5,1:5]
dim(ACC$ACC.pheno)

ACC.log2<- readRDS("ACC_UCSC_RSEM_norm_count_Log2(x+1)_hugo_dataset.rds")
ACC.ACC.RSEM.gene.norm_count.log2 <- ACC.log2$`ACC.RSEM.gene.norm_count_Log2(x+1)`
ACC.ACC.RSEM.gene.norm_count.log2[1:5,1:5]
