#### cluster COAD by CorValueRes ####
head(CorValueRes)
CorValueRes[[1]]
COAD.CorValue <- lapply(CorValueRes, colSums)
COAD.CorValue_df<- data.frame(COAD.CorValue)
saveRDS(COAD.CorValue_df, file = "COAD_referPanel_of_mix_V1.rds")

class(COAD.CorValue[[1]])
