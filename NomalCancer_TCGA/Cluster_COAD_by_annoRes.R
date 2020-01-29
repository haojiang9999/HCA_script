#### Cluster COAD by annoRes ####
# load data
AnnoRes.COAD.datasets <- readRDS("AnnoRes.COAD.datasets.rds")
CorValueRes <- AnnoRes.COAD.datasets$CorValueRes
TopRes <- AnnoRes.COAD.datasets$TopRes
AnnoRes.COAD.merge <- AnnoRes.COAD.datasets$AnnoRes.COAD.merge
#### re-construct data
### Step1: Extract Top 5 correlation value for each cells
TopNcellTypes.corValue <- TopRes$TopNcellTypes.corValue
TopNcellTypes.corValue <- as.data.frame(TopNcellTypes.corValue)
TopNcellTypes.corValue[TopNcellTypes.corValue< 0.5] <- 0
### Step2:AnnoRes parse
AnnoRes.COAD.merge.df<- ToList_df(AnnoRes.COAD.merge)
AnnoRes.COAD.merge.df[[1]]
library(tidyverse)
merged.df <- AnnoRes.COAD.merge.df %>% reduce(full_join, by = "df3")
merged.df <- merged.df[!is.na(merged.df[,1]),]
rownames(merged.df) <- as.character(merged.df[,1])
merged.df <- merged.df[,-1]
colnames(merged.df) <- names(AnnoRes.COAD.merge.df)
merged.df

### Step3:Combine two df
cbind(colnames(TopNcellTypes.corValue), colnames(merged.df))
Combine_df <- rbind(TopNcellTypes.corValue,merged.df)
Combine_df[is.na(Combine_df)] <- 0
head(Combine_df)
saveRDS(Combine_df, file = "COAD_AnnoRes_Combine.rds")
### Step4: Cluster 



