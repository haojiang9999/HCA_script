#### Build AnnoRes Dataframe
#### Cluster COAD by annoRes ####
# load data
AnnoRes.688.merge.datasets <- readRDS("/data8t_4/JH/MyJobs/NormalCancer/AnnoRes.688.merge.datasets.rds")
CorValueRes <- AnnoRes.688.merge.datasets$CorValueRes
TopRes <- AnnoRes.688.merge.datasets$TopRes
AnnoRes.688.merge <- AnnoRes.688.merge.datasets$AnnoRes.688.merge
#### re-construct data
### Step1: Extract Top 5 correlation value for each cells
TopNcellTypes.corValue <- TopRes$TopNcellTypesValue
TopNcellTypes.corValue <- as.data.frame(TopNcellTypes.corValue)
TopNcellTypes.corValue[TopNcellTypes.corValue< 0.5] <- 0
### Step2:AnnoRes parse
source("./Pre_Step2_AnnoRes_parse_to_table.R")
AnnoRes.df<- ToList_df(AnnoRes.688.merge)
AnnoRes.df[[1]]
AnnoRes.merge.df <- MergeData_df(AnnoRes.df)

### Step3:Combine two df
head(cbind(colnames(TopNcellTypes.corValue), colnames(AnnoRes.merge.df)))
Combine_df <- rbind(TopNcellTypes.corValue,AnnoRes.merge.df)
Combine_df[is.na(Combine_df)] <- 0
head(Combine_df)
saveRDS(Combine_df, file = "AnnoRes.688.merge_df.rds")













