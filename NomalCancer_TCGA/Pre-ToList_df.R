## Pre-ToList_df
sample1 <- AnnoRes.COAD.merge$TCGA.QG.A5YX.01
sample1.test<-as.data.frame(sample1)
DescriptionCol <- grep(".Description" , colnames(sample1.test))
sample1.test2 <- sample1.test[,DescriptionCol]
sample1.test3 <- as.character(unlist(sample1.test2))
df2 <- data.frame(sample1.test3, rep(1, length(sample1.test3)))
df1 <- data.frame(sample1.test3, rep(1, length(sample1.test3)))
Merge.test <- merge(df1,df2, by = "sample1.test3",all = T)

### test ToList_df
AnnoRes.COAD.merge.df<- ToList_df(AnnoRes.COAD.merge)
AnnoRes.COAD.merge.df[[1]]

          