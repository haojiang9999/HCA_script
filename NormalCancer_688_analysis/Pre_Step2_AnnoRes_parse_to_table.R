#### Step2_AnnoRes_parse_to_table.R
#### Prepare list samples to list of data frame
ToList_df <- function(AnnoRes){
  lapply(AnnoRes, function(x){
    df1<-as.data.frame(x)
    DescriptionCol <- grep(".Description" , colnames(df1))
    df2 <- df1[,DescriptionCol]
    df3 <- as.character(unlist(df2))
    df4 <- data.frame(df3, rep(1, length(df3)))
    return(df4)
  })
  
}
#### Merge AnnoRes list to dataframe
#df_list = AnnoRes.COAD.merge.df
MergeData_df <- function(df_list){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  if (!require("tidyverse")) 
    BiocManager::install("tidyverse")
  require(tidyverse)
  require(purrr)
  merged.df <- df_list %>% reduce(full_join, by = "df3")
  merged.df <- merged.df[!is.na(merged.df[,1]),]
  rownames(merged.df) <- as.character(merged.df[,1])
  merged.df <- merged.df[,-1]
  colnames(merged.df) <- names(df_list)
  return(merged.df)
}