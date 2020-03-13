#### 4.GGplot_COAD_Oncogenic_Signaling_Alteration_dataset.R
# Oncogenic Signaling Pathways in The Cancer
All_pathway_level <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/TCGA_PanCancerAtlas_Publications/Oncogenic_Signaling_Pathways_in_The_Cancer/All__Oncogenic_Signaling_Alteration_pathway_level.rds")
### 1)Merge table
Plot.df.Alt.pathway <- dplyr::left_join(Cluster.df, All_pathway_level, by = "rownames")
library(ggplot2)
for(i in colnames(Plot.df.Alt.pathway)[5:13]){
  # i = "Cell.Cycle"
  Plot.df.sub <- Plot.df.Alt.pathway[,c("dynamicColors",i)]
    # Remove NA value
  data=subset(Plot.df.sub, !is.na(Plot.df.sub[,i]))
  data[,i] <- as.factor(data[,i])
    p <- ggplot(data, aes(x = dynamicColors, fill = data[,i])) + 
      geom_bar(position = "fill") + theme_minimal()+ scale_fill_discrete(name =i)+
      labs(title =i)
    print(p)
}











