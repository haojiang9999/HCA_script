#### 3.GGplot_Aneuploidy_score.R
# Genomic and Functional Approaches to Understanding Cancer Aneuploidy
COAD.Aneuploidy.score <- COAD_Aneuploidy_score_dataset$COAD.Aneuploidy.score
### 1)Merge table
Plot.df <- dplyr::left_join(Cluster.df, COAD.Aneuploidy.score, by = "rownames")
colnames(Plot.df)
Plot.df$AneuploidyScore.AS.
class(Plot.df$AneuploidyScore.AS.)
class(Plot.df$X1p)
i = "AneuploidyScore.AS."
for(i in colnames(Plot.df)){
  Plot.df.sub <- Plot.df[,c("dynamicColors",i)]
  ##### 1.Discrete variable #####
  if (is.factor(Plot.df.sub[,i])) {
    ######## 1.Discrete variable
    # Remove NA value
    data=subset(Plot.df.sub, !is.na(Plot.df.sub[,i]))
    p <- ggplot(data, aes(x = dynamicColors, fill = data[,i])) + 
      geom_bar(position = "fill") + theme_minimal()+ scale_fill_discrete(name =i)+
      labs(title =i)
    print(p)
  } else if (is.numeric(Plot.df.sub[,i])) {
    ######## 2.Continuous variable
    # Remove NA value
    data=subset(Plot.df.sub, !is.na(Plot.df.sub[,i]))
    p2 <-ggplot(data,aes(x=data[,1], y=data[,2]),color=dynamicColors)  + 
      geom_boxplot(outlier.colour = "red",outlier.shape = 1,colour = c("blue","brown","turquoise")) + 
      labs(title =i, y = i)
    print(p2)
  }  else if (is.integer(Plot.df.sub[,i])) {
      ######## 2.Continuous variable
      # Remove NA value
      data=subset(Plot.df.sub, !is.na(Plot.df.sub[,i]))
      p2 <-ggplot(data,aes(x=data[,1], y=data[,2]))  + 
        geom_boxplot(outlier.colour = "red",outlier.shape = 1,colour = c("blue","brown","turquoise","yellow")) + 
        labs(title =i, y = i)
      print(p2)
  } else {
    #### 3.If something wrong
    print(paste0("Some thing wrong with---",i))
  } 
}





