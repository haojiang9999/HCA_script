#### 2.GGplot_GI_Adenocarcinomas_Characteristics.R
library(ggplot2)
Cluster.df
################## 1.GI_Adenocarcinomas_Characteristics #####################
### 1)Merge table
Plot.df <- dplyr::left_join(Cluster.df, COAD.GI.Adenocarcinomas.Characteristics, by = "rownames")
colnames(Plot.df)

for(i in colnames(Plot.df)[-c(1:5)]){
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
  } else {
    #### 3.If something wrong
    print(paste0("Some thing wrong with---",i))
  } 
}

################## 2.GI_Adenocarcinomas_COAD.Epigenetic.Silencing.Calls #####################
COAD.Epigenetic.Silencing.Calls
COAD.Epigenetic.Silencing.Calls[1:5,1:5]
### 1)Building plot table
Gene = "MLH1"
gene_tb <- COAD.Epigenetic.Silencing.Calls[Gene,]
gene_tb <- as.data.frame(t(gene_tb))
gene_tb$rownames <- rownames(gene_tb)
### 2)merge
Plot.df.gene <- dplyr::left_join(Cluster.df, gene_tb, by = "rownames")
colnames(Plot.df.gene)
data=subset(Plot.df.gene, !is.na(Plot.df.gene[,Gene]))
p <- ggplot(data, aes(x = dynamicColors, fill = data[,Gene])) + 
  geom_bar(position = "fill") + theme_minimal()+ scale_fill_discrete(name =Gene)+
  labs(title =Gene)
print(p)

table(Cluster.df$dynamicColors)





