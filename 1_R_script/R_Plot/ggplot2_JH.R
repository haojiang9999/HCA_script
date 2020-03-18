#### ggplot2_JH.R
require(ggplot2)
### Box plot for continous value
ggplot2_boxPlot_JH <- function(df,groupInfo,character,...){
  #character <- "Leukocyte.Fraction"
  data=subset(df, !is.na(df[,character]))
  xlabs <- paste(levels(data[,groupInfo]),"\n(N=",table(data[,groupInfo]),")",sep="")
  p <- ggplot(data,aes(x = data[,groupInfo], y = data[,character],fill=data[,groupInfo]))+
    geom_boxplot()+scale_x_discrete(labels=xlabs) + theme_minimal()+
    scale_fill_manual(...) +
    labs(title =character, y = character)
  print(p)
}
### Bar plot for categorical data
ggplot2_barPlot_JH <- function(df,groupInfo,character,...){
  #character <- "CDKN2A.methylation"
  data=subset(df, !is.na(df[,character]))
  data[,character] <- as.factor(data[,character])
  xlabs <- paste(levels(data[,groupInfo]),"\n(N=",table(data[,groupInfo]),")",sep="")
  ggplot(data,aes(x = data[,groupInfo],fill=data[,character]))+
    geom_bar(position = "fill") +scale_x_discrete(labels=xlabs) + theme_minimal() +
    scale_fill_manual(...) + #values= alpha(c("#999999", "#E69F00", "#56B4E9"),0.6)
    labs(title =character, y = character)
}

