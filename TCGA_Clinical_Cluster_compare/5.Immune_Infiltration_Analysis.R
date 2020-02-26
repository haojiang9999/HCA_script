#### 5.Immune_Infiltration_Analysis.R
library(ggplot2)
##### TotalLeukocyte
Plot.df$TotalLeukocyte
ggplot(data=subset(Plot.df, !is.na(TotalLeukocyte)), aes(x = dynamicColors, y = TotalLeukocyte,fill=dynamicColors)) + 
  geom_boxplot() + scale_fill_manual(values= c("blue","brown","turquoise","yellow"))
##### Different immune cell type comparison ######
## 1.Extract immune cell info
colnames(Plot.df)
immuneCell.df <- Plot.df[,61:82] 
## But this is relative infiltration of each cell type
## 2.Get the absolute infiltration value for immune cells
# Relative value multiple TotalLeukocyte(ESTIMATE)
immuneCell.Absolute.df <- as.data.frame(apply(immuneCell.df,2,function(x){x*Plot.df$TotalLeukocyte}))
## 3.Construct ggplot dataframe
immuneCell.Absolute.df$dynamicColors <- Plot.df$dynamicColors
library(reshape2)
dat.m <- melt(immuneCell.Absolute.df,id.vars='dynamicColors')
library(ggplot2)
ggplot(dat.m) + geom_boxplot(aes(x=variable, y=value, color=dynamicColors)) + coord_flip()+
scale_color_manual(values= c("blue","brown","turquoise","yellow"))
## 4.Cell type plot ##
colnames(immuneCell.Absolute.df)
# B cells
dat.m.Bcell <- melt(immuneCell.Absolute.df[,c(1:3,23)],id.vars='dynamicColors')
ggplot(dat.m.Bcell) + geom_boxplot(aes(x=variable, y=value, color=dynamicColors)) + 
 scale_color_manual(values= c("blue","brown","turquoise","yellow"))
# T cells
dat.m.Tcell <- melt(immuneCell.Absolute.df[,c(4:10,23)],id.vars='dynamicColors')
ggplot(dat.m.Tcell) + geom_boxplot(aes(x=variable, y=value, color=dynamicColors)) + 
  scale_color_manual(values= c("blue","brown","turquoise","yellow"))
# NK cells
dat.m.NKcell <- melt(immuneCell.Absolute.df[,c(11:12,23)],id.vars='dynamicColors')
ggplot(dat.m.NKcell) + geom_boxplot(aes(x=variable, y=value, color=dynamicColors)) + 
  scale_color_manual(values= c("blue","brown","turquoise","yellow"))
# Monocytopoiesis
dat.m.Monocell <- melt(immuneCell.Absolute.df[,c(13:18,23)],id.vars='dynamicColors')
ggplot(dat.m.Monocell) + geom_boxplot(aes(x=variable, y=value, color=dynamicColors)) +
  scale_color_manual(values= c("blue","brown","turquoise","yellow"))
# Granulopoiesis
dat.m.Grancell <- melt(immuneCell.Absolute.df[,c(19:22,23)],id.vars='dynamicColors')
ggplot(dat.m.Grancell) + geom_boxplot(aes(x=variable, y=value, color=dynamicColors)) + 
  scale_color_manual(values= c("blue","brown","turquoise","yellow"))
########## PD-1 expression ###############
dat.m.Grancell <- melt(Plot.df[,c(2,111:112)],id.vars='dynamicColors')
## Transform the tpm expression to log2(x+1)
ggplot(dat.m.Grancell) + geom_boxplot(aes(x=variable, y=log2(value+1), color=dynamicColors)) +
  scale_color_manual(values= c("blue","brown","turquoise","yellow"))















