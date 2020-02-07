#### 2.Plot
colnames(Plot.df)
table(Plot.df$dynamicColors)
library(ggplot2)
ggplot(data = Plot.df, aes(x = cutree.res, y = gender.x),) + geom_col() 
library(tidyr)
Plot.df %>% gather(key = "var", value = "value", gender.x) %>% head()
test <- Plot.df %>% gather(key = "var", value = "value", gender.x) 
ggplot(data = Plot.df, aes(x = dynamicColors, fill = gender.x)) + geom_bar(position = "fill") 
ggplot(data = Plot.df, aes(x = dynamicColors, fill = msi)) + geom_bar(position = "fill") 
ggplot(data = Plot.df, aes(x = dynamicColors, fill = braf_mut)) + geom_bar() 
### Re-order the fill 
# I need tyo change the level order
Plot.df$msi <- relevel(Plot.df$msi, 'mss')
### Remove NA value
ggplot(data=subset(Plot.df, !is.na(msi)), aes(x = dynamicColors, fill = msi)) + 
  geom_bar(position = "fill") + theme_minimal()+
  scale_fill_manual(values= alpha(c("#999999", "#E69F00", "#56B4E9"),0.6)) # transparency

ggplot(data=subset(Plot.df, !is.na(msi)), aes(x = dynamicColors, fill = msi)) + 
  geom_bar() + theme_minimal()+
  scale_fill_manual(values= alpha(c("#999999", "#E69F00", "#56B4E9"),0.6)) # transparency
### Change the colors



##### TotalLeukocyte
Plot.df$TotalLeukocyte
ggplot(data=subset(Plot.df, !is.na(TotalLeukocyte)), aes(x = dynamicColors, y = TotalLeukocyte,fill=dynamicColors)) + 
  geom_boxplot() + scale_fill_manual(values= alpha(c("#999999", "#E69F00", "#56B4E9","#76EE00","#9A32CD"),0.6)) # transparency
Plot.df$Mast.cells.activated
#### multiple columns
# different immuno cells
library(reshape2)
dat.m <- melt(Plot.df,id.vars='dynamicColors', measure.vars=c('B.cells.memory','NK.cells.resting','Mast.cells.activated'))
library(ggplot2)
p <- ggplot(dat.m) +
  geom_boxplot(aes(x=variable, y=value, color=dynamicColors))
p


data=subset(Plot.df, !is.na(msi))
leu.blue <- data[data$dynamicColors == "blue",]$TotalLeukocyte
leu.brown <- data[data$dynamicColors == "brown",]$TotalLeukocyte
leu.turquoise <- data[data$dynamicColors == "turquoise",]$TotalLeukocyte
t.test(leu.blue,leu.brown)
t.test(leu.blue,leu.turquoise)
wilcox.test(leu.blue,leu.brown)
wilcox.test(leu.blue,leu.turquoise)
data$dynamicColors
