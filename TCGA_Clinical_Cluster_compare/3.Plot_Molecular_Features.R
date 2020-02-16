#### 3.Plot_Seperate_clinical_info
library(ggplot2)
############## Mutation count ####################
Plot.df$totalMut
ggplot(data=subset(Plot.df, !is.na(totalMut)), aes(x = dynamicColors, y = totalMut,fill=dynamicColors)) + 
  geom_boxplot() + scale_fill_manual(values= alpha(c("#999999", "#E69F00", "#56B4E9","#76EE00","#9A32CD"),0.6)) # transparency
############## SCNA count (by genes) #############
Plot.df$SCNA.gene.counts
ggplot(data=subset(Plot.df, !is.na(SCNA.gene.counts)), aes(x = dynamicColors, y = SCNA.gene.counts,fill=dynamicColors)) + 
  geom_boxplot() + scale_fill_manual(values= alpha(c("#999999", "#E69F00", "#56B4E9","#76EE00","#9A32CD"),0.6)) # transparency
############## MSI Porpotion #####################
# Plot stack order
Plot.df$msi <- relevel(Plot.df$msi, 'mss')
# Remove NA value
ggplot(data=subset(Plot.df, !is.na(msi)), aes(x = dynamicColors, fill = msi)) + 
  geom_bar(position = "fill") + theme_minimal()+
  scale_fill_manual(values= alpha(c("#999999", "#E69F00", "#56B4E9"),0.6)) # transparency

############## Gender #############################
ggplot(data = Plot.df, aes(x = dynamicColors, fill = gender.x)) +  
  geom_bar(position = "fill") + theme_minimal()+
  scale_fill_manual(values= alpha(c("#999999", "#E69F00", "#56B4E9"),0.6)) # transparency

############# CIMP ##############################
# remove NA
#Plot.df$msi <- relevel(Plot.df$msi, 'mss')
table(!is.na(Plot.df$cimp))
Plot.df$cimp <- factor(Plot.df$cimp, levels = c("CIMP.Neg", "CIMP.Low", "CIMP.High"))
# Remove NA value
ggplot(data=subset(Plot.df, !is.na(cimp)), aes(x = dynamicColors, fill = cimp)) + 
  geom_bar(position = "fill") + theme_minimal()+
  scale_fill_manual(values= alpha(c("#999999", "#E69F00", "#56B4E9"),0.6)) # transparency

############ kras_mut #####################
Plot.df$kras_mut <- as.factor(Plot.df$kras_mut)
ggplot(data=subset(Plot.df, !is.na(kras_mut)), aes(x = dynamicColors, fill = kras_mut)) + 
  geom_bar(position = "fill") + theme_minimal()+
  scale_fill_manual(values= alpha(c("#999999", "#E69F00", "#56B4E9"),0.6)) # transparency

############ braf_mut #####################
Plot.df$braf_mut <- as.factor(Plot.df$braf_mut)
ggplot(data=subset(Plot.df, !is.na(braf_mut)), aes(x = dynamicColors, fill = braf_mut)) + 
  geom_bar(position = "fill") + theme_minimal()+
  scale_fill_manual(values= alpha(c("#999999", "#E69F00", "#56B4E9"),0.6)) # transparency

########### Gene Mutation analysis #############
## APC
Plot.df$APC <- as.factor(Plot.df$APC)
ggplot(data=subset(Plot.df, !is.na(APC)), aes(x = dynamicColors, fill = APC)) + 
  geom_bar(position = "fill") + theme_minimal()+
  scale_fill_manual(values= alpha(c("#999999", "#E69F00", "#56B4E9"),0.6)) # transparency
## KRAS
Plot.df$KRAS <- as.factor(Plot.df$KRAS)
ggplot(data=subset(Plot.df, !is.na(KRAS)), aes(x = dynamicColors, fill = KRAS)) + 
  geom_bar(position = "fill") + theme_minimal()+
  scale_fill_manual(values= alpha(c("#999999", "#E69F00", "#56B4E9"),0.6)) # transparency
## BRAF
Plot.df$BRAF <- as.factor(Plot.df$BRAF)
ggplot(data=subset(Plot.df, !is.na(BRAF)), aes(x = dynamicColors, fill = BRAF)) + 
  geom_bar(position = "fill") + theme_minimal()+
  scale_fill_manual(values= alpha(c("#999999", "#E69F00", "#56B4E9"),0.6)) # transparency
## TP53
Plot.df$TP53 <- as.factor(Plot.df$TP53)
ggplot(data=subset(Plot.df, !is.na(TP53)), aes(x = dynamicColors, fill = TP53)) + 
  geom_bar(position = "fill") + theme_minimal()+
  scale_fill_manual(values= alpha(c("#999999", "#E69F00", "#56B4E9"),0.6)) # transparency

############### Tumor cell purity ###############

Plot.df$purity
ggplot(Plot.df) + geom_boxplot(aes(x=dynamicColors, y=purity, color=dynamicColors)) 
