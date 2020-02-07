#### 3.Plot_Seperate_clinical_info
library(ggplot2)
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

############ Diagnosis Stage I II III ###############
Plot.df$diagnosis_stage <- gsub("[ABC]","",Plot.df$ajcc_pathologic_tumor_stage)
Plot.df[Plot.df$diagnosis_stage == "",]$diagnosis_stage <-NA
Plot.df[Plot.df$diagnosis_stage %in% "[Discrepancy]",]$diagnosis_stage <-NA
table(Plot.df$ajcc_pathologic_tumor_stage)
Plot.df$diagnosis_stage
ggplot(data=subset(Plot.df, !is.na(diagnosis_stage)), aes(x = dynamicColors, fill = diagnosis_stage)) + 
  geom_bar(position = "fill") + theme_minimal()+
  scale_fill_manual(values= alpha(c("#999999", "#E69F00", "#56B4E9","#76EE00","#9A32CD","#CD6600"),0.6)) # transparency

############ Tumor site location ###############