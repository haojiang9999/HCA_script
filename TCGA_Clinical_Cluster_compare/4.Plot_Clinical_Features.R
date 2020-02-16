#### 4.Plot_Clinical_Features.R
############## Gender #############################
ggplot(data = Plot.df, aes(x = dynamicColors, fill = gender.x)) +  
  geom_bar(position = "fill") + theme_minimal()+
  scale_fill_manual(values= alpha(c("#999999", "#E69F00", "#56B4E9"),0.6)) # transparency

############ Tumor site location ###############
table(Plot.df$tumor_site_location)
ggplot(data=subset(Plot.df, !is.na(tumor_site_location)), aes(x = dynamicColors, fill = tumor_site_location)) + 
  geom_bar(position = "fill") + theme_minimal()+
  scale_fill_manual(values= alpha(c("#999999", "#E69F00", "#56B4E9","#76EE00","#9A32CD","#CD6600"),0.6)) # transparency

table(Plot.df$anatomic_neoplasm_subdivision)
ggplot(data=subset(Plot.df, !is.na(anatomic_neoplasm_subdivision)), aes(x = dynamicColors, fill = anatomic_neoplasm_subdivision)) + 
  geom_bar(position = "fill") + theme_minimal()
  #scale_fill_manual(values= alpha(c("#999999", "#E69F00", "#56B4E9","#76EE00","#9A32CD","#CD6600"),0.6)) # transparency

############ Diagnosis Stage I II III ###############
Plot.df$diagnosis_stage <- gsub("[ABC]","",Plot.df$ajcc_pathologic_tumor_stage)
Plot.df[Plot.df$diagnosis_stage == "",]$diagnosis_stage <-NA
Plot.df[Plot.df$diagnosis_stage %in% "[Discrepancy]",]$diagnosis_stage <-NA
table(Plot.df$ajcc_pathologic_tumor_stage)
Plot.df$diagnosis_stage
ggplot(data=subset(Plot.df, !is.na(diagnosis_stage)), aes(x = dynamicColors, fill = diagnosis_stage)) + 
  geom_bar(position = "fill") + theme_minimal()+
  scale_fill_manual(values= alpha(c("#999999", "#E69F00", "#56B4E9","#76EE00","#9A32CD","#CD6600"),0.6)) # transparency
