#### 3.2.ImmuneCell_Infiltration_Clean.R
## 1) Loading data
### Papaer：Machine Learning Identifies Stemness Features Associated with Oncogenic Dedifferentiation
### https://gdc.cancer.gov/about-data/publications/PanCanStemness-2018
COAD_Machine_Learning_Stemness.ImmuneCell_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/Machine_Learning_Identifies_Stemness_Features/COAD_Machine_Learning_Stemness.ImmuneCell.cibersort.relative_dataset.rds")
#COAD_Machine_Learning_Stemness.ImmuneCell_dataset$Machine_Learning_Stemness.metadata
COAD_Machine_Learning_StemnessimmuneCell.Absolute <- COAD_Machine_Learning_Stemness.ImmuneCell_dataset$COAD_Machine_Learning_StemnessimmuneCell.Absolute

## 2)clean table Remove duplicated one randomly
table(duplicated(COAD_Machine_Learning_StemnessimmuneCell.Absolute$rownames))
COAD_Machine_Learning_StemnessimmuneCell.Absolute.Unique <- COAD_Machine_Learning_StemnessimmuneCell.Absolute[!duplicated(COAD_Machine_Learning_StemnessimmuneCell.Absolute$rownames),]

## 3)Merge table
MergeTable.Machine.Immune.clean <- dplyr::left_join(Cluster.df, COAD_Machine_Learning_StemnessimmuneCell.Absolute.Unique, by = "rownames")
# Remove yellow one
MergeTable.Machine.Immune.clean <- MergeTable.Machine.Immune.clean[MergeTable.Machine.Immune.clean$dynamicColors != "yellow",]
dim(MergeTable.Machine.Immune.clean)
MergeTable.Machine.Immune.clean <- MergeTable.Machine.Immune.clean[complete.cases(MergeTable.Machine.Immune.clean), ]


MergeTable.Machine.Immune.clean<- MergeTable.Machine.Immune.clean[,-c(1,3,4)]

## 4) Plotting
# B cells
table(MergeTable.Machine.Immune.clean$dynamicColors)
sampleSize <- paste(levels(MergeTable.Machine.Immune.clean$dynamicColors),"(N=",table(MergeTable.Machine.Immune.clean$dynamicColors),")",sep="")
sampleSize <- paste(sampleSize,collapse=" ")
dat.m.Bcell <- melt(MergeTable.Machine.Immune.clean[,c(2:4,1)],id.vars='dynamicColors')
ggplot(dat.m.Bcell) + geom_boxplot(aes(x=variable, y=value, color=dynamicColors)) + 
  scale_color_manual(values= c("blue","brown","turquoise","yellow")) +
  theme_classic()+
  theme(axis.line = element_line(colour = 'black', size = 2),axis.ticks = element_line(colour = "black", size = 2),
        axis.ticks.length=unit(.5, "cm"),legend.key.size =unit(3,"line"))+
  labs(title = "B cells",subtitle =sampleSize)
