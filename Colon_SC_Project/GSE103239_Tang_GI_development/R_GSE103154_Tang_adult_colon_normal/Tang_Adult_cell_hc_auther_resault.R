### The authors original cluster resaults
## read file
getwd()
setwd("/data8t_4/JH/MyJobs/Colon_SC_Project/GSE103239_Tang_GI_development/R_GSE103154_Tang_adult_colon_normal/")
hcResault <- read.table("Adult_Binary_hc10_result.txt")
hcResault$hc_10
as.factor(hcResault$hc_10)
library(ggplot2)
p<-ggplot(hcResault,aes(x=tsne1,y=tsne2,color=as.factor(hcResault$hc_10)))
p<-p+geom_point(size = 2)
p

# Scatter plot change color
p + scale_color_brewer(palette="Paired")
#
range(1:10)
seq(10)
hcNames <- cbind(seq(10), c("Goblet_1", "Enter_2","Enter_3", "Goblet_2",
                 "Enter_1", "MKI67_High", "Endocrine", "Mesenchymal",
                 "Immune", "OLFM4_High"))
hcNames <- as.data.frame(hcNames)
colnames(hcNames) <- c("hc_10", "cellType")
hcResault.final <- hcResault[,c(1:7,13)]
hcResault.final$hc_10 <- as.factor(hcResault.final$hc_10)
hcResault.final <- dplyr::left_join(hcResault.final, hcNames, by = "hc_10")
## save data
saveRDS(hcResault.final, file = "Tang_GI_adult_author_cluster_cellType.rds")
write.csv(hcResault.final, file = "Tang_GI_adult_author_cluster_cellType.csv")

table(hcResault.final$cellType)
