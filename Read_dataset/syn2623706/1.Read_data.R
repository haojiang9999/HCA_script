#### 1.Read expressio data 
formatted_crc_data <- read.table("/stor/jianghao/Synapse/classifier/formatted_crc_data.txt")
dim(formatted_crc_data)
formatted_crc_data[1:5,1:5]
rownames(formatted_crc_data)
colnames(formatted_crc_data)
saveRDS(formatted_crc_data, file = "formatted_crc_data.rds")
### 2.Clinical molecular data
clinical_molecular_public_all <- read.table("/stor/jianghao/Synapse/mergedPhenotype/clinical_molecular_public_all.txt",
                                            header = T)
head(clinical_molecular_public_all)
### 3.cms labels
cms_labels_public_all <- read.table("/stor/jianghao/Synapse/mergedPhenotype/cms_labels_public_all.txt",
                                            header = T)
head(cms_labels_public_all)
