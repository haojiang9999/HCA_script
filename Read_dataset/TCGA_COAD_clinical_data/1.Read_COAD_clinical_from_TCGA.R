#### Read clinical data from TCGA portal
# download from https://portal.gdc.cancer.gov/repository?filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.case_id%22%2C%22value%22%3A%5B%22268c01b3-d2ce-44c0-a2fe-ea846a1253cc%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_category%22%2C%22value%22%3A%5B%22Clinical%22%5D%7D%7D%5D%7D&searchTableTab=files
nationwidechildrens.org_clinical_patient_coad <- read.delim("/data8t_4/JH/MyJobs/Read_dataset/TCGA_COAD_clinical_data/COAD_Clinical_data/nationwidechildrens.org_clinical_patient_coad.txt",header = T)
colnames(nationwidechildrens.org_clinical_patient_coad)
clinical_patient_coad <- nationwidechildrens.org_clinical_patient_coad[-(1:2),]
head(clinical_patient_coad)
colnames(clinical_patient_coad)
# What I need:
colsNeed <- c("bcr_patient_barcode","vascular_invasion_indicator","lymphovascular_invasion_indicator",
              "history_colon_polyps","anatomic_neoplasm_subdivision")
### Finally I give up  the infomartion like "microsatellite_instability"
# annotion number less then I expected, 
table(clinical_patient_coad$anatomic_neoplasm_subdivision)
### Well still some important informations
clinical_patient_coad_subset <- clinical_patient_coad[,colsNeed]
### 1.I will re-define colon anatomic_neoplasm_subdivision 
# Into more simple version right or left
table(clinical_patient_coad_subset$anatomic_neoplasm_subdivision)
Tumor.site <- data.frame(anatomic_neoplasm_subdivision = c("Ascending Colon","Cecum","Hepatic Flexure","Transverse Colon", # Right colon
        "Descending Colon","Rectosigmoid Junction","Sigmoid Colon","Splenic Flexure"),# Left colon
        tumor_site_location = c("Right","Right","Right","Right","Left","Left","Left","Left"))
### 2.Merge the two table clinical_patient_coad_subset and Tumor.site
TCGA.COAD.Clinical <- dplyr::left_join(clinical_patient_coad_subset, Tumor.site, by = "anatomic_neoplasm_subdivision")
## Change sample names
sampleID <- as.character(TCGA.COAD.Clinical$bcr_patient_barcode)
TCGA.COAD.Clinical$patient_barcode <- gsub("-",".",sampleID)

table(duplicated(sampleID))












