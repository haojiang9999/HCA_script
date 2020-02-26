## 1.Read data 
# Download from https://gdc.cancer.gov/about-data/publications/PanCanStemness-2018
clinical_PANCAN_patient_with_followup <- read.delim("/stor/jianghao/Paper_data/Machine_Learning_Identifies_Stemness_Features/clinical_PANCAN_patient_with_followup.tsv")
clinical_PANCAN_patient_with_followup[1:10,1:10]
Merged_EPC_reviews <- read.delim("/stor/jianghao/Paper_data/Machine_Learning_Identifies_Stemness_Features/Merged_EPC_reviews.tsv")
merged_sample_quality_annotations <- read.delim("/stor/jianghao/Paper_data/Machine_Learning_Identifies_Stemness_Features/merged_sample_quality_annotations.tsv")
pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16 <- read.csv("/stor/jianghao/Paper_data/Machine_Learning_Identifies_Stemness_Features/pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.csv",
                                                                                    header = T)
Purity_Ploidy_All_Samples_9_28_16 <- read.delim("/stor/jianghao/Paper_data/Machine_Learning_Identifies_Stemness_Features/Purity_Ploidy_All_Samples_9_28_16.tsv")
### Therre were more data down load from 
# https://gdc.cancer.gov/about-data/publications/PanCanStemness-2018
#TCGA-CDR-SupplementalTableS1 <- re("/stor/jianghao/Paper_data/Machine_Learning_Identifies_Stemness_Features/TCGA-CDR-SupplementalTableS1.xlsx")
TCGA.Kallisto.cibersort.relative <- read.delim("/stor/jianghao/Paper_data/Machine_Learning_Identifies_Stemness_Features/TCGA.Kallisto.cibersort.relative.tsv")





