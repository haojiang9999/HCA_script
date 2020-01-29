### combine GI development and embryo reference panel
# GI reference fetal panel (Tang) TPM
referPanel_of_fetal <- readRDS("/data8t_4/JH/MyJobs/Colon_SC_Project/GSE103239_Tang_GI_development/GSE95630_fetal_tissues/referPanel_of_fetal.rds")
# GI reference of Adult GI (Tang) TPM
referPanel_of_adult <- readRDS("/data8t_4/JH/MyJobs/Colon_SC_Project/GSE103239_Tang_GI_development/R_GSE103154_Tang_adult_colon_normal/Tang.GI.Adult.ref.panel.filtered.rds")
# Embryo reference of embryo development (E-MTAB-3929) RPKM
referPanel_of_embryo <- readRDS("/data8t_4/JH/MyJobs/Embryo/E_MTAB_3929/R_E_MTAB_3929/referPanel_of_embryo_E_MTAB_3929_RPKM.rds")
head(referPanel_of_embryo)
# Anno data colnames
colnames(referPanel_of_adult) <- paste("Adult_LIntes", colnames(referPanel_of_adult), sep = "_")
colnames(referPanel_of_fetal) <- paste("Fetal_", colnames(referPanel_of_fetal), sep = "_")
head(referPanel_of_adult)
# Build large reference panel 
geneName <- rownames(referPanel_of_adult[rownames(referPanel_of_adult) %in% rownames(referPanel_of_embryo),])
referPanel_of_mix<- cbind(referPanel_of_fetal[geneName,],
      referPanel_of_adult[geneName,],
      referPanel_of_embryo[geneName,])
head(referPanel_of_mix)
class(referPanel_of_mix)
# save resault
saveRDS(referPanel_of_mix, file = "referPanel_of_mix.rds")
