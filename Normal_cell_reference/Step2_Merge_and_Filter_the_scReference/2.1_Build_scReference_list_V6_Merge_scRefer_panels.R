### 2.1_Build_scReference_list_V6_Merge_scRefer_panels.R
## V6.1 Mearge all reference panel
### Step1 Find common gene across these reference panel
GeneNames.Tang.Adult.colon <- rownames(Tang.Adult.colon.ref)
GeneNames.Tang.Fetal.GI <- rownames(Tang.Fetal.GI.ref)
GeneNames.Tang.Normal.embryo <- rownames(Tang.Normal.embryo.ref)
GeneNames.Lanner.Preim.embryo <- rownames(Lanner.Preim.embryo.ref)
GeneNames.Zemin.CRC.Tcell <- rownames(Zemin.CRC.Tcell.ref)
GeneNames.Zemin.CRC.Tcell.Cluster <- rownames(Zemin.CRC.Tcell.Cluster.ref)
# 19244 comnon genes
GeneNames.common.V6.1 <- Reduce(intersect, list(GeneNames.Tang.Adult.colon,
                                           GeneNames.Tang.Fetal.GI,
                                          GeneNames.Tang.Normal.embryo,
                                          GeneNames.Lanner.Preim.embryo,
                                          GeneNames.Zemin.CRC.Tcell,
                                          GeneNames.Zemin.CRC.Tcell.Cluster))
scRef.merge.V6.1 <- cbind(Tang.Adult.colon.ref[GeneNames.common.V6.1,],
                          Tang.Fetal.GI.ref[GeneNames.common.V6.1,],
                          Tang.Normal.embryo.ref[GeneNames.common.V6.1,],
                          Lanner.Preim.embryo.ref[GeneNames.common.V6.1,],
                          Zemin.CRC.Tcell.ref[GeneNames.common.V6.1,],
                          Zemin.CRC.Tcell.Cluster.ref[GeneNames.common.V6.1,])
head(scRef.merge.V6.1)
dim(scRef.merge.V6.1)
### Step2.Find genes had high coefficience-varience
source("/data8t_4/JH/MyJobs/1_R_script/FUN_TopCV.R")
scRef.merge.V6.1.CV8000 <- TopCV(scRef.merge.V6.1, TopN = 8000, MARGIN = 1)
scRef.merge.V6.1.CV4000 <- TopCV(scRef.merge.V6.1, TopN = 4000, MARGIN = 1)
scRef.merge.V6.1.CV3000 <- TopCV(scRef.merge.V6.1, TopN = 3000, MARGIN = 1)
scRef.merge.V6.1.CV2500 <- TopCV(scRef.merge.V6.1, TopN = 2500, MARGIN = 1)
scRef.merge.V6.1.CV2000 <- TopCV(scRef.merge.V6.1, TopN = 2000, MARGIN = 1)
scRef.merge.V6.1.CV1500 <- TopCV(scRef.merge.V6.1, TopN = 1500, MARGIN = 1)

### Step3.log10(x+1) transformed
scRef.merge.V6.1.CV8000.log10 <- log10(scRef.merge.V6.1.CV8000+1)
scRef.merge.V6.1.CV4000.log10 <- log10(scRef.merge.V6.1.CV4000+1)
scRef.merge.V6.1.CV3000.log10 <- log10(scRef.merge.V6.1.CV3000+1)
scRef.merge.V6.1.CV2500.log10 <- log10(scRef.merge.V6.1.CV2500+1)
scRef.merge.V6.1.CV2000.log10 <- log10(scRef.merge.V6.1.CV2000+1)
scRef.merge.V6.1.CV1500.log10 <- log10(scRef.merge.V6.1.CV1500+1)

### Step4.Save the reference
scRef.merge.V6.1.log10.list <- list(scRef.merge.V6.1.CV8000.log10 = scRef.merge.V6.1.CV8000.log10,
                                    scRef.merge.V6.1.CV4000.log10 = scRef.merge.V6.1.CV4000.log10,
                                    scRef.merge.V6.1.CV3000.log10 = scRef.merge.V6.1.CV3000.log10,
                                    scRef.merge.V6.1.CV2500.log10 = scRef.merge.V6.1.CV2500.log10,
                                    scRef.merge.V6.1.CV2000.log10 = scRef.merge.V6.1.CV2000.log10,
                                    scRef.merge.V6.1.CV1500.log10 = scRef.merge.V6.1.CV1500.log10)
saveRDS(scRef.merge.V6.1.log10.list, file = "scRef.merge.V6.1.log10.CV.rds")
