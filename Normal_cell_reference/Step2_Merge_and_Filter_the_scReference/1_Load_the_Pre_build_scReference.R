#### Read the Step1 build scReference table
Tang.Adult.colon.ref <- readRDS("/data8t_4/JH/MyJobs/Normal_cell_reference/Step1_Normal_cell_reference_panel_construction/Tang.Adult.colon.ref.rds")
Tang.Fetal.GI.ref <- readRDS("/data8t_4/JH/MyJobs/Normal_cell_reference/Step1_Normal_cell_reference_panel_construction/Tang.Fetal.GI.ref.rds")
Tang.Normal.embryo.ref <- readRDS("/data8t_4/JH/MyJobs/Normal_cell_reference/Step1_Normal_cell_reference_panel_construction/Tang.Normal.embryo.ref.rds")
Lanner.Preim.embryo.ref <- readRDS("/data8t_4/JH/MyJobs/Normal_cell_reference/Step1_Normal_cell_reference_panel_construction/Lanner.Preim.embryo.ref.rds")
Zemin.CRC.Tcell.ref <- readRDS("/data8t_4/JH/MyJobs/Normal_cell_reference/Step1_Normal_cell_reference_panel_construction/Zemin.CRC.Tcell.ref.rds")
Zemin.CRC.Tcell.Cluster.ref <- readRDS("/data8t_4/JH/MyJobs/Normal_cell_reference/Step1_Normal_cell_reference_panel_construction/Zemin.CRC.Tcell.Cluster.ref.rds")
RCA.CRC.nonEpi.ref <- readRDS("/data8t_4/JH/MyJobs/Normal_cell_reference/Step1_Normal_cell_reference_panel_construction/RCA.CRC.nonEpi.ref.rds")
#### Gene names
#GeneNames.Tang.Adult.colon <- rownames(Tang.Adult.colon.ref)
#GeneNames.Tang.Fetal.GI <- rownames(Tang.Fetal.GI.ref)
#GeneNames.Tang.Normal.embryo <- rownames(Tang.Normal.embryo.ref)
#GeneNames.Lanner.Preim.embryo <- rownames(Lanner.Preim.embryo.ref)
#### find the common genes
#GeneNames.common <- Reduce(intersect, list(GeneNames.Tang.Adult.colon,GeneNames.Tang.Fetal.GI,
#                       GeneNames.Tang.Normal.embryo,GeneNames.Lanner.Preim.embryo))


### How many non-Zero genes in the reference panel
summary(apply(Tang.Adult.colon.ref, MARGIN = 2 ,FUN = function(x){
sum(x != 0)
}))
summary(apply(Tang.Fetal.GI.ref, MARGIN = 2 ,FUN = function(x){
  sum(x != 0)
}))

summary(apply(Tang.Normal.embryo.ref, MARGIN = 2 ,FUN = function(x){
  sum(x != 0)
}))
summary(apply(Lanner.Preim.embryo.ref, MARGIN = 2 ,FUN = function(x){
  sum(x != 0)
}))

