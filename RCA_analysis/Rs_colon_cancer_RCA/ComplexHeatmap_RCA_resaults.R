### ComplexHeatmap plot
library(ComplexHeatmap)
color_scheme = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                  "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(100);

# scale the data
scaled.mx <- scale(test_for_clust_power,
                   center = T)
## sample annotation
colnames(scaled.mx) <- gsub("_scTrioSeq2Rna_scTrioSeq2Rna_", "_scTrioSeq2Rna_", colnames(scaled.mx))
cellSitesPatient <- sapply(strsplit(as.character(colnames(scaled.mx)), "_"), "[[", 4 )
cellSites <-substr(cellSitesPatient, start = 1, stop = 2)
pateint <- sapply(strsplit(as.character(colnames(scaled.mx)), "_"), "[[", 3 )
GSE97693_Tang_TPM_cells_sample_anno <- cbind(sampleName = colnames(scaled.mx),
                                             cellSitesPatient = cellSitesPatient,
                                             cellSites = cellSites,
                                             pateintID = pateint)

## contruct heatmap column annotation
column_ha = HeatmapAnnotation(cluster = clusterResault[[2]]$dynamicColors, # myclustered by RCA
                              pateints = pateint,
                              sites = cellSites,
                              col = list(cluster = c("red" = "red", "green" = "green", 
                                                     "blue" = "blue","turquoise" = "turquoise",
                                                     "brown" = "brown", "yellow" = "yellow",
                                                     "grey" = "grey"),
                                         pateints = c("CRC03" = "519", "CRC06" = "544",
                                                      "CRC04" = "569", "CRC11" = "594",
                                                      "CRC09" = "619", "CRC10" = "644"),
                                         sites = c("PT" = "burlywood4", "LN" = "darkmagenta",
                                                   "NC" = "darkolivegreen1")
                                         )
                              )



Heatmap(scaled.mx, 
        col = color_scheme, # heatmap color
        show_column_names = F,
        row_dend_reorder = TRUE,
        row_order = order(seq(1,44)),
        cluster_columns = clusterResault[[1]], # my cluster resaults
        top_annotation = column_ha             # my cluster group 
)
class(test_for_clust_power)
