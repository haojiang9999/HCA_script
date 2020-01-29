#### Unscaled plotting of correlation results
#Cell_info
DimeReduPlot(mx = (Cor.Res.CV8000$Cor.merged), color = Pheno.merged$Cell_info, 
             tiltle = "Cor.Res.CV8000_Unscale_Cell_info", print = F)
DimeReduPlot(mx = (Cor.Res.CV4000$Cor.merged), color = Pheno.merged$Cell_info, 
             tiltle = "Cor.Res.CV4000_Unscale_Cell_info", print = F)
DimeReduPlot(mx = (Cor.Res.CV2000$Cor.merged), color = Pheno.merged$Cell_info, 
             tiltle = "Cor.Res.CV2000_Unscale_Cell_info", print = F)
DimeReduPlot(mx = (Cor.Res.CV1500$Cor.merged), color = Pheno.merged$Cell_info, 
             tiltle = "Cor.Res.CV1500_Unscale_Cell_info", print = F)
DimeReduPlot(mx = (Cor.Res.CV1000$Cor.merged), color = Pheno.merged$Cell_info, 
             tiltle = "Cor.Res.CV1000_Unscale_Cell_info", print = F)
# Cell CellType
DimeReduPlot(mx = (Cor.Res.CV8000$Cor.merged), color = Pheno.merged$CellType, 
             tiltle = "Cor.Res.CV8000_Unscale_Cell_Type", print = F)
DimeReduPlot(mx = (Cor.Res.CV4000$Cor.merged), color = Pheno.merged$CellType, 
             tiltle = "Cor.Res.CV4000_Unscale_Cell_Type", print = F)
DimeReduPlot(mx = (Cor.Res.CV2000$Cor.merged), color = Pheno.merged$CellType, 
             tiltle = "Cor.Res.CV2000_Unscale_Cell_Type", print = F)
DimeReduPlot(mx = (Cor.Res.CV1500$Cor.merged), color = Pheno.merged$CellType, 
             tiltle = "Cor.Res.CV1500_Unscale_Cell_Type", print = F)
DimeReduPlot(mx = (Cor.Res.CV1000$Cor.merged), color = Pheno.merged$CellType, 
             tiltle = "Cor.Res.CV1000_Unscale_Cell_Type", print = F)



DimeReduPlot(mx = (Cor.Res.CV8000$Cor.merged)^4, color = Pheno.merged$Cell_info, 
             tiltle = "Cor.Res.CV8000_Unscale_Cell_info", print = T)

DimeReduPlot(mx = (Cor.Res.CV1500$Cor.merged)^4, color = Pheno.merged$Cell_info, 
             tiltle = "Cor.Res.CV1500_Unscale_Cell_info", print = T)
DimeReduPlot(mx =scale((Cor.Res.CV1500$Cor.merged)^4), color = Pheno.merged$Cell_info, 
             tiltle = "Cor.Res.CV1500_Unscale_Cell_info", print = T)

DimeReduPlot(mx =scale((Cor.Res.CV2000$Cor.merged)^4), color = Pheno.merged$Cell_info, 
             tiltle = "Cor.Res.CV2000_Unscale_Cell_info", print = T)
