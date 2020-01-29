### Distance normalization
#### Raw Correlation value transformed into correlation percent for each cell
########## Convert cor value to rang From 1 to 0 max was 1 min was 0
Min_max_norm.1000<- base::apply(Cor.Res.CV1000$Cor.merged, 2, function(x){
  (x-min(x))/(max(x)-min(x))
})
Min_max_norm.1500<- base::apply(Cor.Res.CV1500$Cor.merged, 2, function(x){
  (x-min(x))/(max(x)-min(x))
})
Min_max_norm.2000<- base::apply(Cor.Res.CV2000$Cor.merged, 2, function(x){
  (x-min(x))/(max(x)-min(x))
})
Min_max_norm.4000<- base::apply(Cor.Res.CV4000$Cor.merged, 2, function(x){
  (x-min(x))/(max(x)-min(x))
})
Min_max_norm.8000<- base::apply(Cor.Res.CV8000$Cor.merged, 2, function(x){
  (x-min(x))/(max(x)-min(x))
})

