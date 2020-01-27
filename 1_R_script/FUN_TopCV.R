###Function of select Top cv features
TopCV <- function(df, TopN = 10, MARGIN = 1){
  FeatureCV <- apply(df, MARGIN , function(x) {sd(x) / mean(x)})
  Feature.TopN <- names(head(sort(FeatureCV, decreasing = T), TopN))
  df[Feature.TopN,]
}
