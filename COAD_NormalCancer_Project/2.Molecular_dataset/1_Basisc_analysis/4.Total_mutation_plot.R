#### 4.Total_mutation_plot.R
# Merge table
MergeTable.totalMut <- dplyr::left_join(Cluster.df, COAD.mc3.PUBLIC_sampleSummary, by = "rownames")
# Plot 
library(ggplot2)
# Remove NA value
data=subset(MergeTable.totalMut, !is.na(MergeTable.totalMut$total))
# ggplot
xlabs <- paste(levels(data$dynamicColors),"\n(N=",table(data$dynamicColors),")",sep="")
p<- ggplot(data,aes(x = dynamicColors, y = total,fill=dynamicColors))+
  geom_boxplot()+scale_x_discrete(labels=xlabs) + theme_minimal()+
  scale_fill_manual(values=c("blue","brown","turquoise","yellow")) +
  labs(title =" MC3 totalMut")

#### 
library(scales)
squish_trans <- function(from, to, factor) {
  
  trans <- function(x) {
    
    # get indices for the relevant regions
    isq <- x > from & x < to
    ito <- x >= to
    
    # apply transformation
    x[isq] <- from + (x[isq] - from)/factor
    x[ito] <- from + (to - from)/factor + (x[ito] - to)
    
    return(x)
  }
  
  inv <- function(x) {
    
    # get indices for the relevant regions
    isq <- x > from & x < from + (to - from)/factor
    ito <- x >= from + (to - from)/factor
    
    # apply transformation
    x[isq] <- from + (x[isq] - from) * factor
    x[ito] <- to + (x[ito] - (from + (to - from)/factor))
    
    return(x)
  }
  
  # return the transformation
  return(trans_new("squished", trans, inv))
}
p +   scale_y_continuous(trans = squish_trans(2500, 10000, 20),
                         breaks = c(0, 100, 300,500,1000 ,1500,2000, 2500,5000, 10000,10000))


