#### Heatmap plot function using pheatmap
heatmap.JH <- function(mx,...){
  require(pheatmap)
  ### Set the color scheme
  color_scheme = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                    "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(100);
  
  pheatmap(mx,color = color_scheme,...)
}
