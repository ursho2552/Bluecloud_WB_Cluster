
# Custom viridis-like color palette

viridis_pal <- function(n){
  tmp <- c("#440154FF","#414487FF","#2A788EFF","#22A884FF","#7AD151FF","#FDE725FF")
  pal <- colorRampPalette(tmp)(n)
  
  return(pal)

} # END FUNCTION
