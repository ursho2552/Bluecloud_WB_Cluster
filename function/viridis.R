
# Custom viridis-like color palette(s)

viridis_pal <- function(n){
  tmp <- c("#440154FF","#414487FF","#2A788EFF","#22A884FF","#7AD151FF","#FDE725FF")
  pal <- colorRampPalette(tmp)(n)
  
  return(pal)
} # END FUNCTION

inferno_pal <- function(n){
  tmp <- c("#000004ff","#1b0b40ff","#480e64ff","#761b6bff","#a32d5fff","#cc4546ff","#ea6827ff","#f99913ff","#f9cd3aff","#fcffa4ff")
  pal <- colorRampPalette(tmp)(n)
  
  return(pal)
} # END FUNCTION
