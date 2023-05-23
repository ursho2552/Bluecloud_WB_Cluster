#' Functions for bivariate raster plotting in R
#' 
#' @function colmat :
#' @return : two-dimensional color matrix
#' 
#' @param pal : a color palette used for the first dimension variable
#' @param saturation : the degree of saturation related to the second variable
#' @param xlab @param ylab : second (x) and first (y) variable names
#' 
#' @function colmat_plot
#' @return : a plot of the color matrix
#' 
#' @param colourmatrix : output of colmat
#' @param xlab @param ylab : second (x) and first (y) variable names
#' 
#' @function bivar_map
#' @return a raster with the cell values corresponding to the index of colmat
#' 
#' @param rasterx : raster corresponding to variable 1, same format as rastery
#' @param rastery : raster corresponding to variable 2, same format as rasterx
#' @param colourmatrix : output of colmat
#' @param cutx : user defined interval for the color classes of x
#' @param cuty : user defined interval for the color classes of y
#' 
#' @references
#' Simplified version of :
#' https://gist.github.com/scbrown86/2779137a9378df7b60afd23e0c45c188#file-bivarrasterplot-r

# Function that produces the color matrix
colmat <- function(xmax = "deepskyblue4",
                   ymax = "red4",
                   nbreaks = 100){
  
  # --- 1. Required Packages
  require(RColorBrewer)
  require(colorspace)
  
  # --- 2. Build the two extremes color columns
  left_pal <- colorRampPalette(c("white",ymax))(nbreaks)
  right_pal <- colorRampPalette(c(ymax, xmax))(2*nbreaks-1) %>% 
    .[-c(1:nbreaks-1)] %>% 
    rev()
  
  # --- 3. Interpolate colors in between
  col_matrix <- NULL
  for(i in 1:nbreaks){
    col_matrix <- rbind(col_matrix, colorRampPalette(colors = c(left_pal[i], right_pal[i]))(nbreaks))
  }
  return(col_matrix)
}

# Plotting the color matrix

colmat_plot <- function(colormatrix, xlab, ylab){
  image(x = matrix(seq(1,length(colormatrix)), nrow = nrow(colormatrix), ncol = ncol(colormatrix), byrow = TRUE),
        col = colormatrix,
        main = "",
        xlab = xlab, ylab = ylab,
        axes = FALSE)
  box()
  grid(lty = "solid", col = "black")
}

# Function to assign color-code to raster data
bivar_map <- function(rasterx, rastery, colormatrix, cutx = NULL, cuty = NULL){
  require(raster)
  require(virtualspecies)
  
  colorid <- matrix(seq(1:length(colormatrix)), nrow=nrow(colormatrix), ncol=ncol(colormatrix), byrow=TRUE)
  
  splitx <- cut(rasterx, breaks = nrow(colormatrix))
  splity <- cut(rastery, breaks = ncol(colormatrix))
  
  if(!is.null(cutx) & !is.null(cuty)){
    if(length(cutx) == (ncol(colormatrix)+1) & length(cuty) == (nrow(colormatrix)+1)){
      # Add security if cutx or cuty has a smaller range than the raster
      rasterx[rasterx<min(cutx) | rasterx>max(cutx)] <- NA
      rastery[rastery<min(cuty) | rastery>max(cuty)] <- NA
      
      splitx <- cut(rasterx, breaks = cutx)
      splity <- cut(rastery, breaks = cuty)
    } else {
      stop("cutx and cuty dimension should be respectively nrow(colormatrix)+1 and ncol(colormatrix)+1")
    }
  } # end if
  z <- setValues(rasterx,colorid[cbind(getValues(splitx),getValues(splity))])
  col_plot <- colormatrix[min(getValues(z), na.rm = TRUE):max(getValues(z), na.rm = TRUE)]
  
  return(list(z, col_plot))
}
