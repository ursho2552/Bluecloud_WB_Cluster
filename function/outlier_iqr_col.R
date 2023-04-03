#Outlier exclusion
outlier_iqr_col <- function(x,n) {
  #' Function to identify outliers based on z-score
  #' @param n maximum allowed deviation in number of SD
  
  #Evaluate outliers based on z-score
  # y <- (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
  y <- x - apply(x, 2, function(x)(x = mean(x, na.rm = TRUE)))/ apply(x, 2, function(x)(x = sd(x, na.rm = TRUE)))
  x[y > n] <- NA
  
  return(x)
}
