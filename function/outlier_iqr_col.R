#Outlier exclusion
outlier_iqr_col <- function(y,n) {
  #' @description Function to identify outliers based on z-score. 
  #' Calculated on non-zero values
  #' Adapted from Nielja's code to be used on multiple target matrices
  #' @param n maximum allowed deviation in number of SD
  
  # --- 1. Remove zeros from calculations (i.e. later na.rm = TRUE)
  y[y == 0] <- NA
  
  # --- 2. Calculate score
  score <- apply(y, 2, function(x)(y = (x-mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)))
  
  # --- 3. Identify row numbers over the threshold
  outliers <- apply(score, 2, function(x)(x = which(x > n))) %>% 
    unlist() %>% 
    unique()
  
  return(outliers)
}
