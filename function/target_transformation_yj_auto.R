
#' @name target_transformation
#' @param x vector of target values
#' @param REVERSE boolean, reverse the transformation or not ?
#' @param PARAM parameters used for the transformation if reverse TRUE
#' @return the transformed values of x

target_transformation <- function(x, REVERSE = FALSE, PARAM = NULL){
  library(bestNormalize)
  
  # --- 1. Forward transformation /w. REVERSE = FALSE
  if(REVERSE == FALSE){
    # --- 1.1. Build the yeo-johnson object with optimal parameters
    yeojohnson_obj <- yeojohnson(x)
    
    # --- 1.2. Apply the transformation
    x_out <- predict(yeojohnson_obj)
    
    # --- 1.3. Wrap up
    return(list(out = x_out, yj_obj = yeojohnson_obj))
  } # end reverse false
  
  # --- 2. Reverse transformation
  if(REVERSE == TRUE){
    # --- 2.1. Load transformation parameters from query object
    yeojohnson_obj <- PARAM$yj_obj
    
    # --- 2.2. Apply the transformation
    x_out <- predict(yeojohnson_obj, newdata = x, inverse = TRUE)
    
    # --- 2.3. Wrap up
    return(x_out)
    
  } # end reverse true
  
} # end function
