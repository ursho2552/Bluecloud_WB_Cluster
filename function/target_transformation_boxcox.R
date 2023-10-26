
#' @name target_transformation
#' @param x vector of target values
#' @param REVERSE boolean, reverse the transformation or not ?
#' @param PARAM parameters used for the transformation if reverse TRUE
#' @return the transformed values of x

target_transformation <- function(x, REVERSE = FALSE, PARAM = NULL){
  # --- 1. Forward transformation /w. REVERSE = FALSE
  if(REVERSE == FALSE){
    # --- 1.1. Define transformation parameters
    LAMBDA <- 0.2
    GAMMA <- 1
    # if(length(which(x <= 0)) == 0){GAMMA <- 0}
    
    # --- 1.2. Apply the transformation
    x_out <- bcnPower(x, lambda = LAMBDA, gamma = GAMMA)
    
    # --- 1.3. Wrap up
    return(list(out = x_out, LAMBDA = LAMBDA, GAMMA = GAMMA))
  } # end reverse false
  
  # --- 2. Reverse transformation
  if(REVERSE == TRUE){
    # --- 2.1. Load transformation parameters from query object
    LAMBDA <- PARAM$LAMBDA
    GAMMA <- PARAM$GAMMA
    
    # --- 2.2. Apply the transformation
    x_out <- bcnPowerInverse(x, lambda = LAMBDA, gamma = GAMMA)
    
    # --- 2.3. Wrap up
    return(x_out)
    
  } # end reverse true
  
} # end function
