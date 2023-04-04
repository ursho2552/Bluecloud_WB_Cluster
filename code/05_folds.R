#' =============================================================================
#' @name folds
#' @description this functions adds an ID list to the query. It corresponds to
#' the different folds between train and test splits
#' @param QUERY query object, for tracking
#' @param NFOLD number of folds, used defined integer
#' @param FOLD_METHOD method used to create the folds, integer between "kfold"
#' and "lon"
#' @return CALL$FOLD_METHOD for tracking
#' @return ID : list of line id corresponding to the n defined folds

folds <- function(QUERY = query,
                  NFOLD = 5,
                  FOLD_METHOD = "kfold"){
  
  # ============================ INITIAL SPLIT =================================
  # --- 1. Re-assemble all query tables
  tmp <- cbind(QUERY$Y, QUERY$X, QUERY$S)
  
  # --- 1. Do the initial split
  init_split <- tmp %>% 
    initial_split(prop = 0.8)
  
  QUERY$FOLDS$init_split <- init_split
  
  QUERY$FOLDS$test <- testing(init_split)
  QUERY$FOLDS$train <- training(init_split)
  
  # ====================== TRAIN SET RESAMPLING ================================
  # --- 0. Parameter check
  if(FOLD_METHOD != "kfold" & FOLD_METHOD != "lon"){
    stop("FOLD_METHOD not implemented or incorrect. It should be 'kfold' or 'lon'")
  }
  
  # --- 1. Normal k-fold re sampling
  if(FOLD_METHOD == "kfold"){
    folds <- vfold_cv(data = QUERY$FOLDS$train,
                      v = NFOLD)
  }
  
  
  # --- 2. Longitudinal block re sampling
  if(FOLD_METHOD == "lon"){
    folds <- clustering_cv(data = QUERY$FOLDS$train,
                           vars = c(decimallongitude),
                           v = NFOLD)
  }
  
  # =================== APPEND QUERY OBJECT ====================================
  # S$id keeps track of the initial row numbers within each split/re sample
  QUERY$FOLDS$resample_split <- folds
  
  for(i in 1:nrow(folds)){
    fold_name <- folds$id[i]
    QUERY$FOLDS[["resample_folds"]][[fold_name]][["assessment"]] <- folds$splits[[i]] %>% assessment()
    QUERY$FOLDS[["resample_folds"]][[fold_name]][["analysis"]] <- folds$splits[[i]] %>% analysis()
  }
  
  QUERY$CALL$FOLD_METHOD <- FOLD_METHOD
  
  return(query = QUERY)
  
} # END FUNCTION
