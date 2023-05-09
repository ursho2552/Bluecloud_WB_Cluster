#' =============================================================================
#' @name folds
#' @description this functions adds an ID list to the query. It corresponds to
#' the different folds between train and test splits
#' @param SP_SELECT species to run the analysis for, in form of Aphia ID
#' @param FOLDER_NAME name of the corresponding folder
#' @param NFOLD number of folds, used defined integer
#' @param FOLD_METHOD method used to create the folds, integer between "kfold"
#' and "lon"
#' @return CALL$FOLD_METHOD for tracking
#' @return ID : list of line id corresponding to the n defined folds
#' @return Updates the output in a QUERY.RData and CALL.Rdata files

folds <- function(SP_SELECT = NULL,
                  FOLDER_NAME = NULL,
                  NFOLD = 5,
                  FOLD_METHOD = "kfold"){
  
  # =========================== PARAMETER LOADING ==============================
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/", SP_SELECT, "/QUERY.RData"))
  
  # ============================ INITIAL SPLIT =================================
  # --- 1. Re-assemble all query tables
  tmp <- cbind(QUERY$Y, QUERY$X, QUERY$S)
  
  # --- 1. Do the initial split
  init_split <- tmp %>% 
    initial_split(prop = 0.8,
                  strata = measurementvalue)
  
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
                      strata = measurementvalue,
                      v = NFOLD)
  }
  
  
  # --- 2. Longitudinal block re sampling
  if(FOLD_METHOD == "lon"){
    folds <- clustering_cv(data = QUERY$FOLDS$train,
                           vars = c(decimallongitude),
                           v = NFOLD)
  }
  
  # =================== APPEND QUERY and CALL OBJECT ===========================
  # S$id keeps track of the initial row numbers within each split/re sample
  
  # --- 1. Append object
  QUERY$FOLDS$resample_split <- folds
  
  for(i in 1:nrow(folds)){
    fold_name <- folds$id[i]
    QUERY$FOLDS[["resample_folds"]][[fold_name]][["assessment"]] <- folds$splits[[i]] %>% assessment()
    QUERY$FOLDS[["resample_folds"]][[fold_name]][["analysis"]] <- folds$splits[[i]] %>% analysis()
  }
  
  # --- 2. Save query
  save(QUERY, file = paste0(project_wd, "/output/", FOLDER_NAME,"/", SP_SELECT, "/QUERY.RData"))
  
  # --- 3. Append CALL object
  CALL$FOLD_METHOD <- FOLD_METHOD
  CALL$NFOLD <- NFOLD
  save(CALL, file = paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  
} # END FUNCTION
