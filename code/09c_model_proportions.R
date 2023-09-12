#' =============================================================================
#' @name model_proportions
#' @description sub-pipeline corresponding to the model fitting procedure for
#' proportions data. Called via the model_wrapper.R function, thus no default
#' parameter values are defined.
#' @param CALL the call object from the master pipeline. 
#' @param QUERY the query object from the master pipeline
#' @return returns a path to the python model file

model_proportions <- function(CALL,
                              QUERY){
  
  # --- 1. Initialize function
  # --- 1.1. Storage in MODEL object
  MODEL <- CALL$HP
  
  # --- 1.2. Source the MBTR functions
  source_python(paste0(project_wd,"/function/mbtr_function.py"))
  
  # --- 1.3. Create a clean cache folder in data
  # i.e. to write train and test files to pass to python (necessary for parallel)
  dir.create(paste0(project_wd, "/data/MBTR_cache"))
  to_clean <- list.files(paste0(project_wd, "/data/MBTR_cache"), full.names = TRUE)
  file.remove(to_clean)
  
  # --- 2. Run the model for each fold x hyper parameter
  # --- 2.1. Initialize loss & hp object
  loss <- list()
  hp <- 1:nrow(MODEL$MBTR$model_grid)
  
  # --- 2.2. Initialize the hyperparameter fit function
  # i.e. to use mcmapply across hyperparameters for faster computing
  mbtr_hp_fit <- function(hp, path){
    m <- mbtr_fit(path,
                  loss_type='mse',
                  n_boosts = as.integer(1000),
                  min_leaf= MODEL$MBTR$model_grid$MEAN_LEAF[hp],
                  learning_rate=MODEL$MBTR$model_grid$LEARNING_RATE[hp],
                  lambda_weights=MODEL$MBTR$model_grid$LEARNING_RATE[hp]/100,
                  lambda_leaves=0,
                  n_q= as.integer(MODEL$MBTR$model_grid$N_Q[hp]),
                  early_stopping_rounds = 10)
    tmp <- m[[2]] %>% unlist()
  } # end function
  
  # --- 2.3. Do the resample fit
  for(cv in 1:CALL$NFOLD){
    
    # --- 2.3.1. Prepare the inputs
    message(paste0(Sys.time(), "--- MBTR: writing data to cache for cv = ", cv))
    X_tr <- QUERY$FOLDS$resample_folds[[cv]]$analysis %>%
      dplyr::select(QUERY$SUBFOLDER_INFO$ENV_VAR)
    write_feather(X_tr, paste0(project_wd, "/data/MBTR_cache/", cv, "_X_tr.feather"))
    
    Y_tr <- QUERY$FOLDS$resample_folds[[cv]]$analysis %>%
      dplyr::select(as.character(CALL$SP_SELECT))
    write_feather(Y_tr, paste0(project_wd, "/data/MBTR_cache/", cv, "_Y_tr.feather"))
    
    X_val <- QUERY$FOLDS$resample_folds[[cv]]$assessment %>%
      dplyr::select(QUERY$SUBFOLDER_INFO$ENV_VAR)
    write_feather(X_val, paste0(project_wd, "/data/MBTR_cache/", cv, "_X_val.feather"))
    
    Y_val <- QUERY$FOLDS$resample_folds[[cv]]$assessment %>%
      dplyr::select(as.character(CALL$SP_SELECT))
    write_feather(Y_val, paste0(project_wd, "/data/MBTR_cache/", cv, "_Y_val.feather"))
    
    # --- 2.3.2. Fitting hyperparameters
    message(paste0(Sys.time(), "--- MBTR: fit hyperparameters for cv = ", cv))
    loss[[cv]] <- mcmapply(FUN = mbtr_hp_fit,
                           path = paste0(project_wd, "/data/MBTR_cache/",cv,"_"),
                           hp = hp,
                           mc.cores = nrow(MODEL$MBTR$model_grid))
    message("Algorithm fitting - DONE")
  } # cv loop
  
  # --- 2.4. Compute the minimum loss and corresponding nb. of boost rounds
  min_loss <- nboost <- list()
  for(i in 1:nrow(MODEL$MBTR$model_grid)){
    max_boost <- lapply(loss, function(x)(x = length(x[[i]]))) %>% unlist() %>% min()
    min_loss[[i]] <- lapply(loss, function(x)(x = x[[i]][1:max_boost])) %>% as.data.frame() %>% apply(1, mean) %>% min()
    nboost[[i]] <- which(lapply(loss, function(x)(x = x[[i]][1:max_boost])) %>% as.data.frame() %>% apply(1, mean) == min_loss[[i]])
  }
  best_hp <- which(unlist(min_loss) == min(unlist(min_loss)))
  nboost <- nboost[[best_hp]]
  
  # --- 2.5. Retrieve the corresponding RMSE as well
  MODEL[["MBTR"]][["best_fit"]] <- sqrt(min_loss[[best_hp]])
  
  # --- 3. Final model fit
  # --- 3.1. Extract train and validation fold from initial split
  message(paste0(Sys.time(), "--- MBTR: fit final model"))
  X_tr <- QUERY$FOLDS$train %>% 
    dplyr::select(QUERY$SUBFOLDER_INFO$ENV_VAR)
  write_feather(X_tr, paste0(project_wd, "/data/MBTR_cache/0_X_tr.feather"))
  
  Y_tr <- QUERY$FOLDS$train %>% 
    dplyr::select(as.character(CALL$SP_SELECT))
  write_feather(Y_tr, paste0(project_wd, "/data/MBTR_cache/0_Y_tr.feather"))
  
  X_val <- QUERY$FOLDS$test %>% 
    dplyr::select(QUERY$SUBFOLDER_INFO$ENV_VAR)
  write_feather(X_val, paste0(project_wd, "/data/MBTR_cache/0_X_val.feather"))
  
  Y_val <- QUERY$FOLDS$test %>% 
    dplyr::select(as.character(CALL$SP_SELECT))
  write_feather(Y_val, paste0(project_wd, "/data/MBTR_cache/0_Y_val.feather"))
  
  # --- 3.2. Fit the final model
  final_fit <- mbtr_fit(path = paste0(project_wd, "/data/MBTR_cache/0_"),
                        loss_type='mse',
                        n_boosts = as.integer(100),
                        min_leaf= MODEL$MBTR$model_grid$MEAN_LEAF[best_hp],
                        learning_rate=MODEL$MBTR$model_grid$LEARNING_RATE[best_hp],
                        lambda_weights=MODEL$MBTR$model_grid$LEARNING_RATE[best_hp]/100,
                        lambda_leaves=0,
                        n_q= as.integer(MODEL$MBTR$model_grid$N_Q[best_hp]),
                        early_stopping_rounds = 10)
  
  # --- 3.3. Write model in a file
  py_save_object(final_fit, paste0(project_wd, "/data/MBTR_cache/final_fit"), pickle = "pickle")
  
  # --- 3.4. Pass the file path
  MODEL[["MBTR"]][["final_wf"]] <- MODEL$MBTR$model_grid[best_hp,] %>% mutate(NBOOST = nboost)
  MODEL[["MBTR"]][["final_fit"]] <- paste0(project_wd, "/data/MBTR_cache/final_fit")
  return(MODEL)
  
} # END FUNCTION
