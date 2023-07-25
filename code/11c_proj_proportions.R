#' =============================================================================
#' @name proj_proportions
#' @description computes spatial projections for the proportions data sub-pipeline
#' @param QUERY the query object from the master pipeline
#' @param MODELS the models object from the master pipeline
#' @param N_BOOTSTRAP number of bootstrap to do for the projections
#' @param PROJ_PATH (optional) path to a environmental raster, potentially 
#' different than the one given in the QUERY object. This is the case for 
#' supplementary projections in other time and climate scenarios for example. 
#' To your own risk and only for expert users !
#' @return an updated model list object containing the projections objects
#' embedded in each model sub-list.

proj_proportions <- function(QUERY,
                             MODEL,
                             CALL,
                             N_BOOTSTRAP,
                             PROJ_PATH = NULL){
  
  # --- 1. Initialize function
  # --- 1.1. Source the MBTR functions
  source_python(paste0(project_wd,"/function/mbtr_function.py"))
  
  # --- 1.2. Load environmental data - TO FIX DYNAMICALLY
  features <- stack(CALL$ENV_PATH) %>% 
    readAll() %>% 
    raster::subset(QUERY$SUBFOLDER_INFO$ENV_VAR) %>% 
    rasterToPoints() %>% 
    as.data.frame() %>% 
    dplyr::select(-c(x, y))
  
  # --- 2. Define bootstraps
  # --- 2.1. Re-assemble all query tables
  tmp <- cbind(QUERY$Y, QUERY$X, QUERY$S)
  
  # --- 2.2. Run the bootstrap generation from tidy models
  boot_split <- bootstraps(tmp, times = N_BOOTSTRAP)
  
  # --- 3. Fit model on bootstrap
  # --- 3.1. Extract bootstrap input
  for(b in 1:N_BOOTSTRAP){
    X_tr <- boot_split$splits[[b]] %>% analysis() %>%
      dplyr::select(QUERY$SUBFOLDER_INFO$ENV_VAR)
    write_feather(X_tr, paste0(project_wd, "/data/MBTR_cache/", b, "_X_tr.feather"))
    
    Y_tr <- boot_split$splits[[b]] %>% analysis() %>%
      dplyr::select(as.character(CALL$SP_SELECT))
    write_feather(Y_tr, paste0(project_wd, "/data/MBTR_cache/", b, "_Y_tr.feather"))
    
    X_val <- boot_split$splits[[b]] %>% assessment() %>%
      dplyr::select(QUERY$SUBFOLDER_INFO$ENV_VAR)
    write_feather(X_val, paste0(project_wd, "/data/MBTR_cache/", b, "_X_val.feather"))
    
    Y_val <- boot_split$splits[[b]] %>% assessment() %>%
      dplyr::select(as.character(CALL$SP_SELECT))
    write_feather(Y_val, paste0(project_wd, "/data/MBTR_cache/", b, "_Y_val.feather"))
  } # for b
  
  # --- 3.2. Fit model
  # --- 3.2.1. Computing the model objects
  message(paste0(Sys.time(), "--- MBTR: re-fit model for bootstrap"))
  boot_fit <- mcmapply(FUN = mbtr_fit,
                       path = paste0(project_wd, "/data/MBTR_cache/",1:N_BOOTSTRAP,"_"),
                       loss_type='mse',
                       n_boosts = as.integer(1000),
                       min_leaf= MODEL$MBTR$final_wf$MEAN_LEAF,
                       learning_rate=MODEL$MBTR$final_wf$LEARNING_RATE,
                       lambda_weights=MODEL$MBTR$final_wf$LEARNING_RATE/100,
                       lambda_leaves=0,
                       n_q= as.integer(MODEL$MBTR$final_wf$N_Q),
                       early_stopping_rounds = 10,
                       SIMPLIFY = FALSE,
                       USE.NAMES = FALSE,
                       mc.cores = N_BOOTSTRAP)

  # --- 3.2.2. Reload them in R because of "previous session invalidity"
  for(b in 1:N_BOOTSTRAP){
    boot_fit[[b]] <- py_load_object(paste0(project_wd, "/data/MBTR_cache/",b,"_m"), pickle = "pickle")[[1]]
  }
  message("Bootstrap fitting - DONE")
  
  # --- 4. Computing one prediction per bootstrap
  boot_proj <- mclapply(boot_fit,
                        function(a_boot) mbtr_predict(model = a_boot, X_pred = features, n_boosts = as.integer(MODEL$MBTR$final_wf$NBOOST)),
                        mc.cores = N_BOOTSTRAP) %>% 
    abind(along = 3)
  
  # --- 5. Re-assign to cells
  # Open a raster to have the list of cells
  r_val <- raster(CALL$ENV_PATH) %>% 
    getValues()
  
  # Assign the desired values to the non-NA cells in the list
  y_hat <- apply(boot_proj, c(2,3), function(x){
    r <- r_val
    r[!is.na(r)] <- x
    x <- r
  }) %>% 
    aperm(c(1,3,2))
  
  # --- 6. Compute the average CV across bootstrap runs as a QC
  AVG_CV <- apply(y_hat, c(1,3), function(x)(x = cv(x, na.rm = TRUE))) %>% 
    mean(na.rm = TRUE)
  
  # --- 7. Append the MODEL object
  MODEL[["MBTR"]][["proj"]][["y_hat"]] <- y_hat
  MODEL[["MBTR"]][["eval"]][["AVG_CV"]] <- AVG_CV
  return(MODEL)
  
} # END FUNCTION

