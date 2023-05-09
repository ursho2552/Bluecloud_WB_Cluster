#' =============================================================================
#' @name proj_pres
#' @description computes spatial projections for the presence data sub-pipeline
#' @param QUERY the query object from the master pipeline
#' @param MODELS the models object from the master pipeline
#' @param N_BOOTSTRAP number of bootstrap to do for the projections
#' @param PROJ_PATH (optional) path to a environmental raster, potentially 
#' different than the one given in the QUERY object. This is the case for 
#' supplementary projections in other time and climate scenarios for example. 
#' To your own risk and only for expert users !
#' @return an updated model list object containing the projections objects
#' embedded in each model sub-list.

proj_pres <- function(QUERY,
                      MODEL,
                      CALL,
                      N_BOOTSTRAP,
                      PROJ_PATH){
  
  # --- 1. Load environmental data - TO FIX DYNAMICALLY
  features <- stack(paste0(project_wd, "/data/features_mean_from_monthly")) %>% 
    readAll() %>% 
    raster::subset(CALL$ENV_VAR) %>% 
    rasterToPoints() %>% 
    as.data.frame() %>% 
    dplyr::select(-c(x, y))
  
  # --- 2. Define bootstraps
  # --- 2.1. Re-assemble all query tables
  tmp <- cbind(QUERY$Y, QUERY$X, QUERY$S)
  
  # --- 2.2. Run the bootstrap generation from tidy models
  boot_split <- bootstraps(tmp, times = N_BOOTSTRAP)
  
  # =========================== MODEL LOOP SECTION =============================
  for(i in MODEL$CALL$MODEL_LIST){
    
    # --- 3. Fit model on bootstrap
    # --- 3.1. Register parallel
    # Only if the number of species is less then the number of available clusters
    # Otherwise, the parallel computing is done by species
    if(length(CALL$SP_SELECT) < LOCAL_CLUSTERS){
      cl <- makePSOCKcluster(LOCAL_CLUSTERS)
      doParallel::registerDoParallel(cl)
      message(paste(Sys.time(), "--- Parallel bootstrap for", i, ": START"))
    }
    
    # --- 3.2. Fit
    # fit_resamples() does not save models by default. Thus the control_resamples()
    boot_fit <- MODEL[[i]][["final_wf"]] %>% 
      fit_resamples(resamples = boot_split,
                    control = control_resamples(extract = function (x) extract_fit_parsnip(x),
                                                allow_par = TRUE,
                                                parallel_over = "everything")) %>% 
      unnest(.extracts)
    
    # --- 3.3. Stop parallel backend
    stopCluster(cl)
    message(paste(Sys.time(), "--- Parallel grid tuning for", i, ": DONE"))
    
    # --- 4. Compute one prediction per bootstrap
    # As we extracted the model information in a supplementary column, we can
    # directly compute the bootstrap within the synthetic resample object.
    boot_proj <- boot_fit %>% 
      mutate(proj = map(.extracts, function(x)(x = predict(x, features))))
    
    # --- 5. Compute average and CV across bootstraps
    # First transform the object into a cell x bootstrap matrix
    # /!\ Need to create a unique row identifier for pivot_wider to work...
    tmp <- boot_proj %>% 
      dplyr::select(id, proj) %>% 
      unnest(c(id, proj)) %>% 
      as.data.frame() %>% 
      group_by(id) %>%
      mutate(row = row_number()) %>%
      pivot_wider(names_from = id, values_from = .pred) %>%
      dplyr::select(-row)
    
    # Open a raster to have the list of cells
    r_val <- raster(paste0(project_wd, "/data/features_mean_from_monthly")) %>% 
      getValues()
    
    # Assign the desired values to the non-NA cells in the list
    y_hat <- apply(tmp, 2, function(x){
      r <- r_val
      r[!is.na(r)] <- x
      x <- r
    })
    
    # --- 6. Append the MODEL object
    MODEL[[i]][["proj"]][["y_hat"]] <- y_hat
    
  } # for i model loop
  
  return(MODEL)
  
} # END FUNCTION


