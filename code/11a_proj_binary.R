#' =============================================================================
#' @name proj_binary
#' @description computes spatial projections for the binary data sub-pipeline
#' @param QUERY the query object from the master pipeline
#' @param MODELS the models object from the master pipeline
#' @param N_BOOTSTRAP number of bootstrap to do for the projections
#' @param PROJ_PATH (optional) path to a environmental raster, potentially 
#' different than the one given in the QUERY object. This is the case for 
#' supplementary projections in other time and climate scenarios for example. 
#' To your own risk and only for expert users !
#' @return an updated model list object containing the projections objects
#' embedded in each model sub-list.

proj_binary <- function(QUERY,
                        MODEL,
                        CALL,
                        N_BOOTSTRAP,
                        PROJ_PATH = NULL,
                        CUT = NULL){
  
  # --- 1. Load environmental data - TO FIX DYNAMICALLY
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
  
  # =========================== MODEL LOOP SECTION =============================
  for(i in MODEL$CALL$MODEL_LIST){
    
    # --- 3. Fit model on bootstrap
    # fit_resamples() does not save models by default. Thus the control_resamples()
    boot_fit <- MODEL[[i]][["final_wf"]] %>% 
      fit_resamples(resamples = boot_split,
                    control = control_resamples(extract = function (x) extract_fit_parsnip(x))) %>% 
      unnest(.extracts)

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
    r_val <- raster(CALL$ENV_PATH) %>% 
      getValues()
    
    # Assign the desired values to the non-NA cells in the list
    y_hat <- apply(tmp, 2, function(x){
      r <- r_val
      r[!is.na(r)] <- x
      x <- r
    })
    
    # --- 6. Cut spatial discontinuities
    if(!is.null(CUT)){
      r0 <- raster(CALL$ENV_PATH)
      tmp <- apply(y_hat, 2, function(x){
        # --- 6.1. Open presence data
        xy <- QUERY$S %>% 
          dplyr::select(decimallongitude, decimallatitude)
        xy <- xy[which(QUERY$Y == 1),] # specific to presence data
        
        # --- 6.2. Cut y_hat
        x[x < CUT] <- 0
        r <- setValues(r0, x)
        
        # --- 6.3. Define patches and overlap with presence points
        r_patch <- clump(r)
        id_patch <- r_patch %>% 
          raster::extract(xy) %>% unique() %>% 
          .[!is.na(.)]

        # --- 6.4. Subset values from a patch overlapping with presences
        r_patch[!(r_patch %in% id_patch)] <- 0
        r_patch <- getValues(r_patch)
        r_patch[r_patch > 0] <- 1
        x <- x*r_patch
        return(x)
      })
      y_hat <- tmp
      
    } # if CUT
    
    # --- 7. Compute the average CV across bootstrap runs as a QC
    AVG_CV <- apply(y_hat, 1, function(x)(x = cv(x, na.rm = TRUE))) %>% 
      mean(na.rm = TRUE)
    
    # --- 8. Append the MODEL object
    MODEL[[i]][["proj"]][["y_hat"]] <- y_hat
    MODEL[[i]][["eval"]][["AVG_CV"]] <- AVG_CV
    
  } # for i model loop
  
  return(MODEL)
  
} # END FUNCTION


