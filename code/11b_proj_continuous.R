#' =============================================================================
#' @name proj_continuous
#' @description computes spatial projections for the continuous data sub-pipeline
#' @param CALL the call object from the master pipeline
#' @param QUERY the query object from the master pipeline
#' @param MODEL the models object from the master pipeline
#' @return an updated model list object containing the projections objects
#' embedded in each model sub-list.

proj_continuous <- function(QUERY,
                            MODEL,
                            CALL){
  
  # --- 1. Load environmental data
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
  boot_split <- bootstraps(tmp, times = CALL$N_BOOTSTRAP)
  
  # --- 3. Define the projections to compute
  # All algorithms if FAST == FALSE; only the ones that passed QC otherwise
  if(CALL$FAST == FALSE){
    loop_over <- CALL$HP$MODEL_LIST
  } else {
    loop_over <- MODEL$MODEL_LIST
  }
  
  for(i in loop_over){
    # --- 4. Fit model on bootstrap
    # fit_resamples() does not save models by default. Thus the control_resamples()
    boot_fit <- MODEL[[i]][["final_wf"]] %>% 
      fit_resamples(resamples = boot_split,
                    control = control_resamples(extract = function (x) extract_fit_parsnip(x))) %>% 
      unnest(.extracts)
    
    # --- 5. Compute one prediction per bootstrap
    # As we extracted the model information in a supplementary column, we can
    # directly compute the bootstrap within the synthetic resample object.
    boot_proj <- boot_fit %>% 
      mutate(proj = map(.extracts, function(x)(x = predict(x, features))))
    
    # --- 6. Compute average and CV across bootstraps
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
    
    
    # --- 7. Cut spatial discontinuities
    if(!is.null(CALL$CUT)){
      r0 <- raster(CALL$ENV_PATH)
      tmp <- apply(y_hat, 2, function(x){
        # --- 7.1. Open presence data
        xy <- QUERY$S %>% 
          dplyr::select(decimallongitude, decimallatitude)
        xy <- xy[which(QUERY$Y != 0),] # specific to presence data
        
        # --- 7.2. Cut y_hat
        x[x < CALL$CUT] <- 0
        r <- setValues(r0, x)
        
        # --- 7.3. Define patches and overlap with presence points
        r_patch <- clump(r)
        id_patch <- r_patch %>% 
          raster::extract(xy) %>% unique() %>% 
          .[!is.na(.)]
        
        # --- 7.4. Subset values from a patch overlapping with presences
        r_patch[!(r_patch %in% id_patch)] <- 0
        r_patch <- getValues(r_patch)
        r_patch[r_patch > 0] <- 1
        x <- x*r_patch
        return(x)
      })
      y_hat <- tmp
      
    } # if CUT
    
    # --- 8. Compute the average CV across bootstrap runs as a QC
    AVG_CV <- apply(y_hat, 1, function(x)(x = cv(x, na.rm = TRUE))) %>% 
      mean(na.rm = TRUE)
    
    # --- 9. Append the MODEL object
    MODEL[[i]][["proj"]][["y_hat"]] <- y_hat
    MODEL[[i]][["eval"]][["AVG_CV"]] <- AVG_CV
    
  } # for i model loop
  
  return(MODEL)
  
} # END FUNCTION