#' =============================================================================
#' @name standard_maps
#' @description Simple function for building standard mean and CV maps from
#' projection data
#' @param FOLDER_NAME name of the corresponding folder
#' @param SUBFOLDER_NAME list of sub_folders to parallelize on.
#' @param ENSEMBLE if TRUE, computes the ensemble map ? 
#' @param MESS if TRUE, considers the mess analysis in the maps
#' @return plots mean and uncertainty maps per model or ensemble

standard_maps <- function(FOLDER_NAME = NULL,
                          SUBFOLDER_NAME = NULL,
                          ENSEMBLE = FALSE,
                          MESS = FALSE){
  
  # =========================== PARAMETER LOADING ==============================
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/QUERY.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/MODEL.RData"))
  
  # ============================= BUILDING MAPS ================================
  # With PDF saving
  pdf(paste0(project_wd,"/output/",FOLDER_NAME,"/",SUBFOLDER_NAME,"/standard_maps.pdf"))
  
  # --- 1. Build model-level outputs
  # --- 1.1. Set initial plot layout & requirements
  par(mfrow = c(4,2), mar = c(2,8,3,3))
  r <- raster(paste0(project_wd, "/data/features_mean_from_monthly"))
  
  for(i in MODEL$CALL$MODEL_LIST){
    # --- 1.2. Compute mean and CV
    # --- Mean value
    val <- MODEL[[i]][["proj"]][["y_hat"]] %>% 
      apply(1, function(x)(x = mean(x, na.rm = TRUE)))
    r_m <- r %>% 
      setValues(val)
    
    # --- Coefficient of variation
    val <- MODEL[[i]][["proj"]][["y_hat"]] %>% 
      apply(1, function(x)(x = cv(x, na.rm = TRUE)))
    r_cv <- r %>% 
      setValues(val)
    r_cv[r_cv > 100] <- 100
    
    # --- 1.3. Plot the corresponding maps
    plot(r_m, col = viridis_pal(100),
         main = paste("Average proj. for", i, "\n", 
                      names(MODEL[[i]][["eval"]]), "=", MODEL[[i]][["eval"]])
         )
    
    plot(r_cv, col = viridis_pal(100),
         main = paste("Uncertainty (CV) proj. for", i))
    
  } # End i model loop
  
  # --- 2. Build ensemble outputs
  if(ENSEMBLE == TRUE){
    # --- 2.1. Set initial plot layout & requirements
    par(mfrow = c(2,1), mar = c(2,8,3,3))
    r <- raster(paste0(project_wd, "/data/features_mean_from_monthly"))
    
    # --- 2.2. Build ensemble array, weighted by evaluation values, re-scale max=1
    y_ens <- NULL
    for(i in MODEL$CALL$MODEL_LIST){
      y_ens <- cbind(y_ens,
                     MODEL[[i]][["proj"]][["y_hat"]]*MODEL[[i]][["eval"]][[1]])
    }
    y_ens <- y_ens/max(y_ens, na.rm = TRUE)
    
    # --- 2.3. Compute mean and CV
    # --- Mean value
    val <- apply(y_ens, 1, function(x)(x = mean(x, na.rm = TRUE)))
    r_m <- r %>% 
      setValues(val)
    
    # --- Coefficient of variation
    val <- apply(y_ens, 1, function(x)(x = cv(x, na.rm = TRUE)))
    r_cv <- r %>% 
      setValues(val)
    r_cv[r_cv > 100] <- 100
    
    # --- 2.4. Plot the corresponding maps
    plot(r_m, col = viridis_pal(100),
         main = "Average ensemble proj.")
    
    plot(r_cv, col = viridis_pal(100),
         main = "Average ensemble uncertainty (CV)")
  } # End if ENSEMBLE = TRUE
  
  # Stop pdf saving
  dev.off()
  
} # END FUNCTION

