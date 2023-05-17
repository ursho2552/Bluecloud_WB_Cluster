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
  
  # --- 1. Initialize function
  # --- 1.1. Parameter loading
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/QUERY.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/MODEL.RData"))
  
  # --- 1.2. Create PDF saving
  pdf(paste0(project_wd,"/output/",FOLDER_NAME,"/",SUBFOLDER_NAME,"/standard_maps.pdf"))
  
  # --- 2. Build model-level outputs
  # --- 2.1. Set initial plot layout & requirements
  par(mfrow = c(4,2), mar = c(2,8,3,3))
  r <- raster(paste0(project_wd, "/data/features_mean_from_monthly"))
  
  for(i in MODEL$CALL$MODEL_LIST){
    # --- 2.2. Compute mean and CV
    # --- 2.2.1. Mean value
    val <- MODEL[[i]][["proj"]][["y_hat"]] %>% 
      apply(1, function(x)(x = mean(x, na.rm = TRUE)))
    r_m <- r %>% 
      setValues(val)
    
    # --- 2.2.2. Coefficient of variation
    val <- MODEL[[i]][["proj"]][["y_hat"]] %>% 
      apply(1, function(x)(x = cv(x, na.rm = TRUE)))
    r_cv <- r %>% 
      setValues(val)
    r_cv[r_cv > 100] <- 100
    
    # --- 2.3. Plot the corresponding maps
    plot(r_m, col = viridis_pal(100),
         main = paste("Average proj. for", i, "\n", 
                      names(MODEL[[i]][["eval"]]), "=", MODEL[[i]][["eval"]])
         )
    
    plot(r_cv, col = viridis_pal(100),
         main = paste("Uncertainty (CV) proj. for", i))
    
  } # End i model loop
  
  # --- 3. Build ensemble outputs
  if(ENSEMBLE == TRUE){
    # --- 3.1. Set initial plot layout & requirements
    par(mfrow = c(2,1), mar = c(2,8,3,3))
    r <- raster(paste0(project_wd, "/data/features_mean_from_monthly"))
    
    # --- 3.2. Build ensemble array, weighted by evaluation values, re-scale max=1
    y_ens <- NULL
    for(i in MODEL$CALL$MODEL_LIST){
      y_ens <- cbind(y_ens,
                     MODEL[[i]][["proj"]][["y_hat"]]*MODEL[[i]][["eval"]][[1]])
    }
    y_ens <- y_ens/max(y_ens, na.rm = TRUE)
    
    # --- 3.3. Compute mean and CV
    # --- 3.3.1. Mean value
    val <- apply(y_ens, 1, function(x)(x = mean(x, na.rm = TRUE)))
    r_m <- r %>% 
      setValues(val)
    
    # --- 3.3.2. Coefficient of variation
    val <- apply(y_ens, 1, function(x)(x = cv(x, na.rm = TRUE)))
    r_cv <- r %>% 
      setValues(val)
    r_cv[r_cv > 100] <- 100
    
    # --- 3.4. Plot the corresponding maps
    plot(r_m, col = viridis_pal(100),
         main = "Average ensemble proj.")
    
    plot(r_cv, col = viridis_pal(100),
         main = "Average ensemble uncertainty (CV)")
  } # End if ENSEMBLE = TRUE
  
  # --- 4. Wrap up and save
  # --- 4.1. Stop PDF saving
  dev.off()
  
} # END FUNCTION

