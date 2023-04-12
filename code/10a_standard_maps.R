#' =============================================================================
#' @name standard_maps
#' @description Simple function for building standard mean and CV maps from
#' projection data
#' @param QUERY the query object from the master pipeline
#' @param MODELS the models object from the master pipeline
#' @param ENSEMBLE if TRUE, computes the ensemble map ? 
#' @param MESS if TRUE, considers the mess analysis in the maps
#' @return plots mean and uncertainty maps per model or ensemble

standard_maps <- function(QUERY = query,
                          MODELS = models,
                          ENSEMBLE = FALSE,
                          MESS = FALSE){
  
  # --- 1. Build model-level outputs
  # --- 1.1. Set initial plot layout & requirements
  par(mfrow = c(4,2), mar = c(2,8,3,3))
  r <- raster(paste0(project_wd, "/data/features_mean_from_monthly"))
  
  for(i in MODELS$CALL$MODEL_LIST){
    # --- 1.2. Compute mean and CV
    # --- Mean value
    val <- MODELS[[i]][["proj"]][["y_hat"]] %>% 
      apply(1, function(x)(x = mean(x, na.rm = TRUE)))
    r_m <- r %>% 
      setValues(val)
    
    # --- Coefficient of variation
    val <- MODELS[[i]][["proj"]][["y_hat"]] %>% 
      apply(1, function(x)(x = cv(x, na.rm = TRUE)))
    r_cv <- r %>% 
      setValues(val)
    
    # --- 1.3. Plot the corresponding maps
    plot(r_m, col = viridis_pal(100),
         main = paste("Average proj. for", i, "\n", 
                      names(MODELS[[i]][["eval"]]), "=", MODELS[[i]][["eval"]])
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
    for(i in MODELS$CALL$MODEL_LIST){
      y_ens <- cbind(y_ens,
                     MODELS[[i]][["proj"]][["y_hat"]]*MODELS[[i]][["eval"]][[1]])
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
    
    # --- 2.4. Plot the corresponding maps
    plot(r_m, col = viridis_pal(100),
         main = "Average ensemble proj.")
    
    plot(r_cv, col = viridis_pal(100),
         main = "Average ensemble uncertainty (CV)")
  } # End if ENSEMBLE = TRUE
  
} # END FUNCTION


# ==================
# PLOT PROTOTYPE
i = "MLP"
MODELS = models

par(mfrow = c(2,1), mar = c(2,2,3,5))
# Mean value
val <- MODELS[[i]][["proj"]][["y_hat_m"]]
r <- raster(paste0(project_wd, "/data/features_mean_from_monthly")) %>% 
  setValues(val)
plot(r, main = "mean", col = viridis_pal(100))

# CV value
val <- MODELS[[i]][["proj"]][["y_hat_cv"]]
r <- raster(paste0(project_wd, "/data/features_mean_from_monthly")) %>% 
  setValues(val)
plot(r, main = "cv", col = viridis_pal(100))

