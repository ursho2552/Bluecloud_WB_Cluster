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
  pdf(paste0(project_wd,"/output/",FOLDER_NAME,"/",SUBFOLDER_NAME,"/standard_maps2.pdf"))
  
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
    
    # --- 2.2.4. MESS
    r_mess <- QUERY$MESS*-1
    r_mess[r_mess<0] <- NA
    r_mess[r_mess>100] <- 100
    
    # --- 2.2.5. Landmask
    land <- r_cv
    land[is.na(land)] <- 9999
    land[land != 9999] <- NA
    
    # --- 2.3. Construct color palettes
    cv_pal <- alpha(colorRampPalette(brewer.pal(9,"Blues"))(100), seq(0,0.5,length.out = 100))
    mess_pal <- alpha(colorRampPalette(brewer.pal(9,"Reds"))(100), seq(0,0.5,length.out = 100))
    
    # --- 2.4. Plot the corresponding maps
    # --- 2.4.1. Plot the abundance
    plot(r_m, col = viridis_pal(100),
         main = paste("Average proj. for", i, "\n", 
                      names(MODEL[[i]][["eval"]]), "=", MODEL[[i]][["eval"]])
    )
    
    # --- 2.4.2. Plot the CV
    plot(r_cv, col = cv_pal, main = paste("Uncertainties (CV, MESS, 80%) proj. for", i))
    plot(r_mess, col = mess_pal, add = TRUE)
    
    contour(r_m, levels =quantile(r_m, 0.8), labels = "", lwd = 2, col = "black", add = TRUE)
    plot(land, col = "antiquewhite3", add = TRUE)
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
    # --- 3.4.1. Plot the abundance
    plot(r_m, col = viridis_pal(100),
         main = paste("Average proj. for", i, "\n", 
                      names(MODEL[[i]][["eval"]]), "=", MODEL[[i]][["eval"]])
    )
    
    # --- 3.4.2. Plot the CV
    plot(r_cv, col = cv_pal, main = paste("Uncertainties (CV, MESS, 80%) proj. for", i))
    plot(r_mess, col = mess_pal, add = TRUE)
    
    contour(r_m, levels =quantile(r_m, 0.8), labels = "", lwd = 2, col = "black", add = TRUE)
    plot(land, col = "antiquewhite3", add = TRUE)
  } # End if ENSEMBLE = TRUE
  
  # --- 4. Wrap up and save
  # --- 4.1. Stop PDF saving
  dev.off()
  
} # END FUNCTION

