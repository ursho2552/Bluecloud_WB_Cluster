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
  
  # --- 1.3. Set initial plot layout & requirements
  par(mfrow = c(4,3), mar = c(2,2,4,1))
  r0 <- raster(paste0(project_wd, "/data/features_mean_from_monthly"))
  
  # --- 1.4 Land mask
  land <- r0
  land[is.na(land)] <- 9999
  land[land != 9999] <- NA
  
  # --- 2. Plot the legends
  # --- 2.1. Abundance or habitat suitability legend
  hsi_pal <- inferno_pal(100)
  plot.new()
  colorbar.plot(x = 0.5, y = 0, strip = seq(0,1,length.out = 100),
                strip.width = 0.3, strip.length = 2.7,
                col = hsi_pal, border = "black")
  axis(side = 1)
  text(x = 0.5, y = 0.3, "Habitat Suitability Index", adj = 0.5)
  # --- 2.2. Observation vs 75% quartile
  plot.new()
  points(x = 0.1, y = 0.4, pch = 22, col = "black", bg = "gray80", cex = 5)
  text(x = 0.2, y = 0.4, "75% Habitat Suitability Index Quantile", pos = 4)
  points(x = 0.1, y = 0.1, pch = 3, col = "red", cex = 2)
  text(x = 0.2, y = 0.1, "Observation", pos = 4)
  # --- 2.3. MESS x CV 2D color scale
  par(mar = c(5,5,3,2))
  bivar_pal <- colmat(xmax = "deepskyblue4", ymax = "darkgoldenrod2", nbreaks = 100)
  colmat_plot(bivar_pal, xlab = "Coefficient of variation", ylab = "MESS value")
  axis(side = 1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c(0, 0.2, 0.4, 0.6, 0.8, 1))
  axis(side = 2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c(0, -20, -40, -60, -80, -100), las = 2)
  par(mar = c(2,2,4,1))
  
  # --- 2. Build model-level outputs
  for(i in MODEL$CALL$MODEL_LIST){
    # --- 2.1. Compute the different layers
    # --- 2.1.1. Mean value
    val <- MODEL[[i]][["proj"]][["y_hat"]] %>% 
      apply(1, function(x)(x = mean(x, na.rm = TRUE)))
    r_m <- r0 %>% 
      setValues(val)
    
    # --- 2.1.2. Coefficient of variation
    val <- MODEL[[i]][["proj"]][["y_hat"]] %>% 
      apply(1, function(x)(x = cv(x, na.rm = TRUE)))
    r_cv <- r0 %>% 
      setValues(val)
    r_cv[r_cv > 100] <- 100
    
    # --- 2.1.3. MESS
    r_mess <- QUERY$MESS*-1
    r_mess[r_mess<0] <- NA
    r_mess[r_mess>100] <- 100
    
    # --- 2.2. Plot the corresponding maps
    # --- 2.2.1. Plot the abundance
    # Abundance or habitat suitability values
    plot(r_m, col = hsi_pal, legend=FALSE,
         main = paste("Average proj. for", i, names(MODEL[[i]][["eval"]][[1]])))
    # Land mask
    plot(land, col = "antiquewhite4", legend=FALSE, add = TRUE)
    
    # --- 2.2.2. Plot the observations
    # Top abundance quartile as contour
    plot(r_m > quantile(r_m, 0.75), col = c("white","gray80"), legend=FALSE, main = "Observations")
    # Land mask
    plot(land, col = "antiquewhite4", legend=FALSE, add = TRUE)
    # Observations
    tmp <- QUERY$S[which(QUERY$Y$measurementvalue > 0),]
    points(tmp$decimallongitude, tmp$decimallatitude,
           col = "red", pch = 3)
    
    # --- 2.2.3. Plot the uncertainties
    # Land mask
    plot(land, col = "antiquewhite4", legend=FALSE, main = "Uncertainties")
    # CV x MESS plot
    r <- bivar_map(rasterx = r_cv, rastery = r_mess, colormatrix = bivar_pal,
                   cutx = 0:100, cuty = 0:100)
    plot(r[[1]], col = r[[2]], legend=FALSE, add = TRUE)
  
  } # End i model loop
  
  # --- 3. Build ensemble outputs
  if(ENSEMBLE == TRUE){
    # --- 3.1. Build ensemble array, weighted by evaluation values, re-scale max=1
    y_ens <- NULL
    for(i in MODEL$CALL$MODEL_LIST){
      y_ens <- cbind(y_ens,
                     MODEL[[i]][["proj"]][["y_hat"]]*MODEL[[i]][["eval"]][[1]])
    }
    y_ens <- y_ens/max(y_ens, na.rm = TRUE)
    # --- 3.2. Compute the different layers
    # --- 3.2.1. Mean value
    val <- apply(y_ens, 1, function(x)(x = mean(x, na.rm = TRUE)))
    r_m <- r0 %>% 
      setValues(val)
    
    # --- 3.2.2. Coefficient of variation
    val <- apply(y_ens, 1, function(x)(x = cv(x, na.rm = TRUE)))
    r_cv <- r0 %>% 
      setValues(val)
    r_cv[r_cv > 100] <- 100
    
    # --- 3.3. Plot the corresponding maps
    # --- 3.3.1. Plot the abundance
    # Abundance or habitat suitability values
    plot(r_m, col = hsi_pal, legend=FALSE,
         main = "Average Ensemble proj.")
    # Land mask
    plot(land, col = "antiquewhite4", legend=FALSE, add = TRUE)
    
    # --- 3.3.2. Plot the observations
    # Top abundance quartile as contour
    plot(r_m > quantile(r_m, 0.75), col = c("white","gray80"), legend=FALSE, main = "Observations")
    # Land mask
    plot(land, col = "antiquewhite4", legend=FALSE, add = TRUE)
    # Observations
    tmp <- QUERY$S[which(QUERY$Y$measurementvalue > 0),]
    points(tmp$decimallongitude, tmp$decimallatitude,
           col = "red", pch = 3)
    
    # --- 3.3.3. Plot the uncertainties
    # Land mask
    plot(land, col = "antiquewhite4", legend=FALSE, main = "Uncertainties")
    # CV x MESS plot
    r <- bivar_map(rasterx = r_cv, rastery = r_mess, colormatrix = bivar_pal,
                   cutx = 0:100, cuty = 0:100)
    plot(r[[1]], col = r[[2]], legend=FALSE, add = TRUE)
  } # End if ENSEMBLE = TRUE
  
  # --- 4. Wrap up and save
  # --- 4.1. Stop PDF saving
  dev.off()
  
} # END FUNCTION

