#' =============================================================================
#' @name standard_maps
#' @description Simple function for building standard mean and CV maps from
#' projection data
#' @param FOLDER_NAME name of the corresponding folder
#' @param SUBFOLDER_NAME list of sub_folders to parallelize on.
#' @param ENSEMBLE if TRUE, computes the ensemble map ?
#' @return plots mean and uncertainty maps per model or ensemble

standard_maps <- function(FOLDER_NAME = NULL,
                          SUBFOLDER_NAME = NULL,
                          ENSEMBLE = FALSE){
  
  # --- 1. Initialize function
  # --- 1.1. Parameter loading
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/QUERY.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/MODEL.RData"))
  
  # --- 1.2. Check for projections
  if(length(MODEL$CALL$MODEL_LIST) == 0){
    message("No validated algorithms to display projections from")
    return(NULL)
  }
  
  # --- 1.3. Create PDF saving
  pdf(paste0(project_wd,"/output/",FOLDER_NAME,"/",SUBFOLDER_NAME,"/05_standard_maps.pdf"))
  
  # --- 1.4. Set initial plot layout & requirements
  par(mfrow = c(4,3), mar = c(2,2,4,1))
  r0 <- raster(CALL$ENV_PATH)
  
  # --- 1.5 Land mask
  land <- r0
  land[is.na(land)] <- 9999
  land[land != 9999] <- NA
  
  # --- 2. Plot the quality checks
  # --- 2.1. Compute the recommendation table
  rec <- qc_recommandations(MODEL = MODEL)
  traffic_col <- rep(rec$COL, each = 3)
  traffic_val <- rec[,1:3] %>% as.matrix() %>% t() %>% c()
  traffic_col[which(traffic_val == 0)] <- "white"
  # --- 2.2. Plot the traffic lights and recommendations
  plot.new()
  if(CALL$DATA_TYPE != "proportions"){mtext(paste("QC for", QUERY$annotations$scientificname))}
  par(mar = c(1,5,3,1), xpd = NA)
  plot(x = rep(1:3, nrow(rec)), y = rep(nrow(rec):1, each = 3), axes = FALSE, cex = 4,
       xlim = c(0,4), ylim = c(0,nrow(rec)+1), ylab = "", xlab = "",
       pch = 21, col = "black", bg = traffic_col)
  axis(side = 3, at = 1:3, labels = c("FIT","VIP","DEV"), tick = FALSE, line = NA, cex.axis = 1)
  axis(side = 2, at = nrow(rec):1, labels = rownames(rec), tick = FALSE, line = NA, las = 2, cex.axis = 1)
  axis(side = 4, at = nrow(rec):1, labels = rec$Recommandation, tick = FALSE, line = NA, las = 2, cex.axis = 1)
  abline(h = 0) # PDF wide separator
  plot.new()
  par(mar = c(2,2,4,1))

  # --- 3. Plot the legends
  # --- 3.1. Abundance or habitat suitability legend
  hsi_pal <- inferno_pal(100)
  plot.new()
  # --- 3.1.1. Extract the plot true scale (average max across bootstrap)
  plot_scale <- lapply(MODEL, FUN = function(z)(z = apply(MODEL$GLM$proj$y_hat, 1, function(x)(x = mean(x, na.rm = TRUE))) %>% max(na.rm = TRUE))) %>% 
    unlist() %>% 
    max(na.rm = TRUE)
  # --- 3.1.2. Plot the colorbar
  colorbar.plot(x = 0.5, y = 0, strip = seq(0,1,length.out = 100),
                strip.width = 0.3, strip.length = 2.7,
                col = hsi_pal, border = "black")
  # --- 3.1.3. Apply the scale to the axis caption
  # This is only informative and the raster will be rescaled by the maximum
  axis(side = 1, at = seq(0, 1, length.out = 5), labels = round(seq(0, 1 * plot_scale, length.out = 5), 2))
  text(x = 0.5, y = 0.3, "Habitat Suitability Index", adj = 0.5)
  abline(h = -0.5) # PDF wide separator
  # --- 3.2. Observation vs 75% quartile
  plot.new()
  points(x = 0.1, y = 0.4, pch = 22, col = "black", bg = "gray80", cex = 5)
  text(x = 0.2, y = 0.4, "Q75 Habitat Suitability Index", pos = 4)
  points(x = 0.1, y = 0.1, pch = 20, col = "black", cex = 2)
  text(x = 0.2, y = 0.1, "Observation", pos = 4)
  # --- 3.3. MESS x CV 2D color scale
  par(mar = c(5,5,3,2), xpd = FALSE)
  bivar_pal <- colmat(nbreaks = 100)
  colmat_plot(bivar_pal, xlab = "Coefficient of variation", ylab = "MESS value")
  axis(side = 1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c(0, 0.2, 0.4, 0.6, 0.8, 1))
  axis(side = 2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c(0, -20, -40, -60, -80, -100), las = 2)
  par(mar = c(2,2,4,1))

  # --- 4. Build model-level outputs
  # --- 4.1. Initialize loop
  if(CALL$DATA_TYPE == "proportions"){loop_over <- 1:ncol(QUERY$Y)
  }else{loop_over <- MODEL$CALL$MODEL_LIST}
  
  for(i in loop_over){
    # --- 4.2. Compute the different layers
    # --- 4.2.1. Extract projection data for the iteration
    if(CALL$DATA_TYPE == "proportions"){val_raw <- MODEL[["MBTR"]][["proj"]][["y_hat"]][,,i]
    }else{val_raw <- MODEL[[i]][["proj"]][["y_hat"]]}
    
    # --- 4.2.2. Mean value
    # Rescaled by the maximum to match the colorbar
    val <- apply(val_raw, 1, function(x)(x = mean(x, na.rm = TRUE)))
    r_m <- r0 %>% setValues(val / plot_scale)

    # --- 4.2.3. Coefficient of variation
    val <- apply(val_raw, 1, function(x)(x = cv(x, na.rm = TRUE)))
    r_cv <- r0 %>% setValues(val)
    r_cv[r_cv > 100] <- 100
    r_cv[r_cv <= 0] <- 1e-10 # temporary fix : modify bivarmap so that 0 is included

    # --- 4.2.4. MESS
    r_mess <- QUERY$MESS*-1
    r_mess[r_mess<0] <- 1e-10 # temporary fix : modify bivarmap so that 0 is included
    r_mess[r_mess>100] <- 100

    # --- 4.3. Plot the corresponding maps
    # --- 4.3.1. Plot the abundance
    if(CALL$DATA_TYPE == "proportions"){
      # Proportions
      tmp <- which(QUERY$annotations$worms_id == colnames(QUERY$Y)[i])
      plot(r_m, col = hsi_pal[max(1, floor(r_m@data@min*100)):min(100, ceiling(r_m@data@max*100))], legend=FALSE,
           main = paste("Average proj. for", QUERY$annotations$scientificname[tmp]))
    } else {
      # Abundance or habitat suitability values
      plot(r_m, col = hsi_pal[max(1, floor(r_m@data@min*100)):min(100, ceiling(r_m@data@max*100))], legend=FALSE,
           main = paste("Average proj. for", i, names(MODEL[[i]][["eval"]][[1]]), "\n", 
                        names(MODEL[[i]][["eval"]])[1], "=", MODEL[[i]][["eval"]][[1]]))
    }
    
    # Land mask
    plot(land, col = "antiquewhite4", legend=FALSE, add = TRUE)

    # --- 4.3.2. Plot the observations
    # --- 4.3.2.1. Top abundance quartile as contour
    plot(r_m > quantile(r_m, 0.75), col = c("white","gray80"), legend=FALSE, main = "Observations")
    # --- 4.3.2.2. Land mask
    plot(land, col = "antiquewhite4", legend=FALSE, add = TRUE)
    # --- 4.3.2.3. Observations location
    if(CALL$DATA_TYPE == "proportions"){tmp <- QUERY$S
    } else {tmp <- QUERY$S[which(QUERY$Y$measurementvalue > 0),]}
    # --- 4.3.2.4. Observation colors
    if(CALL$DATA_TYPE == "continuous"){
      points(tmp$decimallongitude, tmp$decimallatitude,
             col = col_numeric("inferno", domain = range(QUERY$Y$measurementvalue))(QUERY$Y$measurementvalue), pch = 20)
    } else {
      points(tmp$decimallongitude, tmp$decimallatitude,
             col = "black", pch = 20)
    }
    
    # --- 4.3.3. Plot the uncertainties
    # Land mask
    plot(land, col = "antiquewhite4", legend=FALSE, main = "Uncertainties")
    # CV x MESS plot
    r <- bivar_map(rasterx = r_cv, rastery = r_mess, colormatrix = bivar_pal,
                   cutx = 0:100, cuty = 0:100)
    plot(r[[1]], col = r[[2]], legend=FALSE, add = TRUE)

  } # End i model loop

  # --- 5. Build ensemble outputs
  if(ENSEMBLE == TRUE & CALL$DATA_TYPE != "proportions"){
    # --- 5.1. Build ensemble array, weighted by evaluation values, re-scale max=1
    y_ens <- NULL
    for(i in MODEL$CALL$MODEL_LIST){
      y_ens <- cbind(y_ens,
                     MODEL[[i]][["proj"]][["y_hat"]])
    }
    # --- 5.2. Compute the different layers
    # --- 5.2.1. Mean value
    # Rescaled by the maximum to match the colorbar
    val <- apply(y_ens, 1, function(x)(x = mean(x, na.rm = TRUE)))
    r_m <- r0 %>%
      setValues(val / plot_scale)

    # --- 5.2.2. Coefficient of variation
    val <- apply(y_ens, 1, function(x)(x = cv(x, na.rm = TRUE)))
    r_cv <- r0 %>%
      setValues(val)
    r_cv[r_cv > 100] <- 100
    r_cv[r_cv <= 0] <- 1e-10 # temporary fix : modify bivarmap so that 0 is included

    # --- 5.3. Plot the corresponding maps
    # --- 5.3.1. Plot the abundance
    # Abundance or habitat suitability values
    plot(r_m, col = hsi_pal[max(1, floor(r_m@data@min*100)):min(100, ceiling(r_m@data@max*100))], legend=FALSE,
         main = "Average Ensemble proj.")
    # Land mask
    plot(land, col = "antiquewhite4", legend=FALSE, add = TRUE)

    # --- 5.3.2. Plot the observations
    # --- 5.3.2.1. Top abundance quartile as contour
    plot(r_m > quantile(r_m, 0.75), col = c("white","gray80"), legend=FALSE, main = "Observations")
    # --- 5.3.2.2. Land mask
    plot(land, col = "antiquewhite4", legend=FALSE, add = TRUE)
    # --- 5.3.2.3. Observations 
    tmp <- QUERY$S[which(QUERY$Y$measurementvalue > 0),]
    if(CALL$DATA_TYPE == "continuous"){
      points(tmp$decimallongitude, tmp$decimallatitude,
             col = col_numeric("inferno", domain = range(QUERY$Y$measurementvalue))(QUERY$Y$measurementvalue), pch = 20)
    } else {
      points(tmp$decimallongitude, tmp$decimallatitude,
             col = "black", pch = 20)
    }

    # --- 5.3.3. Plot the uncertainties
    # Land mask
    plot(land, col = "antiquewhite4", legend=FALSE, main = "Uncertainties")
    # CV x MESS plot
    r <- bivar_map(rasterx = r_cv, rastery = r_mess, colormatrix = bivar_pal,
                   cutx = 0:100, cuty = 0:100)
    plot(r[[1]], col = r[[2]], legend=FALSE, add = TRUE)
  } # End if ENSEMBLE = TRUE

  # --- 6. Wrap up and save
  # --- 6.1. Stop PDF saving
  dev.off()
  
} # END FUNCTION

