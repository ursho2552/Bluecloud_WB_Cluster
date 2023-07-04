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
    stop("No validated algorithms to display projections from")
  }
  
  # --- 1.3. Create PDF saving
  pdf(paste0(project_wd,"/output/",FOLDER_NAME,"/",SUBFOLDER_NAME,"/05_standard_maps.pdf"))
  
  # --- 1.4. Set initial plot layout & requirements
  par(mfrow = c(4,3), mar = c(2,2,4,1))
  r0 <- raster(paste0(project_wd, "/data/features_mean_from_monthly"))
  
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
  if(CALL$DATA_SOURCE != "omic"){mtext(paste("QC for", QUERY$annotations$scientificname))}
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
  colorbar.plot(x = 0.5, y = 0, strip = seq(0,1,length.out = 100),
                strip.width = 0.3, strip.length = 2.7,
                col = hsi_pal, border = "black")
  axis(side = 1)
  text(x = 0.5, y = 0.3, "Habitat Suitability Index", adj = 0.5)
  abline(h = -0.5) # PDF wide separator
  # --- 3.2. Observation vs 75% quartile
  plot.new()
  points(x = 0.1, y = 0.4, pch = 22, col = "black", bg = "gray80", cex = 5)
  text(x = 0.2, y = 0.4, "Q75 Habitat Suitability Index", pos = 4)
  points(x = 0.1, y = 0.1, pch = 3, col = "red", cex = 2)
  text(x = 0.2, y = 0.1, "Observation", pos = 4)
  # --- 3.3. MESS x CV 2D color scale
  par(mar = c(5,5,3,2), xpd = FALSE)
  bivar_pal <- colmat(xmax = "deepskyblue4", ymax = "darkgoldenrod2", nbreaks = 100)
  colmat_plot(bivar_pal, xlab = "Coefficient of variation", ylab = "MESS value")
  axis(side = 1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c(0, 0.2, 0.4, 0.6, 0.8, 1))
  axis(side = 2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c(0, -20, -40, -60, -80, -100), las = 2)
  par(mar = c(2,2,4,1))

  # --- 4. Build model-level outputs
  for(i in MODEL$CALL$MODEL_LIST){
    # --- 4.1. Compute the different layers
    # --- 4.1.1. Mean value
    val <- MODEL[[i]][["proj"]][["y_hat"]] %>%
      apply(1, function(x)(x = mean(x, na.rm = TRUE)))
    r_m <- r0 %>%
      setValues(val)

    # --- 4.1.2. Coefficient of variation
    val <- MODEL[[i]][["proj"]][["y_hat"]] %>%
      apply(1, function(x)(x = cv(x, na.rm = TRUE)))
    r_cv <- r0 %>%
      setValues(val)
    r_cv[r_cv > 100] <- 100
    r_cv[r_cv <= 0] <- 1e-10 # temporary fix : modify bivarmap so that 0 is included

    # --- 4.1.3. MESS
    r_mess <- QUERY$MESS*-1
    r_mess[r_mess<0] <- 1e-10 # temporary fix : modify bivarmap so that 0 is included
    r_mess[r_mess>100] <- 100

    # --- 4.2. Plot the corresponding maps
    # --- 4.2.1. Plot the abundance
    # Abundance or habitat suitability values
    plot(r_m, col = hsi_pal[max(1, floor(r_m@data@min*100)):min(100, ceiling(r_m@data@max*100))], legend=FALSE,
         main = paste("Average proj. for", i, names(MODEL[[i]][["eval"]][[1]]), "\n", 
                      names(MODEL[[i]][["eval"]])[1], "=", MODEL[[i]][["eval"]][[1]]))
    
    # Land mask
    plot(land, col = "antiquewhite4", legend=FALSE, add = TRUE)

    # --- 4.2.2. Plot the observations
    # Top abundance quartile as contour
    plot(r_m > quantile(r_m, 0.75), col = c("white","gray80"), legend=FALSE, main = "Observations")
    # Land mask
    plot(land, col = "antiquewhite4", legend=FALSE, add = TRUE)
    # Observations
    tmp <- QUERY$S[which(QUERY$Y$measurementvalue > 0),]
    points(tmp$decimallongitude, tmp$decimallatitude,
           col = "red", pch = 3)

    # --- 4.2.3. Plot the uncertainties
    # Land mask
    plot(land, col = "antiquewhite4", legend=FALSE, main = "Uncertainties")
    # CV x MESS plot
    r <- bivar_map(rasterx = r_cv, rastery = r_mess, colormatrix = bivar_pal,
                   cutx = 0:100, cuty = 0:100)
    plot(r[[1]], col = r[[2]], legend=FALSE, add = TRUE)

  } # End i model loop

  # --- 5. Build ensemble outputs
  if(ENSEMBLE == TRUE){
    # --- 5.1. Build ensemble array, weighted by evaluation values, re-scale max=1
    y_ens <- NULL
    for(i in MODEL$CALL$MODEL_LIST){
      y_ens <- cbind(y_ens,
                     MODEL[[i]][["proj"]][["y_hat"]])
    }
    # --- 5.2. Compute the different layers
    # --- 5.2.1. Mean value
    val <- apply(y_ens, 1, function(x)(x = mean(x, na.rm = TRUE)))
    r_m <- r0 %>%
      setValues(val)

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
    # Top abundance quartile as contour
    plot(r_m > quantile(r_m, 0.75), col = c("white","gray80"), legend=FALSE, main = "Observations")
    # Land mask
    plot(land, col = "antiquewhite4", legend=FALSE, add = TRUE)
    # Observations
    tmp <- QUERY$S[which(QUERY$Y$measurementvalue > 0),]
    points(tmp$decimallongitude, tmp$decimallatitude,
           col = "red", pch = 3)

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

